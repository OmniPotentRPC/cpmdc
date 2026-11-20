#define _GNU_SOURCE
/**
 * CPMDC_ORBITAL_SLOT_LIMIT bounds parked orbital keys.
 *
 * The reference evaluator has no Fortran store. A keyed evaluation still
 * parks the key it leaves, and a park past the limit fails the call. The
 * result message and cpmdc_last_error() name the variable. Keys already
 * parked stay usable. Each case runs in a child so the process-global
 * ledger and the environment do not leak into the next case.
 */
#include "cpmdc.h"
#include "setup_child.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

static const char *g_params = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
/* Distinct from every other step in this file, and far enough that the
 * session cache does not reuse an earlier geometry. */
static const double k_geometry_step = 0.25;

static const char k_limit_name[] = "CPMDC_ORBITAL_SLOT_LIMIT";
static const char k_exceeded[] = "exceeded";
static const char k_unchanged[] = "unchanged";

enum { k_msg = 512, k_keys = 4 };

struct slot_result {
  int ok[k_keys];
  char message[k_keys][k_msg];
  char last_error[k_keys][k_msg];
};

struct slot_ctx {
  const unsigned char *params;
  size_t params_n;
  const char *limit;
  int n_keys;
  long long keys[k_keys];
};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  long n;
  unsigned char *buf;
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  buf = (unsigned char *)malloc((size_t)n);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void positions_for(int step, double out[9]) {
  memcpy(out, k_pos, sizeof(k_pos));
  out[0] = k_pos[0] + k_geometry_step * (double)step;
}

static void copy_text(char *dst, size_t dst_n, const char *src) {
  if (!src)
    src = "";
  snprintf(dst, dst_n, "%s", src);
}

static void eval_keys(void *ctx_void, void *out_void) {
  struct slot_ctx *ctx = ctx_void;
  struct slot_result *out = out_void;
  CPMDCSession *session;
  int step;
  if (ctx->limit)
    setenv(k_limit_name, ctx->limit, 1);
  else
    unsetenv(k_limit_name);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session) {
    copy_text(out->last_error[0], k_msg, cpmdc_last_error());
    return;
  }
  for (step = 0; step < ctx->n_keys; ++step) {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(step, pos);
    if (cpmdc_session_select_orbitals(session, ctx->keys[step]) != 0) {
      copy_text(out->message[step], k_msg, "select failed");
      copy_text(out->last_error[step], k_msg, cpmdc_last_error());
      break;
    }
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->ok[step] = r.ok;
    copy_text(out->message[step], k_msg, r.message);
    copy_text(out->last_error[step], k_msg, cpmdc_last_error());
  }
  cpmdc_session_destroy(session);
}

static int run_case(const unsigned char *params, size_t params_n, const char *limit,
                    const long long *keys, int n_keys, struct slot_result *out) {
  struct slot_ctx ctx;
  char why[256];
  memset(out, 0, sizeof(*out));
  memset(&ctx, 0, sizeof(ctx));
  ctx.params = params;
  ctx.params_n = params_n;
  ctx.limit = limit;
  ctx.n_keys = n_keys;
  memcpy(ctx.keys, keys, (size_t)n_keys * sizeof(*keys));
  return run_setup_child(eval_keys, &ctx, out, sizeof(*out), why, sizeof(why));
}

static void test_limit_rejects_a_new_park(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  const long long keys[4] = {11, 12, 13, 11};
  struct slot_result out;
  (void)state;
  assert_non_null(params);
  assert_int_equal(run_case(params, n, "1", keys, 4, &out), 0);
  assert_int_equal(out.ok[0], 1);
  assert_int_equal(out.ok[1], 1);
  assert_int_equal(out.ok[2], 0);
  assert_non_null(strstr(out.message[2], k_limit_name));
  assert_non_null(strstr(out.message[2], k_exceeded));
  assert_non_null(strstr(out.message[2], k_unchanged));
  assert_non_null(strstr(out.last_error[2], k_limit_name));
  assert_non_null(strstr(out.last_error[2], k_exceeded));
  assert_non_null(strstr(out.last_error[2], k_unchanged));
  assert_int_equal(out.ok[3], 1);
  assert_int_equal(out.last_error[3][0], '\0');
  free(params);
}

/* One child evaluates three keys under a wide limit, then lowers it. */
struct lower_ctx {
  const unsigned char *params;
  size_t params_n;
};

struct lower_result {
  int wide_ok;
  int kept_ok;
  int fresh_ok;
  char fresh_message[k_msg];
  char fresh_error[k_msg];
  char kept_error[k_msg];
};

static void eval_lower(void *ctx_void, void *out_void) {
  struct lower_ctx *ctx = ctx_void;
  struct lower_result *out = out_void;
  const long long wide_keys[3] = {21, 22, 23};
  CPMDCSession *session;
  int step;
  setenv(k_limit_name, "3", 1);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session)
    return;
  out->wide_ok = 1;
  for (step = 0; step < 3; ++step) {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(10 + step, pos);
    if (cpmdc_session_select_orbitals(session, wide_keys[step]) != 0) {
      out->wide_ok = 0;
      break;
    }
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    if (!r.ok)
      out->wide_ok = 0;
  }
  setenv(k_limit_name, "1", 1);
  {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(20, pos);
    cpmdc_session_select_orbitals(session, 22);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->kept_ok = r.ok;
    copy_text(out->kept_error, k_msg, cpmdc_last_error());
  }
  {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(21, pos);
    cpmdc_session_select_orbitals(session, 24);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->fresh_ok = r.ok;
    copy_text(out->fresh_message, k_msg, r.message);
    copy_text(out->fresh_error, k_msg, cpmdc_last_error());
  }
  cpmdc_session_destroy(session);
}

static void test_lowering_keeps_stored_keys(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct lower_ctx ctx;
  struct lower_result out;
  char why[256];
  (void)state;
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  assert_int_equal(
      run_setup_child(eval_lower, &ctx, &out, sizeof(out), why, sizeof(why)), 0);
  assert_int_equal(out.wide_ok, 1);
  assert_int_equal(out.kept_ok, 1);
  assert_int_equal(out.kept_error[0], '\0');
  assert_int_equal(out.fresh_ok, 0);
  assert_non_null(strstr(out.fresh_message, k_limit_name));
  assert_non_null(strstr(out.fresh_message, k_unchanged));
  assert_non_null(strstr(out.fresh_error, k_limit_name));
  free(params);
}

static void test_invalid_limit_names_the_variable(void **state) {
  static const char *values[] = {"0", "-1", "", "abc", "999999999999999999999",
                                 " "};
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  size_t i;
  (void)state;
  assert_non_null(params);
  for (i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
    const long long key = 31;
    struct slot_result out;
    assert_int_equal(run_case(params, n, values[i], &key, 1, &out), 0);
    assert_int_equal(out.ok[0], 0);
    assert_non_null(strstr(out.message[0], k_limit_name));
    assert_non_null(strstr(out.last_error[0], k_limit_name));
    assert_null(strstr(out.message[0], k_exceeded));
  }
  free(params);
}

static void test_unset_limit_parks_every_key(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  const long long keys[3] = {41, 42, 43};
  struct slot_result out;
  (void)state;
  assert_non_null(params);
  assert_int_equal(run_case(params, n, NULL, keys, 3, &out), 0);
  assert_int_equal(out.ok[0], 1);
  assert_int_equal(out.ok[1], 1);
  assert_int_equal(out.ok[2], 1);
  assert_int_equal(out.last_error[2][0], '\0');
  free(params);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_limit_rejects_a_new_park),
      cmocka_unit_test(test_lowering_keeps_stored_keys),
      cmocka_unit_test(test_invalid_limit_names_the_variable),
      cmocka_unit_test(test_unset_limit_parks_every_key),
  };
  if (argc < 2)
    return 1;
  g_params = argv[1];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
