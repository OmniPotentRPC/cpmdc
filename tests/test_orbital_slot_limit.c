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
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

#include <cmocka.h>

#if defined(__SANITIZE_ADDRESS__)
#define CPMDC_TEST_ASAN 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define CPMDC_TEST_ASAN 1
#endif
#endif

/* Slots the reference ledger holds for parked keys. Zero on a CPMD build,
 * whose store lives in Fortran. Exported by the embed layer for tests. */
size_t cpmdc_orbital_ledger_capacity(void);

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

/* Two keys that alternate, as the images of a band do between iterations.
 * Each switch parks the key it leaves and restores the other, so at most
 * one key is parked and the table never has to grow past its first block. */
enum {
  k_switches = 1000,
  k_first_block = 16,
  k_address_space_margin = 64 * 1024 * 1024,
  k_bytes_per_kib = 1024
};

struct alternate_ctx {
  const unsigned char *params;
  size_t params_n;
  int limit_address_space;
};

struct alternate_result {
  int switches_ok;
  int rlimit_set;
  size_t capacity;
  char first_error[k_msg];
};

/* Virtual size of this process, in bytes. Zero when it cannot be read. */
static size_t virtual_size_bytes(void) {
  unsigned long pages = 0;
  FILE *fp = fopen("/proc/self/statm", "r");
  if (!fp)
    return 0;
  if (fscanf(fp, "%lu", &pages) != 1)
    pages = 0;
  fclose(fp);
  return (size_t)pages * (size_t)sysconf(_SC_PAGESIZE);
}

static void alternate_keys(void *ctx_void, void *out_void) {
  struct alternate_ctx *ctx = ctx_void;
  struct alternate_result *out = out_void;
  CPMDCSession *session;
  int step;
  unsetenv(k_limit_name);
  if (ctx->limit_address_space) {
    size_t base = virtual_size_bytes();
    struct rlimit lim;
    if (base != 0) {
      lim.rlim_cur = (rlim_t)(base + (size_t)k_address_space_margin);
      lim.rlim_max = lim.rlim_cur;
      out->rlimit_set = setrlimit(RLIMIT_AS, &lim) == 0;
    }
  }
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session) {
    copy_text(out->first_error, k_msg, cpmdc_last_error());
    return;
  }
  for (step = 0; step < k_switches; ++step) {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(0, pos);
    pos[0] += 0.001 * (double)(step + 1);
    cpmdc_session_select_orbitals(session, step % 2 == 0 ? 1 : 2);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    if (!r.ok) {
      if (out->first_error[0] == '\0')
        copy_text(out->first_error, k_msg, r.message);
      break;
    }
    out->switches_ok++;
  }
  out->capacity = cpmdc_orbital_ledger_capacity();
  cpmdc_session_destroy(session);
}

static void test_alternating_keys_keep_bounded_storage(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct alternate_ctx ctx;
  struct alternate_result out;
  char why[256];
  (void)state;
#if defined(CPMDC_HAS_CPMD)
  free(params);
  skip();
#endif
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  ctx.limit_address_space = 1;
#if defined(CPMDC_TEST_ASAN)
  /* The sanitizer maps address space the limit would refuse. The capacity
   * check below still catches a table that grows on every switch. */
  ctx.limit_address_space = 0;
#endif
  assert_int_equal(run_setup_child(alternate_keys, &ctx, &out, sizeof(out), why,
                                   sizeof(why)),
                   0);
  if (ctx.limit_address_space)
    assert_int_equal(out.rlimit_set, 1);
  assert_string_equal(out.first_error, "");
  assert_int_equal(out.switches_ok, k_switches);
  assert_true(out.capacity > 0);
  assert_true(out.capacity <= (size_t)k_first_block);
  free(params);
}

/* An invalid limit does not leak into the call a bad argument fails, and a
 * limit that cannot be used fails the copy that would have parked. */
struct order_ctx {
  const unsigned char *params;
  size_t params_n;
  const char *limit;
};

struct order_result {
  int null_session_ok;
  char null_session[k_msg];
  char null_session_last[k_msg];
  char null_buffer[k_msg];
  char bad_atomic_number[k_msg];
  char first_valid[k_msg];
  char topology[k_msg];
  char valid_args[k_msg];
  int copy_null_rc;
  char copy_null_last[k_msg];
  int copy_unevaluated_rc;
  char copy_unevaluated_last[k_msg];
  int copy_evaluated_rc;
  char copy_evaluated_last[k_msg];
  int select_rc;
};

static void order_checks(void *ctx_void, void *out_void) {
  struct order_ctx *ctx = ctx_void;
  struct order_result *out = out_void;
  static const int k_bad_z[3] = {8, 0, 1};
  static const int k_other_z[3] = {8, 8, 1};
  CPMDCSession *session;
  CPMDCSession *fresh;
  double pos[9];
  double forces[9];
  CPMDCResult r;
  positions_for(0, pos);
  setenv(k_limit_name, ctx->limit, 1);
  r = cpmdc_session_energy_forces(NULL, 3, pos, k_z, forces);
  out->null_session_ok = r.ok;
  copy_text(out->null_session, k_msg, r.message);
  copy_text(out->null_session_last, k_msg, cpmdc_last_error());
  out->copy_null_rc = cpmdc_session_copy_orbitals(NULL, 1, 2);
  copy_text(out->copy_null_last, k_msg, cpmdc_last_error());

  unsetenv(k_limit_name);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  fresh = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session || !fresh) {
    copy_text(out->first_valid, k_msg, cpmdc_last_error());
    return;
  }
  out->select_rc = cpmdc_session_select_orbitals(session, 1);

  r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  copy_text(out->first_valid, k_msg, r.message);
  setenv(k_limit_name, ctx->limit, 1);
  r = cpmdc_session_energy_forces(fresh, 3, NULL, k_z, forces);
  copy_text(out->null_buffer, k_msg, r.message);
  r = cpmdc_session_energy_forces(fresh, 3, pos, k_bad_z, forces);
  copy_text(out->bad_atomic_number, k_msg, r.message);
  out->copy_unevaluated_rc = cpmdc_session_copy_orbitals(fresh, 1, 2);
  copy_text(out->copy_unevaluated_last, k_msg, cpmdc_last_error());
  r = cpmdc_session_energy_forces(fresh, 3, pos, k_z, forces);
  copy_text(out->valid_args, k_msg, r.message);
  r = cpmdc_session_energy_forces(fresh, 3, pos, k_other_z, forces);
  copy_text(out->topology, k_msg, r.message);

  out->copy_evaluated_rc = cpmdc_session_copy_orbitals(session, 1, 2);
  copy_text(out->copy_evaluated_last, k_msg, cpmdc_last_error());
  cpmdc_session_destroy(fresh);
  cpmdc_session_destroy(session);
}

static void check_error_order(const char *limit) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct order_ctx ctx;
  struct order_result out;
  char why[256];
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  ctx.limit = limit;
  assert_int_equal(
      run_setup_child(order_checks, &ctx, &out, sizeof(out), why, sizeof(why)), 0);
  /* A bad call is reported as that, not as the limit. */
  assert_int_equal(out.null_session_ok, 0);
  assert_string_equal(out.null_session, "null session");
  assert_string_equal(out.null_session_last, "null session");
  assert_string_equal(out.null_buffer, "invalid arguments");
  assert_non_null(strstr(out.bad_atomic_number, "atomic numbers"));
  assert_null(strstr(out.bad_atomic_number, k_limit_name));
  assert_int_equal(out.copy_null_rc, -1);
  assert_string_equal(out.copy_null_last, "null session");
  assert_int_equal(out.copy_unevaluated_rc, -1);
  assert_non_null(strstr(out.copy_unevaluated_last, "no stored orbitals"));
  assert_null(strstr(out.copy_unevaluated_last, k_limit_name));
  /* A well formed call names the limit. */
  assert_int_equal(out.select_rc, 0);
  assert_string_equal(out.first_valid, "ok");
  assert_non_null(strstr(out.valid_args, k_limit_name));
  assert_int_equal(out.copy_evaluated_rc, -1);
  assert_non_null(strstr(out.copy_evaluated_last, k_limit_name));
  assert_null(strstr(out.copy_evaluated_last, k_exceeded));
  /* The first evaluation fixed the topology, so a change of it is the
   * reason for the failure while the limit stays unusable. */
  assert_non_null(strstr(out.topology, "topology change"));
  free(params);
}

static void test_bad_call_is_not_reported_as_the_limit(void **state) {
  (void)state;
  check_error_order("abc");
}

static void test_empty_limit_orders_errors_the_same_way(void **state) {
  (void)state;
  check_error_order("");
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_limit_rejects_a_new_park),
      cmocka_unit_test(test_lowering_keeps_stored_keys),
      cmocka_unit_test(test_invalid_limit_names_the_variable),
      cmocka_unit_test(test_unset_limit_parks_every_key),
      cmocka_unit_test(test_alternating_keys_keep_bounded_storage),
      cmocka_unit_test(test_bad_call_is_not_reported_as_the_limit),
      cmocka_unit_test(test_empty_limit_orders_errors_the_same_way),
  };
  if (argc < 2)
    return 1;
  g_params = argv[1];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
