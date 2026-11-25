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
/* Keys the reference ledger holds parked, and whether it holds a selected
 * key. Both stay empty on a CPMD build. */
size_t cpmdc_orbital_ledger_parked(void);
int cpmdc_orbital_ledger_selected(long long *key);
void cpmdc_orbital_ledger_reset(void);
/* Keys with parked orbitals in the OpenCPMD table, the selected key
 * excluded. Zero on the reference evaluator. */
int cpmdc_embed_orbital_slots(void);
int cpmdc_embed_init(void);
/* Converged c0 in the selected extrapolation ring and the count the last
 * SCF combined. Both stay 0 on the reference evaluator. */
void cpmdc_embed_extrapolation_state(int *count, int *used);
/* The embed copy. The link stub of a reference build refuses every copy,
 * a copy of a key onto itself included. */
int cpmdc_embed_copy_orbitals(long long from_key, long long to_key);

static const char *g_params = NULL;
/* The water of g_params with EXTRAPOLATE WFN ASPC, for the ring cases. */
static const char *g_aspc_params = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
/* Small enough that the real SCF converges from a cold start at every step,
 * and still a new geometry for each step. */
static const double k_geometry_step = 0.25;
/* Close enough that each step starts from the stored c0 and extends the
 * ring instead of dropping it. */
static const double k_ring_step = 0.02;

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

static void ring_positions_for(int step, double out[9]) {
  memcpy(out, k_pos, sizeof(k_pos));
  out[0] = k_pos[0] + k_ring_step * (double)step;
}

static int ring_count(void) {
  int count = -1, used = -1;
  cpmdc_embed_extrapolation_state(&count, &used);
  return count;
}

static int ring_used(void) {
  int count = -1, used = -1;
  cpmdc_embed_extrapolation_state(&count, &used);
  return used;
}

/* Parked keys in whichever store this build counts against the limit. */
static int orbital_slots(void) {
#if defined(CPMDC_HAS_CPMD)
  return cpmdc_embed_orbital_slots();
#else
  return (int)cpmdc_orbital_ledger_parked();
#endif
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
    positions_for(step, pos);
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
    positions_for(3, pos);
    cpmdc_session_select_orbitals(session, 22);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->kept_ok = r.ok;
    copy_text(out->kept_error, k_msg, cpmdc_last_error());
  }
  {
    double pos[9];
    double forces[9];
    CPMDCResult r;
    positions_for(4, pos);
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

#if !defined(CPMDC_HAS_CPMD)
static void test_finalize_frees_and_reinitializes_the_ledger(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  CPMDCSession *session;
  CPMDCSession *fresh;
  double forces[9];
  double pos[9];
  CPMDCResult result;
  int step;
  (void)state;
  assert_non_null(params);
  setenv(k_limit_name, "2", 1);
  assert_int_equal(cpmdc_embed_init(), 1);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  for (step = 0; step < 2; ++step) {
    positions_for(step, pos);
    assert_int_equal(cpmdc_session_select_orbitals(session, 61 + step), 0);
    result = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    assert_int_equal(result.ok, 1);
  }
  assert_true(cpmdc_orbital_ledger_capacity() > 0);
  cpmdc_session_destroy(session);

  /* A fresh session starts from an empty ledger after destroy. */
  fresh = cpmdc_session_create(params, n);
  assert_non_null(fresh);
  positions_for(2, pos);
  assert_int_equal(cpmdc_session_select_orbitals(fresh, 63), 0);
  result = cpmdc_session_energy_forces(fresh, 3, pos, k_z, forces);
  assert_int_equal(result.ok, 1);
  cpmdc_session_destroy(fresh);

  cpmdc_finalize();
  assert_int_equal(cpmdc_orbital_ledger_capacity(), 0);

  assert_int_equal(cpmdc_embed_init(), 1);
  assert_int_equal(cpmdc_orbital_ledger_capacity(), 0);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  positions_for(2, pos);
  assert_int_equal(cpmdc_session_select_orbitals(session, 63), 0);
  result = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(result.ok, 1);
  assert_int_equal(cpmdc_orbital_ledger_capacity(), 0);
  cpmdc_session_destroy(session);
  cpmdc_finalize();
  assert_int_equal(cpmdc_orbital_ledger_capacity(), 0);
  unsetenv(k_limit_name);
  free(params);
}
#endif

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
  char after_rejected_gate[k_msg];
  char valid_args[k_msg];
  int fresh_reconfig_rc;
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
  out->fresh_reconfig_rc =
      cpmdc_session_set_params(fresh, ctx->params, ctx->params_n);
  r = cpmdc_session_energy_forces(session, 3, pos, k_other_z, forces);
  copy_text(out->topology, k_msg, r.message);
  unsetenv(k_limit_name);
  r = cpmdc_session_energy_forces(fresh, 3, pos, k_other_z, forces);
  copy_text(out->after_rejected_gate, k_msg, r.message);

  setenv(k_limit_name, ctx->limit, 1);
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
  /* A well formed call names the limit and leaves the topology open. */
  assert_int_equal(out.select_rc, 0);
  assert_string_equal(out.first_valid, "ok");
  assert_non_null(strstr(out.valid_args, k_limit_name));
  assert_int_equal(out.fresh_reconfig_rc, 0);
  assert_int_equal(out.copy_evaluated_rc, -1);
  assert_non_null(strstr(out.copy_evaluated_last, k_limit_name));
  assert_null(strstr(out.copy_evaluated_last, k_exceeded));
  /* A gate rejection does not latch the fresh session's topology. */
  assert_string_equal(out.after_rejected_gate, "ok");
  /* The first evaluation on `session` fixed its topology, so a change of it
   * is the reason for the failure while the limit stays unusable. */
  assert_non_null(strstr(out.topology, "topology change"));
  assert_null(strstr(out.topology, k_limit_name));
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

struct preflight_result {
  char topology_oom[k_msg];
  int configure_after_oom;
  char stateless_invalid[k_msg];
};

static void preflight_checks(void *ctx_void, void *out_void) {
  struct order_ctx *ctx = ctx_void;
  struct preflight_result *out = out_void;
  CPMDCSession *session;
  double pos[9];
  double grad[9];
  CPMDCResult r;
  positions_for(0, pos);
  unsetenv(k_limit_name);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session)
    return;
  setenv("CPMDC_TEST_FAULT", "topology-allocation", 1);
  r = cpmdc_session_energy_gradient(session, 3, pos, k_z, grad);
  copy_text(out->topology_oom, k_msg, r.message);
  unsetenv("CPMDC_TEST_FAULT");
  out->configure_after_oom =
      cpmdc_session_set_params(session, ctx->params, ctx->params_n);
  cpmdc_session_destroy(session);

  setenv(k_limit_name, "abc", 1);
  r = cpmdc_energy_gradient(3, NULL, k_z, ctx->params, ctx->params_n, grad);
  copy_text(out->stateless_invalid, k_msg, r.message);
}

static void test_preflight_errors_do_not_latch_topology(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct order_ctx ctx;
  struct preflight_result out;
  char why[256];
  (void)state;
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  ctx.limit = NULL;
  assert_int_equal(run_setup_child(preflight_checks, &ctx, &out, sizeof(out),
                                   why, sizeof(why)),
                   0);
  assert_string_equal(out.topology_oom, "out of memory");
  assert_int_equal(out.configure_after_oom, 0);
  assert_string_equal(out.stateless_invalid, "invalid arguments");
  free(params);
}

struct unkeyed_ctx {
  const unsigned char *params;
  size_t params_n;
};

struct destroy_ctx {
  const unsigned char *params;
  size_t params_n;
  int destroy_evaluating;
};

struct unkeyed_result {
  int keyed_ok;
  int ring_before_unkeyed;
  int selected_after_keyed;
  long long selected_key_after_keyed;
  int unkeyed_ok;
  int ring_after_unkeyed;
  int used_unkeyed;
  int selected_after_unkeyed;
  int keyed_after_unkeyed_ok;
  char keyed_after_unkeyed_message[k_msg];
  char keyed_after_unkeyed_error[k_msg];
  int slots_after_keyed;
  int next_key_ok;
  char next_key_message[k_msg];
  int slots_after_next_key;
  int selected_after_next_key;
  long long selected_key_after_next_key;
  int self_copy_rc;
  int stub_self_copy_rc;
  int ring_before_self_copy;
  int ring_after_self_copy;
  int slots_after_self_copy;
};

/* Evaluates `session` under key 302 at two nearby geometries, so an
 * OpenCPMD build holds one converged c0 in the key's ring. */
static void keyed_twice(CPMDCSession *session, struct unkeyed_result *out) {
  int step;
  out->keyed_ok = 1;
  cpmdc_session_select_orbitals(session, 302);
  for (step = -1; step <= 0; ++step) {
    CPMDCResult r;
    double pos[9];
    double forces[9];
    ring_positions_for(step, pos);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->keyed_ok = out->keyed_ok && r.ok;
  }
  out->ring_before_unkeyed = ring_count();
  out->selected_after_keyed =
      cpmdc_orbital_ledger_selected(&out->selected_key_after_keyed);
}

/* Evaluates `session` under key 303 and then 304. Key 303 starts with no
 * key selected, so it parks nothing, and the switch to 304 parks 303 alone
 * under a limit of 1. A copy of 304 onto itself then changes nothing. */
static void keyed_after_unkeyed(CPMDCSession *session,
                                struct unkeyed_result *out) {
  CPMDCResult r;
  double pos[9];
  double forces[9];
  ring_positions_for(2, pos);
  cpmdc_session_select_orbitals(session, 303);
  r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  out->keyed_after_unkeyed_ok = r.ok;
  copy_text(out->keyed_after_unkeyed_message, k_msg, r.message);
  copy_text(out->keyed_after_unkeyed_error, k_msg, cpmdc_last_error());
  out->slots_after_keyed = orbital_slots();

  ring_positions_for(3, pos);
  cpmdc_session_select_orbitals(session, 304);
  r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  out->next_key_ok = r.ok;
  copy_text(out->next_key_message, k_msg, r.message);
  out->slots_after_next_key = orbital_slots();
  out->selected_after_next_key =
      cpmdc_orbital_ledger_selected(&out->selected_key_after_next_key);

  out->ring_before_self_copy = ring_count();
  out->self_copy_rc = cpmdc_session_copy_orbitals(session, 304, 304);
  out->ring_after_self_copy = ring_count();
  out->slots_after_self_copy = orbital_slots();
#if defined(CPMDC_HAS_CPMD)
  out->stub_self_copy_rc = -1;
#else
  out->stub_self_copy_rc = cpmdc_embed_copy_orbitals(304, 304);
#endif
}

/* The OpenCPMD build keeps the key and ring in Fortran, and the reference
 * ledger stays empty. The reference build keeps the key in its ledger, and
 * the link stub reports an empty ring and refuses every embed copy. */
static void assert_keyed_after_unkeyed(const struct unkeyed_result *out) {
  assert_int_equal(out->keyed_after_unkeyed_ok, 1);
  assert_string_equal(out->keyed_after_unkeyed_error, "");
  assert_int_equal(out->slots_after_keyed, 0);
  assert_int_equal(out->next_key_ok, 1);
  assert_string_equal(out->next_key_message, "ok");
  assert_int_equal(out->slots_after_next_key, 1);
  assert_int_equal(out->self_copy_rc, 0);
  assert_int_equal(out->ring_after_self_copy, out->ring_before_self_copy);
  assert_int_equal(out->slots_after_self_copy, 1);
  assert_int_equal(out->stub_self_copy_rc, -1);
#if defined(CPMDC_HAS_CPMD)
  assert_int_equal(out->selected_after_next_key, 0);
#else
  assert_int_equal(out->selected_after_next_key, 1);
  assert_int_equal(out->selected_key_after_next_key, 304);
  assert_int_equal(out->ring_before_self_copy, 0);
#endif
}

static void assert_unkeyed_result(const struct unkeyed_result *out) {
  assert_int_equal(out->keyed_ok, 1);
#if defined(CPMDC_HAS_CPMD)
  assert_true(out->ring_before_unkeyed > 0);
  assert_int_equal(out->selected_after_keyed, 0);
#else
  assert_int_equal(out->ring_before_unkeyed, 0);
  assert_int_equal(out->selected_after_keyed, 1);
  assert_int_equal(out->selected_key_after_keyed, 302);
#endif
  assert_int_equal(out->unkeyed_ok, 1);
  assert_int_equal(out->ring_after_unkeyed, 0);
  assert_int_equal(out->used_unkeyed, 0);
  assert_int_equal(out->selected_after_unkeyed, 0);
  assert_keyed_after_unkeyed(out);
}

static void unkeyed_session_switch_case(void *ctx_void, void *out_void) {
  struct unkeyed_ctx *ctx = ctx_void;
  struct unkeyed_result *out = out_void;
  CPMDCSession *session_a = NULL;
  CPMDCSession *session_b = NULL;
  CPMDCResult r;
  double pos[9];
  double forces[9];

  setenv(k_limit_name, "1", 1);
  session_a = cpmdc_session_create(ctx->params, ctx->params_n);
  session_b = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session_a || !session_b) {
    copy_text(out->keyed_after_unkeyed_error, k_msg, cpmdc_last_error());
    goto done;
  }

  keyed_twice(session_a, out);

  ring_positions_for(1, pos);
  r = cpmdc_session_energy_forces(session_b, 3, pos, k_z, forces);
  out->unkeyed_ok = r.ok;
  out->ring_after_unkeyed = ring_count();
  out->used_unkeyed = ring_used();
  out->selected_after_unkeyed = cpmdc_orbital_ledger_selected(NULL);

  keyed_after_unkeyed(session_b, out);

done:
  cpmdc_session_destroy(session_b);
  cpmdc_session_destroy(session_a);
}

static void global_unkeyed_switch_case(void *ctx_void, void *out_void) {
  struct unkeyed_ctx *ctx = ctx_void;
  struct unkeyed_result *out = out_void;
  CPMDCSession *session = NULL;
  CPMDCResult r;
  double pos[9];
  double forces[9];

  setenv(k_limit_name, "1", 1);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session) {
    copy_text(out->keyed_after_unkeyed_error, k_msg, cpmdc_last_error());
    goto done;
  }

  keyed_twice(session, out);

  ring_positions_for(1, pos);
  r = cpmdc_energy_forces(3, pos, k_z, ctx->params, ctx->params_n, forces);
  out->unkeyed_ok = r.ok;
  out->ring_after_unkeyed = ring_count();
  out->used_unkeyed = ring_used();
  out->selected_after_unkeyed = cpmdc_orbital_ledger_selected(NULL);

  keyed_after_unkeyed(session, out);

done:
  cpmdc_session_destroy(session);
}

static void test_unkeyed_context_clears_selected_key(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_aspc_params, &n);
  struct unkeyed_ctx ctx;
  struct unkeyed_result session_out;
  struct unkeyed_result global_out;
  char why[256];
  (void)state;
  assert_non_null(params);
  memset(&session_out, 0, sizeof(session_out));
  memset(&global_out, 0, sizeof(global_out));
  ctx.params = params;
  ctx.params_n = n;

  assert_int_equal(
      run_setup_child(unkeyed_session_switch_case, &ctx, &session_out,
                      sizeof(session_out), why, sizeof(why)),
      0);
  assert_unkeyed_result(&session_out);

  assert_int_equal(
      run_setup_child(global_unkeyed_switch_case, &ctx, &global_out,
                      sizeof(global_out), why, sizeof(why)),
      0);
  assert_unkeyed_result(&global_out);
  free(params);
}

/* Session B is created after keyed session A evaluated, so B is the
 * configured session when A is destroyed. A still owns the selected key and
 * its history, and destroying A clears them. Destroying B instead leaves
 * them in place. */
struct destroy_result {
  int keyed_ok[2];
  int ring_before_destroy;
  int selected_before_destroy;
  long long selected_key_before_destroy;
  int ring_after_destroy;
  int selected_after_destroy;
  long long selected_key;
  struct unkeyed_result after;
};

static void destroy_after_later_create_case(void *ctx_void, void *out_void) {
  struct destroy_ctx *ctx = ctx_void;
  struct destroy_result *out = out_void;
  CPMDCSession *session_a = NULL;
  CPMDCSession *session_b = NULL;
  int step;

  setenv(k_limit_name, "1", 1);
  session_a = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session_a)
    goto done;
  cpmdc_session_select_orbitals(session_a, 302);
  for (step = 0; step < 2; ++step) {
    CPMDCResult r;
    double pos[9];
    double forces[9];
    ring_positions_for(step, pos);
    r = cpmdc_session_energy_forces(session_a, 3, pos, k_z, forces);
    out->keyed_ok[step] = r.ok;
  }
  session_b = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session_b)
    goto done;
  out->ring_before_destroy = ring_count();
  out->selected_before_destroy =
      cpmdc_orbital_ledger_selected(&out->selected_key_before_destroy);
  if (ctx->destroy_evaluating) {
    cpmdc_session_destroy(session_a);
    session_a = NULL;
  } else {
    cpmdc_session_destroy(session_b);
    session_b = NULL;
  }
  out->ring_after_destroy = ring_count();
  out->selected_after_destroy =
      cpmdc_orbital_ledger_selected(&out->selected_key);
  if (ctx->destroy_evaluating)
    keyed_after_unkeyed(session_b, &out->after);

done:
  cpmdc_session_destroy(session_b);
  cpmdc_session_destroy(session_a);
}

/* Before either destroy, key 302 holds a ring on an OpenCPMD build and is
 * the reference ledger's selected key on the reference build, so the
 * checks after the destroy do not pass on state that was never there. */
static void assert_destroy_starting_state(const struct destroy_result *out) {
#if defined(CPMDC_HAS_CPMD)
  assert_true(out->ring_before_destroy > 0);
  assert_int_equal(out->selected_before_destroy, 0);
#else
  assert_int_equal(out->ring_before_destroy, 0);
  assert_int_equal(out->selected_before_destroy, 1);
  assert_int_equal(out->selected_key_before_destroy, 302);
#endif
}

static void test_destroy_clears_the_evaluating_session_key(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_aspc_params, &n);
  struct destroy_ctx ctx;
  struct destroy_result out;
  char why[256];
  (void)state;
  assert_non_null(params);
  ctx.params = params;
  ctx.params_n = n;

  memset(&out, 0, sizeof(out));
  ctx.destroy_evaluating = 1;
  assert_int_equal(run_setup_child(destroy_after_later_create_case, &ctx, &out,
                                   sizeof(out), why, sizeof(why)),
                   0);
  assert_int_equal(out.keyed_ok[0], 1);
  assert_int_equal(out.keyed_ok[1], 1);
  assert_destroy_starting_state(&out);
  assert_int_equal(out.ring_after_destroy, 0);
  assert_int_equal(out.selected_after_destroy, 0);
  assert_keyed_after_unkeyed(&out.after);

  memset(&out, 0, sizeof(out));
  ctx.destroy_evaluating = 0;
  assert_int_equal(run_setup_child(destroy_after_later_create_case, &ctx, &out,
                                   sizeof(out), why, sizeof(why)),
                   0);
  assert_int_equal(out.keyed_ok[0], 1);
  assert_int_equal(out.keyed_ok[1], 1);
  assert_destroy_starting_state(&out);
  assert_int_equal(out.ring_after_destroy, out.ring_before_destroy);
#if !defined(CPMDC_HAS_CPMD)
  assert_int_equal(out.selected_after_destroy, 1);
  assert_int_equal(out.selected_key, 302);
#endif
  free(params);
}

#if !defined(CPMDC_HAS_CPMD)
/* A reset empties the ledger, the selected key included, so the store the
 * next evaluation leaves is not credited to the key selected before it. */
static void ledger_reset_case(void *ctx_void, void *out_void) {
  struct unkeyed_ctx *ctx = ctx_void;
  struct unkeyed_result *out = out_void;
  CPMDCSession *session;
  CPMDCResult r;
  double pos[9];
  double forces[9];
  unsetenv(k_limit_name);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session)
    return;
  positions_for(0, pos);
  cpmdc_session_select_orbitals(session, 71);
  r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  out->keyed_ok = r.ok;
  copy_text(out->keyed_after_unkeyed_message, k_msg, r.message);
  out->unkeyed_ok = cpmdc_orbital_ledger_selected(NULL);
  cpmdc_orbital_ledger_reset();
  out->selected_after_unkeyed = cpmdc_orbital_ledger_selected(NULL);
  out->slots_after_keyed = orbital_slots();
  cpmdc_session_destroy(session);
}

static void test_ledger_reset_clears_the_selected_key(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct unkeyed_ctx ctx;
  struct unkeyed_result out;
  char why[256];
  (void)state;
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  assert_int_equal(run_setup_child(ledger_reset_case, &ctx, &out, sizeof(out),
                                   why, sizeof(why)),
                   0);
  assert_string_equal(out.keyed_after_unkeyed_message, "ok");
  assert_int_equal(out.keyed_ok, 1);
  assert_int_equal(out.unkeyed_ok, 1);
  assert_int_equal(out.selected_after_unkeyed, 0);
  assert_int_equal(out.slots_after_keyed, 0);
  free(params);
}
#endif

#if defined(CPMDC_HAS_CPMD)
/* A basis change frees every key's orbitals, so a limit the old basis
 * filled does not refuse the first key of the new one. */
struct basis_ctx {
  const unsigned char *params;
  size_t params_n;
};

struct basis_result {
  int old_ok[3];
  int new_ok;
  char new_message[k_msg];
  char new_error[k_msg];
};

static void new_basis_after_saturation(void *ctx_void, void *out_void) {
  struct basis_ctx *ctx = ctx_void;
  struct basis_result *out = out_void;
  static const long long old_keys[3] = {51, 52, 53};
  static const int h2_z[2] = {1, 1};
  static const double h2_pos[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.74};
  CPMDCSession *session;
  double forces[9];
  double h2_forces[6];
  int step;
  setenv(k_limit_name, "1", 1);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session)
    return;
  for (step = 0; step < 2; ++step) {
    double pos[9];
    CPMDCResult r;
    positions_for(step, pos);
    cpmdc_session_select_orbitals(session, old_keys[step]);
    r = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
    out->old_ok[step] = r.ok;
  }
  cpmdc_session_destroy(session);
  session = cpmdc_session_create(ctx->params, ctx->params_n);
  if (!session)
    return;
  {
    CPMDCResult r;
    cpmdc_session_select_orbitals(session, old_keys[2]);
    r = cpmdc_session_energy_forces(session, 2, h2_pos, h2_z, h2_forces);
    out->new_ok = r.ok;
    copy_text(out->new_message, k_msg, r.message);
    copy_text(out->new_error, k_msg, cpmdc_last_error());
  }
  cpmdc_session_destroy(session);
}

static void test_new_basis_clears_a_saturated_limit(void **state) {
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  struct basis_ctx ctx;
  struct basis_result out;
  char why[256];
  (void)state;
  assert_non_null(params);
  memset(&out, 0, sizeof(out));
  ctx.params = params;
  ctx.params_n = n;
  assert_int_equal(run_setup_child(new_basis_after_saturation, &ctx, &out,
                                   sizeof(out), why, sizeof(why)),
                   0);
  assert_int_equal(out.old_ok[0], 1);
  assert_int_equal(out.old_ok[1], 1);
  assert_int_equal(out.new_ok, 1);
  assert_string_equal(out.new_error, "");
  free(params);
}
#endif

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_limit_rejects_a_new_park),
      cmocka_unit_test(test_lowering_keeps_stored_keys),
      cmocka_unit_test(test_invalid_limit_names_the_variable),
      cmocka_unit_test(test_unset_limit_parks_every_key),
      cmocka_unit_test(test_alternating_keys_keep_bounded_storage),
      cmocka_unit_test(test_bad_call_is_not_reported_as_the_limit),
      cmocka_unit_test(test_empty_limit_orders_errors_the_same_way),
      cmocka_unit_test(test_preflight_errors_do_not_latch_topology),
      cmocka_unit_test(test_unkeyed_context_clears_selected_key),
      cmocka_unit_test(test_destroy_clears_the_evaluating_session_key),
#if defined(CPMDC_HAS_CPMD)
      cmocka_unit_test(test_new_basis_clears_a_saturated_limit),
#else
      cmocka_unit_test(test_ledger_reset_clears_the_selected_key),
      cmocka_unit_test(test_finalize_frees_and_reinitializes_the_ledger),
#endif
  };
  if (argc < 3)
    return 1;
  g_params = argv[1];
  g_aspc_params = argv[2];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
