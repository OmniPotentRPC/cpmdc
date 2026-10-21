/**
 * A positive cell volume does not make the stress snapshot valid.
 * An isolated (symmetry 0 / Hockney) cell leaves it unset. A periodic cell
 * marks it valid only when the tensor was computed (CPMDC_STRESS is not 0).
 */
#define _GNU_SOURCE
#include "cpmdc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

#include "setup_child.h"

#include "schema/Potentials.capnp.h"
#include <capn.h>

static const char *g_isolated = NULL;
static const char *g_periodic = NULL;
static const char *g_step = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
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

struct stress_ctx {
  CPMDCSession *session;
  const unsigned char *step;
  size_t step_size;
  size_t need;
};

struct stress_out {
  CPMDCResult r;
  size_t wrote;
  unsigned char out[4096];
  int rc;
  CPMDCStressTensor stress;
  int session_rc;
  CPMDCStressTensor session_stress;
};

static void stress_child(void *vctx, void *vout) {
  struct stress_ctx *ctx = (struct stress_ctx *)vctx;
  struct stress_out *o = (struct stress_out *)vout;
  o->r = cpmdc_session_calculate_result(ctx->session, ctx->step,
                                        ctx->step_size, o->out, ctx->need,
                                        &o->wrote);
  o->rc = cpmdc_last_stress(&o->stress);
  o->session_rc = cpmdc_session_last_stress(ctx->session, &o->session_stress);
}

/* Returns the cpmdc_last_stress status. *wire_stress is 1 when the
 * PotentialResult carries a 9-vector. Each deck is its own OpenCPMD setup,
 * so the evaluation runs in a child of its own. */
static int eval_stress(const char *params_path, int *wire_stress,
                       CPMDCStressTensor *stress) {
  size_t ps = 0, ss = 0;
  unsigned char *params = read_file(params_path, &ps);
  unsigned char *step = read_file(g_step, &ss);
  assert_non_null(params);
  assert_non_null(step);
  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  assert_int_equal(cpmdc_last_error()[0], '\0');
  size_t need = cpmdc_potential_result_size_for_force_input(step, ss);
  assert_true(need > 0);
  static struct stress_out o;
  assert_true(need <= sizeof(o.out));
  struct stress_ctx ctx = {session, step, ss, need};
  RUN_SETUP_CHILD(stress_child, &ctx, &o);
  CPMDCResult r = o.r;
  size_t wrote = o.wrote;
  assert_int_equal(r.ok, 1);

  *stress = o.stress;
  int rc = o.rc;
  CPMDCStressTensor session_stress = o.session_stress;
  int session_rc = o.session_rc;
  assert_int_equal(session_rc, rc);
  assert_int_equal(session_stress.valid, stress->valid);

  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, o.out, (int)wrote, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);
  /* An omitted tensor is an empty list, not a null pointer. Nine values
   * means the snapshot was written. */
  int nstress = 0;
  if (view.stress.p.type != CAPN_NULL)
    nstress = capn_len(view.stress);
  *wire_stress = nstress == 9;
  capn_free(&arena);
  cpmdc_session_destroy(session);
  free(step);
  free(params);
  return rc;
}

static void test_isolated_stress_invalid(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "1", 1);
  CPMDCStressTensor stress;
  int wire = 1;
  int rc = eval_stress(g_isolated, &wire, &stress);
  assert_int_equal(rc, -1);
  assert_int_equal(stress.valid, 0);
  assert_int_equal(wire, 0);
}

static void test_periodic_stress_valid(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "1", 1);
  CPMDCStressTensor stress;
  int wire = 0;
  int rc = eval_stress(g_periodic, &wire, &stress);
  assert_int_equal(rc, 0);
  assert_int_equal(stress.valid, 1);
  assert_int_equal(wire, 1);
}

static void test_periodic_stress_skipped_when_disabled(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "0", 1);
  CPMDCStressTensor stress;
  int wire = 1;
  int rc = eval_stress(g_periodic, &wire, &stress);
  assert_int_equal(rc, -1);
  assert_int_equal(stress.valid, 0);
  assert_int_equal(wire, 0);
  setenv("CPMDC_STRESS", "1", 1);
}

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "usage: %s ISOLATED_PARAMS PERIODIC_PARAMS FORCE_INPUT\n",
            argv[0]);
    return 2;
  }
  g_isolated = argv[1];
  g_periodic = argv[2];
  g_step = argv[3];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_isolated_stress_invalid),
      cmocka_unit_test(test_periodic_stress_valid),
      cmocka_unit_test(test_periodic_stress_skipped_when_disabled),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
