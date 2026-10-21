/**
 * End-to-end single-point evaluation through Cap'n Proto carriers.
 * Uses the embed-shell reference PEF (no OpenCPMD archives required).
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <errno.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

#include "setup_child.h"

static const char *g_params = NULL;
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

/* The ForceInput result and the session forces share one basis (the step
 * cell). The C-array entry point has no cell, which is a new basis and
 * so a second OpenCPMD setup: it runs in a child of its own. */
struct step_ctx {
  const unsigned char *params;
  size_t params_size;
  const unsigned char *step;
  size_t step_size;
  size_t need;
  CPMDCSession *session;
};

struct step_out {
  CPMDCResult r;
  size_t out_size;
  unsigned char out[4096];
  CPMDCResult f;
  double forces[6];
};

static void step_child(void *vctx, void *vout) {
  struct step_ctx *ctx = (struct step_ctx *)vctx;
  struct step_out *o = (struct step_out *)vout;
  o->r = cpmdc_calculate_result(ctx->params, ctx->params_size, ctx->step,
                                ctx->step_size, o->out, ctx->need,
                                &o->out_size);
  o->f = cpmdc_session_calculate_forces(ctx->session, ctx->step,
                                        ctx->step_size, o->forces, 6);
}

struct nocell_out {
  CPMDCResult fnc;
  double forces_nc[6];
};

static void nocell_child(void *vctx, void *vout) {
  struct step_ctx *ctx = (struct step_ctx *)vctx;
  struct nocell_out *o = (struct nocell_out *)vout;
  double positions[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.7414};
  int atmnrs[2] = {1, 8};
  o->fnc = cpmdc_energy_forces(2, positions, atmnrs, ctx->params,
                               ctx->params_size, o->forces_nc);
}

static void test_single_point_calculate_result(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked "
                  "(set meson -Dwith_cpmd=true -Dcpmd_root=/path/to/CPMD with lib/libcpmd.a)\n");
    skip();
  }
  assert_int_equal(cpmdc_available(), 1);

  size_t params_size = 0, step_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(params);
  assert_non_null(step);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need >= 32u + 6u * 8u);

  static struct step_out so;
  assert_true(need <= sizeof(so.out));
  /* Session forces path on the same ForceInput (includes cell term). */
  CPMDCSession *session = cpmdc_session_create(params, params_size);
  assert_non_null(session);
  struct step_ctx ctx = {params, params_size, step, step_size, need, session};
  RUN_SETUP_CHILD(step_child, &ctx, &so);

  CPMDCResult r = so.r;
  size_t out_size = so.out_size;
  unsigned char *out = so.out;
  assert_int_equal(r.ok, 1);
  assert_true(out_size > 0);
  assert_true(out_size <= need);
  assert_true(isfinite(r.energy_h));
  assert_true(fabs(r.energy_h) > 1e-8);

  /* Decode PotentialResult */
  struct capn arena;
  memset(&arena, 0, sizeof(arena));
  assert_int_equal(capn_init_mem(&arena, out, out_size, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  assert_int_equal(pr.p.type, CAPN_STRUCT);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);
  assert_true(isfinite(view.energy));
  assert_float_equal(view.energy, r.energy_h, 1e-9);
  assert_int_equal(view.forces.p.type == CAPN_FAR_POINTER ||
                       view.forces.p.type == CAPN_LIST,
                   1);
  capn_resolve(&view.forces.p);
  assert_int_equal(view.forces.p.len, 6);

  CPMDCResult f = so.f;
  assert_int_equal(f.ok, 1);
  assert_true(isfinite(f.energy_h));
  assert_float_equal(f.energy_h, r.energy_h, 1e-9);
  for (int i = 0; i < 6; ++i)
    assert_true(isfinite(so.forces[i]));
  /* C-array entry point has no cell; energy is slightly lower than ForceInput. */
  static struct nocell_out no;
  RUN_SETUP_CHILD(nocell_child, &ctx, &no);
  CPMDCResult fnc = no.fnc;
  assert_int_equal(fnc.ok, 1);
  assert_true(isfinite(fnc.energy_h));
  cpmdc_session_destroy(session);

  cpmdc_params_release(&arena);
  free(params);
  free(step);
  cpmdc_finalize();
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin force_input.bin\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_single_point_calculate_result),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
