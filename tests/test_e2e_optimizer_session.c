/**
 * Sustained multi-step evaluation (optimizer / MD style): one CPMDCSession,
 * many ForceInput geometries, fixed topology, changing coordinates/cell.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

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
static const char *g_step_a = NULL;
static const char *g_step_b = NULL;
static const char *g_step_species = NULL;

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

static CPMDCResult eval_step(CPMDCSession *session, const unsigned char *step,
                             size_t step_size, unsigned char *out, size_t need,
                             size_t *out_size) {
  *out_size = 0;
  return cpmdc_session_calculate_result(session, step, step_size, out, need,
                                        out_size);
}

/* Geometry A and geometry B have different cells, and the alternate
 * topology has a different composition: each is a new basis and so an
 * OpenCPMD setup of its own. Every setup runs in a child process on the
 * parent's session; the parent makes the assertions. */
struct opt_ctx {
  CPMDCSession *session;
  const unsigned char *step;
  size_t step_size;
  size_t need;
  /* Second call on the same session and basis, or NULL. */
  const unsigned char *then_step;
  size_t then_size;
};

struct opt_out {
  CPMDCResult r;
  size_t out_size;
  CPMDCResult then;
  size_t then_out_size;
};

static void eval_child(void *vctx, void *vout) {
  struct opt_ctx *ctx = (struct opt_ctx *)vctx;
  struct opt_out *o = (struct opt_out *)vout;
  unsigned char *out = (unsigned char *)malloc(ctx->need);
  if (!out) {
    snprintf(o->r.message, sizeof(o->r.message), "out of memory");
    return;
  }
  o->r = eval_step(ctx->session, ctx->step, ctx->step_size, out, ctx->need,
                   &o->out_size);
  if (ctx->then_step)
    o->then = eval_step(ctx->session, ctx->then_step, ctx->then_size, out,
                        ctx->need, &o->then_out_size);
  free(out);
}

struct forces_out {
  CPMDCResult rf;
  double forces[6];
};

static void forces_child(void *vctx, void *vout) {
  struct opt_ctx *ctx = (struct opt_ctx *)vctx;
  struct forces_out *o = (struct forces_out *)vout;
  o->rf = cpmdc_session_calculate_forces(ctx->session, ctx->step,
                                         ctx->step_size, o->forces, 6);
}

static void test_optimizer_style_session_loop(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked "
                  "(set meson -Dwith_cpmd=true -Dcpmd_root=/path/to/CPMD with lib/libcpmd.a)\n");
    skip();
  }
  assert_int_equal(cpmdc_available(), 1);

  size_t params_size = 0, a_size = 0, b_size = 0, sp_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step_a = read_file(g_step_a, &a_size);
  unsigned char *step_b = read_file(g_step_b, &b_size);
  unsigned char *step_sp = read_file(g_step_species, &sp_size);
  assert_non_null(params);
  assert_non_null(step_a);
  assert_non_null(step_b);
  assert_non_null(step_sp);

  CPMDCSession *session = cpmdc_session_create(params, params_size);
  assert_non_null(session);

  size_t need_a = cpmdc_potential_result_size_for_force_input(step_a, a_size);
  size_t need_b = cpmdc_potential_result_size_for_force_input(step_b, b_size);
  assert_int_equal(need_a, need_b);
  assert_true(need_a > 0);

  struct opt_out o0, o1, o2;
  struct opt_ctx c0 = {session, step_a, a_size, need_a, NULL, 0};
  RUN_SETUP_CHILD(eval_child, &c0, &o0);
  CPMDCResult r0 = o0.r;
  assert_int_equal(r0.ok, 1);
  assert_true(isfinite(r0.energy_h));
  double e0 = r0.energy_h;

  /* Second optimizer step: stretched O-H / different cell (step_ev fixture). */
  struct opt_ctx c1 = {session, step_b, b_size, need_a, NULL, 0};
  RUN_SETUP_CHILD(eval_child, &c1, &o1);
  CPMDCResult r1 = o1.r;
  assert_int_equal(r1.ok, 1);
  assert_true(isfinite(r1.energy_h));
  assert_true(fabs(r1.energy_h - e0) > 1.0e-8);

  /* Third step returns to geometry A; energy matches first step (same PEF).
   * The topology change follows on the same live session. */
  struct opt_ctx c2 = {session, step_a, a_size, need_a, step_sp, sp_size};
  RUN_SETUP_CHILD(eval_child, &c2, &o2);
  CPMDCResult r2 = o2.r;
  assert_int_equal(r2.ok, 1);
  assert_float_equal(r2.energy_h, e0, 1e-6);

  /* Forces path on the live session (optimizer gradient-style). */
  struct forces_out of;
  struct opt_ctx cf = {session, step_b, b_size, need_a, NULL, 0};
  RUN_SETUP_CHILD(forces_child, &cf, &of);
  CPMDCResult rf = of.rf;
  assert_int_equal(rf.ok, 1);
  assert_float_equal(rf.energy_h, r1.energy_h, 1e-6);
  /* Force on oxygen z should be non-zero (atom at z=0.96 angstrom). */
  assert_true(isfinite(of.forces[5]));

  /* Topology change must fail without creating a new session. */
  CPMDCResult bad = o2.then;
  assert_int_equal(bad.ok, 0);
  assert_non_null(strstr(bad.message, "topology"));

  /* New session accepts the alternate topology. */
  CPMDCSession *session2 = cpmdc_session_create(params, params_size);
  assert_non_null(session2);
  size_t need_sp =
      cpmdc_potential_result_size_for_force_input(step_sp, sp_size);
  assert_true(need_sp > 0);
  struct opt_out osp;
  struct opt_ctx csp = {session2, step_sp, sp_size, need_sp, NULL, 0};
  RUN_SETUP_CHILD(eval_child, &csp, &osp);
  CPMDCResult ok_sp = osp.r;
  assert_int_equal(ok_sp.ok, 1);

  cpmdc_session_destroy(session2);
  cpmdc_session_destroy(session);
  free(params);
  free(step_a);
  free(step_b);
  free(step_sp);
  cpmdc_finalize();
}

int main(int argc, char **argv) {
  if (argc < 5) {
    fprintf(stderr,
            "usage: %s params.bin step_a.bin step_b.bin step_species.bin\n",
            argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step_a = argv[2];
  g_step_b = argv[3];
  g_step_species = argv[4];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_optimizer_style_session_loop),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
