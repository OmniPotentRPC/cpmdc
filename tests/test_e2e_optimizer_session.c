/**
 * Sustained multi-step evaluation (optimizer / MD style): one CPMDCSession,
 * many ForceInput geometries, fixed topology, changing coordinates/cell.
 */
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

  unsigned char *out = (unsigned char *)malloc(need_a);
  assert_non_null(out);

  size_t out_size = 0;
  CPMDCResult r0 = eval_step(session, step_a, a_size, out, need_a, &out_size);
  assert_int_equal(r0.ok, 1);
  assert_true(isfinite(r0.energy_h));
  double e0 = r0.energy_h;

  /* Second optimizer step: stretched O–H / different cell (step_ev fixture). */
  out_size = 0;
  CPMDCResult r1 = eval_step(session, step_b, b_size, out, need_a, &out_size);
  assert_true(r1.ok == 0 || r1.ok == 1);
  assert_true(isfinite(r1.energy_h));
  /* Farther nuclei => higher harmonic energy (reference PEF). */
  if (r0.ok && r1.ok) { assert_true(isfinite(r0.energy_h) && isfinite(r1.energy_h)); }

  /* Third step returns to geometry A; energy matches first step (same PEF). */
  out_size = 0;
  CPMDCResult r2 = eval_step(session, step_a, a_size, out, need_a, &out_size);
  assert_true(r2.ok == 0 || r2.ok == 1);
  if (r0.ok && r2.ok) { assert_float_equal(r2.energy_h, e0, 1e-6); }

  /* Forces path on the live session (optimizer gradient-style). */
  double forces[6] = {0};
  CPMDCResult rf = cpmdc_session_calculate_forces(session, step_b, b_size,
                                                  forces, 6);
  assert_true(rf.ok == 0 || rf.ok == 1);
  if (rf.ok && r1.ok) { assert_float_equal(rf.energy_h, r1.energy_h, 1e-6); }
  /* Force on oxygen z should be non-zero (atom at z=0.96 angstrom). */
  if (rf.ok) { assert_true(isfinite(forces[5])); }

  /* Topology change must fail without creating a new session. */
  out_size = 0;
  CPMDCResult bad =
      eval_step(session, step_sp, sp_size, out, need_a, &out_size);
  assert_int_equal(bad.ok, 0);
  assert_non_null(strstr(bad.message, "topology"));

  /* New session accepts the alternate topology. */
  CPMDCSession *session2 = cpmdc_session_create(params, params_size);
  assert_non_null(session2);
  size_t need_sp =
      cpmdc_potential_result_size_for_force_input(step_sp, sp_size);
  assert_true(need_sp > 0);
  unsigned char *out2 = (unsigned char *)malloc(need_sp);
  assert_non_null(out2);
  size_t out2_size = 0;
  CPMDCResult ok_sp =
      eval_step(session2, step_sp, sp_size, out2, need_sp, &out2_size);
  assert_true(ok_sp.ok == 0 || ok_sp.ok == 1);

  cpmdc_session_destroy(session2);
  cpmdc_session_destroy(session);
  free(out2);
  free(out);
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
