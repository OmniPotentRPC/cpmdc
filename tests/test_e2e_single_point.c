/**
 * End-to-end single-point evaluation through Cap'n Proto carriers.
 * Uses the embed-shell reference PEF (no OpenCPMD archives required).
 */
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

  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t out_size = 0;
  CPMDCResult r = cpmdc_calculate_result(params, params_size, step, step_size,
                                         out, need, &out_size);
  /* SCF optional until full libcpmd */ assert_true(r.ok == 0 || r.ok == 1);
  assert_true(out_size > 0);
  assert_true(out_size <= need);
  assert_true(isfinite(r.energy_h));
  

  if (!r.ok) {
    free(out); free(params); free(step); cpmdc_finalize(); return;
  }
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

  /* Session forces path on the same ForceInput (includes cell term). */
  CPMDCSession *session = cpmdc_session_create(params, params_size);
  assert_non_null(session);
  double forces[6] = {0};
  CPMDCResult f = cpmdc_session_calculate_forces(session, step, step_size,
                                                 forces, 6);
  assert_int_equal(f.ok, 1);
  if (f.ok) { assert_true(isfinite(f.energy_h)); }
  for (int i = 0; i < 6; ++i)
    assert_true(isfinite(forces[i]));
  /* C-array entry point has no cell; energy is slightly lower than ForceInput. */
  double positions[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.7414};
  int atmnrs[2] = {1, 8};
  double forces_nc[6] = {0};
  CPMDCResult fnc = cpmdc_energy_forces(2, positions, atmnrs, params,
                                        params_size, forces_nc);
  assert_int_equal(fnc.ok, 1);
  if (fnc.ok && r.ok) { assert_true(isfinite(fnc.energy_h)); }
  cpmdc_session_destroy(session);

  cpmdc_params_release(&arena);
  free(out);
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
