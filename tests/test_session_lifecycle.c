/**
 * Session lifecycle on the reference evaluator and on a linked OpenCPMD.
 *
 * A session accepts a constant cutoff until the first evaluation fixes the
 * topology. After that, another parameter buffer and another cutoff are
 * refused, and the energy already computed is unchanged. cpmdc_finalize()
 * is process shutdown: a later evaluation fails, and the diagnostic says
 * the embed is not available.
 */
#include "cpmdc.h"

#include <limits.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};

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

static void test_topology_freeze_keeps_the_energy(void **state) {
  (void)state;
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, n);
  assert_non_null(session);
  assert_int_equal(cpmdc_last_error()[0], '\0');

  assert_int_equal(cpmdc_session_set_params(session, NULL, 0), -1);
  assert_non_null(strstr(cpmdc_last_error(), "CPMDParams buffer is empty"));
  assert_int_equal(cpmdc_session_set_constant_cutoff(session, 30.0, 0.0, 0.0),
                   0);
  assert_int_equal(cpmdc_last_error()[0], '\0');

  double forces[9];
  CPMDCResult first =
      cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(first.ok, 1);
  assert_int_equal(cpmdc_last_error()[0], '\0');

  CPMDCEnergyComponents components;
  memset(&components, 0xab, sizeof(components));
  assert_int_equal(cpmdc_session_last_energy_components(session, &components),
                   0);
  assert_int_equal(components.valid, 1);
  assert_true(components.etot == first.energy_h);

  CPMDCStressTensor stress;
  memset(&stress, 0xab, sizeof(stress));
  assert_int_equal(cpmdc_session_last_stress(session, &stress), -1);
  assert_int_equal(stress.valid, 0);
  assert_int_equal(cpmdc_session_last_stress(NULL, &stress), -1);
  assert_int_equal(cpmdc_session_last_energy_components(session, NULL), -1);
  assert_int_equal(cpmdc_last_charge_integrals(NULL), -1);
  assert_int_equal(cpmdc_last_rank_spread(NULL), -1);

  double energy = 0.0;
  double au_forces[9];
  assert_int_equal(cpmdc_session_calculate_au(session, INT_MAX / 3 + 1, k_pos,
                                              k_z, NULL, &energy, au_forces,
                                              NULL),
                   -1);
  assert_string_equal(cpmdc_last_error(), "invalid arguments");

  assert_int_equal(cpmdc_session_set_constant_cutoff(session, 20.0, 0.0, 0.0),
                   -1);
  assert_string_equal(cpmdc_last_error(), "session topology is already fixed");
  assert_int_equal(cpmdc_session_set_params(session, params, n), -1);
  assert_string_equal(cpmdc_last_error(), "session topology is already fixed");

  double again_forces[9];
  CPMDCResult again =
      cpmdc_session_energy_forces(session, 3, k_pos, k_z, again_forces);
  assert_int_equal(again.ok, 1);
  assert_true(again.energy_h == first.energy_h);
  assert_memory_equal(again_forces, forces, sizeof(forces));

  cpmdc_session_destroy(session);
  cpmdc_session_destroy(NULL);
  free(params);
}

static void test_finalize_refuses_later_evaluations(void **state) {
  (void)state;
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, n);
  assert_non_null(session);
  double forces[9];
  CPMDCResult ready =
      cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(ready.ok, 1);
  cpmdc_session_destroy(session);

  cpmdc_finalize();

  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  CPMDCResult after =
      cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(after.ok, 0);
  assert_non_null(strstr(after.message, "CPMD embed not available"));
  assert_string_equal(cpmdc_last_error(), after.message);
  cpmdc_session_destroy(session);
  free(params);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_topology_freeze_keeps_the_energy),
      cmocka_unit_test(test_finalize_refuses_later_evaluations),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
