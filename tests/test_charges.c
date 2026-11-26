#define _POSIX_C_SOURCE 200809L

/**
 * Hirshfeld charges through the C ABI.
 *
 * The reference evaluator leaves the snapshot unset. A linked OpenCPMD
 * stores q = Z_v - n from atomc, in the caller's atom order. Neutral
 * water sums to the total charge, 0, within 2e-3 e. This 50 Ry Hockney
 * deck leaves a real-space partition residual of about 1.2e-3 e, so the
 * test bound is 2e-3 e. A permutation
 * of the atoms carries the same charges. K points, broken symmetry, and a
 * species without a weight function stay invalid, as does
 * CPMDC_CHARGE_SNAPSHOT=0.
 */
#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
static const double k_perm[9] = {
    0.0, 0.7572, -0.4692, 0.0, 0.0, 0.1173, 0.0, -0.7572, -0.4692,
};
static const int k_perm_z[3] = {1, 8, 1};

/* Sum of a neutral molecule, electrons. */
static const double k_sum_tol = 2.0e-3;
/* Two cold evaluations of the same atoms. */
static const double k_perm_tol = 1.0e-8;

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

static int buffer_untouched(const double *buf, size_t n) {
  const unsigned char *bytes = (const unsigned char *)buf;
  size_t i;
  for (i = 0; i < n * sizeof(double); ++i) {
    if (bytes[i] != 0x5a)
      return 0;
  }
  return 1;
}

static void test_unset_charges_do_not_write_the_buffer(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  CPMDCResult result;
  CPMDCPropertySnapshot prop;
  double held[4];
  double forces[9];
  (void)state;
  memset(held, 0x5a, sizeof(held));
  assert_int_equal(cpmdc_last_charges(NULL, 4), -1);
  assert_int_equal(cpmdc_last_charges(held, 0), -1);
  assert_int_equal(cpmdc_session_last_charges(NULL, held, 4), -1);
  assert_true(buffer_untouched(held, 4));

  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  free(params);
  assert_non_null(session);
  assert_int_equal(cpmdc_session_last_charges(session, held, 4), -1);
  assert_true(buffer_untouched(held, 4));

  if (cpmdc_available()) {
    cpmdc_session_destroy(session);
    return;
  }
  result = cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(result.ok, 1);
  assert_int_equal(cpmdc_last_property_snapshot(&prop), 0);
  assert_int_equal(prop.charges_valid, 0);
  assert_int_equal(prop.charge_count, 0);
  assert_int_equal(cpmdc_last_charges(held, 4), -1);
  assert_int_equal(cpmdc_session_last_charges(session, held, 3), -1);
  assert_true(buffer_untouched(held, 4));
  cpmdc_session_destroy(session);
}

static int pseudopotentials_present(void) {
  const char *dir = getenv("CPMDC_PSEUDO_DIR");
  if (!dir || !dir[0])
    dir = getenv("CPMD_PP_LIBRARY_PATH");
  return dir && dir[0];
}

static void assert_permutation(const double *original, const double *permuted) {
  assert_true(fabs(permuted[0] - original[1]) < k_perm_tol);
  assert_true(fabs(permuted[1] - original[0]) < k_perm_tol);
  assert_true(fabs(permuted[2] - original[2]) < k_perm_tol);
}

static void test_opencpmd_water_charges_sum_and_follow_atoms(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  CPMDCSession *other;
  CPMDCResult result;
  CPMDCPropertySnapshot prop;
  double forces[9];
  double charges[4];
  double permuted[3];
  double cached[3];
  double off[3];
  double energy;
  double sum;
  (void)state;
  if (!cpmdc_available() || !pseudopotentials_present())
    skip();

  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  result = cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(result.ok, 1);
  energy = result.energy_h;
  assert_int_equal(cpmdc_session_last_property_snapshot(session, &prop), 0);
  assert_int_equal(prop.charges_valid, 1);
  assert_int_equal(prop.charge_count, 3);
  memset(charges, 0x5a, sizeof(charges));
  assert_int_equal(cpmdc_session_last_charges(session, charges, 2), -1);
  assert_true(buffer_untouched(charges, 4));
  assert_int_equal(cpmdc_session_last_charges(session, charges, 4), 0);
  assert_true(((unsigned char *)&charges[3])[0] == 0x5a);
  sum = charges[0] + charges[1] + charges[2];
  assert_true(fabs(sum) < k_sum_tol);

  result = cpmdc_session_energy_forces(session, 3, k_perm, k_perm_z, forces);
  assert_int_equal(result.ok, 1);
  assert_int_equal(cpmdc_session_last_charges(session, cached, 3), 0);
  assert_permutation(charges, cached);

  other = cpmdc_session_create(params, n);
  assert_non_null(other);
  result = cpmdc_session_energy_forces(other, 3, k_perm, k_perm_z, forces);
  assert_int_equal(result.ok, 1);
  assert_int_equal(cpmdc_session_last_charges(other, permuted, 3), 0);
  assert_permutation(charges, permuted);
  cpmdc_session_destroy(other);
  cpmdc_session_destroy(session);

  assert_int_equal(setenv("CPMDC_CHARGE_SNAPSHOT", "0", 1), 0);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  result = cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  unsetenv("CPMDC_CHARGE_SNAPSHOT");
  assert_int_equal(result.ok, 1);
  assert_true(fabs(result.energy_h - energy) < 1.0e-8);
  assert_int_equal(cpmdc_session_last_property_snapshot(session, &prop), 0);
  assert_int_equal(prop.charges_valid, 0);
  memset(off, 0x5a, sizeof(off));
  assert_int_equal(cpmdc_session_last_charges(session, off, 3), -1);
  assert_true(buffer_untouched(off, 3));
  cpmdc_session_destroy(session);
  free(params);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_unset_charges_do_not_write_the_buffer),
      cmocka_unit_test(test_opencpmd_water_charges_sum_and_follow_atoms),
  };
  if (argc < 2)
    return 1;
  g_params = argv[1];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
