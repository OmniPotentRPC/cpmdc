/**
 * A non-finite QM position is an error. It must not hit the geometry
 * cache, and a non-finite result must not be stored for a later step.
 */
#include "cpmdc.h"

#include <math.h>
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
static const char k_msg[] = "positions contain a non-finite component";

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

static void copy_pos(double out[9]) { memcpy(out, k_pos, sizeof(k_pos)); }

static void test_nan_after_finite_is_rejected(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  double pos[9];
  double forces[9];
  CPMDCResult finite;
  CPMDCResult bad;
  CPMDCResult again;
  (void)state;

  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  copy_pos(pos);
  finite = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(finite.ok, 1);
  assert_true(isfinite(finite.energy_h));

  pos[0] = nan("");
  bad = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(bad.ok, 0);
  assert_string_equal(bad.message, k_msg);
  assert_string_equal(cpmdc_last_error(), k_msg);
  assert_false(isfinite(bad.energy_h) && bad.ok);

  copy_pos(pos);
  again = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(again.ok, 1);
  assert_true(isfinite(again.energy_h));
  assert_true(again.energy_h == finite.energy_h);

  cpmdc_session_destroy(session);
  free(params);
}

static void test_nan_first_is_not_stored(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *clean;
  CPMDCSession *session;
  double pos[9];
  double forces[9];
  CPMDCResult reference;
  CPMDCResult bad;
  CPMDCResult finite;
  (void)state;

  params = read_file(g_params, &n);
  assert_non_null(params);
  clean = cpmdc_session_create(params, n);
  assert_non_null(clean);
  reference = cpmdc_session_energy_forces(clean, 3, k_pos, k_z, forces);
  assert_int_equal(reference.ok, 1);
  assert_true(isfinite(reference.energy_h));
  cpmdc_session_destroy(clean);

  session = cpmdc_session_create(params, n);
  assert_non_null(session);
  copy_pos(pos);
  pos[4] = nan("");
  bad = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(bad.ok, 0);
  assert_string_equal(bad.message, k_msg);
  assert_false(isnan(bad.energy_h) && bad.ok);

  copy_pos(pos);
  pos[8] = INFINITY;
  bad = cpmdc_session_energy_forces(session, 3, pos, k_z, forces);
  assert_int_equal(bad.ok, 0);
  assert_string_equal(bad.message, k_msg);

  finite = cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
  assert_int_equal(finite.ok, 1);
  assert_true(isfinite(finite.energy_h));
  assert_true(finite.energy_h == reference.energy_h);

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
      cmocka_unit_test(test_nan_after_finite_is_rejected),
      cmocka_unit_test(test_nan_first_is_not_stored),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
