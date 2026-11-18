/**
 * Atomic numbers outside 1..118 are refused on every evaluation path with
 * the same message, before a session fixes its topology. A session that
 * refused them stays usable.
 */
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static unsigned char *g_params = NULL;
static size_t g_params_size = 0;

static const double kPositions[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572,
                                     -0.4692, 0.0, -0.7572, -0.4692};

static int load_params(const char *path) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return -1;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return -1;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return -1;
  }
  rewind(fp);
  g_params = (unsigned char *)malloc((size_t)n);
  if (!g_params || fread(g_params, 1, (size_t)n, fp) != (size_t)n) {
    fclose(fp);
    return -1;
  }
  fclose(fp);
  g_params_size = (size_t)n;
  return 0;
}

static void expect_range_error(CPMDCResult r) {
  assert_int_equal(r.ok, 0);
  assert_non_null(strstr(cpmdc_last_error(), "1..118"));
  assert_string_equal(cpmdc_last_error(), r.message);
}

static void test_session_refuses_out_of_range(void **state) {
  (void)state;
  static const int bad[][3] = {{0, 1, 1}, {8, 1, 119}, {-1, 1, 1}, {8, 1, 1 << 20}};
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  assert_non_null(session);
  double forces[9];
  for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i)
    expect_range_error(
        cpmdc_session_energy_forces(session, 3, kPositions, bad[i], forces));
  cpmdc_session_destroy(session);
}

static void test_global_path_refuses_out_of_range(void **state) {
  (void)state;
  static const int bad[3] = {8, 1, 0};
  double forces[9];
  expect_range_error(cpmdc_energy_forces(3, kPositions, bad, g_params,
                                         g_params_size, forces));
  expect_range_error(
      cpmdc_energy_gradient(3, kPositions, bad, g_params, g_params_size, forces));
}

#if !defined(CPMDC_HAS_CPMD)
static void test_refusal_leaves_the_topology_open(void **state) {
  (void)state;
  static const int bad[3] = {0, 1, 1};
  static const int good[3] = {8, 1, 1};
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  assert_non_null(session);
  double forces[9];
  expect_range_error(
      cpmdc_session_energy_forces(session, 3, kPositions, bad, forces));
  CPMDCResult first =
      cpmdc_session_energy_forces(session, 3, kPositions, good, forces);
  assert_int_equal(first.ok, 1);
  CPMDCResult again =
      cpmdc_session_energy_forces(session, 3, kPositions, good, forces);
  assert_int_equal(again.ok, 1);
  expect_range_error(
      cpmdc_session_energy_forces(session, 3, kPositions, bad, forces));
  cpmdc_session_destroy(session);
}
#endif

int main(int argc, char **argv) {
  if (argc != 2 || load_params(argv[1]) != 0) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_session_refuses_out_of_range),
      cmocka_unit_test(test_global_path_refuses_out_of_range),
#if !defined(CPMDC_HAS_CPMD)
      cmocka_unit_test(test_refusal_leaves_the_topology_open),
#endif
  };
  int status = cmocka_run_group_tests(tests, NULL, NULL);
  free(g_params);
  return status;
}
