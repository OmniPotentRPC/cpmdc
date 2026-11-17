/**
 * The message of an SCF that stops short of the orbital threshold names the
 * minimizer of the last pass, the iteration it stopped at, MAXITER, and the
 * final gemax, and the first pass when the embed continued after it.
 */
#include "cpmdc_scf_report.h"

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static void test_one_pass(void **state) {
  (void)state;
  char buf[512];
  int n = cpmdc_format_scf_unconverged(buf, sizeof buf, CPMDC_SCF_MINIMIZER_ODIIS, 40,
                                       40, 3.25e-4, 1.0e-5, 0, 0);
  assert_true(n > 0 && (size_t)n < sizeof buf);
  assert_string_equal(buf, "CPMD SCF did not converge: ODIIS stopped at iteration 40 of "
                           "MAXITER 40 with gemax 3.250000E-04 above CONVERGENCE "
                           "ORBITALS 1.000000E-05");
}

static void test_fallback(void **state) {
  (void)state;
  char buf[512];
  int n = cpmdc_format_scf_unconverged(buf, sizeof buf, CPMDC_SCF_MINIMIZER_PCG_MINIMIZE,
                                       7, 50, 1.5e-3, 1.0e-6, CPMDC_SCF_MINIMIZER_ODIIS,
                                       50);
  assert_true(n > 0 && (size_t)n < sizeof buf);
  assert_string_equal(buf, "CPMD SCF did not converge: PCG MINIMIZE stopped at "
                           "iteration 7 of MAXITER 50 with gemax 1.500000E-03 above "
                           "CONVERGENCE ORBITALS 1.000000E-06, after ODIIS stopped at "
                           "iteration 50 of MAXITER 50");
}

static void test_gemax_reads_back(void **state) {
  (void)state;
  char buf[512];
  double gemax = 0.0;
  int maxiter = 0;
  cpmdc_format_scf_unconverged(buf, sizeof buf, CPMDC_SCF_MINIMIZER_PCG, 3, 3,
                               2.718281828e-2, 1.0e-5, 0, 0);
  const char *at = strstr(buf, "MAXITER ");
  assert_non_null(at);
  assert_int_equal(sscanf(at, "MAXITER %d with gemax %lf", &maxiter, &gemax), 2);
  assert_int_equal(maxiter, 3);
  assert_true(gemax > 2.7182815e-2 && gemax < 2.7182821e-2);
}

static void test_minimizer_names(void **state) {
  (void)state;
  assert_string_equal(cpmdc_scf_minimizer_name(CPMDC_SCF_MINIMIZER_ODIIS), "ODIIS");
  assert_string_equal(cpmdc_scf_minimizer_name(CPMDC_SCF_MINIMIZER_PCG), "PCG");
  assert_string_equal(cpmdc_scf_minimizer_name(CPMDC_SCF_MINIMIZER_PCG_MINIMIZE),
                      "PCG MINIMIZE");
  assert_string_equal(cpmdc_scf_minimizer_name(CPMDC_SCF_MINIMIZER_STEEPEST_DESCENT),
                      "STEEPEST DESCENT");
  assert_string_equal(cpmdc_scf_minimizer_name(CPMDC_SCF_MINIMIZER_UNKNOWN),
                      "an unknown minimizer");
  assert_string_equal(cpmdc_scf_minimizer_name(99), "an unknown minimizer");
}

/* The widest message, every number at its widest, fits the 512 bytes of
 * CPMDCResult.message and of the embed's detail buffer. */
static void test_fits_result_message(void **state) {
  (void)state;
  char buf[512];
  int n = cpmdc_format_scf_unconverged(
      buf, sizeof buf, CPMDC_SCF_MINIMIZER_STEEPEST_DESCENT, -2147483647 - 1,
      -2147483647 - 1, -1.0e-300, -1.0e-300, CPMDC_SCF_MINIMIZER_UNKNOWN + 99,
      -2147483647 - 1);
  assert_true(n > 0 && (size_t)n < sizeof buf);
}

static void test_truncates_safely(void **state) {
  (void)state;
  char small[16];
  memset(small, 'x', sizeof small);
  int n = cpmdc_format_scf_unconverged(small, sizeof small, CPMDC_SCF_MINIMIZER_ODIIS, 1,
                                       1, 1.0, 1.0, 0, 0);
  assert_true((size_t)n >= sizeof small);
  assert_int_equal(small[sizeof small - 1], '\0');
  assert_int_equal(strncmp(small, "CPMD SCF did no", sizeof small - 1), 0);
  assert_int_equal(cpmdc_format_scf_unconverged(NULL, 0, 1, 1, 1, 1.0, 1.0, 0, 0), 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_one_pass),
      cmocka_unit_test(test_fallback),
      cmocka_unit_test(test_gemax_reads_back),
      cmocka_unit_test(test_minimizer_names),
      cmocka_unit_test(test_fits_result_message),
      cmocka_unit_test(test_truncates_safely),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
