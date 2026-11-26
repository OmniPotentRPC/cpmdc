#include "ipi_extra.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <cmocka.h>

static void test_absent_dipole_is_empty(void **state) {
  char buf[CPMDC_IPI_DIPOLE_EXTRA_CAP];
  int32_t length = 7;
  (void)state;
  memset(buf, 'x', sizeof(buf));
  assert_int_equal(cpmdc_ipi_format_dipole(0, NULL, buf, sizeof(buf), &length),
                   0);
  assert_int_equal(length, 0);
  assert_int_equal(buf[0], 'x');
}

static void test_finite_dipole_is_json(void **state) {
  char buf[CPMDC_IPI_DIPOLE_EXTRA_CAP];
  int32_t length = 0;
  const double values[3] = {1.0, -2.5, 0.0};
  const char *expected =
      "{\"dipole\":[1.0000000000000000e+00,-2.5000000000000000e+00,"
      "0.0000000000000000e+00]}";
  (void)state;
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), 0);
  assert_string_equal(buf, expected);
  assert_int_equal(length, (int32_t)strlen(expected));
}

static void test_negative_zero_keeps_its_sign(void **state) {
  char buf[CPMDC_IPI_DIPOLE_EXTRA_CAP];
  int32_t length = 0;
  const double values[3] = {-0.0, 0.0, 1.0};
  (void)state;
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), 0);
  assert_string_equal(buf, "{\"dipole\":[-0.0000000000000000e+00,"
                           "0.0000000000000000e+00,1.0000000000000000e+00]}");
}

static void test_nonfinite_dipole_is_rejected(void **state) {
  char buf[CPMDC_IPI_DIPOLE_EXTRA_CAP];
  int32_t length = 4;
  double values[3] = {NAN, 0.0, 0.0};
  (void)state;
  memset(buf, 'x', sizeof(buf));
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), -1);
  assert_int_equal(length, 0);
  values[0] = INFINITY;
  length = 4;
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), -1);
  assert_int_equal(length, 0);
  values[0] = 0.0;
  values[2] = -INFINITY;
  length = 4;
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), -1);
  assert_int_equal(length, 0);
}

static void test_short_buffer_is_rejected(void **state) {
  char buf[8];
  int32_t length = 3;
  const double values[3] = {1.0, 2.0, 3.0};
  (void)state;
  assert_int_equal(
      cpmdc_ipi_format_dipole(1, values, buf, sizeof(buf), &length), -1);
  assert_int_equal(length, 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_absent_dipole_is_empty),
      cmocka_unit_test(test_finite_dipole_is_json),
      cmocka_unit_test(test_negative_zero_keeps_its_sign),
      cmocka_unit_test(test_nonfinite_dipole_is_rejected),
      cmocka_unit_test(test_short_buffer_is_rejected),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
