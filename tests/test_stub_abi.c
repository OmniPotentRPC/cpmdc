#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <cmocka.h>

static void test_stub_reports_unavailable(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 0);
  const char *version = cpmdc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "stub"));
  assert_int_not_equal(cpmdc_set_params(NULL, 0), 0);
  CPMDCResult energy_result = cpmdc_energy(0, NULL, NULL, NULL, 0);
  assert_int_equal(energy_result.ok, 0);
  CPMDCResult result = cpmdc_energy_gradient(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(result.ok, 0);
  CPMDCResult forces_result = cpmdc_energy_forces(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(forces_result.ok, 0);
  assert_null(cpmdc_session_create(NULL, 0));
  assert_int_not_equal(cpmdc_session_set_params(NULL, NULL, 0), 0);
  cpmdc_session_destroy(NULL);
  CPMDCResult session_energy = cpmdc_session_energy(NULL, 0, NULL, NULL);
  assert_int_equal(session_energy.ok, 0);
  CPMDCResult session_gradient =
      cpmdc_session_energy_gradient(NULL, 0, NULL, NULL, NULL);
  assert_int_equal(session_gradient.ok, 0);
  CPMDCResult session_forces =
      cpmdc_session_energy_forces(NULL, 0, NULL, NULL, NULL);
  assert_int_equal(session_forces.ok, 0);
  CPMDCResult session_step =
      cpmdc_session_calculate_forces(NULL, NULL, 0, NULL, 0);
  assert_int_equal(session_step.ok, 0);
  CPMDCResult session_result =
      cpmdc_session_calculate_result(NULL, NULL, 0, NULL, 0, NULL);
  assert_int_equal(session_result.ok, 0);
  CPMDCResult one_shot =
      cpmdc_calculate_result(NULL, 0, NULL, 0, NULL, 0, NULL);
  assert_int_equal(one_shot.ok, 0);
  assert_int_equal(cpmdc_potential_result_size_for_force_input(NULL, 0), 0);
  cpmdc_finalize();
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_stub_reports_unavailable),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
