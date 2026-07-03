#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <cmocka.h>

static const char *const required_abi_features[] = {
    "abi.cpmdc_set_params",
    "abi.cpmdc_energy_gradient",
    "abi.cpmdc_energy",
    "abi.cpmdc_energy_forces",
    "abi.cpmdc_session_create",
    "abi.cpmdc_session_set_params",
    "abi.cpmdc_session_destroy",
    "abi.cpmdc_session_energy_gradient",
    "abi.cpmdc_session_energy",
    "abi.cpmdc_session_energy_forces",
    "abi.cpmdc_last_energy_components",
    "abi.cpmdc_last_property_snapshot",
    "abi.cpmdc_last_md_trajectory_row",
    "abi.cpmdc_last_multi_state_energies",
    "abi.cpmdc_last_charge_integrals",
    "abi.cpmdc_session_calculate_forces",
    "abi.cpmdc_session_calculate_result",
    "abi.cpmdc_calculate_result",
    "abi.cpmdc_potential_result_size_for_force_input",
    "abi.cpmdc_version",
    "abi.cpmdc_abi_version",
    "abi.cpmdc_configure",
    "abi.cpmdc_session_create_from_config",
    "abi.cpmdc_session_configure",
    "abi.cpmdc_last_error",
    "abi.cpmdc_available",
    "abi.cpmdc_finalize",
    "abi.cpmdc_feature_count",
    "abi.cpmdc_feature_table",
    "abi.cpmdc_feature_find",
};

static void test_stub_reports_unavailable(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 0);
  const char *version = cpmdc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "stub"));
  assert_int_equal(cpmdc_abi_version(), CPMDC_ABI_VERSION);
  assert_int_not_equal(cpmdc_configure(NULL, 0), 0);
  assert_non_null(cpmdc_last_error());
  assert_non_null(strstr(cpmdc_last_error(), "stub"));
  assert_null(cpmdc_session_create_from_config(NULL, 0));
  assert_int_not_equal(cpmdc_session_configure(NULL, NULL, 0), 0);
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
  CPMDCEnergyComponents components;
  memset(&components, 0xab, sizeof(components));
  assert_int_equal(cpmdc_last_energy_components(&components), -1);
  assert_int_equal(components.valid, 0);
  assert_true(components.etot == 0.0);
  const CPMDCFeatureEntry *features = cpmdc_feature_table();
  assert_non_null(features);
  assert_true(cpmdc_feature_count() > 0);
  for (size_t i = 0;
       i < sizeof(required_abi_features) / sizeof(required_abi_features[0]);
       ++i) {
    const CPMDCFeatureEntry *feature =
        cpmdc_feature_find(required_abi_features[i]);
    assert_non_null(feature);
    assert_string_equal(feature->feature_id, required_abi_features[i]);
    assert_int_equal(feature->kind, CPMDC_FEATURE_ABI);
    assert_int_equal(feature->stub_applicable, 1);
  }
  cpmdc_finalize();
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_stub_reports_unavailable),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
