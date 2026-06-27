#include "cpmdc_features.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

static void test_feature_table_nonempty(void **state) {
  (void)state;
  assert_true(cpmdc_feature_count() > 50);
}

static void test_inscan_sections(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("catalog.section.EAM"));
  assert_non_null(cpmdc_feature_find("catalog.section.MOLSTATES"));
  assert_non_null(cpmdc_feature_find("catalog.section.NLCC"));
  assert_non_null(cpmdc_feature_find("catalog.section.VECTORS"));
  assert_non_null(cpmdc_feature_find("catalog.section.PIMD"));
  assert_non_null(cpmdc_feature_find("catalog.section.PROP"));
  /* keyword-only must not be catalog.section.* */
  assert_null(cpmdc_feature_find("catalog.section.QMMM"));
  assert_null(cpmdc_feature_find("catalog.section.PROPERTIES"));
  assert_null(cpmdc_feature_find("catalog.section.BICANONICAL"));
  assert_null(cpmdc_feature_find("catalog.section.CDFT"));
}

static void test_cp_keywords_not_sections(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("catalog.cpmd.QMMM"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.EMASS"));
  assert_non_null(cpmdc_feature_find("catalog.dft.FUNCTIONAL_BLYP"));
}

static void test_abi(void **state) {
  (void)state;
  static const char *const abi_features[] = {
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
      "abi.cpmdc_session_calculate_forces",
      "abi.cpmdc_session_calculate_result",
      "abi.cpmdc_calculate_result",
      "abi.cpmdc_potential_result_size_for_force_input",
      "abi.cpmdc_version",
      "abi.cpmdc_available",
      "abi.cpmdc_finalize",
      "abi.cpmdc_feature_count",
      "abi.cpmdc_feature_table",
      "abi.cpmdc_feature_find",
  };

  const CPMDCFeatureEntry *table = cpmdc_feature_table();
  assert_non_null(table);
  assert_true(cpmdc_feature_count() >=
              sizeof(abi_features) / sizeof(abi_features[0]));

  for (size_t i = 0; i < sizeof(abi_features) / sizeof(abi_features[0]); i++) {
    const CPMDCFeatureEntry *entry = cpmdc_feature_find(abi_features[i]);
    assert_non_null(entry);
    assert_string_equal(entry->feature_id, abi_features[i]);
    assert_int_equal(entry->kind, CPMDC_FEATURE_ABI);
    assert_int_equal(entry->stub_applicable, 1);
    assert_int_equal(entry->embed_applicable, 1);
  }

  assert_null(cpmdc_feature_find("abi.cpmdc_missing"));
  assert_non_null(cpmdc_feature_find("params.inputBlocks"));
}

static void test_param_applicability(void **state) {
  (void)state;
  const CPMDCFeatureEntry *memory = cpmdc_feature_find("params.memoryMb");
  const CPMDCFeatureEntry *engine = cpmdc_feature_find("params.enginePath");
  const CPMDCFeatureEntry *scratch = cpmdc_feature_find("params.scratchDir");
  const CPMDCFeatureEntry *permanent = cpmdc_feature_find("params.permanentDir");
  assert_non_null(memory);
  assert_non_null(engine);
  assert_non_null(scratch);
  assert_non_null(permanent);
  assert_int_equal(memory->embed_applicable, 0);
  assert_int_equal(engine->embed_applicable, 0);
  assert_int_equal(scratch->embed_applicable, 1);
  assert_int_equal(permanent->embed_applicable, 1);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_feature_table_nonempty),
      cmocka_unit_test(test_inscan_sections),
      cmocka_unit_test(test_cp_keywords_not_sections),
      cmocka_unit_test(test_abi),
      cmocka_unit_test(test_param_applicability),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
