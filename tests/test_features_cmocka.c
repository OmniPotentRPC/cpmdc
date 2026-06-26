#include "cpmdc_features.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

static void test_feature_table_nonempty(void **state) {
  (void)state;
  assert_true(cpmdc_feature_count() > 50);
  assert_non_null(cpmdc_feature_table());
}

static void test_section_catalog(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("section.generic"));
  assert_non_null(cpmdc_feature_find("catalog.section.CPMD"));
  assert_non_null(cpmdc_feature_find("catalog.section.PIMD"));
  assert_non_null(cpmdc_feature_find("catalog.section.TDDFT"));
  assert_non_null(cpmdc_feature_find("catalog.section.QMMM"));
  assert_non_null(cpmdc_feature_find("catalog.section.VDW"));
}

static void test_cp_and_dft_keywords(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("catalog.cpmd.EMASS"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.MOLECULAR_DYNAMICS_CP"));
  assert_non_null(cpmdc_feature_find("catalog.dft.FUNCTIONAL_BLYP"));
  assert_non_null(cpmdc_feature_find("catalog.dft.FUNCTIONAL_PBE0"));
  assert_non_null(cpmdc_feature_find("catalog.dft.LSD"));
  assert_non_null(cpmdc_feature_find("catalog.dft.GC_CUTOFF"));
}

static void test_abi_and_params(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("params.cutOffRy"));
  assert_non_null(cpmdc_feature_find("params.inputBlocks"));
  assert_non_null(cpmdc_feature_find("abi.cpmdc_calculate_result"));
  assert_non_null(cpmdc_feature_find("abi.cpmdc_feature_find"));
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_feature_table_nonempty),
      cmocka_unit_test(test_section_catalog),
      cmocka_unit_test(test_cp_and_dft_keywords),
      cmocka_unit_test(test_abi_and_params),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
