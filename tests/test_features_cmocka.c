#include "cpmdc_features.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

static void test_feature_table_nonempty(void **state) {
  (void)state;
  assert_true(cpmdc_feature_count() > 10);
  assert_non_null(cpmdc_feature_table());
}

static void test_section_kinds_present(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("section.generic"));
  assert_non_null(cpmdc_feature_find("section.raw"));
  assert_non_null(cpmdc_feature_find("section.atoms"));
  const CPMDCFeatureEntry *set = cpmdc_feature_find("section.set");
  assert_non_null(set);
  assert_int_equal(set->embed_applicable, 0);
}

static void test_abi_and_params(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("params.cutOffRy"));
  assert_non_null(cpmdc_feature_find("params.inputBlocks"));
  assert_non_null(cpmdc_feature_find("abi.cpmdc_calculate_result"));
  assert_non_null(cpmdc_feature_find("abi.cpmdc_available"));
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_feature_table_nonempty),
      cmocka_unit_test(test_section_kinds_present),
      cmocka_unit_test(test_abi_and_params),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
