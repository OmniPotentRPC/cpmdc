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
  assert_non_null(cpmdc_feature_find("abi.cpmdc_feature_find"));
  assert_non_null(cpmdc_feature_find("params.inputBlocks"));
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_feature_table_nonempty),
      cmocka_unit_test(test_inscan_sections),
      cmocka_unit_test(test_cp_keywords_not_sections),
      cmocka_unit_test(test_abi),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
