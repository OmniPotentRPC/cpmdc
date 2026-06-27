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
