/**
 * Capabilities discovery: size query, round-trip decode, and stub semantics
 * (same operation surface as the embed build with available = false).
 */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const capn_text no_text = {0, "", 0};

static int text_equals(capn_text text, const char *expected) {
  size_t expected_len = strlen(expected);
  return text.len == (int)expected_len && text.str &&
         memcmp(text.str, expected, expected_len) == 0;
}

static void test_capabilities_round_trip(void **state) {
  (void)state;
  size_t required = 0;
  assert_int_equal(cpmdc_capabilities_result(NULL, 0, &required), -1);
  assert_true(required > 0);

  unsigned char *buf = (unsigned char *)malloc(required);
  assert_non_null(buf);
  size_t wrote = 0;
  assert_int_equal(cpmdc_capabilities_result(buf, required, &wrote), 0);
  assert_int_equal(wrote, required);

  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, buf, wrote, 0), 0);
  Capabilities_ptr caps;
  caps.p = capn_getp(capn_root(&arena), 0, 1);
  assert_int_equal(caps.p.type, CAPN_STRUCT);
  struct Capabilities view;
  read_Capabilities(&view, caps);

  assert_true(text_equals(view.backendName, "cpmdc"));
  assert_true(view.backendVersion.len > 0);
  assert_int_equal(view.abiVersion, cpmdc_abi_version());
  assert_int_equal(view.available ? 1 : 0, cpmdc_available() ? 1 : 0);
  assert_true(view.schemaVersion.len > 0);

  capn_resolve(&view.operations.p);
  assert_int_equal(view.operations.p.len, 3);
  assert_int_equal(capn_get16(view.operations, 0),
                   Capabilities_Operation_energy);
  assert_int_equal(capn_get16(view.operations, 1),
                   Capabilities_Operation_forces);
  assert_int_equal(capn_get16(view.operations, 2),
                   Capabilities_Operation_gradient);

  capn_resolve(&view.loweredCommonFields);
  assert_int_equal(view.loweredCommonFields.len, 11);
  assert_true(text_equals(capn_get_text(view.loweredCommonFields, 0, no_text),
                          "xcFunctionals"));
  assert_true(text_equals(capn_get_text(view.loweredCommonFields, 10, no_text),
                          "vanDerWaalsS6"));
  /* relativityMethod stays out: the overlay rejects it with last_error. */
  for (int i = 0; i < view.loweredCommonFields.len; ++i)
    assert_false(text_equals(
        capn_get_text(view.loweredCommonFields, i, no_text),
        "relativityMethod"));

  capn_resolve(&view.configKinds);
  assert_int_equal(view.configKinds.len, 1);
  assert_true(
      text_equals(capn_get_text(view.configKinds, 0, no_text), "cpmd"));

  capn_free(&arena);
  free(buf);
}

static void test_capabilities_rejects_small_buffer(void **state) {
  (void)state;
  size_t required = 0;
  assert_int_equal(cpmdc_capabilities_result(NULL, 0, &required), -1);
  unsigned char small[8];
  size_t wrote = 0;
  assert_int_equal(cpmdc_capabilities_result(small, sizeof(small), &wrote),
                   -1);
  assert_int_equal(wrote, required);
  assert_int_equal(cpmdc_capabilities_result(NULL, 0, NULL), -1);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_capabilities_round_trip),
      cmocka_unit_test(test_capabilities_rejects_small_buffer),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
