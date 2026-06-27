#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include <cmocka.h>

static void test_embed_library_surface(void **state) {
  (void)state;
  const char *version = cpmdc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "cpmdc/0.1.0"));
  assert_int_equal(cpmdc_available(), 1);
  assert_null(cpmdc_session_create(NULL, 0));
  cpmdc_finalize();
  assert_int_equal(cpmdc_available(), 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_embed_library_surface),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
