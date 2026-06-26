#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <cmocka.h>

static void test_embed_shell_unavailable_without_cpmd(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 0);
  const char *version = cpmdc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "cpmdc/0.1.0"));
  /* Session may be created from params only when Cap'n Proto parses; without
   * a params buffer it must fail. */
  assert_null(cpmdc_session_create(NULL, 0));
  cpmdc_finalize();
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_embed_shell_unavailable_without_cpmd),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
