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
  /* available() is 1 only when built with -Dwith_cpmd and libcpmd.a linked. */
  int av = cpmdc_available();
  if (!av) {
    print_message(
        "[  INFO   ] embed library built without libcpmd.a "
        "(meson -Dwith_cpmd=true -Dcpmd_root=...)\n");
  }
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
