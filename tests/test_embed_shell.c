#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include <cmocka.h>

/* 1 only when this binary is the OpenCPMD link. The reference build is 0. */
#if defined(CPMDC_HAS_CPMD)
#define CPMDC_AVAILABLE_WHEN_LINKED 1
#else
#define CPMDC_AVAILABLE_WHEN_LINKED 0
#endif

static void test_embed_library_surface(void **state) {
  (void)state;
  const char *version = cpmdc_version();
  assert_non_null(version);
  assert_non_null(strstr(version, "cpmdc/"));
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
