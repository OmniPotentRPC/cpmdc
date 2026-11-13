/**
 * A session that already holds a result must reject a missing geometry.
 * The cached-gradient lookup used to read the position pointer before
 * checking it, so a second call with a null buffer aborted the process.
 */
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void test_cached_null_geometry(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  double positions[3] = {0.2, -0.4, 0.7};
  int atomic_numbers[1] = {8};
  double forces[3] = {0.0, 0.0, 0.0};
  CPMDCResult first;
  CPMDCResult missing_positions;
  CPMDCResult missing_numbers;
  (void)state;

  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);

  first = cpmdc_session_energy_forces(session, 1, positions, atomic_numbers,
                                      forces);
  assert_int_equal(first.ok, 1);

  missing_positions = cpmdc_session_energy_forces(session, 1, NULL,
                                                  atomic_numbers, forces);
  assert_int_equal(missing_positions.ok, 0);
  assert_non_null(strstr(missing_positions.message, "invalid arguments"));
  assert_string_equal(cpmdc_last_error(), missing_positions.message);

  missing_numbers =
      cpmdc_session_energy_forces(session, 1, positions, NULL, forces);
  assert_int_equal(missing_numbers.ok, 0);
  assert_non_null(strstr(missing_numbers.message, "invalid arguments"));

  cpmdc_session_destroy(session);
  free(params);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_cached_null_geometry),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
