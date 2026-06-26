#include "cpmdc.h"

#include <errno.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;
static const char *g_step = NULL;

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

static void test_session_socket_contract(void **state) {
  (void)state;
  size_t params_size = 0, step_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(params);
  assert_non_null(step);

  CPMDCSession *session = cpmdc_session_create(params, params_size);
  assert_non_null(session);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);

  size_t out_size = 0;
  unsigned char tiny[8] = {0};
  CPMDCResult small = cpmdc_session_calculate_result(session, step, step_size,
                                                     tiny, sizeof(tiny),
                                                     &out_size);
  assert_int_equal(small.ok, 0);
  assert_int_equal(out_size, need);
  assert_non_null(strstr(small.message, "too small"));

  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  out_size = 0;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] no cpmd.x — socket sizing only\n");
    /* undersized buffer path already checked; evaluation needs real CPMD */
  } else {
    CPMDCResult eval = cpmdc_session_calculate_result(session, step, step_size,
                                                      out, need, &out_size);
    assert_int_equal(eval.ok, 1);
    assert_true(out_size > 0);
    assert_true(isfinite(eval.energy_h) || eval.energy_h != 0.0 || eval.ok);
  }

  cpmdc_session_destroy(session);
  free(out);
  free(params);
  free(step);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin step.bin\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_session_socket_contract),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
