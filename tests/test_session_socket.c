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
static const char *g_step_b = NULL;
static const char *g_step_species = NULL;

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
  size_t params_size = 0, step_size = 0, step_b_size = 0, species_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step = read_file(g_step, &step_size);
  unsigned char *step_b = read_file(g_step_b, &step_b_size);
  unsigned char *step_species = read_file(g_step_species, &species_size);
  assert_non_null(params);
  assert_non_null(step);
  assert_non_null(step_b);
  assert_non_null(step_species);

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
  assert_int_equal(cpmdc_available(), 1);
  CPMDCResult eval = cpmdc_session_calculate_result(session, step, step_size,
                                                    out, need, &out_size);
  assert_int_equal(eval.ok, 1);
  assert_true(out_size > 0);
  assert_true(isfinite(eval.energy_h));
  double step_energy = eval.energy_h;

  size_t need_b =
      cpmdc_potential_result_size_for_force_input(step_b, step_b_size);
  assert_int_equal(need_b, need);
  out_size = 0;
  CPMDCResult eval_b = cpmdc_session_calculate_result(
      session, step_b, step_b_size, out, need, &out_size);
  assert_int_equal(eval_b.ok, 1);
  assert_true(isfinite(eval_b.energy_h));
  assert_true(fabs(eval_b.energy_h - step_energy) > 1.0e-8);

  out_size = 0;
  CPMDCResult repeated = cpmdc_session_calculate_result(
      session, step, step_size, out, need, &out_size);
  assert_int_equal(repeated.ok, 1);
  assert_float_equal(repeated.energy_h, step_energy, 1.0e-12);

  out_size = 0;
  CPMDCResult bad = cpmdc_session_calculate_result(
      session, step_species, species_size, out, need, &out_size);
  assert_int_equal(bad.ok, 0);
  assert_non_null(strstr(bad.message, "topology"));

  cpmdc_session_destroy(session);
  free(out);
  free(params);
  free(step);
  free(step_b);
  free(step_species);
}

int main(int argc, char **argv) {
  if (argc < 5) {
    fprintf(stderr, "usage: %s params.bin step_a.bin step_b.bin species.bin\n",
            argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  g_step_b = argv[3];
  g_step_species = argv[4];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_session_socket_contract),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
