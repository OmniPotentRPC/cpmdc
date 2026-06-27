/**
 * Structured CPMD sections define the effective embed scalar configuration.
 */
#include "cpmdc.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmocka.h>

static const char *g_top = NULL;
static const char *g_sections = NULL;
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

static CPMDCResult eval_result(const unsigned char *params, size_t params_size,
                               const unsigned char *step, size_t step_size) {
  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t out_size = 0;
  CPMDCResult result = cpmdc_calculate_result(
      params, params_size, step, step_size, out, need, &out_size);
  assert_true(out_size == 0 || out_size <= need);
  free(out);
  return result;
}

static void test_sections_override_top_level_scalars(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }

  size_t top_size = 0, section_size = 0, step_size = 0;
  unsigned char *top = read_file(g_top, &top_size);
  unsigned char *sections = read_file(g_sections, &section_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(top);
  assert_non_null(sections);
  assert_non_null(step);

  CPMDCResult r_top = eval_result(top, top_size, step, step_size);
  CPMDCResult r_sections = eval_result(sections, section_size, step, step_size);
  assert_int_equal(r_top.ok, 1);
  assert_int_equal(r_sections.ok, 1);
  assert_true(isfinite(r_top.energy_h));
  assert_true(isfinite(r_sections.energy_h));
  assert_true(fabs(r_sections.energy_h - r_top.energy_h) <= 1e-10);

  free(top);
  free(sections);
  free(step);
  cpmdc_finalize();
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s top.bin sections.bin force_input.bin\n", argv[0]);
    return 2;
  }
  g_top = argv[1];
  g_sections = argv[2];
  g_step = argv[3];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_sections_override_top_level_scalars),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
