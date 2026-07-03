/**
 * PotentialConfig common-overlay lowering: session creation, evaluation, and
 * loud rejection of unsound combinations.
 */
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_overlay = NULL;
static const char *g_reject = NULL;
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
  if (n <= 0 || n > (16 << 20)) {
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

static void test_overlay_session_evaluates(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] embed shell not available\n");
    skip();
  }
  size_t config_size = 0, step_size = 0;
  unsigned char *config = read_file(g_overlay, &config_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(config);
  assert_non_null(step);

  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(session);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  CPMDCResult r =
      cpmdc_session_calculate_result(session, step, step_size, out, need,
                                     &wrote);
  assert_int_equal(r.ok, 1);
  assert_true(wrote > 0);

  cpmdc_session_destroy(session);
  free(out);
  free(step);
  free(config);
}

static void test_overlay_rejects_energy_tolerance(void **state) {
  (void)state;
  size_t config_size = 0;
  unsigned char *config = read_file(g_reject, &config_size);
  assert_non_null(config);
  assert_null(cpmdc_session_create_from_config(config, config_size));
  assert_non_null(strstr(cpmdc_last_error(), "scfEnergyToleranceEv"));
  assert_int_not_equal(cpmdc_configure(config, config_size), 0);
  free(config);
}

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "usage: %s OVERLAY_BIN REJECT_BIN STEP_BIN\n", argv[0]);
    return 2;
  }
  g_overlay = argv[1];
  g_reject = argv[2];
  g_step = argv[3];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_overlay_session_evaluates),
      cmocka_unit_test(test_overlay_rejects_energy_tolerance),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
