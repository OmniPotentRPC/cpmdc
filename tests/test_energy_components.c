/**
 * In-process ener_com snapshot via cpmdc_last_energy_components().
 * Reference PEF (default build): valid + etot == energy_h.
 * Live OpenCPMD (CPMDC_HAS_CPMD): also requires a non-total DFT/PP term non-zero.
 */
#include "cpmdc.h"

#include <math.h>
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

static void test_last_energy_components_after_eval(void **state) {
  (void)state;
  size_t params_size = 0;
  size_t step_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(params);
  assert_non_null(step);

  CPMDCEnergyComponents before;
  memset(&before, 0, sizeof(before));
  /* May be -1 before any successful eval in a fresh process. */
  (void)cpmdc_last_energy_components(&before);

  CPMDCSession *session = cpmdc_session_create(params, params_size);
  assert_non_null(session);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *result_buf = (unsigned char *)malloc(need);
  assert_non_null(result_buf);
  size_t wrote = 0;
  CPMDCResult r = cpmdc_session_calculate_result(session, step, step_size,
                                                 result_buf, need, &wrote);
  assert_int_equal(r.ok, 1);
  assert_true(isfinite(r.energy_h));

  CPMDCEnergyComponents c;
  memset(&c, 0, sizeof(c));
  assert_int_equal(cpmdc_last_energy_components(&c), 0);
  assert_int_equal(c.valid, 1);
  assert_true(isfinite(c.etot));
  assert_true(fabs(c.etot - r.energy_h) < 1e-12);

  /* At least one non-total field is finite (always); live embed also non-zero PP/XC. */
  assert_true(isfinite(c.ekin) && isfinite(c.exc) && isfinite(c.eht) &&
              isfinite(c.epseu) && isfinite(c.enl));
#if defined(CPMDC_HAS_CPMD)
  double non_total = fabs(c.ekin) + fabs(c.exc) + fabs(c.eht) + fabs(c.epseu) +
                     fabs(c.enl) + fabs(c.ehii) + fabs(c.esr);
  assert_true(non_total > 1e-8);
  printf("live_ener_com etot=%.12f ekin=%.12f exc=%.12f eht=%.12f epseu=%.12f "
         "enl=%.12f\n",
         c.etot, c.ekin, c.exc, c.eht, c.epseu, c.enl);
#else
  printf("ref_pef_ener_com etot=%.12f (components zero except etot)\n", c.etot);
#endif

  /* Second eval is consistent (same geometry). */
  CPMDCResult r2 = cpmdc_session_calculate_result(session, step, step_size,
                                                  result_buf, need, &wrote);
  assert_int_equal(r2.ok, 1);
  CPMDCEnergyComponents c2;
  assert_int_equal(cpmdc_last_energy_components(&c2), 0);
  assert_true(fabs(c2.etot - c.etot) < 1e-10);

  free(result_buf);
  cpmdc_session_destroy(session);
  free(params);
  free(step);
  cpmdc_finalize();
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin force_input.bin\n", argv[0]);
    return 1;
  }
  g_params = argv[1];
  g_step = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_last_energy_components_after_eval),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
