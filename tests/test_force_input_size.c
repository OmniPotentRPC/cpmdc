#include "cpmdc.h"
#include "cpmdc_params.h"

#include <errno.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_step_a = NULL;
static const char *g_step_ev = NULL;

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

static void test_sizes_match_atom_count(void **state) {
  (void)state;
  size_t a_size = 0, ev_size = 0;
  unsigned char *a = read_file(g_step_a, &a_size);
  unsigned char *ev = read_file(g_step_ev, &ev_size);
  assert_non_null(a);
  assert_non_null(ev);

  size_t need_a = cpmdc_potential_result_size_for_force_input(a, a_size);
  size_t need_ev = cpmdc_potential_result_size_for_force_input(ev, ev_size);
  /* Extended PotentialResult is larger than legacy 32+8*nforces estimate (80 for 2 atoms). */
  assert_true(need_a >= 80);
  assert_true(need_ev >= 80);
  assert_int_equal(cpmdc_potential_result_size_for_force_input(NULL, 0), 0);

  free(a);
  free(ev);
}

static void test_flat_size_fits_large_system(void **state) {
  (void)state;
  /* 100 atoms: the writer emits forces AND gradient (both natoms*3), so the
   * advertised size must cover both plus the fixed property lists. */
  const size_t force_count = 300;
  size_t need = cpmdc_potential_result_flat_size(force_count);
  assert_true(need > 0);
  double *forces = (double *)calloc(force_count, sizeof(double));
  assert_non_null(forces);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  assert_int_equal(
      cpmdc_potential_result_write(-1.25, forces, force_count, out, need,
                                   &wrote),
      0);
  assert_true(wrote > 0);
  assert_true(wrote <= need);
  free(out);
  free(forces);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s step_a.bin step_ev.bin\n", argv[0]);
    return 2;
  }
  g_step_a = argv[1];
  g_step_ev = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_sizes_match_atom_count),
      cmocka_unit_test(test_flat_size_fits_large_system),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
