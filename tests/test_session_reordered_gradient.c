/**
 * A cache hit that reorders identical atoms must reorder the gradient
 * snapshot with the forces. cpmdc_session_last_gradient reads that
 * snapshot, and it has to stay the negative of the returned forces.
 */
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};

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

static void assert_gradient_is_neg_forces(const double *gradient,
                                          const double *forces) {
  int i;
  for (i = 0; i < 9; ++i)
    assert_true(gradient[i] == -forces[i]);
}

static void test_reordered_identical_atoms(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  double first_forces[9];
  double swapped_forces[9];
  double gradient[9];
  double swapped_pos[9];
  size_t count = 0;
  CPMDCResult first;
  CPMDCResult swapped;
  (void)state;

  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  assert_non_null(session);

  first = cpmdc_session_energy_forces(session, 3, k_pos, k_z, first_forces);
  assert_int_equal(first.ok, 1);
  assert_string_equal(first.message, "ok");
  assert_int_equal(cpmdc_session_last_gradient(session, gradient, 9, &count),
                   0);
  assert_int_equal(count, 9);
  assert_gradient_is_neg_forces(gradient, first_forces);
  /* The two hydrogens are not the same force, so a missed swap is visible. */
  assert_int_not_equal(memcmp(first_forces + 3, first_forces + 6,
                              3 * sizeof(double)),
                       0);

  memcpy(swapped_pos, k_pos, 3 * sizeof(double));
  memcpy(swapped_pos + 3, k_pos + 6, 3 * sizeof(double));
  memcpy(swapped_pos + 6, k_pos + 3, 3 * sizeof(double));
  swapped = cpmdc_session_energy_forces(session, 3, swapped_pos, k_z,
                                        swapped_forces);
  assert_int_equal(swapped.ok, 1);
  assert_true(swapped.energy_h == first.energy_h);
  /* An empty message is the cache hit. A recompute would say "ok". */
  assert_int_equal(swapped.message[0], '\0');
  assert_memory_equal(swapped_forces + 3, first_forces + 6, 3 * sizeof(double));
  assert_memory_equal(swapped_forces + 6, first_forces + 3, 3 * sizeof(double));
  count = 0;
  assert_int_equal(cpmdc_session_last_gradient(session, gradient, 9, &count),
                   0);
  assert_int_equal(count, 9);
  assert_gradient_is_neg_forces(gradient, swapped_forces);

  /* Back to the stored order. The snapshot is rewritten from the stored
   * gradient, so this hit is not a second permutation of the last one. */
  swapped = cpmdc_session_energy_forces(session, 3, k_pos, k_z, swapped_forces);
  assert_int_equal(swapped.ok, 1);
  assert_int_equal(swapped.message[0], '\0');
  assert_memory_equal(swapped_forces, first_forces, sizeof(first_forces));
  count = 0;
  assert_int_equal(cpmdc_session_last_gradient(session, gradient, 9, &count),
                   0);
  assert_gradient_is_neg_forces(gradient, swapped_forces);

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
      cmocka_unit_test(test_reordered_identical_atoms),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
