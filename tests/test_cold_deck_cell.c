/**
 * The cold deck carries the caller's box to CPMD with twelve decimals. The
 * CELL line holds a and the ratios b/a and c/a; six decimals on a ratio
 * moved b by up to 7e-6 A in a 14 A cell.
 */
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);

static void test_orthorhombic_cell_keeps_twelve_decimals(void **state) {
  (void)state;
  const double pos[6] = {1.0, 1.0, 1.0, 1.74, 1.0, 1.0};
  const int z[2] = {1, 1};
  const double cell[9] = {13.37, 0.0, 0.0, 0.0, 14.123456789, 0.0,
                          0.0,   0.0, 15.5};
  char deck[8192];
  int len = 0;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   1);
  assert_non_null(strstr(deck, " CELL\n"));
  assert_non_null(strstr(deck, "13.370000000000"));
  assert_non_null(strstr(deck, "1.056354284892"));
  assert_non_null(strstr(deck, "1.159311892296"));
}

static void test_cell_vectors_keep_twelve_decimals(void **state) {
  (void)state;
  const double pos[6] = {1.0, 1.0, 1.0, 1.74, 1.0, 1.0};
  const int z[2] = {1, 1};
  const double cell[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                          0.0,  0.0, 12.0};
  char deck[8192];
  int len = 0;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   1);
  assert_non_null(strstr(deck, " CELL VECTORS\n"));
  assert_non_null(strstr(deck, "0.123456789012"));
  /* Every line stays within the 80 columns CPMD reads. */
  for (const char *line = deck; *line;) {
    const char *end = strchr(line, '\n');
    size_t n = end ? (size_t)(end - line) : strlen(line);
    assert_true(n <= 80);
    line += n + (end ? 1 : 0);
  }
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_orthorhombic_cell_keeps_twelve_decimals),
      cmocka_unit_test(test_cell_vectors_keep_twelve_decimals),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
