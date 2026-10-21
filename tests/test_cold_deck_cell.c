/**
 * The caller's box is the cell. The cold deck carries the ForceInput box
 * to CPMD with twelve decimals and replaces any CELL the deck states:
 * SYMMETRY 1 for a cubic deck and a cubic box, SYMMETRY 8 for another
 * orthorhombic box, CELL VECTORS without SYMMETRY for a box with
 * off-diagonal components, Bohr when the deck has no ANGSTROM. An isolated
 * deck (SYMMETRY 0, the built-in minimal deck) takes only an orthorhombic
 * box. An all-zero box is no cell, and the deck's CELL stays.
 */
#include "cpmdc.h"
#include "cpmdc_embed_image.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

int cpmdc_embed_init(void);
int cpmdc_embed_set_deck(const char *deck, int deck_len,
                         CPMDCEmbedImage *image);
int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);
int cpmdc_embed_compose_cold_deck_image(int n_atoms,
                                        const double *positions_ang,
                                        const int *atomic_numbers,
                                        const double *cell_ang, int has_cell,
                                        char *deck_out, int deck_cap,
                                        int *deck_len, CPMDCEmbedImage *image);

static const double kPos[6] = {1.0, 1.0, 1.0, 1.74, 1.0, 1.0};
static const int kZ[2] = {1, 1};

/* A periodic deck with real atoms and its own 14 A cubic cell. */
static const char kPeriodicDeck[] =
    "&CPMD\n OPTIMIZE WAVEFUNCTION\n&END\n"
    "&SYSTEM\n SYMMETRY\n  1\n ANGSTROM\n CELL\n  14.0 1.0 1.0 0.0 0.0 0.0\n"
    " CHECK SYMMETRY OFF\n CUTOFF\n  70.0\n&END\n"
    "&DFT\n FUNCTIONAL BLYP\n&END\n"
    "&ATOMS\n*H_CVB_BLYP.psp KLEINMAN-BYLANDER\n LMAX=S\n 2\n"
    " 1.0 1.0 1.0\n 1.74 1.0 1.0\n&END\n";

static int compose_with_deck(const char *deck_in, const double *cell,
                             char *out, int cap) {
  CPMDCEmbedImage *image = (CPMDCEmbedImage *)calloc(1, sizeof(*image));
  int len = 0;
  assert_non_null(image);
  assert_int_equal(cpmdc_embed_set_deck(deck_in, (int)strlen(deck_in), image),
                   1);
  int ok = cpmdc_embed_compose_cold_deck_image(2, kPos, kZ, cell, 1, out, cap,
                                               &len, image);
  free(image->input_deck);
  free(image);
  return ok;
}

static void lines_within_80(const char *deck) {
  for (const char *line = deck; *line;) {
    const char *end = strchr(line, '\n');
    size_t n = end ? (size_t)(end - line) : strlen(line);
    assert_true(n <= 80);
    line += n + (end ? 1 : 0);
  }
}

static int setup(void **state) {
  (void)state;
  return cpmdc_embed_init() != 0 ? 0 : -1;
}

static void test_minimal_deck_orthorhombic_box(void **state) {
  (void)state;
  const double cell[9] = {13.37, 0.0, 0.0, 0.0, 14.123456789, 0.0,
                          0.0,   0.0, 15.5};
  char deck[8192];
  int len = 0;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, kPos, kZ, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   1);
  assert_non_null(strstr(deck, " SYMMETRY\n0\n CELL\n"));
  assert_non_null(strstr(deck, "13.370000000000"));
  assert_non_null(strstr(deck, "1.056354284892"));
  assert_non_null(strstr(deck, "1.159311892296"));
}

static void test_minimal_deck_refuses_a_skewed_box(void **state) {
  (void)state;
  const double cell[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                          0.0,  0.0, 12.0};
  char deck[8192];
  int len = 0;
  /* The minimal deck is an isolated cluster with a Hockney solver. */
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, kPos, kZ, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   0);
}

static void test_box_replaces_the_deck_cell(void **state) {
  (void)state;
  char deck[8192];
  const double cubic[9] = {14.5, 0.0, 0.0, 0.0, 14.5, 0.0, 0.0, 0.0, 14.5};
  assert_int_equal(compose_with_deck(kPeriodicDeck, cubic, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " SYMMETRY\n1\n CELL\n"));
  assert_non_null(strstr(deck, "14.500000000000"));
  assert_null(strstr(deck, "14.0 1.0 1.0"));
  /* CHECK SYMMETRY is not SYMMETRY, and the rest of the deck stays. */
  assert_non_null(strstr(deck, " CHECK SYMMETRY OFF\n"));
  assert_non_null(strstr(deck, " CUTOFF\n  70.0\n"));
  assert_non_null(strstr(deck, " 1.74 1.0 1.0\n"));

  const double ortho[9] = {14.5, 0.0, 0.0, 0.0, 15.0, 0.0, 0.0, 0.0, 16.0};
  assert_int_equal(compose_with_deck(kPeriodicDeck, ortho, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " SYMMETRY\n8\n CELL\n"));
  assert_non_null(strstr(deck, "1.034482758621"));

  const double skew[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                          0.0,  0.0, 12.0};
  assert_int_equal(compose_with_deck(kPeriodicDeck, skew, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " CELL VECTORS\n"));
  assert_non_null(strstr(deck, "0.123456789012"));
  assert_null(strstr(deck, " SYMMETRY\n"));
  lines_within_80(deck);
}

static void test_box_in_bohr_without_angstrom(void **state) {
  (void)state;
  char bohr_deck[sizeof(kPeriodicDeck)];
  char deck[8192];
  const char *a = strstr(kPeriodicDeck, " ANGSTROM\n");
  size_t head = (size_t)(a - kPeriodicDeck);
  memcpy(bohr_deck, kPeriodicDeck, head);
  strcpy(bohr_deck + head, a + strlen(" ANGSTROM\n"));
  const double cubic[9] = {14.5, 0.0, 0.0, 0.0, 14.5, 0.0, 0.0, 0.0, 14.5};
  assert_int_equal(compose_with_deck(bohr_deck, cubic, deck, (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, "27.401028809352"));
}

static void test_zero_box_keeps_the_deck_cell(void **state) {
  (void)state;
  char deck[8192];
  const double zero[9] = {0.0};
  assert_int_equal(compose_with_deck(kPeriodicDeck, zero, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " CELL\n  14.0 1.0 1.0 0.0 0.0 0.0\n"));
  assert_non_null(strstr(deck, " SYMMETRY\n  1\n"));
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_minimal_deck_orthorhombic_box),
      cmocka_unit_test(test_minimal_deck_refuses_a_skewed_box),
      cmocka_unit_test(test_box_replaces_the_deck_cell),
      cmocka_unit_test(test_box_in_bohr_without_angstrom),
      cmocka_unit_test(test_zero_box_keeps_the_deck_cell),
  };
  return cmocka_run_group_tests(tests, setup, NULL);
}
