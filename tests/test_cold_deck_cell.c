/**
 * The caller's box is the cell. The cold deck carries the ForceInput box
 * to CPMD with twelve decimals and replaces any CELL the deck states:
 * SYMMETRY 1 for a cubic deck and a cubic box, SYMMETRY 8 for another
 * orthorhombic box, CELL VECTORS without SYMMETRY for a box with
 * off-diagonal components, Bohr when the deck has no ANGSTROM. An isolated
 * deck (SYMMETRY 0, the built-in minimal deck, or CLUSTER, POISSON SOLVER,
 * ISOLATED MOLECULE, or MOLECULE ISOLATED in &SYSTEM) uses CLUSTER with
 * CELL VECTORS for a non-orthorhombic box. A box that is not rank 3, has a
 * non-finite component, or has a volume below binary64 is rejected. An
 * all-zero box is no cell, and the deck's CELL stays.
 */
#include "cpmdc.h"
#include "cpmdc_embed_image.h"

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
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

/* Number of lines of deck that are exactly line, without the newline. */
static int count_lines(const char *deck, const char *line) {
  size_t n = strlen(line);
  int count = 0;
  for (const char *p = deck; *p;) {
    const char *end = strchr(p, '\n');
    size_t len = end ? (size_t)(end - p) : strlen(p);
    if (len == n && memcmp(p, line, n) == 0)
      ++count;
    p += len + (end ? 1 : 0);
  }
  return count;
}

/* Number of lines whose first token is CLUSTER in any case. */
static int count_cluster_lines(const char *deck) {
  int count = 0;
  for (const char *p = deck; *p;) {
    const char *end = strchr(p, '\n');
    size_t len = end ? (size_t)(end - p) : strlen(p);
    size_t i = 0;
    while (i < len && (p[i] == ' ' || p[i] == '\t'))
      ++i;
    size_t k = 0;
    while (k < 7 && i + k < len &&
           toupper((unsigned char)p[i + k]) == "CLUSTER"[k])
      ++k;
    if (k == 7 && (i + 7 == len || p[i + 7] == ' ' || p[i + 7] == '\t'))
      ++count;
    p += len + (end ? 1 : 0);
  }
  return count;
}

/* kPeriodicDeck with its SYMMETRY lines replaced by symmetry and with
 * before ahead of the deck and after behind it. */
static void periodic_deck_with(char *out, size_t cap, const char *before,
                               const char *symmetry, const char *after) {
  const char *sym = strstr(kPeriodicDeck, " SYMMETRY\n  1\n");
  assert_non_null(sym);
  int n = snprintf(out, cap, "%s%.*s%s%s%s", before,
                   (int)(sym - kPeriodicDeck), kPeriodicDeck, symmetry,
                   sym + strlen(" SYMMETRY\n  1\n"), after);
  assert_true(n > 0 && (size_t)n < cap);
}

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

static void assert_cell_vectors_equal(const char *deck,
                                      const double expected[9]) {
  const char *vectors = strstr(deck, " CELL VECTORS\n");
  double parsed[9];
  assert_non_null(vectors);
  vectors += strlen(" CELL VECTORS\n");
  assert_int_equal(
      sscanf(vectors, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &parsed[0],
             &parsed[1], &parsed[2], &parsed[3], &parsed[4], &parsed[5],
             &parsed[6], &parsed[7], &parsed[8]),
      9);
  assert_memory_equal(parsed, expected, sizeof(parsed));
}

static void assert_cell_equal(const char *deck, const double expected[3]) {
  const char *cell = strstr(deck, " CELL\n");
  double parsed[6];
  assert_non_null(cell);
  cell += strlen(" CELL\n");
  assert_int_equal(sscanf(cell, "%lf %lf %lf %lf %lf %lf", &parsed[0],
                          &parsed[1], &parsed[2], &parsed[3], &parsed[4],
                          &parsed[5]),
                   6);
  assert_memory_equal(parsed, expected, 3 * sizeof(double));
  assert_true(parsed[3] == 0.0 && parsed[4] == 0.0 && parsed[5] == 0.0);
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
  const double expected[3] = {cell[0], cell[4] / cell[0],
                              cell[8] / cell[0]};
  assert_cell_equal(deck, expected);
  lines_within_80(deck);
}

static void test_minimal_deck_uses_cluster_for_a_skewed_box(void **state) {
  (void)state;
  const double cell[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                          0.0,  0.0, 12.0};
  char deck[8192];
  int len = 0;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, kPos, kZ, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   1);
  assert_non_null(strstr(deck, " CLUSTER\n CELL VECTORS\n"));
  assert_cell_vectors_equal(deck, cell);
  assert_non_null(strstr(deck, " POISSON SOLVER HOCKNEY\n"));
  assert_null(strstr(deck, " SYMMETRY\n"));
  lines_within_80(deck);
}

static void assert_minimal_deck_rejects(const double cell[9],
                                        const char *message) {
  char deck[8192];
  int len = -1;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, kPos, kZ, cell, 1, deck,
                                                 (int)sizeof(deck), &len),
                   0);
  assert_int_equal(len, 0);
  assert_string_equal(cpmdc_last_error(), message);
}

static void test_minimal_deck_rejects_invalid_skewed_boxes(void **state) {
  (void)state;
  const double singular[9] = {10.0, 1.0, 0.0, 20.0, 2.0, 0.0,
                              0.0,  0.0, 12.0};
  const double nonfinite[9] = {10.0, NAN, 0.0, 0.0, 11.0, 0.0,
                               0.0,  0.0, 12.0};
  assert_minimal_deck_rejects(singular,
                              "cell does not have three independent vectors");
  assert_minimal_deck_rejects(nonfinite, "cell contains a non-finite component");
}

/* Positive diagonals pass the diagonal check, so the rank-3 state decides. */
static void test_minimal_deck_rejects_invalid_diagonal_boxes(void **state) {
  (void)state;
  const double underflow[9] = {1.0e-110, 0.0, 0.0, 0.0, 1.0e-110, 0.0,
                               0.0,      0.0, 1.0e-110};
  const double flat[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0e-20};
  const double nan_only[9] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
  const double nan_diagonal[9] = {NAN, 0.0, 0.0, 0.0, 11.0, 0.0,
                                  0.0, 0.0, 12.0};
  assert_minimal_deck_rejects(underflow, "cell volume is too small for binary64");
  assert_minimal_deck_rejects(nan_diagonal,
                              "cell contains a non-finite component");
  assert_minimal_deck_rejects(flat,
                              "cell does not have three independent vectors");
  assert_minimal_deck_rejects(nan_only, "cell contains a non-finite component");
}

static const double kOrtho[9] = {14.5, 0.0, 0.0, 0.0, 15.0,
                                 0.0,  0.0, 0.0, 16.0};
static const double kSkew[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                                0.0,  0.0, 12.0};

#ifndef CPMDC_HAS_CPMD
int cpmdc_embed_energy_grad(int n_atoms, const double *positions_ang,
                            const int *atomic_numbers, const double *cell_ang,
                            int has_cell, double *energy_h, double *grad_h_bohr,
                            int n_external, const double *external_positions_ang,
                            const double *external_charges_e,
                            const double *external_widths_ang,
                            double *external_forces_h_bohr,
                            CPMDCEmbedImage *image);
#endif

/* The reference evaluator computes a stress tensor for a periodic deck
 * and none for an isolated one. OpenCPMD builds have no reference
 * evaluator, so there is nothing to check. */
static void assert_reference_periodic(const char *deck_in, int periodic) {
#ifndef CPMDC_HAS_CPMD
  CPMDCEmbedImage *image = (CPMDCEmbedImage *)calloc(1, sizeof(*image));
  double energy = 0.0;
  double grad[6] = {0.0};
  assert_non_null(image);
  assert_int_equal(cpmdc_embed_set_deck(deck_in, (int)strlen(deck_in), image),
                   1);
  assert_int_equal(cpmdc_embed_energy_grad(2, kPos, kZ, kOrtho, 1, &energy, grad,
                                          0, NULL, NULL, NULL, NULL, image),
                   1);
  assert_int_equal(image->stress.valid, periodic);
  free(image->input_deck);
  free(image);
#else
  (void)deck_in;
  (void)periodic;
#endif
}

/* A deck without SYMMETRY 0 is isolated when &SYSTEM names CLUSTER, a
 * Poisson solver, or an isolated molecule, as the periodicity injection
 * and the reference evaluator read it. CLUSTER is a keyword in any case
 * and at any indent; a deck that has it keeps its one line. */
static void test_isolated_keyword_decks_stay_isolated(void **state) {
  (void)state;
  static const struct {
    const char *keyword;
    const char *cluster_line;
  } cases[] = {
      {" CLUSTER\n", " CLUSTER"},
      {" cluster\n", " cluster"},
      {"CLUSTER\n", "CLUSTER"},
      {"\tCluster\n", "\tCluster"},
      {" POISSON SOLVER TUCKERMAN\n", " CLUSTER"},
      {" ISOLATED MOLECULE\n", " CLUSTER"},
      {" MOLECULE ISOLATED\n", " CLUSTER"},
  };
  const double singular[9] = {10.0, 1.0, 0.0, 20.0, 2.0, 0.0,
                              0.0,  0.0, 12.0};
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    char deck_in[sizeof(kPeriodicDeck) + 64];
    char deck[8192];
    periodic_deck_with(deck_in, sizeof(deck_in), "", cases[i].keyword, "");

    assert_int_equal(
        compose_with_deck(deck_in, kOrtho, deck, (int)sizeof(deck)), 1);
    assert_non_null(strstr(deck, " SYMMETRY\n  0\n CELL\n"));
    assert_non_null(strstr(deck, cases[i].keyword));

    assert_int_equal(compose_with_deck(deck_in, kSkew, deck, (int)sizeof(deck)),
                     1);
    assert_int_equal(count_lines(deck, cases[i].cluster_line), 1);
    assert_int_equal(count_cluster_lines(deck), 1);
    assert_int_equal(count_lines(deck, " CELL VECTORS"), 1);
    assert_cell_vectors_equal(deck, kSkew);
    assert_null(strstr(deck, " SYMMETRY\n"));
    lines_within_80(deck);

    assert_int_equal(
        compose_with_deck(deck_in, singular, deck, (int)sizeof(deck)), 0);

    assert_reference_periodic(deck_in, 0);
  }
}

/* The cell a periodic deck gets from an orthorhombic and a skewed box:
 * SYMMETRY 8 and CELL, then CELL VECTORS with no CLUSTER line before it. */
static void assert_box_stays_periodic(const char *deck_in) {
  char deck[8192];
  assert_int_equal(compose_with_deck(deck_in, kOrtho, deck, (int)sizeof(deck)),
                   1);
  const char *sym = strstr(deck, " SYMMETRY\n");
  int code = -1;
  assert_non_null(sym);
  assert_int_equal(sscanf(sym + strlen(" SYMMETRY\n"), "%d", &code), 1);
  assert_int_equal(code, 8);
  assert_null(strstr(sym + 1, " SYMMETRY\n"));
  const double expected[3] = {kOrtho[0], kOrtho[4] / kOrtho[0],
                              kOrtho[8] / kOrtho[0]};
  assert_cell_equal(deck, expected);

  assert_int_equal(compose_with_deck(deck_in, kSkew, deck, (int)sizeof(deck)),
                   1);
  assert_int_equal(count_lines(deck, " CELL VECTORS"), 1);
  assert_null(strstr(deck, "\n CLUSTER\n CELL VECTORS\n"));
  assert_cell_vectors_equal(deck, kSkew);
  lines_within_80(deck);

  assert_reference_periodic(deck_in, 1);
}

/* CLUSTER names an isolated system only as a keyword of &SYSTEM. In the
 * comment text around the sections or in another section it is not one. */
static void test_cluster_outside_system_stays_periodic(void **state) {
  (void)state;
  static const struct {
    const char *before;
    const char *after;
  } cases[] = {
      {"Water CLUSTER, periodic run\n", ""},
      {"CLUSTER\n", ""},
      {" CLUSTER\n", ""},
      {"", " CLUSTER\n"},
      {"", "&PROP\n CLUSTER\n&END\n"},
  };
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    char deck_in[sizeof(kPeriodicDeck) + 128];
    periodic_deck_with(deck_in, sizeof(deck_in), cases[i].before,
                       " SYMMETRY\n  1\n", cases[i].after);
    assert_box_stays_periodic(deck_in);
  }

  /* A CLUSTER line in &CPMD, ahead of &SYSTEM. */
  char deck_in[sizeof(kPeriodicDeck) + 64];
  const char *opt = strstr(kPeriodicDeck, " OPTIMIZE WAVEFUNCTION\n");
  assert_non_null(opt);
  snprintf(deck_in, sizeof(deck_in), "%.*s CLUSTER\n%s",
           (int)(opt - kPeriodicDeck), kPeriodicDeck, opt);
  assert_box_stays_periodic(deck_in);
}

/* CHECK CLUSTER is not the CLUSTER keyword, with or without a SYMMETRY
 * line, as CHECK SYMMETRY is not SYMMETRY. */
static void test_check_cluster_is_not_cluster(void **state) {
  (void)state;
  static const char *const symmetry[] = {
      " SYMMETRY\n  1\n CHECK CLUSTER\n", " CHECK CLUSTER\n",
      " check cluster\n"};
  for (size_t i = 0; i < sizeof(symmetry) / sizeof(symmetry[0]); ++i) {
    char deck_in[sizeof(kPeriodicDeck) + 64];
    char deck[8192];
    periodic_deck_with(deck_in, sizeof(deck_in), "", symmetry[i], "");
    assert_box_stays_periodic(deck_in);
    assert_int_equal(
        compose_with_deck(deck_in, kSkew, deck, (int)sizeof(deck)), 1);
    assert_int_equal(count_cluster_lines(deck), 0);
  }
}

static void test_box_replaces_the_deck_cell(void **state) {
  (void)state;
  char deck[8192];
  const double cubic[9] = {14.5, 0.0, 0.0, 0.0, 14.5, 0.0, 0.0, 0.0, 14.5};
  assert_int_equal(compose_with_deck(kPeriodicDeck, cubic, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " SYMMETRY\n1\n CELL\n"));
  const double expected_cubic[3] = {cubic[0], 1.0, 1.0};
  assert_cell_equal(deck, expected_cubic);
  lines_within_80(deck);
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
  const double expected_ortho[3] = {ortho[0], ortho[4] / ortho[0],
                                    ortho[8] / ortho[0]};
  assert_cell_equal(deck, expected_ortho);
  lines_within_80(deck);

  const double skew[9] = {10.0, 0.0, 0.0, 0.123456789012, 11.0, 0.0,
                          0.0,  0.0, 12.0};
  assert_int_equal(compose_with_deck(kPeriodicDeck, skew, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_non_null(strstr(deck, " CELL VECTORS\n"));
  assert_cell_vectors_equal(deck, skew);
  assert_null(strstr(deck, " SYMMETRY\n"));
  lines_within_80(deck);
}

static void test_small_triclinic_components_are_preserved(void **state) {
  (void)state;
  char deck[8192];
  const double cell[9] = {1.0,      0.0, 0.0, 2.58e-12, 1.0,
                          0.0,      1.0e-13, 0.0, 1.0e-13};
  assert_int_equal(compose_with_deck(kPeriodicDeck, cell, deck,
                                     (int)sizeof(deck)),
                   1);
  assert_cell_vectors_equal(deck, cell);
  assert_null(strstr(deck, " SYMMETRY\n"));
  lines_within_80(deck);
}

static void test_triclinic_tilt_deck_rejects_transpose(void **state) {
  (void)state;
  char deck[8192];
  const double tilted[9] = {6.5, 0.0, 0.0, 1.2, 6.2, 0.0, -0.9, 0.8, 6.4};
  const double transposed[9] = {6.5, 1.2, -0.9, 0.0, 6.2, 0.8, 0.0, 0.0, 6.4};
  const double straight[9] = {6.5, 0.0, 0.0, 0.0, 6.2, 0.0, 0.0, 0.0, 6.4};
  const double expected_straight[3] = {6.5, 6.2 / 6.5, 6.4 / 6.5};
  double parsed[9];
  const char *vectors;

  assert_int_equal(
      compose_with_deck(kPeriodicDeck, tilted, deck, (int)sizeof(deck)), 1);
  assert_cell_vectors_equal(deck, tilted);
  assert_null(strstr(deck, " SYMMETRY\n"));
  lines_within_80(deck);

  assert_int_equal(
      compose_with_deck(kPeriodicDeck, transposed, deck, (int)sizeof(deck)), 1);
  vectors = strstr(deck, " CELL VECTORS\n");
  assert_non_null(vectors);
  vectors += strlen(" CELL VECTORS\n");
  assert_int_equal(
      sscanf(vectors, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &parsed[0],
             &parsed[1], &parsed[2], &parsed[3], &parsed[4], &parsed[5],
             &parsed[6], &parsed[7], &parsed[8]),
      9);
  assert_memory_not_equal(parsed, tilted, sizeof(parsed));

  assert_int_equal(
      compose_with_deck(kPeriodicDeck, straight, deck, (int)sizeof(deck)), 1);
  assert_null(strstr(deck, " CELL VECTORS\n"));
  assert_non_null(strstr(deck, " SYMMETRY\n8\n CELL\n"));
  assert_cell_equal(deck, expected_straight);
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
  assert_non_null(strstr(deck, "2.740102880935"));
  lines_within_80(deck);
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
      cmocka_unit_test(test_minimal_deck_uses_cluster_for_a_skewed_box),
      cmocka_unit_test(test_minimal_deck_rejects_invalid_skewed_boxes),
      cmocka_unit_test(test_minimal_deck_rejects_invalid_diagonal_boxes),
      cmocka_unit_test(test_isolated_keyword_decks_stay_isolated),
      cmocka_unit_test(test_cluster_outside_system_stays_periodic),
      cmocka_unit_test(test_check_cluster_is_not_cluster),
      cmocka_unit_test(test_box_replaces_the_deck_cell),
      cmocka_unit_test(test_small_triclinic_components_are_preserved),
      cmocka_unit_test(test_triclinic_tilt_deck_rejects_transpose),
      cmocka_unit_test(test_box_in_bohr_without_angstrom),
      cmocka_unit_test(test_zero_box_keeps_the_deck_cell),
  };
  return cmocka_run_group_tests(tests, setup, NULL);
}
