#include "cpmdc_features.h"
#include "cpmdc_params.h"

#include <errno.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params_path = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp) {
    fprintf(stderr, "open failed for %s: %s\n", path, strerror(errno));
    return NULL;
  }
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

static void assert_deck_has(const char *deck, const char *needle) {
  if (!strstr(deck, needle))
    fail_msg("missing token [%s] in deck:\n%s", needle, deck);
}

static int count_occurrences(const char *haystack, const char *needle) {
  int count = 0;
  size_t needle_len = strlen(needle);
  const char *cursor = haystack;
  while ((cursor = strstr(cursor, needle)) != NULL) {
    ++count;
    cursor += needle_len;
  }
  return count;
}

static void assert_known_sections_are_unique(const char *deck) {
  assert_int_equal(count_occurrences(deck, "&CPMD\n"), 1);
  assert_int_equal(count_occurrences(deck, "&SYSTEM\n"), 1);
  assert_int_equal(count_occurrences(deck, "&DFT\n"), 1);
}

static void test_set_sections_render_into_target_sections(void **state) {
  (void)state;
  const CPMDCFeatureEntry *set = cpmdc_feature_find("section.set");
  assert_non_null(set);
  assert_int_equal(set->embed_applicable, 1);

  size_t message_size = 0;
  unsigned char *message = read_file(g_params_path, &message_size);
  assert_non_null(message);

  struct capn arena;
  CPMDParams_ptr params_root;
  assert_int_equal(
      cpmdc_params_root(message, message_size, &arena, &params_root), 0);

  char deck[CPMDC_BLOCKS];
  assert_int_equal(
      cpmdc_params_render_input_deck(params_root, deck, sizeof(deck)), 0);
  assert_deck_has(deck, "&CPMD");
  assert_deck_has(deck, "PRINT FORCES ON");
  assert_deck_has(deck, "MAXSTEP");
  assert_deck_has(deck, "  3");
  assert_deck_has(deck, "&SYSTEM");
  assert_deck_has(deck, "POISSON SOLVER");
  assert_deck_has(deck, "  TUCKERMAN");
  assert_deck_has(deck, "&DFT");
  assert_deck_has(deck, "GC-CUTOFF");
  assert_deck_has(deck, "  1.0e-7");
  assert_known_sections_are_unique(deck);

  double positions[6] = {0.0, 0.0, 0.7414, 0.0, 0.0, 0.0};
  int atomic_numbers[2] = {8, 1};
  double cell[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  char geometry_deck[CPMDC_BLOCKS];
  assert_int_equal(cpmdc_params_render_deck_with_geometry(
                       params_root, 2, positions, atomic_numbers, cell, 1,
                       geometry_deck, sizeof(geometry_deck)),
                   0);
  assert_deck_has(geometry_deck, "OPTIMIZE WAVEFUNCTION");
  assert_deck_has(geometry_deck, "PRINT FORCES ON");
  assert_deck_has(geometry_deck, "MAXSTEP");
  assert_deck_has(geometry_deck, "POISSON SOLVER");
  assert_deck_has(geometry_deck, "GC-CUTOFF");
  assert_known_sections_are_unique(geometry_deck);

  cpmdc_params_release(&arena);
  free(message);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s CPMD_PARAMS.bin\n", argv[0]);
    return 2;
  }
  g_params_path = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_set_sections_render_into_target_sections),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
