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
static const char *g_default_directives_path = NULL;

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

static void test_geometry_atoms_accept_silicon_symbol(void **state) {
  (void)state;
  size_t message_size = 0;
  unsigned char *message = read_file(g_params_path, &message_size);
  assert_non_null(message);

  struct capn arena;
  CPMDParams_ptr params_root;
  assert_int_equal(
      cpmdc_params_root(message, message_size, &arena, &params_root), 0);

  const double positions[3] = {0.25, 0.5, 0.75};
  const int atomic_numbers[1] = {14};
  const double cell[9] = {5.43, 0.0, 0.0, 0.0, 5.43, 0.0, 0.0, 0.0, 5.43};
  char deck[CPMDC_BLOCKS];
  assert_int_equal(cpmdc_params_render_deck_with_geometry(
                       params_root, 1, positions, atomic_numbers, cell, 1,
                       deck, sizeof(deck)),
                   0);
  assert_deck_has(deck, "&ATOMS");
  assert_deck_has(deck, "*Si_MT_PBE.psp");
  assert_deck_has(deck, "LMAX=D");
  assert_deck_has(deck, "   1");
  assert_deck_has(deck, "0.250000  0.500000  0.750000");

  const int missing_atomic_numbers[1] = {1};
  assert_int_equal(cpmdc_params_render_deck_with_geometry(
                       params_root, 1, positions, missing_atomic_numbers, cell,
                       1, deck, sizeof(deck)),
                   -1);

  cpmdc_params_release(&arena);
  free(message);
}

static void test_default_atoms_keep_directives(void **state) {
  (void)state;
  size_t message_size = 0;
  unsigned char *message = read_file(g_default_directives_path, &message_size);
  assert_non_null(message);

  struct capn arena;
  CPMDParams_ptr params_root;
  assert_int_equal(
      cpmdc_params_root(message, message_size, &arena, &params_root), 0);

  const double positions[9] = {
      0.000000, 0.000000, 0.000000,
      0.757160, 0.586260, 0.000000,
      -0.757160, 0.586260, 0.000000,
  };
  const int atomic_numbers[3] = {8, 1, 1};
  const double cell[9] = {10.0, 0.0, 0.0, 0.0, 10.0,
                          0.0,  0.0, 0.0, 10.0};
  char deck[CPMDC_BLOCKS];
  assert_int_equal(cpmdc_params_render_deck_with_geometry(
                       params_root, 3, positions, atomic_numbers, cell, 1,
                       deck, sizeof(deck)),
                   0);
  assert_deck_has(deck, "&ATOMS");
  assert_deck_has(deck, "*O_MT_BLYP.psp");
  assert_deck_has(deck, "*H_CVB_BLYP.psp");
  assert_deck_has(deck, "ISOLATED MOLECULE");

  cpmdc_params_release(&arena);
  free(message);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s silicon_params.bin default_atoms_params.bin\n",
            argv[0]);
    return 2;
  }
  g_params_path = argv[1];
  g_default_directives_path = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_geometry_atoms_accept_silicon_symbol),
      cmocka_unit_test(test_default_atoms_keep_directives),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
