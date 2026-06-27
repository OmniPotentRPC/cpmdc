/**
 * Drive shipped cpmdc_params_* render paths on multi-section fixtures.
 * Wire types come from generated Potentials.capnp.h (capnpc-c); asserts use
 * read_CPMDParams / deck substrings only — no parallel config structs.
 */
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

static int text_equals(capn_text text, const char *expected) {
  size_t n = strlen(expected);
  return text.len == (int)n && text.str && memcmp(text.str, expected, n) == 0;
}

static void test_parser_renders_structured_input(void **state) {
  (void)state;

  size_t message_size = 0;
  unsigned char *message = read_file(g_params_path, &message_size);
  assert_non_null(message);

  struct capn arena;
  CPMDParams_ptr params_root;
  assert_int_equal(
      cpmdc_params_root(message, message_size, &arena, &params_root), 0);

  struct CPMDParams params;
  read_CPMDParams(&params, params_root);
  assert_true(text_equals(params.functional, "PBE"));
  assert_float_equal(params.cutOffRy, 90.0, 1e-12);
  assert_int_equal(params.charge, 0);
  assert_int_equal(params.multiplicity, 1);
  assert_int_equal(params.memoryMb, 2048);
  assert_true(text_equals(params.cpmdRoot, "/opt/cpmd"));
  assert_true(text_equals(params.scratchDir, "/scratch/cpmd"));
  assert_true(text_equals(params.permanentDir, "/perm/cpmd"));

  char deck[CPMDC_BLOCKS];
  assert_int_equal(
      cpmdc_params_render_input_deck(params_root, deck, sizeof(deck)), 0);
  assert_non_null(strstr(deck, "&CPMD"));
  assert_non_null(strstr(deck, "OPTIMIZE WAVEFUNCTION"));
  assert_non_null(strstr(deck, "CONVERGENCE ORBITALS"));
  assert_non_null(strstr(deck, "FILEPATH"));
  assert_non_null(strstr(deck, "  /perm/cpmd/"));
  assert_non_null(strstr(deck, "&SYSTEM"));
  assert_non_null(strstr(deck, "CUTOFF"));
  assert_non_null(strstr(deck, "90"));
  assert_non_null(strstr(deck, "&DFT"));
  assert_non_null(strstr(deck, "FUNCTIONAL PBE"));
  assert_non_null(strstr(deck, "&ATOMS"));
  assert_non_null(strstr(deck, "*H_MT_PBE.psp"));
  assert_non_null(strstr(deck, "*O_MT_PBE.psp"));
  assert_non_null(strstr(deck, "LMAX=S"));
  assert_non_null(strstr(deck, "&PIMD"));
  assert_non_null(strstr(deck, "PROC_NPROC"));
  assert_non_null(strstr(deck, "TEST BLOCK"));

  cpmdc_params_release(&arena);
  free(message);
}

static void test_parser_rejects_empty(void **state) {
  (void)state;
  struct capn arena;
  CPMDParams_ptr params_root;
  assert_int_equal(cpmdc_params_root(NULL, 0, &arena, &params_root), -1);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s CPMD_PARAMS.bin\n", argv[0]);
    return 2;
  }
  g_params_path = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_parser_renders_structured_input),
      cmocka_unit_test(test_parser_rejects_empty),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
