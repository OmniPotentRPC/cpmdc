#include "cpmdc_params.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmocka.h>

static const char *g_path;
static unsigned char *slurp(const char *p, size_t *n) {
  FILE *f = fopen(p, "rb");
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  long L = ftell(f);
  rewind(f);
  unsigned char *b = malloc((size_t)L);
  *n = fread(b, 1, (size_t)L, f);
  fclose(f);
  return b;
}
static void test_escape_tokens_in_deck(void **state) {
  (void)state;
  size_t sz = 0;
  unsigned char *msg = slurp(g_path, &sz);
  assert_non_null(msg);
  struct capn arena;
  CPMDParams_ptr root;
  assert_int_equal(cpmdc_params_root(msg, sz, &arena, &root), 0);
  char deck[CPMDC_BLOCKS];
  assert_int_equal(cpmdc_params_render_input_deck(root, deck, sizeof(deck)), 0);
  assert_non_null(strstr(deck, "&INFO"));
  assert_non_null(strstr(deck, "DEMO KEYWORD"));
  assert_non_null(strstr(deck, "&PROPERTIES"));
  assert_non_null(strstr(deck, "DIPOLE"));
  assert_non_null(strstr(deck, "&EXTE"));
  assert_non_null(strstr(deck, "O_MT_BLYP.psp"));
  assert_non_null(strstr(deck, "H_CVB_BLYP.psp"));
  assert_non_null(strstr(deck, "FUNCTIONAL BLYP"));
  cpmdc_params_release(&arena);
  free(msg);
}
int main(int argc, char **argv) {
  g_path = argc > 1 ? argv[1] : "params_escape_generic.bin";
  const struct CMUnitTest tests[] = {cmocka_unit_test(test_escape_tokens_in_deck)};
  return cmocka_run_group_tests(tests, NULL, NULL);
}
