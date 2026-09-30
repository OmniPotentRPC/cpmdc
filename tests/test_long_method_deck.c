/**
 * A rendered method deck longer than the old fixed buffers must come back
 * whole. A caller buffer that cannot hold it fails and is left untouched.
 */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

/* 1 only when this binary is the OpenCPMD link. The reference build is 0. */
#if defined(CPMDC_HAS_CPMD)
#define CPMDC_AVAILABLE_WHEN_LINKED 1
#else
#define CPMDC_AVAILABLE_WHEN_LINKED 0
#endif

int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len, char *cpmd_root,
                           int cpmd_root_len);
int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);

enum { k_raw_len = 20000, k_mark_at = 18000, k_out_cap = 65536 };

static const char k_mark[] = "LONGDECKTAIL";

static char *make_raw(void) {
  char *raw = (char *)malloc((size_t)k_raw_len + 1u);
  size_t mark_len;
  size_t tail;
  if (!raw)
    return NULL;
  memset(raw, 'A', (size_t)k_raw_len);
  memcpy(raw, "&INFO\n", 6);
  mark_len = strlen(k_mark);
  memcpy(raw + k_mark_at, k_mark, mark_len);
  raw[k_mark_at + mark_len] = '\n';
  tail = 6;
  memcpy(raw + (k_raw_len - (int)tail), "\n&END\n", tail);
  raw[k_raw_len] = '\0';
  return raw;
}

static int pack_raw_params(const char *raw, unsigned char **out, size_t *out_size) {
  struct capn arena;
  capn_ptr root;
  CPMDParams_ptr params;
  struct CPMDParams view;
  CPMDInputSection_list sections;
  struct CPMDInputSection sec;
  size_t cap = 4096u;
  unsigned char *buffer = NULL;
  int written = -1;
  *out = NULL;
  *out_size = 0;
  capn_init_malloc(&arena);
  root = capn_root(&arena);
  params = new_CPMDParams(root.seg);
  memset(&view, 0, sizeof(view));
  view.functional.str = "BLYP";
  view.functional.len = 4;
  view.cutOffRy = 70.0;
  view.multiplicity = 1;
  view.task.str = "gradient";
  view.task.len = 8;
  sections = new_CPMDInputSection_list(root.seg, 1);
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_raw;
  sec.raw.str = raw;
  sec.raw.len = (int)strlen(raw);
  set_CPMDInputSection(&sec, sections, 0);
  view.inputSections = sections;
  write_CPMDParams(&view, params);
  if (capn_setp(root, 0, params.p) != 0) {
    capn_free(&arena);
    return -1;
  }
  for (int attempt = 0; attempt < 16 && written < 0; ++attempt) {
    unsigned char *next = (unsigned char *)realloc(buffer, cap);
    if (!next) {
      free(buffer);
      capn_free(&arena);
      return -1;
    }
    buffer = next;
    written = capn_write_mem(&arena, (uint8_t *)buffer, cap, 0);
    if (written < 0)
      cap *= 2u;
  }
  capn_free(&arena);
  if (written < 0) {
    free(buffer);
    return -1;
  }
  *out = buffer;
  *out_size = (size_t)written;
  return 0;
}

static void test_long_method_deck_is_kept(void **state) {
  char *raw;
  unsigned char *msg = NULL;
  size_t msg_n = 0;
  char *stored;
  char *cold;
  char tiny[8];
  char small[4];
  char functional[64];
  char root[128];
  double pos[6] = {0.0, 0.0, 0.0, 0.74, 0.0, 0.0};
  int z[2] = {1, 1};
  double cell[9] = {0};
  double cutoff = 0.0;
  int charge = 0;
  int mult = 0;
  int cold_len = 0;
  (void)state;
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
  raw = make_raw();
  assert_non_null(raw);
  assert_true(k_mark_at > 16384);
  assert_non_null(strstr(raw, k_mark));
  assert_int_equal(pack_raw_params(raw, &msg, &msg_n), 0);
  assert_int_equal(cpmdc_set_params(msg, msg_n), 0);

  memset(small, 0x5a, sizeof(small));
  assert_int_equal(cpmdc_embed_get_config(functional, (int)sizeof(functional),
                                          &cutoff, &charge, &mult, small,
                                          (int)sizeof(small), root,
                                          (int)sizeof(root)),
                   -1);
  assert_int_equal((unsigned char)small[0], 0x5a);

  stored = (char *)calloc(k_out_cap, 1);
  assert_non_null(stored);
  assert_int_equal(cpmdc_embed_get_config(functional, (int)sizeof(functional),
                                          &cutoff, &charge, &mult, stored,
                                          k_out_cap, root, (int)sizeof(root)),
                   0);
  assert_true(strlen(stored) > 16384);
  assert_non_null(strstr(stored, k_mark));
  assert_non_null(strstr(stored, "FUNCTIONAL"));

  memset(tiny, 0x5a, sizeof(tiny));
  cold_len = 7;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 0, tiny,
                                                 (int)sizeof(tiny), &cold_len),
                   0);
  assert_int_equal(cold_len, 0);
  assert_int_equal((unsigned char)tiny[0], 0x5a);

  cold = (char *)calloc(k_out_cap, 1);
  assert_non_null(cold);
  cold_len = 0;
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 0, cold,
                                                 k_out_cap, &cold_len),
                   1);
  assert_true(cold_len > 16384);
  assert_non_null(strstr(cold, k_mark));
  assert_true(cold_len < k_out_cap);
  assert_int_equal((unsigned char)cold[cold_len], 0);

  free(cold);
  free(stored);
  free(msg);
  free(raw);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_long_method_deck_is_kept),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
