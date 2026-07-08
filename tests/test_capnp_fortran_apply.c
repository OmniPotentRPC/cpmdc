/**
 * Drive the shipped set_params / session_create path and observe embed knobs
 * via cpmdc_embed_get_config (capnp-fortran apply-from-bytes).
 */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

/* Real bind(C) getter from cpmdc_embed_apply_params.f90 */
int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len, char *cpmd_root,
                           int cpmd_root_len);
/* Cold-deck compose (method merge / real atoms / minimal) without SCF. */
int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);

static const char *g_top = NULL;
static const char *g_sections = NULL;
static const char *g_parser = NULL;
static const char *g_atoms_extras = NULL;
static const char *g_method_only = NULL;

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

static void read_applied(char *functional, size_t fsz, double *cutoff,
                         int *charge, int *mult, char *deck, size_t dsz,
                         char *root, size_t rsz) {
  memset(functional, 0, fsz);
  memset(deck, 0, dsz);
  memset(root, 0, rsz);
  *cutoff = 0.0;
  *charge = 0;
  *mult = 0;
  assert_int_equal(cpmdc_embed_get_config(functional, (int)fsz, cutoff, charge,
                                          mult, deck, (int)dsz, root, (int)rsz),
                   0);
}

static void test_set_params_applies_top_level_via_fortran(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 1);
  size_t n = 0;
  unsigned char *msg = read_file(g_top, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);

  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  assert_string_equal(functional, "PBE0");
  assert_true(fabs(cutoff - 80.0) < 1e-12);
  assert_int_equal(charge, 2);
  assert_int_equal(mult, 3);
  assert_non_null(strstr(deck, "CUTOFF"));
  assert_non_null(strstr(deck, "80"));
  free(msg);
}

static void test_set_params_applies_section_overrides_via_fortran(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 1);
  size_t n = 0;
  unsigned char *msg = read_file(g_sections, &n);
  assert_non_null(msg);
  /* Top-level is LDA/10/-1/1; system+dft sections override to PBE0/80/2/3. */
  assert_int_equal(cpmdc_set_params(msg, n), 0);

  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  assert_string_equal(functional, "PBE0");
  assert_true(fabs(cutoff - 80.0) < 1e-12);
  assert_int_equal(charge, 2);
  assert_int_equal(mult, 3);
  free(msg);
}

static void test_session_create_applies_parser_fixture(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 1);
  size_t n = 0;
  unsigned char *msg = read_file(g_parser, &n);
  assert_non_null(msg);
  CPMDCSession *session = cpmdc_session_create(msg, n);
  assert_non_null(session);

  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  assert_string_equal(functional, "PBE");
  assert_true(fabs(cutoff - 90.0) < 1e-12);
  assert_string_equal(root, "/opt/cpmd");
  assert_non_null(strstr(deck, "FUNCTIONAL"));
  assert_non_null(strstr(deck, "PBE"));

  cpmdc_session_destroy(session);
  free(msg);
}


/* wdwj: configure path stores Cap'n-rendered deck with typed atoms/DFT extras
 * so cold OpenCPMD path can consume applied_input_deck (&ATOMS present).
 * 3ba9: long-tail typed sections (&VDW/&PROP/&LINRES/&PIMD/&PATH/&TDDFT) land
 * in applied deck via shipped render. */
static void test_set_params_stores_typed_section_deck(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 1);
  assert_non_null(g_atoms_extras);
  size_t n = 0;
  unsigned char *msg = read_file(g_atoms_extras, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);

  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  assert_string_equal(functional, "PBE");
  /* Typed section tokens from shipped render (params_atoms_extras fixture). */
  assert_non_null(strstr(deck, "CONSTRAINTS"));
  assert_non_null(strstr(deck, "ISOTOPE"));
  assert_non_null(strstr(deck, "VELOCITIES"));
  assert_non_null(strstr(deck, "DUMMY ATOMS"));
  assert_non_null(strstr(deck, "HUBBARD U") || strstr(deck, "HUBBARD"));
  assert_non_null(strstr(deck, "&ATOMS") || strstr(deck, "&atoms"));
  assert_non_null(strstr(deck, "FUNCTIONAL"));
  /* Real PP line present → cold compose uses deck as-is (not method merge). */
  assert_non_null(strstr(deck, "O_MT_BLYP.psp") || strstr(deck, "*O"));
  /* 3ba9 long-tail tokens from shipped cpmdc_params.c render. */
  assert_non_null(strstr(deck, "EMPIRICAL CORRECTION"));
  assert_non_null(strstr(deck, "GRIMME"));
  assert_non_null(strstr(deck, "DIPOLE MOMENT"));
  assert_non_null(strstr(deck, "LOCALIZE"));
  assert_non_null(strstr(deck, "HTHRS"));
  assert_true(strstr(deck, "TROTTER DIMENSION") != NULL ||
              strstr(deck, "REPLICA NUMBER") != NULL);
  assert_non_null(strstr(deck, "TAMM-DANCOFF"));

  /* wdwj: real-PP applied deck survives compose (geometry merge must not run). */
  {
    double pos[6] = {0.0, 0.0, 0.0, 0.74, 0.0, 0.0};
    int z[2] = {1, 1};
    double cell[9] = {0};
    char cold[16384];
    int cold_len = 0;
    memset(cold, 0, sizeof(cold));
    assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 0, cold,
                                                   (int)sizeof(cold), &cold_len),
                     1);
    assert_true(cold_len > 0);
    assert_non_null(strstr(cold, "O_MT_BLYP.psp") || strstr(cold, "*O"));
    assert_non_null(strstr(cold, "EMPIRICAL CORRECTION"));
    assert_non_null(strstr(cold, "DIPOLE MOMENT"));
    /* Must not replace real PP deck with H-only geometry merge. */
    assert_null(strstr(cold, "H_CVB_BLYP.psp"));
  }
  free(msg);
}


/* wdwj: method-only Cap'n deck (no real PP under &ATOMS) keeps typed method
 * text; cold compose strips empty &ATOMS placeholder and merges geometry. */
static void test_set_params_method_only_keeps_dft_section(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), 1);
  assert_non_null(g_method_only);
  size_t n = 0;
  unsigned char *msg = read_file(g_method_only, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);
  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  assert_string_equal(functional, "PBE");
  assert_true(strstr(deck, "&DFT") != NULL || strstr(deck, "&dft") != NULL);
  assert_non_null(strstr(deck, "FUNCTIONAL"));
  assert_non_null(strstr(deck, "PBE"));
  /* C render may emit empty &ATOMS/&END; no real PP (*) until geometry merge. */
  assert_null(strstr(deck, "*H_"));
  assert_null(strstr(deck, "*.psp"));
  assert_null(strstr(deck, ".psp"));

  /* Shipped cold compose: keep method FUNCTIONAL PBE, append geometry atoms. */
  {
    double pos[6] = {0.0, 0.0, 0.0, 0.74, 0.0, 0.0};
    int z[2] = {1, 1};
    double cell[9] = {0};
    char cold[16384];
    int cold_len = 0;
    memset(cold, 0, sizeof(cold));
    assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 0, cold,
                                                   (int)sizeof(cold), &cold_len),
                     1);
    assert_true(cold_len > 0);
    assert_non_null(strstr(cold, "FUNCTIONAL"));
    assert_non_null(strstr(cold, "PBE"));
    /* Method merge keeps Cap'n FUNCTIONAL PBE (not rebuilt minimal BLYP).
     * PP filenames may still contain BLYP (H_CVB_BLYP.psp library names). */
    assert_null(strstr(cold, "FUNCTIONAL BLYP"));
    assert_null(strstr(cold, "FUNCTIONAL\n  BLYP"));
    assert_true(strstr(cold, "&ATOMS") != NULL || strstr(cold, "&atoms") != NULL);
    assert_non_null(strstr(cold, "H_CVB_BLYP.psp") || strstr(cold, "*H_"));
    assert_non_null(strstr(cold, "0.740000") || strstr(cold, "0.74"));
    /* SYSTEM cell/cutoff from Cap'n method sections survive geometry merge. */
    assert_true(strstr(cold, "&SYSTEM") != NULL || strstr(cold, "&system") != NULL);
    assert_non_null(strstr(cold, "CELL"));
    assert_non_null(strstr(cold, "12"));
    assert_non_null(strstr(cold, "CUTOFF"));
    assert_non_null(strstr(cold, "70"));
    assert_non_null(strstr(cold, "MAXSTEP") || strstr(cold, "MAX STEP") ||
                    strstr(cold, "50"));
  }
  free(msg);
}

int main(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr, "usage: %s top.bin sections.bin parser.bin atoms_extras.bin method_only.bin\n", argv[0]);
    return 2;
  }
  g_top = argv[1];
  g_sections = argv[2];
  g_parser = argv[3];
  g_atoms_extras = argv[4];
  g_method_only = argv[5];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_set_params_applies_top_level_via_fortran),
      cmocka_unit_test(test_set_params_applies_section_overrides_via_fortran),
      cmocka_unit_test(test_session_create_applies_parser_fixture),
      cmocka_unit_test(test_set_params_stores_typed_section_deck),
      cmocka_unit_test(test_set_params_method_only_keeps_dft_section),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
