/* Method-deck rebuild keeps non-coordinate &ATOMS directives. */
#include "cpmdc.h"
#include "cpmdc_atoms_extras.h"
#include "cpmdc_params.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len, char *cpmd_root,
                           int cpmd_root_len);
int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);

static const char *g_rebuild = NULL;
static const char *g_bad = NULL;

static const char k_deck[] =
    "&DFT\n"
    " CONSTRAINTS\n"
    "&END\n"
    "&ATOMS\n"
    "! CONSTRAINTS in a comment\n"
    "!SPECIES H\n"
    "*H_CVB_BLYP.psp\n"
    " LMAX=S\n"
    "   0\n"
    "  1.250000  2.250000  3.250000\n"
    " CONSTRAINTS\n"
    "  FIX\n"
    "  1 2\n"
    "  1 2 1.800000 GROWTH 0.050000\n"
    " END CONSTRAINTS\n"
    " ISOTOPE\n"
    "  18\n"
    " VELOCITIES\n"
    "  1\n"
    "  1 0.100000 0.200000 0.300000\n"
    " END VELOCITIES\n"
    " DUMMY ATOMS\n"
    "  2\n"
    "  TYPE4 1 2 0.500000 0.500000\n"
    " CHANGE BONDS\n"
    "  1\n"
    "  3 2 +1\n"
    " CONFINEMENT POTENTIAL\n"
    "  1.5 2.0\n"
    " ATOMIC CHARGES\n"
    "  0.0\n"
    " MOVIE TYPE\n"
    "  2\n"
    " GENERATE COORDINATES\n"
    "  1 1\n"
    " META DYNAMICS\n"
    "  COLVAR\n"
    " END METADYNAMICS\n"
    "&END\n";

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  long n;
  unsigned char *buf;
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  buf = (unsigned char *)malloc((size_t)n);
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

static int keep(const char *deck, const int *z, int n_atoms, char *out,
                int cap, int *out_len, char *err, int err_cap) {
  return cpmdc_atoms_extras_keep(deck, (int)strlen(deck), z, n_atoms, out, cap,
                                 out_len, err, err_cap);
}

static void test_dropped_directive_constant(void **state) {
  (void)state;
  assert_string_equal(cpmdc_atoms_dropped_directives[0], "VELOCITIES");
  assert_string_equal(cpmdc_atoms_dropped_directives[1], "METADYNAMICS");
  assert_null(cpmdc_atoms_dropped_directives[2]);
}

static void test_keep_remaps_and_drops(void **state) {
  int z[] = {1, 8, 1};
  char out[4096];
  char again[4096];
  char err[192];
  int out_len = -1;
  int again_len = -1;
  const char *constraints;
  const char *movie;
  int hits = 0;
  const char *p;
  (void)state;

  assert_int_equal(keep(k_deck, z, 3, out, (int)sizeof(out), &out_len, err,
                        (int)sizeof(err)),
                   0);
  assert_true(out_len > 0);
  assert_int_equal(keep(k_deck, z, 3, again, (int)sizeof(again), &again_len,
                        err, (int)sizeof(err)),
                   0);
  assert_int_equal(out_len, again_len);
  assert_memory_equal(out, again, (size_t)out_len);

  assert_non_null(strstr(out, "CONSTRAINTS"));
  assert_non_null(strstr(out, "  1 3\n"));
  assert_non_null(strstr(out, "  1 3 1.800000 GROWTH 0.050000\n"));
  assert_null(strstr(out, "  1 2\n"));
  assert_non_null(strstr(out, "END CONSTRAINTS"));
  assert_non_null(strstr(out, "ISOTOPE"));
  assert_non_null(strstr(out, "  18\n"));
  assert_non_null(strstr(out, "DUMMY ATOMS"));
  assert_non_null(strstr(out, "  2\n"));
  assert_non_null(strstr(out, "TYPE4 1 3 0.500000 0.500000"));
  assert_non_null(strstr(out, "CHANGE BONDS"));
  assert_non_null(strstr(out, "  2 3 +1\n"));
  assert_non_null(strstr(out, "CONFINEMENT POTENTIAL"));
  assert_non_null(strstr(out, "  1.5 2.0\n"));
  assert_non_null(strstr(out, "ATOMIC CHARGES"));
  assert_non_null(strstr(out, "  0.0\n"));
  assert_non_null(strstr(out, "MOVIE TYPE"));
  movie = strstr(out, "MOVIE TYPE");
  assert_non_null(movie);
  assert_non_null(strstr(movie, "  2\n"));

  assert_null(strstr(out, "VELOCITIES"));
  assert_null(strstr(out, "METADYNAMICS"));
  assert_null(strstr(out, "META DYNAMICS"));
  assert_null(strstr(out, "COLVAR"));
  assert_null(strstr(out, "GENERATE"));
  assert_null(strstr(out, "0.100000"));
  assert_null(strstr(out, "1.250000"));
  assert_null(strstr(out, "*H_CVB"));
  assert_null(strstr(out, "LMAX"));

  constraints = strstr(out, "CONSTRAINTS");
  assert_true(constraints < strstr(out, "ISOTOPE"));
  assert_true(strstr(out, "ISOTOPE") < strstr(out, "DUMMY ATOMS"));
  assert_true(strstr(out, "DUMMY ATOMS") < strstr(out, "CHANGE BONDS"));
  assert_true(strstr(out, "CHANGE BONDS") < strstr(out, "CONFINEMENT POTENTIAL"));
  assert_true(strstr(out, "CONFINEMENT POTENTIAL") < strstr(out, "ATOMIC CHARGES"));
  assert_true(strstr(out, "ATOMIC CHARGES") < movie);

  for (p = out; (p = strstr(p, "CONSTRAINTS")) != NULL; p++)
    hits++;
  assert_int_equal(hits, 2);
}

static void test_identity_order_and_rejection(void **state) {
  int grouped[] = {8, 1};
  int z[] = {1, 8, 1};
  char out[1024];
  char err[192];
  int out_len = 0;
  const char *same =
      "&ATOMS\n"
      " CONSTRAINTS\n"
      "  1 2\n"
      " END CONSTRAINTS\n"
      "&END\n";
  const char *bad =
      "&ATOMS\n"
      " CONSTRAINTS\n"
      "  9 1\n"
      " END CONSTRAINTS\n"
      "&END\n";
  const char *lone =
      "&atoms\n"
      " constraints\n"
      "  9\n"
      " end constraints\n"
      "&end\n";
  char tiny[4];
  (void)state;

  assert_int_equal(keep(same, grouped, 2, out, (int)sizeof(out), &out_len, err,
                        (int)sizeof(err)),
                   0);
  assert_non_null(strstr(out, "  1 2\n"));

  err[0] = '\0';
  assert_int_equal(keep(bad, z, 3, out, (int)sizeof(out), &out_len, err,
                        (int)sizeof(err)),
                   -1);
  assert_non_null(strstr(err, "CONSTRAINTS"));

  assert_int_equal(keep(lone, z, 3, out, (int)sizeof(out), &out_len, err,
                        (int)sizeof(err)),
                   0);
  assert_non_null(strstr(out, "  9\n"));
  assert_non_null(strstr(out, "constraints"));

  assert_int_equal(keep(same, grouped, 2, tiny, (int)sizeof(tiny), &out_len,
                        err, (int)sizeof(err)),
                   -1);
  assert_non_null(strstr(err, "do not fit"));
}

static void test_compose_keeps_rendered_directives(void **state) {
  size_t n = 0;
  unsigned char *msg;
  char functional[64];
  char applied[CPMDC_BLOCKS];
  char root[64];
  char cold[CPMDC_BLOCKS];
  char cold2[CPMDC_BLOCKS];
  double cutoff = 0.0;
  int charge = 0;
  int mult = 0;
  int cold_len = 0;
  int cold_len2 = 0;
  double pos[9] = {0.0, 0.0, 0.0, 0.75, 0.58, 0.12, 0.0, 0.6, 0.0};
  int z[3] = {1, 8, 1};
  double cell[9] = {0};
  (void)state;

  assert_non_null(g_rebuild);
  msg = read_file(g_rebuild, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);
  memset(applied, 0, sizeof(applied));
  assert_int_equal(cpmdc_embed_get_config(functional, (int)sizeof(functional),
                                          &cutoff, &charge, &mult, applied,
                                          (int)sizeof(applied), root,
                                          (int)sizeof(root)),
                   0);
  assert_non_null(strstr(applied, "CONSTRAINTS"));
  assert_non_null(strstr(applied, "1 2"));
  assert_non_null(strstr(applied, "ISOTOPE"));
  assert_non_null(strstr(applied, "DUMMY ATOMS"));
  assert_non_null(strstr(applied, "CHANGE BONDS"));
  assert_non_null(strstr(applied, "GENERATE COORDINATES"));
  assert_non_null(strstr(applied, "CONFINEMENT POTENTIAL"));
  assert_non_null(strstr(applied, "ATOMIC CHARGES"));
  assert_non_null(strstr(applied, "MOVIE TYPE"));
  assert_non_null(strstr(applied, "VELOCITIES"));
  assert_non_null(strstr(applied, "0.100000"));
  assert_non_null(strstr(applied, "METADYNAMICS"));
  assert_non_null(strstr(applied, "COLVAR"));

  memset(cold, 0, sizeof(cold));
  memset(cold2, 0, sizeof(cold2));
  assert_int_equal(cpmdc_embed_compose_cold_deck(3, pos, z, cell, 0, cold,
                                                 (int)sizeof(cold), &cold_len),
                   1);
  assert_int_equal(cpmdc_embed_compose_cold_deck(
                       3, pos, z, cell, 0, cold2, (int)sizeof(cold2), &cold_len2),
                   1);
  assert_int_equal(cold_len, cold_len2);
  assert_true(cold_len > 0);
  assert_memory_equal(cold, cold2, (size_t)cold_len);

  assert_non_null(strstr(cold, "CONSTRAINTS"));
  assert_non_null(strstr(cold, "1 3"));
  assert_null(strstr(cold, "1 2"));
  assert_non_null(strstr(cold, "ISOTOPE"));
  assert_non_null(strstr(cold, "18.000000"));
  assert_non_null(strstr(cold, "DUMMY ATOMS"));
  assert_non_null(strstr(cold, "TYPE4 1 1.000000"));
  assert_non_null(strstr(cold, "CHANGE BONDS"));
  assert_non_null(strstr(cold, "1 +1"));
  assert_non_null(strstr(cold, "CONFINEMENT POTENTIAL"));
  assert_non_null(strstr(cold, "1.5 2.0"));
  assert_non_null(strstr(cold, "ATOMIC CHARGES"));
  assert_non_null(strstr(cold, "MOVIE TYPE"));
  assert_non_null(strstr(cold, "*H_CVB_BLYP.psp"));
  assert_non_null(strstr(cold, "*O_MT_BLYP.psp"));
  assert_non_null(strstr(cold, "0.580000"));
  assert_non_null(strstr(cold, "0.600000"));
  assert_non_null(strstr(cold, "0.750000"));
  assert_non_null(strstr(cold, "0.120000"));
  assert_null(strstr(cold, "VELOCITIES"));
  assert_null(strstr(cold, "0.100000"));
  assert_null(strstr(cold, "METADYNAMICS"));
  assert_null(strstr(cold, "COLVAR"));
  assert_null(strstr(cold, "GENERATE"));
  free(msg);
}

static void test_compose_rejects_index(void **state) {
  size_t n = 0;
  unsigned char *msg;
  double pos[6] = {0.0, 0.0, 0.0, 0.1, 0.0, 0.0};
  int z[2] = {1, 8};
  double cell[9] = {0};
  char cold[CPMDC_BLOCKS];
  int cold_len = 0;
  (void)state;

  assert_non_null(g_bad);
  msg = read_file(g_bad, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);
  memset(cold, 0, sizeof(cold));
  assert_int_equal(cpmdc_embed_compose_cold_deck(2, pos, z, cell, 0, cold,
                                                 (int)sizeof(cold), &cold_len),
                   0);
  free(msg);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_dropped_directive_constant),
      cmocka_unit_test(test_keep_remaps_and_drops),
      cmocka_unit_test(test_identity_order_and_rejection),
      cmocka_unit_test(test_compose_keeps_rendered_directives),
      cmocka_unit_test(test_compose_rejects_index),
  };
  if (argc >= 2)
    g_rebuild = argv[1];
  if (argc >= 3)
    g_bad = argv[2];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
