/**
 * Drive the shipped set_params / session_create path and observe embed knobs
 * via cpmdc_embed_get_config (capnp-fortran apply-from-bytes).
 */
#define _POSIX_C_SOURCE 200809L

#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

/* 1 only when this binary is the OpenCPMD link. The reference build is 0. */
#if defined(CPMDC_HAS_CPMD)
#define CPMDC_AVAILABLE_WHEN_LINKED 1
#else
#define CPMDC_AVAILABLE_WHEN_LINKED 0
#endif

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
static const char *g_atoms_message = NULL;

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
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
  assert_string_equal(root, "cpmd-root");
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
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
  /* CPMD reads RAGGIO and NLCC from the *file line, not from the LMAX line. */
  assert_non_null(strstr(deck, "*O_MT_BLYP.psp KLEINMAN-BYLANDER RAGGIO=1.200000 NLCC\n"
                               " LMAX=D LOC=S SKIP=P\n"));
  /* 3ba9 long-tail tokens from shipped cpmdc_params.c render. */
  assert_non_null(strstr(deck, "EMPIRICAL CORRECTION"));
  assert_non_null(strstr(deck, "GRIMME"));
  assert_non_null(strstr(deck, "DIPOLE MOMENT"));
  assert_non_null(strstr(deck, "LOCALIZE"));
  assert_non_null(strstr(deck, "HTHRS"));
  assert_true(strstr(deck, "TROTTER DIMENSION") != NULL ||
              strstr(deck, "REPLICA NUMBER") != NULL);
  assert_non_null(strstr(deck, "TAMM-DANCOFF"));
  /* 3ba9 RESP/EXTE/VECTORS (atoms_extras fixture, shipped render). */
  assert_non_null(strstr(deck, "&RESP") || strstr(deck, "&resp"));
  assert_non_null(strstr(deck, "HYPERBOLIC"));
  assert_non_null(strstr(deck, "WEIGHT"));
  assert_non_null(strstr(deck, "&EXTE") || strstr(deck, "&exte"));
  assert_non_null(strstr(deck, "EFIELD"));
  assert_non_null(strstr(deck, "&VECTORS") || strstr(deck, "&vectors"));
  assert_non_null(strstr(deck, "NEWORTHO"));
  assert_non_null(strstr(deck, "OVERLAP"));

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
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
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
  assert_int_equal(charge, 2);
  assert_int_equal(mult, 3);
  assert_true(strstr(deck, "&DFT") != NULL || strstr(deck, "&dft") != NULL);
  assert_non_null(strstr(deck, "FUNCTIONAL"));
  assert_non_null(strstr(deck, "PBE"));
  /* C render may emit empty &ATOMS/&END; no real PP (*) until geometry merge. */
  assert_null(strstr(deck, "*H_"));
  assert_null(strstr(deck, "*.psp"));
  assert_null(strstr(deck, ".psp"));
  /* CPMD method knobs in applied deck before compose. */
  assert_non_null(strstr(deck, "CONVERGENCE ORBITALS"));
  assert_non_null(strstr(deck, "MAXSTEP"));
  assert_non_null(strstr(deck, "CHARGE"));

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
    /* CPMD maxStep + convergence + SYSTEM CHARGE survive merge (wdwj). */
    assert_non_null(strstr(cold, "MAXSTEP"));
    assert_non_null(strstr(cold, "50"));
    assert_non_null(strstr(cold, "CONVERGENCE ORBITALS"));
    assert_non_null(strstr(cold, "CHARGE"));
    assert_non_null(strstr(cold, "2"));
    /* CPMD OPTIMIZE WAVEFUNCTION + DFT LSD survive method merge (wdwj). */
    assert_non_null(strstr(cold, "OPTIMIZE WAVEFUNCTION"));
    assert_non_null(strstr(cold, "LSD"));
  }
  free(msg);
}

static int line_equals(const char *deck, const char *exact) {
  size_t n;
  const char *p;
  if (!deck || !exact)
    return 0;
  n = strlen(exact);
  p = deck;
  while (*p) {
    if (strncmp(p, exact, n) == 0 && (p[n] == '\n' || p[n] == '\0'))
      return 1;
    p = strchr(p, '\n');
    if (!p)
      break;
    p++;
  }
  return 0;
}

#define require_line(deck, exact)                                              \
  do {                                                                         \
    if (!line_equals((deck), (exact)))                                         \
      fail_msg("missing line [%s] in deck:\n%s", (exact), (deck));             \
  } while (0)

/* Compose while Fortran's error unit is a temporary file. */
static int compose_capture_stderr(int n_atoms, const double *pos, const int *z,
                                  const double *cell, char *deck, int cap,
                                  int *dlen, char *err, size_t errcap) {
  char path[] = "cpmdc-pp-errXXXXXX";
  int fd;
  int saved;
  int rc;
  FILE *fp;
  size_t nread;
  fd = mkstemp(path);
  if (fd < 0)
    return -1;
  unlink(path);
  fflush(stderr);
  saved = dup(fileno(stderr));
  if (saved < 0) {
    close(fd);
    return -1;
  }
  if (dup2(fd, fileno(stderr)) < 0) {
    close(saved);
    close(fd);
    return -1;
  }
  rc = cpmdc_embed_compose_cold_deck(n_atoms, pos, z, cell, 0, deck, cap, dlen);
  fflush(stderr);
  if (dup2(saved, fileno(stderr)) < 0) {
    close(saved);
    close(fd);
    return -1;
  }
  close(saved);
  if (lseek(fd, 0, SEEK_SET) < 0) {
    close(fd);
    return -1;
  }
  fp = fdopen(fd, "r");
  if (!fp) {
    close(fd);
    return -1;
  }
  nread = fread(err, 1, errcap - 1, fp);
  err[nread] = '\0';
  fclose(fp);
  return rc;
}

/* Listed elements keep file, LMAX, LOC, and the options of the *file line:
 * per-species KLEINMAN-BYLANDER, RAGGIO, and UPF for a .UPF path. Unlisted
 * table elements fall back without inheriting them. */
static void test_compose_uses_message_pseudopotentials(void **state) {
  (void)state;
  double pos[15] = {0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2, 4.0, 4.1, 4.2};
  int z[5] = {14, 7, 1, 6, 8};
  int missing[1] = {26};
  double missing_pos[3] = {0.0, 0.0, 0.0};
  double cell[9] = {0};
  char cold[16384];
  char err[1024];
  int cold_len = 0;
  int rc;
  size_t n = 0;
  unsigned char *msg;
  char functional[64], deck[CPMDC_BLOCKS], root[1024];
  double cutoff = 0.0;
  int charge = 0, mult = 0;
  const char *si;
  const char *nitrogen;
  const char *hydrogen;
  const char *carbon;
  const char *oxygen;

  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);
  assert_non_null(g_atoms_message);
  msg = read_file(g_atoms_message, &n);
  assert_non_null(msg);
  assert_int_equal(cpmdc_set_params(msg, n), 0);
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root, sizeof(root));
  require_line(deck, "!SPECIES Si");
  require_line(deck, "!SPECIES N");
  require_line(deck, "*custom_silicon.psp RAGGIO=1.500000");
  require_line(deck, "*custom_nitrogen.psp KLEINMAN-BYLANDER");
  require_line(deck, "*O.pbe-mt_fhi.UPF UPF");
  require_line(deck, " LMAX=D LOC=S");
  require_line(deck, " LMAX=P LOC=P SKIP=S");
  assert_null(strstr(deck, "H_CVB_BLYP.psp"));
  assert_non_null(strstr(deck, "KLEINMAN-BYLANDER"));

  memset(cold, 0, sizeof(cold));
  assert_int_equal(cpmdc_embed_compose_cold_deck(5, pos, z, cell, 0, cold,
                                                 (int)sizeof(cold), &cold_len),
                   1);
  assert_true(cold_len > 0);
  require_line(cold, "*custom_silicon.psp RAGGIO=1.500000");
  require_line(cold, " LMAX=D LOC=S");
  require_line(cold, "*custom_nitrogen.psp KLEINMAN-BYLANDER");
  require_line(cold, " LMAX=P LOC=P SKIP=S");
  require_line(cold, "*O.pbe-mt_fhi.UPF UPF");
  require_line(cold, "*H_CVB_BLYP.psp");
  require_line(cold, " LMAX=S");
  require_line(cold, "*C_MT_BLYP.psp");
  require_line(cold, " LMAX=P");
  assert_null(strstr(cold, "*H_CVB_BLYP.psp KLEINMAN-BYLANDER"));
  assert_null(strstr(cold, "*C_MT_BLYP.psp KLEINMAN-BYLANDER"));
  assert_null(strstr(cold, "Si_MT_BLYP.psp"));
  assert_null(strstr(cold, "N_MT_BLYP.psp"));
  assert_non_null(strstr(cold, "NEWCODE"));
  assert_null(strstr(cold, "OLDCODE"));
  assert_non_null(strstr(cold, "0.100000"));
  si = strstr(cold, "*custom_silicon.psp");
  nitrogen = strstr(cold, "*custom_nitrogen.psp");
  hydrogen = strstr(cold, "*H_CVB_BLYP.psp");
  carbon = strstr(cold, "*C_MT_BLYP.psp");
  oxygen = strstr(cold, "*O.pbe-mt_fhi.UPF");
  assert_non_null(si);
  assert_non_null(nitrogen);
  assert_non_null(hydrogen);
  assert_non_null(carbon);
  assert_non_null(oxygen);
  assert_true(si < nitrogen && nitrogen < hydrogen && hydrogen < carbon && carbon < oxygen);

  memset(cold, 0, sizeof(cold));
  memset(err, 0, sizeof(err));
  cold_len = 0;
  rc = compose_capture_stderr(1, missing_pos, missing, cell, cold,
                              (int)sizeof(cold), &cold_len, err, sizeof(err));
  assert_int_equal(rc, 0);
  if (!strstr(err, "atomic number 26"))
    fail_msg("missing atomic number in error [%s]", err);
  free(msg);
}

static void
test_set_params_system_has_charge_zero_overrides_top_level(void **state) {
  (void)state;
  assert_int_equal(cpmdc_available(), CPMDC_AVAILABLE_WHEN_LINKED);

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  assert_int_not_equal(root.type, CAPN_NULL);

  CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
  struct CPMDSystemSection sys;
  read_CPMDSystemSection(&sys, sys_ptr);
  sys.symmetry = 0;
  sys.angstrom = 1;
  sys.charge = 0;
  sys.hasCharge = 1;
  write_CPMDSystemSection(&sys, sys_ptr);

  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = sys_ptr;
  set_CPMDInputSection(&sec, secs, 0);

  CPMDParams_ptr params_ptr = new_CPMDParams(root.seg);
  struct CPMDParams params;
  memset(&params, 0, sizeof(params));
  params.functional = (capn_text){3, "PBE", NULL};
  params.cutOffRy = 70.0;
  params.charge = 2;
  params.multiplicity = 1;
  params.task = (capn_text){8, "gradient", NULL};
  params.inputSections = secs;
  write_CPMDParams(&params, params_ptr);
  assert_int_equal(capn_setp(root, 0, params_ptr.p), 0);

  uint8_t buf[2048];
  int written = capn_write_mem(&arena, buf, sizeof(buf), 0);
  capn_free(&arena);
  assert_true(written > 0);

  assert_int_equal(cpmdc_set_params(buf, (size_t)written), 0);

  char functional[64], deck[CPMDC_BLOCKS], root_dir[1024];
  double cutoff = 0.0;
  int charge = 99, mult = 0;
  read_applied(functional, sizeof(functional), &cutoff, &charge, &mult, deck,
               sizeof(deck), root_dir, sizeof(root_dir));
  assert_string_equal(functional, "PBE");
  assert_int_equal(charge, 0);
  assert_null(strstr(deck, "CHARGE"));
}

int main(int argc, char **argv) {
  if (argc != 7) {
    fprintf(stderr,
            "usage: %s top.bin sections.bin parser.bin atoms_extras.bin "
            "method_only.bin atoms_message.bin\n",
            argv[0]);
    return 2;
  }
  g_top = argv[1];
  g_sections = argv[2];
  g_parser = argv[3];
  g_atoms_extras = argv[4];
  g_method_only = argv[5];
  g_atoms_message = argv[6];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_set_params_applies_top_level_via_fortran),
      cmocka_unit_test(test_set_params_applies_section_overrides_via_fortran),
      cmocka_unit_test(test_session_create_applies_parser_fixture),
      cmocka_unit_test(test_set_params_stores_typed_section_deck),
      cmocka_unit_test(test_set_params_method_only_keeps_dft_section),
      cmocka_unit_test(test_compose_uses_message_pseudopotentials),
      cmocka_unit_test(
          test_set_params_system_has_charge_zero_overrides_top_level),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
