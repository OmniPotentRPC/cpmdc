/**
 * PotentialConfig common-overlay lowering: session creation, evaluation, and
 * loud rejection of unsound combinations.
 */
#define _GNU_SOURCE
#include "cpmdc.h"
#include "Potentials.capnp.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

static const char *g_overlay = NULL;
static const char *g_reject = NULL;
static const char *g_step = NULL;
static const char *g_bad_params = NULL;
static const char *g_arm_overlay = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0 || n > (16 << 20)) {
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

static void test_overlay_session_evaluates(void **state) {
  (void)state;
  size_t config_size = 0, step_size = 0;
  unsigned char *config = read_file(g_overlay, &config_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(config);
  assert_non_null(step);

  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(session);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  CPMDCResult r =
      cpmdc_session_calculate_result(session, step, step_size, out, need,
                                     &wrote);
  assert_int_equal(r.ok, 1);
  assert_true(wrote > 0);

  size_t wrote_one_shot = 0;
  CPMDCResult one_shot = cpmdc_calculate_result_from_config(
      config, config_size, step, step_size, out, need, &wrote_one_shot);
  assert_int_equal(one_shot.ok, 1);
  assert_int_equal(wrote_one_shot, wrote);

  cpmdc_session_destroy(session);
  free(out);
  free(step);
  free(config);
}

static void test_overlay_rejects_unlowerable_field(void **state) {
  (void)state;
  size_t config_size = 0;
  unsigned char *config = read_file(g_reject, &config_size);
  assert_non_null(config);
  assert_null(cpmdc_session_create_from_config(config, config_size));
  assert_non_null(strstr(cpmdc_last_error(), "relativityMethod"));
  assert_int_not_equal(cpmdc_configure(config, config_size), 0);
  free(config);
}

static void test_overlay_kmesh_reaches_deck(void **state) {
  (void)state;
  size_t config_size = 0;
  unsigned char *config = read_file(g_overlay, &config_size);
  assert_non_null(config);
  char deck_path[] = "overlay_kmesh_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  unsetenv("CPMDC_DECK_OUT");
  assert_non_null(session);
  size_t deck_size = 0;
  unsigned char *deck = read_file(deck_path, &deck_size);
  remove(deck_path);
  assert_non_null(deck);
  unsigned char *text = (unsigned char *)realloc(deck, deck_size + 1);
  assert_non_null(text);
  text[deck_size] = '\0';
  /* potential_config_common_cpmd.capnp.txt sets kMesh = [2, 2, 2]. */
  assert_non_null(strstr((const char *)text, "KPOINTS MONKHORST-PACK\n  2 2 2"));
  free(text);
  cpmdc_session_destroy(session);
  free(config);
}

/* Renders the deck of a constant-cutoff call on the session into a text
 * file, which shows the overlay the session currently renders with. */
static char *deck_after_constant_cutoff(CPMDCSession *session) {
  char deck_path[] = "overlay_atomic_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  int rc = cpmdc_session_set_constant_cutoff(session, 10.0, 0.0, 0.0);
  unsetenv("CPMDC_DECK_OUT");
  assert_int_equal(rc, 0);
  size_t size = 0;
  unsigned char *deck = read_file(deck_path, &size);
  remove(deck_path);
  assert_non_null(deck);
  char *text = (char *)realloc(deck, size + 1);
  assert_non_null(text);
  text[size] = '\0';
  return text;
}

static void test_failed_set_params_keeps_overlay(void **state) {
  (void)state;
  size_t config_size = 0, bad_size = 0;
  unsigned char *config = read_file(g_arm_overlay, &config_size);
  unsigned char *bad = read_file(g_bad_params, &bad_size);
  assert_non_null(config);
  assert_non_null(bad);

  CPMDCSession *control =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(control);
  char *expected = deck_after_constant_cutoff(control);
  assert_non_null(strstr(expected, "FUNCTIONAL PBE"));
  assert_non_null(strstr(expected, "CHARGE\n  -1\n"));
  assert_non_null(strstr(expected, "MULTIPLICITY\n  3\n"));

  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const unsigned char garbage[] = {0xde, 0xad, 0xbe, 0xef, 0x01, 0x02,
                                   0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
  assert_int_equal(cpmdc_session_set_params(session, garbage, sizeof(garbage)),
                   -1);
  assert_string_equal(cpmdc_last_error(), "CPMDParams parse failed");

  assert_int_equal(cpmdc_session_set_params(session, bad, bad_size), -1);
  assert_non_null(strstr(cpmdc_last_error(), "deck render failed"));

  char *actual = deck_after_constant_cutoff(session);
  assert_string_equal(actual, expected);

  free(actual);
  free(expected);
  cpmdc_session_destroy(session);
  cpmdc_session_destroy(control);
  free(bad);
  free(config);
}

static capn_text capn_str(const char *s) {
  capn_text t = {s ? (int)strlen(s) : 0, s ? s : "", NULL};
  return t;
}

static unsigned char *serialize_config(struct capn *arena,
                                       PotentialConfig_ptr cfg, size_t *size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, cfg.p), 0);
  size_t cap = 8192u;
  unsigned char *buf = (unsigned char *)malloc(cap);
  assert_non_null(buf);
  int n = capn_write_mem(arena, buf, cap, 0);
  assert_true(n > 0);
  *size = (size_t)n;
  return buf;
}

static char *render_deck_from_config(const unsigned char *config,
                                     size_t config_size) {
  char deck_path[] = "overlay_merge_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  unsetenv("CPMDC_DECK_OUT");
  assert_non_null(session);
  cpmdc_session_destroy(session);
  size_t deck_size = 0;
  unsigned char *deck = read_file(deck_path, &deck_size);
  remove(deck_path);
  assert_non_null(deck);
  char *text = (char *)realloc(deck, deck_size + 1);
  assert_non_null(text);
  text[deck_size] = '\0';
  return text;
}

static PotentialConfig_ptr make_cpmd_plus_common(
    struct capn *arena, struct CPMDParams *arm_view,
    struct CommonMethodSpec *common_view, int with_arm) {
  capn_ptr root = capn_root(arena);
  PotentialConfig_ptr cfg_ptr = new_PotentialConfig(root.seg);
  struct PotentialConfig cfg;
  memset(&cfg, 0, sizeof(cfg));
  if (with_arm) {
    CPMDParams_ptr arm_ptr = new_CPMDParams(root.seg);
    write_CPMDParams(arm_view, arm_ptr);
    cfg.which = PotentialConfig_cpmd;
    cfg.cpmd = arm_ptr;
  } else {
    cfg.which = PotentialConfig_none;
  }
  CommonMethodSpec_ptr com_ptr = new_CommonMethodSpec(root.seg);
  write_CommonMethodSpec(common_view, com_ptr);
  cfg.common = com_ptr;
  write_PotentialConfig(&cfg, cfg_ptr);
  return cfg_ptr;
}

static void init_default_arm(struct CPMDParams *arm) {
  memset(arm, 0, sizeof(*arm));
  arm->functional = capn_str("BLYP");
  arm->cutOffRy = 70.0;
  arm->charge = 0;
  arm->multiplicity = 1;
  arm->task = capn_str("gradient");
}

static void expect_config_rejected(const unsigned char *config,
                                   size_t config_size,
                                   const char *expected_substr) {
  assert_int_equal(cpmdc_configure(config, config_size), -1);
  assert_non_null(strstr(cpmdc_last_error(), expected_substr));
  assert_null(cpmdc_session_create_from_config(config, config_size));
  assert_non_null(strstr(cpmdc_last_error(), expected_substr));
}

static void test_merge_xc_functionals(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    capn_ptr xc = capn_new_ptr_list(root.seg, 1);
    capn_set_text(xc, 0, capn_str("PBE"));
    com.xcFunctionals = xc;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, "FUNCTIONAL PBE\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.functional = capn_str("LDA");
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    capn_ptr xc = capn_new_ptr_list(root.seg, 1);
    capn_set_text(xc, 0, capn_str("PBE"));
    com.xcFunctionals = xc;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, "FUNCTIONAL LDA\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    capn_ptr xc = capn_new_ptr_list(root.seg, 1);
    capn_set_text(xc, 0, capn_str("BOGUS"));
    com.xcFunctionals = xc;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "xcFunctionals");
    free(buf);
  }
}

static void test_merge_basis_set(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.basisSet = capn_str("PlaneWave");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    assert_int_equal(cpmdc_configure(buf, sz), 0);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.basisSet = capn_str("6-31G*");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "basisSet");
    free(buf);
  }
}

static void test_merge_cutoff_charge_multiplicity(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = 1088.45544983952;
    com.charge = -1;
    com.spinMultiplicity = 3;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, "CUTOFF\n  80\n"));
    assert_non_null(strstr(deck, "CHARGE\n  -1\n"));
    assert_non_null(strstr(deck, "MULTIPLICITY\n  3\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.cutOffRy = 50.0;
    arm.charge = 2;
    arm.multiplicity = 2;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = 1088.45544983952;
    com.charge = -1;
    com.spinMultiplicity = 3;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, "CUTOFF\n  50\n"));
    assert_non_null(strstr(deck, "CHARGE\n  2\n"));
    assert_non_null(strstr(deck, "MULTIPLICITY\n  2\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = -5.0;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "planewaveCutoffEv");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = -1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
}

static void test_merge_rejects_kmesh_scf_vdw_method_relativity(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.kMesh = capn_new_list32(root.seg, 3);
    capn_set32(com.kMesh, 0, 2);
    capn_set32(com.kMesh, 1, 2);
    capn_set32(com.kMesh, 2, 2);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "kMesh");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfMaxIterations = 50;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "SCF");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfEnergyToleranceEv = 1e-4;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "SCF");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsMethod = capn_str("DFT-D3");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaals");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.relativityMethod = capn_str("ZORA");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "relativity");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfMaxIterations = -5;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "scfMaxIterations");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfEnergyToleranceEv = -1e-4;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "scfEnergyToleranceEv");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfEnergyToleranceEv = NAN;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "scfEnergyToleranceEv");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.scfEnergyToleranceEv = INFINITY;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "scfEnergyToleranceEv");
    free(buf);
  }
}

static void test_merge_smearing(void **state) {
  (void)state;
  {
    /* Default arm + fermi smearing applies FREE-ENERGY + TEMPERATURE ELECTRON. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.025851999786;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, " FREE-ENERGY\n"));
    assert_non_null(strstr(deck, " TEMPERATURE ELECTRON\n  300\n"));
    free(deck);
    free(buf);
  }
  {
    /* Arm has &CPMD section without freeEnergy/temperatureElectron -> merged. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDCpmdSection cpmd_sec;
    memset(&cpmd_sec, 0, sizeof(cpmd_sec));
    cpmd_sec.optimizeWavefunction = 1;
    cpmd_sec.convergenceOrbitals = 1.0e-6;
    cpmd_sec.maxIter = 50;
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.025851999786;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, " MAXITER\n  50\n"));
    assert_non_null(strstr(deck, " FREE-ENERGY\n"));
    assert_non_null(strstr(deck, " TEMPERATURE ELECTRON\n  300\n"));
    free(deck);
    free(buf);
  }
  {
    /* Arm already sets freeEnergy/temperatureElectron -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDCpmdSection cpmd_sec;
    memset(&cpmd_sec, 0, sizeof(cpmd_sec));
    cpmd_sec.optimizeWavefunction = 1;
    cpmd_sec.convergenceOrbitals = 1.0e-6;
    cpmd_sec.freeEnergy = 1;
    cpmd_sec.temperatureElectron = 500.0;
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.0258519976;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    /* Gaussian smearing has no CPMD lowering -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_gaussian;
    sm.widthEv = 0.05;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    /* Fermi smearing with widthEv <= 0 -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.0;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "widthEv");
    free(buf);
  }
  {
    /* kind = none with nonzero widthEv -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_none;
    sm.widthEv = 0.05;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
}

static void test_merge_vdw_s6(void **state) {
  (void)state;
  {
    /* Arm has vdwParams with s6 == 0 (default) -> overlay s6 copied. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDVdwSection vdw;
    memset(&vdw, 0, sizeof(vdw));
    vdw.empiricalCorrection = 1;
    vdw.grimme = capn_str("D3");
    CPMDVdwSection_ptr vdw_ptr = new_CPMDVdwSection(root.seg);
    write_CPMDVdwSection(&vdw, vdw_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_vdwParams;
    sec.vdwParams = vdw_ptr;
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, " GRIMME D3\n"));
    assert_non_null(strstr(deck, " S6\n  1.050000\n"));
    free(deck);
    free(buf);
  }
  {
    /* Arm already sets vdwParams.s6 -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDVdwSection vdw;
    memset(&vdw, 0, sizeof(vdw));
    vdw.empiricalCorrection = 1;
    vdw.grimme = capn_str("D3");
    vdw.s6 = 1.0;
    CPMDVdwSection_ptr vdw_ptr = new_CPMDVdwSection(root.seg);
    write_CPMDVdwSection(&vdw, vdw_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_vdwParams;
    sec.vdwParams = vdw_ptr;
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
  {
    /* Arm has no &VDW section -> vanDerWaalsS6 rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
  {
    /* Overlay-only with vanDerWaalsS6 and no vanDerWaalsMethod -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, NULL, &com, 0), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
  {
    /* Arm has raw &VDW block -> vanDerWaalsS6 rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_raw;
    sec.raw = capn_str("&VDW\n GRIMME D3\n&END");
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
  {
    /* Arm has directive &VDW section -> vanDerWaalsS6 rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDDirectiveSection vdw_dir;
    memset(&vdw_dir, 0, sizeof(vdw_dir));
    CPMDDirectiveSection_ptr vdw_dir_ptr = new_CPMDDirectiveSection(root.seg);
    write_CPMDDirectiveSection(&vdw_dir, vdw_dir_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_vdw;
    sec.vdw = vdw_dir_ptr;
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = 1.05;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
  {
    /* Negative vanDerWaalsS6 -> rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.vanDerWaalsS6 = -0.5;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "vanDerWaalsS6");
    free(buf);
  }
}

static void test_merge_raw_cpmd_and_session_configure(void **state) {
  (void)state;
  {
    /* Arm has raw &CPMD block -> smearing rejected. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_raw;
    sec.raw = capn_str("&CPMD\n OPTIMIZE WAVEFUNCTION\n&END");
    set_CPMDInputSection(&sec, secs, 0);
    arm.inputSections = secs;
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.025851999786;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    /* cpmdc_session_configure and cpmdc_calculate_result_from_config apply merge. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.basisSet = capn_str("planewave");
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.025851999786;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
    capn_free(&arena);

    size_t step_size = 0;
    unsigned char *step = read_file(g_step, &step_size);
    assert_non_null(step);
    size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
    assert_true(need > 0);
    unsigned char *out = (unsigned char *)malloc(need);
    assert_non_null(out);
    size_t wrote = 0;
    CPMDCResult r = cpmdc_calculate_result_from_config(buf, sz, step, step_size,
                                                       out, need, &wrote);
    assert_int_equal(r.ok, 1);
    assert_true(wrote > 0);

    size_t overlay_size = 0;
    unsigned char *overlay = read_file(g_overlay, &overlay_size);
    assert_non_null(overlay);
    CPMDCSession *session =
        cpmdc_session_create_from_config(overlay, overlay_size);
    assert_non_null(session);
    assert_int_equal(cpmdc_session_configure(session, buf, sz), 0);
    cpmdc_session_destroy(session);
    free(overlay);
    free(out);
    free(step);
    free(buf);
  }
}

int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len,
                           char *cpmd_root, int cpmd_root_len);

static void set_xc(struct capn *arena, struct CommonMethodSpec *com,
                   const char *name) {
  capn_ptr root = capn_root(arena);
  capn_ptr xc = capn_new_ptr_list(root.seg, 1);
  capn_set_text(xc, 0, capn_str(name));
  com->xcFunctionals = xc;
}

static void attach_section(struct capn *arena, struct CPMDParams *arm,
                           struct CPMDInputSection *sec) {
  capn_ptr root = capn_root(arena);
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
  set_CPMDInputSection(sec, secs, 0);
  arm->inputSections = secs;
}

static char *deck_keeping_session_all(const unsigned char *config,
                                      size_t config_size, char **functional_out,
                                      double *cutoff_out, int *charge_out,
                                      int *mult_out) {
  char deck_path[] = "overlay_merge_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  unsetenv("CPMDC_DECK_OUT");
  assert_non_null(session);
  char functional[64];
  char root[64];
  char *embed_deck = (char *)malloc(1 << 16);
  double cutoff = 0.0;
  int charge = 0;
  int mult = 0;
  assert_non_null(embed_deck);
  memset(functional, 0, sizeof(functional));
  memset(root, 0, sizeof(root));
  assert_int_equal(cpmdc_embed_get_config(functional, (int)sizeof(functional),
                                          &cutoff, &charge, &mult, embed_deck,
                                          1 << 16, root, (int)sizeof(root)),
                   0);
  free(embed_deck);
  if (functional_out) {
    *functional_out = (char *)malloc(64);
    assert_non_null(*functional_out);
    memcpy(*functional_out, functional, 64);
  }
  if (cutoff_out)
    *cutoff_out = cutoff;
  if (charge_out)
    *charge_out = charge;
  if (mult_out)
    *mult_out = mult;
  cpmdc_session_destroy(session);
  size_t deck_size = 0;
  unsigned char *deck = read_file(deck_path, &deck_size);
  remove(deck_path);
  assert_non_null(deck);
  char *text = (char *)realloc(deck, deck_size + 1);
  assert_non_null(text);
  text[deck_size] = '\0';
  return text;
}

static char *deck_keeping_session(const unsigned char *config,
                                  size_t config_size, int *charge_out,
                                  char **functional_out) {
  return deck_keeping_session_all(config, config_size, functional_out, NULL,
                                  charge_out, NULL);
}

static PotentialConfig_ptr make_cpmd_only(struct capn *arena,
                                          struct CPMDParams *arm_view) {
  capn_ptr root = capn_root(arena);
  PotentialConfig_ptr cfg_ptr = new_PotentialConfig(root.seg);
  struct PotentialConfig cfg;
  memset(&cfg, 0, sizeof(cfg));
  CPMDParams_ptr arm_ptr = new_CPMDParams(root.seg);
  write_CPMDParams(arm_view, arm_ptr);
  cfg.which = PotentialConfig_cpmd;
  cfg.cpmd = arm_ptr;
  write_PotentialConfig(&cfg, cfg_ptr);
  return cfg_ptr;
}

static CPMDSystemSection_ptr make_wire_default_system(struct capn *arena) {
  capn_ptr root = capn_root(arena);
  CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
  struct CPMDSystemSection sys;
  read_CPMDSystemSection(&sys, sys_ptr);
  write_CPMDSystemSection(&sys, sys_ptr);
  return sys_ptr;
}

static CPMDDftSection_ptr make_wire_default_dft(struct capn *arena) {
  capn_ptr root = capn_root(arena);
  CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
  struct CPMDDftSection dft;
  read_CPMDDftSection(&dft, dft_ptr);
  dft.oldCode = 1;
  write_CPMDDftSection(&dft, dft_ptr);
  return dft_ptr;
}

static void attach_two_sections(struct capn *arena, struct CPMDParams *arm,
                                struct CPMDInputSection *s0,
                                struct CPMDInputSection *s1) {
  capn_ptr root = capn_root(arena);
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 2);
  set_CPMDInputSection(s0, secs, 0);
  set_CPMDInputSection(s1, secs, 1);
  arm->inputSections = secs;
}

static CPMDDirective_list make_single_directive(struct capn *arena,
                                                const char *keyword,
                                                const char *const *args,
                                                int nargs) {
  capn_ptr root = capn_root(arena);
  struct CPMDDirective d;
  memset(&d, 0, sizeof(d));
  d.keyword = capn_str(keyword);
  if (nargs > 0) {
    capn_ptr list = capn_new_ptr_list(root.seg, nargs);
    for (int i = 0; i < nargs; ++i)
      capn_set_text(list, i, capn_str(args[i]));
    d.args = list;
  }
  CPMDDirective_list dirs = new_CPMDDirective_list(root.seg, 1);
  set_CPMDDirective(&d, dirs, 0);
  return dirs;
}

static void set_fermi_smearing(struct capn *arena,
                               struct CommonMethodSpec *com) {
  capn_ptr root = capn_root(arena);
  CommonMethodSpec_Smearing_ptr sm_ptr =
      new_CommonMethodSpec_Smearing(root.seg);
  struct CommonMethodSpec_Smearing sm;
  memset(&sm, 0, sizeof(sm));
  sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
  sm.widthEv = 0.025851999786;
  write_CommonMethodSpec_Smearing(&sm, sm_ptr);
  com->smearing = sm_ptr;
}

/* An unmapped token is rejected before arm-wins, including when the arm
 * functional is already not the schema default. */
static void test_merge_rejects_unmapped_xc_when_arm_wins(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  arm.functional = capn_str("LDA");
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "BOGUS");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

/* system=() leaves charge at the schema default 0. That is unset, so the
 * root charge stays (both with common unset and with common set), and an
 * explicit non-zero system charge wins over both. */
static void test_merge_unset_system_charge_keeps_root(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.charge = 2;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.charge = -1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int charge = 0;
    char *deck = deck_keeping_session(buf, sz, &charge, NULL);
    assert_int_equal(charge, 2);
    assert_non_null(strstr(deck, "CHARGE\n  2\n"));
    assert_null(strstr(deck, "CHARGE\n  -1\n"));
    free(deck);
    free(buf);
  }
  {
    /* Same with common omitted: root charge=2 is kept. */
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.charge = 2;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_only(&arena, &arm), &sz);
    capn_free(&arena);
    int charge = 0;
    char *deck = deck_keeping_session(buf, sz, &charge, NULL);
    assert_int_equal(charge, 2);
    assert_non_null(strstr(deck, "CHARGE\n  2\n"));
    free(deck);
    free(buf);
  }
  {
    /* An explicit non-zero system charge still wins over the root and the
     * overlay. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.charge = 2;
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    sys.charge = 4;
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.charge = -1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int charge = 0;
    char *deck = deck_keeping_session(buf, sz, &charge, NULL);
    assert_int_equal(charge, 4);
    assert_non_null(strstr(deck, "CHARGE\n  4\n"));
    assert_null(strstr(deck, "CHARGE\n  -1\n"));
    free(deck);
    free(buf);
  }
}

/* These messages deliberately leave the typed field absent on the wire.
 * Populating a zeroed generated C view is not equivalent for fields whose
 * schema default is nonzero. */
static void test_merge_omitted_dft_functional_uses_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_dft;
  sec.dft = dft_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  char *deck = render_deck_from_config(buf, sz);
  assert_non_null(strstr(deck, "FUNCTIONAL PBE\n"));
  assert_null(strstr(deck, "FUNCTIONAL BLYP\n"));
  free(deck);
  free(buf);
}

static void test_merge_omitted_system_cutoff_uses_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = sys_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.planewaveCutoffEv = 1088.45544983952;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  char *deck = render_deck_from_config(buf, sz);
  assert_non_null(strstr(deck, "CUTOFF\n  80\n"));
  assert_null(strstr(deck, "CUTOFF\n  70\n"));
  free(deck);
  free(buf);
}

static void test_merge_omitted_system_multiplicity_keeps_arm(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  arm.multiplicity = 2;
  CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = sys_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.spinMultiplicity = 3;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  char *deck = render_deck_from_config(buf, sz);
  assert_non_null(strstr(deck, "MULTIPLICITY\n  2\n"));
  assert_null(strstr(deck, "MULTIPLICITY\n  3\n"));
  free(deck);
  free(buf);
}

static void test_merge_rejects_typed_dft_directive_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDDirective_list directives = new_CPMDDirective_list(root.seg, 1);
  struct CPMDDirective directive;
  memset(&directive, 0, sizeof(directive));
  directive.keyword = capn_str("FUNCTIONAL");
  set_CPMDDirective(&directive, directives, 0);
  struct CPMDDftSection dft;
  memset(&dft, 0, sizeof(dft));
  dft.directives = directives;
  CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
  write_CPMDDftSection(&dft, dft_ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_dft;
  sec.dft = dft_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

static void test_merge_rejects_raw_dft_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, 1);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_raw;
  sec.raw = capn_str("&DFT\n FUNCTIONAL LDA\n&END");
  set_CPMDInputSection(&sec, secs, 0);
  arm.inputSections = secs;
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

static void test_merge_rejects_inputblocks_dft_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  capn_ptr blocks = capn_new_ptr_list(root.seg, 1);
  capn_set_text(blocks, 0, capn_str("&DFT\n FUNCTIONAL LDA\n&END"));
  arm.inputBlocks = blocks;
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

static void test_merge_rejects_generic_dft_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDGenericSection gen;
  memset(&gen, 0, sizeof(gen));
  gen.name = capn_str("DFT");
  CPMDGenericSection_ptr gen_ptr = new_CPMDGenericSection(root.seg);
  write_CPMDGenericSection(&gen, gen_ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_generic;
  sec.generic = gen_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

static void test_merge_rejects_set_functional_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDSetDirective setd;
  memset(&setd, 0, sizeof(setd));
  setd.key = capn_str("DFT.FUNCTIONAL");
  setd.value = capn_str("LDA");
  CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
  write_CPMDSetDirective(&setd, set_ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_set;
  sec.set = set_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "xcFunctionals");
  free(buf);
}

static void test_merge_rejects_raw_system_charge_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_raw;
  sec.raw = capn_str("&SYSTEM\n SYMMETRY\n  0\n&END");
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.charge = -1;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "charge");
  free(buf);
}

static void test_merge_rejects_set_charge_overlay(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDSetDirective setd;
  memset(&setd, 0, sizeof(setd));
  setd.key = capn_str("SYSTEM.CHARGE");
  setd.value = capn_str("3");
  CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
  write_CPMDSetDirective(&setd, set_ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_set;
  sec.set = set_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.charge = -1;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  expect_config_rejected(buf, sz, "charge");
  free(buf);
}

/* A set directive that does not name the functional still accepts the
 * overlay, and the rendered deck carries that functional. */
static void test_merge_set_other_dft_keeps_overlay_functional(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDSetDirective setd;
  memset(&setd, 0, sizeof(setd));
  setd.key = capn_str("DFT.OLDCODE");
  CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
  write_CPMDSetDirective(&setd, set_ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_set;
  sec.set = set_ptr;
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "PBE");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  char *functional = NULL;
  char *deck = deck_keeping_session(buf, sz, NULL, &functional);
  assert_string_equal(functional, "PBE");
  assert_non_null(strstr(deck, "FUNCTIONAL PBE\n"));
  assert_non_null(strstr(deck, "OLDCODE\n"));
  free(functional);
  free(deck);
  free(buf);
}

/* Wire-default typed system cutOffRy=70.0 takes common planewaveCutoffEv. */
static void test_merge_wire_default_system_cutoff_takes_common(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = make_wire_default_system(&arena);
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.planewaveCutoffEv = 1224.51238106946; /* 90 Ry */
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  double cutoff = 0.0;
  char *deck = deck_keeping_session_all(buf, sz, NULL, &cutoff, NULL, NULL);
  assert_float_equal(cutoff, 90.0, 1e-6);
  assert_non_null(strstr(deck, "CUTOFF\n  90\n"));
  assert_null(strstr(deck, "CUTOFF\n  70\n"));
  free(deck);
  free(buf);
}

/* Wire-default typed system multiplicity=1 takes common spinMultiplicity. */
static void
test_merge_wire_default_system_multiplicity_takes_common(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = make_wire_default_system(&arena);
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.spinMultiplicity = 3;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  int mult = 0;
  char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, NULL, &mult);
  assert_int_equal(mult, 3);
  assert_non_null(strstr(deck, "MULTIPLICITY\n  3\n"));
  assert_non_null(strstr(deck, " LSD\n"));
  free(deck);
  free(buf);
}

/* Wire-default typed system charge=0 takes common charge. */
static void test_merge_wire_default_system_charge_takes_common(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = make_wire_default_system(&arena);
  attach_section(&arena, &arm, &sec);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.charge = -1;
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  int charge = 0;
  char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, &charge, NULL);
  assert_int_equal(charge, -1);
  assert_non_null(strstr(deck, "CHARGE\n  -1\n"));
  free(deck);
  free(buf);
}

/* Wire-default typed dft functional="BLYP" takes common hybrid xcFunctionals
 * in embed, rendered deck, and deck_needs_cluster_hfx (TUCKERMAN solver). */
static void test_merge_wire_default_dft_functional_takes_common(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection s0, s1;
  memset(&s0, 0, sizeof(s0));
  s0.which = CPMDInputSection_system;
  s0.system = make_wire_default_system(&arena);
  memset(&s1, 0, sizeof(s1));
  s1.which = CPMDInputSection_dft;
  s1.dft = make_wire_default_dft(&arena);
  attach_two_sections(&arena, &arm, &s0, &s1);
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  set_xc(&arena, &com, "HYB_GGA_XC_PBEH");
  size_t sz = 0;
  unsigned char *buf =
      serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                       &sz);
  capn_free(&arena);
  char *functional = NULL;
  char *deck = deck_keeping_session_all(buf, sz, &functional, NULL, NULL, NULL);
  assert_string_equal(functional, "PBE0");
  assert_non_null(strstr(deck, "FUNCTIONAL PBE0\n"));
  assert_null(strstr(deck, "FUNCTIONAL BLYP\n"));
  assert_non_null(strstr(deck, "POISSON SOLVER TUCKERMAN\n"));
  assert_null(strstr(deck, "POISSON SOLVER HOCKNEY\n"));
  free(functional);
  free(deck);
  free(buf);
}

/* Wire-default system cutOffRy=70 keeps explicit root cutOffRy (with or
 * without common); an explicit non-default system.cutOffRy wins over both. */
static void test_merge_unset_system_cutoff_keeps_root(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.cutOffRy = 90.0;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = 1088.45544983952; /* 80 Ry */
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    double cutoff = 0.0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, &cutoff, NULL, NULL);
    assert_float_equal(cutoff, 90.0, 1e-6);
    assert_non_null(strstr(deck, "CUTOFF\n  90\n"));
    assert_null(strstr(deck, "CUTOFF\n  70\n"));
    assert_null(strstr(deck, "CUTOFF\n  80\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.cutOffRy = 90.0;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_only(&arena, &arm), &sz);
    capn_free(&arena);
    double cutoff = 0.0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, &cutoff, NULL, NULL);
    assert_float_equal(cutoff, 90.0, 1e-6);
    assert_non_null(strstr(deck, "CUTOFF\n  90\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.cutOffRy = 90.0;
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    sys.cutOffRy = 100.0;
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = 1088.45544983952; /* 80 Ry */
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    double cutoff = 0.0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, &cutoff, NULL, NULL);
    assert_float_equal(cutoff, 100.0, 1e-6);
    assert_non_null(strstr(deck, "CUTOFF\n  100\n"));
    free(deck);
    free(buf);
  }
}

/* Wire-default system multiplicity=1 keeps explicit root multiplicity (with or
 * without common); an explicit system.multiplicity > 1 wins over both. */
static void test_merge_unset_system_multiplicity_keeps_root(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.multiplicity = 3;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 2;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int mult = 0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, NULL, &mult);
    assert_int_equal(mult, 3);
    assert_non_null(strstr(deck, "MULTIPLICITY\n  3\n"));
    assert_null(strstr(deck, "MULTIPLICITY\n  2\n"));
    assert_non_null(strstr(deck, " LSD\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.multiplicity = 3;
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_only(&arena, &arm), &sz);
    capn_free(&arena);
    int mult = 0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, NULL, &mult);
    assert_int_equal(mult, 3);
    assert_non_null(strstr(deck, "MULTIPLICITY\n  3\n"));
    assert_non_null(strstr(deck, " LSD\n"));
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.multiplicity = 3;
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    sys.multiplicity = 5;
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 2;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int mult = 0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, NULL, &mult);
    assert_int_equal(mult, 5);
    assert_non_null(strstr(deck, "MULTIPLICITY\n  5\n"));
    assert_non_null(strstr(deck, " LSD\n"));
    free(deck);
    free(buf);
  }
}

/* Wire-default dft functional="BLYP" keeps explicit root functional (with or
 * without common, including hybrid HFX solver selection); an explicit
 * non-default dft.functional wins over both. */
static void test_merge_unset_dft_functional_keeps_root(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.functional = capn_str("LDA");
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_dft;
    sec.dft = make_wire_default_dft(&arena);
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_xc(&arena, &com, "PBE");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    char *functional = NULL;
    char *deck =
        deck_keeping_session_all(buf, sz, &functional, NULL, NULL, NULL);
    assert_string_equal(functional, "LDA");
    assert_non_null(strstr(deck, "FUNCTIONAL LDA\n"));
    assert_null(strstr(deck, "FUNCTIONAL BLYP\n"));
    assert_null(strstr(deck, "FUNCTIONAL PBE\n"));
    free(functional);
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.functional = capn_str("PBE0");
    struct CPMDInputSection s0, s1;
    memset(&s0, 0, sizeof(s0));
    s0.which = CPMDInputSection_system;
    s0.system = make_wire_default_system(&arena);
    memset(&s1, 0, sizeof(s1));
    s1.which = CPMDInputSection_dft;
    s1.dft = make_wire_default_dft(&arena);
    attach_two_sections(&arena, &arm, &s0, &s1);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_only(&arena, &arm), &sz);
    capn_free(&arena);
    char *functional = NULL;
    char *deck =
        deck_keeping_session_all(buf, sz, &functional, NULL, NULL, NULL);
    assert_string_equal(functional, "PBE0");
    assert_non_null(strstr(deck, "FUNCTIONAL PBE0\n"));
    assert_non_null(strstr(deck, "POISSON SOLVER TUCKERMAN\n"));
    free(functional);
    free(deck);
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    arm.functional = capn_str("LDA");
    CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
    struct CPMDDftSection dft;
    read_CPMDDftSection(&dft, dft_ptr);
    dft.functional = capn_str("PBE");
    write_CPMDDftSection(&dft, dft_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_dft;
    sec.dft = dft_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_xc(&arena, &com, "PBE0");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    char *functional = NULL;
    char *deck =
        deck_keeping_session_all(buf, sz, &functional, NULL, NULL, NULL);
    assert_string_equal(functional, "PBE");
    assert_non_null(strstr(deck, "FUNCTIONAL PBE\n"));
    free(functional);
    free(deck);
    free(buf);
  }
}

/* Typed dft.directives carrying FUNCTIONAL conflicts with common xcFunctionals;
 * an unrelated dft.directives entry still accepts the overlay. */
static void test_merge_rejects_dft_directive_functional_overlay(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
    struct CPMDDftSection dft;
    read_CPMDDftSection(&dft, dft_ptr);
    dft.directives = make_single_directive(&arena, "FUNCTIONAL LDA", NULL, 0);
    write_CPMDDftSection(&dft, dft_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_dft;
    sec.dft = dft_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_xc(&arena, &com, "PBE");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "xcFunctionals");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
    struct CPMDDftSection dft;
    read_CPMDDftSection(&dft, dft_ptr);
    const char *args[] = {"LDA"};
    dft.directives = make_single_directive(&arena, "FUNCTIONAL", args, 1);
    write_CPMDDftSection(&dft, dft_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_dft;
    sec.dft = dft_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_xc(&arena, &com, "PBE");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "xcFunctionals");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDDftSection_ptr dft_ptr = new_CPMDDftSection(root.seg);
    struct CPMDDftSection dft;
    read_CPMDDftSection(&dft, dft_ptr);
    dft.directives = make_single_directive(&arena, "OLDCODE", NULL, 0);
    write_CPMDDftSection(&dft, dft_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_dft;
    sec.dft = dft_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_xc(&arena, &com, "PBE");
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    char *functional = NULL;
    char *deck =
        deck_keeping_session_all(buf, sz, &functional, NULL, NULL, NULL);
    assert_string_equal(functional, "PBE");
    assert_non_null(strstr(deck, "FUNCTIONAL PBE\n"));
    assert_non_null(strstr(deck, "OLDCODE\n"));
    free(functional);
    free(deck);
    free(buf);
  }
}

/* Typed system.directives carrying CUTOFF, CHARGE, or MULTIPLICITY conflicts
 * with the matching common overlay scalar; unrelated directives still merge. */
static void test_merge_rejects_system_directive_scalar_overlays(void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    const char *args[] = {"90"};
    sys.directives = make_single_directive(&arena, "CUTOFF", args, 1);
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.planewaveCutoffEv = 1088.45544983952;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "planewaveCutoffEv");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    const char *args[] = {"3"};
    sys.directives = make_single_directive(&arena, "CHARGE", args, 1);
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.charge = -1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "charge");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    const char *args[] = {"3"};
    sys.directives = make_single_directive(&arena, "MULTIPLICITY", args, 1);
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 2;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    const char *args[] = {"HOCKNEY"};
    sys.directives = make_single_directive(&arena, "POISSON SOLVER", args, 1);
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.charge = -1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int charge = 0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, &charge, NULL);
    assert_int_equal(charge, -1);
    assert_non_null(strstr(deck, "CHARGE\n  -1\n"));
    free(deck);
    free(buf);
  }
}

/* Explicit singlet spinMultiplicity=1 (where 0 means unset) still conflicts
 * with opaque &SYSTEM carriers (raw, set, system.directives). */
static void test_merge_rejects_explicit_singlet_against_opaque_system(
    void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_raw;
    sec.raw = capn_str("&SYSTEM\n SYMMETRY\n  0\n MULTIPLICITY\n  3\n&END");
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = capn_str("SYSTEM.MULTIPLICITY");
    setd.value = capn_str("3");
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
    struct CPMDSystemSection sys;
    read_CPMDSystemSection(&sys, sys_ptr);
    const char *args[] = {"3"};
    sys.directives = make_single_directive(&arena, "MULTIPLICITY", args, 1);
    write_CPMDSystemSection(&sys, sys_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = sys_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
  {
    /* Clean typed system without MULTIPLICITY directive accepts singlet. */
    struct capn arena;
    capn_init_malloc(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_system;
    sec.system = make_wire_default_system(&arena);
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    com.spinMultiplicity = 1;
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    int mult = 0;
    char *deck = deck_keeping_session_all(buf, sz, NULL, NULL, NULL, &mult);
    assert_int_equal(mult, 1);
    assert_null(strstr(deck, "MULTIPLICITY\n"));
    free(deck);
    free(buf);
  }
}

/* Split-argument TEMPERATURE + ELECTRON in cpmd.directives and CPMD.TEMPERATURE
 * set directives conflict with common Fermi smearing; ionic TEMPERATURE still
 * accepts Fermi smearing. */
static void test_merge_rejects_split_temperature_electron_smearing(
    void **state) {
  (void)state;
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    struct CPMDCpmdSection cpmd_sec;
    read_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    const char *args[] = {"ELECTRON", "500"};
    cpmd_sec.directives = make_single_directive(&arena, "TEMPERATURE", args, 2);
    write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_fermi_smearing(&arena, &com);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    struct CPMDCpmdSection cpmd_sec;
    read_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    const char *args[] = {"ELECTRON 500"};
    cpmd_sec.directives = make_single_directive(&arena, "TEMPERATURE", args, 1);
    write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_fermi_smearing(&arena, &com);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = capn_str("CPMD.TEMPERATURE");
    setd.value = capn_str("ELECTRON 500");
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_fermi_smearing(&arena, &com);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    expect_config_rejected(buf, sz, "smearing");
    free(buf);
  }
  {
    /* Ionic TEMPERATURE directive does not conflict with Fermi smearing. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    struct CPMDCpmdSection cpmd_sec;
    read_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    const char *args[] = {"300"};
    cpmd_sec.directives = make_single_directive(&arena, "TEMPERATURE", args, 1);
    write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_fermi_smearing(&arena, &com);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, " FREE-ENERGY\n"));
    assert_non_null(strstr(deck, " TEMPERATURE ELECTRON\n  300\n"));
    free(deck);
    free(buf);
  }
  {
    /* Ionic CPMD.TEMPERATURE set directive does not conflict with smearing. */
    struct capn arena;
    capn_init_malloc(&arena);
    capn_ptr root = capn_root(&arena);
    struct CPMDParams arm;
    init_default_arm(&arm);
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = capn_str("CPMD.TEMPERATURE");
    setd.value = capn_str("300");
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    struct CPMDInputSection sec;
    memset(&sec, 0, sizeof(sec));
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    attach_section(&arena, &arm, &sec);
    struct CommonMethodSpec com;
    memset(&com, 0, sizeof(com));
    set_fermi_smearing(&arena, &com);
    size_t sz = 0;
    unsigned char *buf =
        serialize_config(&arena, make_cpmd_plus_common(&arena, &arm, &com, 1),
                         &sz);
    capn_free(&arena);
    char *deck = render_deck_from_config(buf, sz);
    assert_non_null(strstr(deck, " FREE-ENERGY\n"));
    assert_non_null(strstr(deck, " TEMPERATURE ELECTRON\n  300\n"));
    free(deck);
    free(buf);
  }
}

/* One config for the precedence matrix: top-level fields, typed &DFT and
 * &SYSTEM sections, and the common overlay. A section field of NULL, 0, or
 * the schema default is left unset. */
struct PrecedenceCase {
  const char *name;
  const char *root_functional;
  double root_cutoff_ry;
  int root_charge;
  int root_mult;
  int has_dft;
  const char *dft_functional;
  int has_sys;
  double sys_cutoff_ry;
  int sys_charge;
  int sys_mult;
  const char *xc;
  double overlay_cutoff_ry;
  int overlay_charge;
  int overlay_mult;
  const char *want_functional;
  double want_cutoff_ry;
  int want_charge;
  int want_mult;
};

#define RY_EV 13.605693122994

struct Recorded {
  char functional[64];
  double cutoff_ry;
  int charge;
  int mult;
};

static void read_recorded(struct Recorded *rec) {
  char root[64];
  char *deck = (char *)malloc(1 << 16);
  assert_non_null(deck);
  memset(rec, 0, sizeof(*rec));
  memset(root, 0, sizeof(root));
  assert_int_equal(
      cpmdc_embed_get_config(rec->functional, (int)sizeof(rec->functional),
                             &rec->cutoff_ry, &rec->charge, &rec->mult, deck,
                             1 << 16, root, (int)sizeof(root)),
      0);
  free(deck);
}

static int count_text(const char *hay, const char *needle) {
  int n = 0;
  size_t len = strlen(needle);
  for (const char *p = hay; (p = strstr(p, needle)) != NULL; p += len)
    n++;
  return n;
}

/* The scalars a deck states. A deck that omits CHARGE or MULTIPLICITY states
 * the CPMD defaults 0 and 1. Each keyword appears at most once so a
 * duplicate cannot hide behind the first. */
static void read_deck_scalars(const char *deck, struct Recorded *out) {
  memset(out, 0, sizeof(*out));
  assert_int_equal(count_text(deck, " FUNCTIONAL "), 1);
  assert_int_equal(count_text(deck, " CUTOFF\n  "), 1);
  assert_true(count_text(deck, " CHARGE\n  ") <= 1);
  assert_true(count_text(deck, " MULTIPLICITY\n  ") <= 1);
  const char *fn = strstr(deck, " FUNCTIONAL ") + strlen(" FUNCTIONAL ");
  size_t n = strcspn(fn, "\n");
  assert_true(n < sizeof(out->functional));
  memcpy(out->functional, fn, n);
  out->functional[n] = '\0';
  out->cutoff_ry = strtod(strstr(deck, " CUTOFF\n  ") + strlen(" CUTOFF\n  "),
                          NULL);
  const char *ch = strstr(deck, " CHARGE\n  ");
  out->charge = ch ? (int)strtol(ch + strlen(" CHARGE\n  "), NULL, 10) : 0;
  const char *mu = strstr(deck, " MULTIPLICITY\n  ");
  out->mult = mu ? (int)strtol(mu + strlen(" MULTIPLICITY\n  "), NULL, 10) : 1;
}

static CPMDInputSection_list new_sections(struct capn *arena, int n) {
  capn_ptr root = capn_root(arena);
  return new_CPMDInputSection_list(root.seg, n);
}

static void put_dft(struct capn *arena, CPMDInputSection_list secs, int idx,
                    const char *functional, CPMDDirective_list directives,
                    int has_directives) {
  capn_ptr root = capn_root(arena);
  struct CPMDDftSection dft;
  memset(&dft, 0, sizeof(dft));
  if (functional)
    dft.functional = capn_str(functional);
  if (has_directives)
    dft.directives = directives;
  CPMDDftSection_ptr ptr = new_CPMDDftSection(root.seg);
  write_CPMDDftSection(&dft, ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_dft;
  sec.dft = ptr;
  set_CPMDInputSection(&sec, secs, idx);
}

static void put_system(struct capn *arena, CPMDInputSection_list secs, int idx,
                       double cutoff_ry, int charge, int mult,
                       CPMDDirective_list directives, int has_directives) {
  capn_ptr root = capn_root(arena);
  struct CPMDSystemSection sys;
  memset(&sys, 0, sizeof(sys));
  sys.cutOffRy = cutoff_ry;
  sys.charge = charge;
  sys.multiplicity = mult;
  if (has_directives)
    sys.directives = directives;
  CPMDSystemSection_ptr ptr = new_CPMDSystemSection(root.seg);
  write_CPMDSystemSection(&sys, ptr);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  sec.which = CPMDInputSection_system;
  sec.system = ptr;
  set_CPMDInputSection(&sec, secs, idx);
}

static unsigned char *build_precedence_config(const struct PrecedenceCase *pc,
                                              size_t *size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  if (pc->root_functional)
    arm.functional = capn_str(pc->root_functional);
  if (pc->root_cutoff_ry > 0.0)
    arm.cutOffRy = pc->root_cutoff_ry;
  arm.charge = pc->root_charge;
  if (pc->root_mult > 0)
    arm.multiplicity = pc->root_mult;
  int nsec = (pc->has_dft ? 1 : 0) + (pc->has_sys ? 1 : 0);
  if (nsec > 0) {
    CPMDInputSection_list secs = new_sections(&arena, nsec);
    int idx = 0;
    CPMDDirective_list none;
    memset(&none, 0, sizeof(none));
    if (pc->has_dft)
      put_dft(&arena, secs, idx++, pc->dft_functional, none, 0);
    if (pc->has_sys)
      put_system(&arena, secs, idx++, pc->sys_cutoff_ry, pc->sys_charge,
                 pc->sys_mult, none, 0);
    arm.inputSections = secs;
  }
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  if (pc->xc)
    set_xc(&arena, &com, pc->xc);
  com.planewaveCutoffEv = pc->overlay_cutoff_ry * RY_EV;
  com.charge = pc->overlay_charge;
  com.spinMultiplicity = pc->overlay_mult;
  unsigned char *buf = serialize_config(
      &arena, make_cpmd_plus_common(&arena, &arm, &com, 1), size);
  capn_free(&arena);
  return buf;
}

/* The deck text, the configuration the embed records, and the expected
 * values all agree, through a session and through cpmdc_configure. */
static int check_precedence_case(const struct PrecedenceCase *pc) {
  int failures = 0;
  size_t sz = 0;
  unsigned char *buf = build_precedence_config(pc, &sz);
  struct Recorded want;
  memset(&want, 0, sizeof(want));
  snprintf(want.functional, sizeof(want.functional), "%s",
           pc->want_functional);
  want.cutoff_ry = pc->want_cutoff_ry;
  want.charge = pc->want_charge;
  want.mult = pc->want_mult;

  for (int via_session = 0; via_session < 2; ++via_session) {
    char deck_path[] = "overlay_precedence_deck_XXXXXX";
    int fd = mkstemp(deck_path);
    assert_true(fd >= 0);
    close(fd);
    setenv("CPMDC_DECK_OUT", deck_path, 1);
    CPMDCSession *session = NULL;
    if (via_session) {
      session = cpmdc_session_create_from_config(buf, sz);
      assert_non_null(session);
    } else {
      assert_int_equal(cpmdc_configure(buf, sz), 0);
    }
    unsetenv("CPMDC_DECK_OUT");
    struct Recorded rec;
    read_recorded(&rec);
    if (session)
      cpmdc_session_destroy(session);

    size_t deck_size = 0;
    unsigned char *raw = read_file(deck_path, &deck_size);
    remove(deck_path);
    assert_non_null(raw);
    char *deck = (char *)realloc(raw, deck_size + 1);
    assert_non_null(deck);
    deck[deck_size] = '\0';
    struct Recorded text;
    read_deck_scalars(deck, &text);

    if (strcmp(text.functional, want.functional) != 0 ||
        strcmp(rec.functional, want.functional) != 0 ||
        fabs(text.cutoff_ry - want.cutoff_ry) > 1e-9 ||
        fabs(rec.cutoff_ry - want.cutoff_ry) > 1e-9 ||
        text.charge != want.charge || rec.charge != want.charge ||
        text.mult != want.mult || rec.mult != want.mult) {
      print_error(
          "%s (%s): want %s %.10g %d %d, deck %s %.10g %d %d, recorded %s "
          "%.10g %d %d\n",
          pc->name, via_session ? "session" : "configure", want.functional,
          want.cutoff_ry, want.charge, want.mult, text.functional,
          text.cutoff_ry, text.charge, text.mult, rec.functional,
          rec.cutoff_ry, rec.charge, rec.mult);
      failures++;
    }
    free(deck);
  }
  free(buf);
  return failures;
}

/* Typed sections at their schema defaults do not hide the overlay. */
static void test_precedence_overlay_fills_typed_defaults(void **state) {
  (void)state;
  const struct PrecedenceCase cases[] = {
      {.name = "dft at default text, overlay xc",
       .has_dft = 1, .dft_functional = "BLYP", .xc = "PBE",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "dft unset text, overlay xc",
       .has_dft = 1, .xc = "PBE",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "system at default 70, overlay cutoff",
       .has_sys = 1, .sys_cutoff_ry = 70.0, .overlay_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "system at 0, overlay cutoff",
       .has_sys = 1, .overlay_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "dft and system at defaults, overlay xc and cutoff",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .xc = "PBE", .overlay_cutoff_ry = 80.0,
       .want_functional = "PBE", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "dft and system at defaults, overlay all four",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .sys_mult = 1, .xc = "PBE",
       .overlay_cutoff_ry = 50.0, .overlay_charge = -1, .overlay_mult = 3,
       .want_functional = "PBE", .want_cutoff_ry = 50.0, .want_charge = -1,
       .want_mult = 3},
      {.name = "system default charge and multiplicity, overlay both",
       .has_sys = 1, .overlay_charge = 2, .overlay_mult = 2,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_charge = 2,
       .want_mult = 2},
      {.name = "no overlay leaves the defaults",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 70.0,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_mult = 1},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_precedence_case(&cases[i]);
  assert_int_equal(failures, 0);
}

/* A value a typed field sets beats the overlay. */
static void test_precedence_typed_value_beats_overlay(void **state) {
  (void)state;
  const struct PrecedenceCase cases[] = {
      {.name = "typed dft functional beats overlay xc",
       .has_dft = 1, .dft_functional = "LDA", .xc = "PBE",
       .want_functional = "LDA", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "typed system cutoff beats overlay cutoff",
       .has_sys = 1, .sys_cutoff_ry = 50.0, .overlay_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "typed system charge and multiplicity beat the overlay",
       .has_sys = 1, .sys_charge = 4, .sys_mult = 2, .overlay_charge = -1,
       .overlay_mult = 3,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_charge = 4,
       .want_mult = 2},
      {.name = "typed dft, default system, overlay both",
       .has_dft = 1, .dft_functional = "LDA", .has_sys = 1,
       .sys_cutoff_ry = 70.0, .xc = "PBE", .overlay_cutoff_ry = 80.0,
       .want_functional = "LDA", .want_cutoff_ry = 80.0, .want_mult = 1},
      {.name = "default dft, typed system cutoff, overlay both",
       .has_dft = 1, .dft_functional = "BLYP", .has_sys = 1,
       .sys_cutoff_ry = 50.0, .xc = "PBE", .overlay_cutoff_ry = 80.0,
       .want_functional = "PBE", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "typed dft and system beat the overlay",
       .has_dft = 1, .dft_functional = "LDA", .has_sys = 1,
       .sys_cutoff_ry = 50.0, .sys_charge = 1, .sys_mult = 2, .xc = "PBE",
       .overlay_cutoff_ry = 80.0, .overlay_charge = -1, .overlay_mult = 3,
       .want_functional = "LDA", .want_cutoff_ry = 50.0, .want_charge = 1,
       .want_mult = 2},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_precedence_case(&cases[i]);
  assert_int_equal(failures, 0);
}

/* A typed section at its default does not erase an explicit top-level
 * field, and that field still beats the overlay. */
static void test_precedence_root_value_with_default_sections(void **state) {
  (void)state;
  const struct PrecedenceCase cases[] = {
      {.name = "root functional, default dft, no overlay",
       .root_functional = "PBE", .has_dft = 1, .dft_functional = "BLYP",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "root functional beats overlay xc",
       .root_functional = "PBE", .has_dft = 1, .dft_functional = "BLYP",
       .xc = "LDA",
       .want_functional = "PBE", .want_cutoff_ry = 70.0, .want_mult = 1},
      {.name = "root cutoff beats overlay cutoff",
       .root_cutoff_ry = 50.0, .has_sys = 1, .sys_cutoff_ry = 70.0,
       .overlay_cutoff_ry = 80.0,
       .want_functional = "BLYP", .want_cutoff_ry = 50.0, .want_mult = 1},
      {.name = "root charge and multiplicity, default system",
       .root_charge = 2, .root_mult = 3, .has_sys = 1, .overlay_charge = -1,
       .overlay_mult = 2,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_charge = 2,
       .want_mult = 3},
      {.name = "root multiplicity, default system, no overlay",
       .root_mult = 3, .has_sys = 1,
       .want_functional = "BLYP", .want_cutoff_ry = 70.0, .want_mult = 3},
      {.name = "typed section beats the root field",
       .root_functional = "PBE", .root_cutoff_ry = 50.0, .has_dft = 1,
       .dft_functional = "LDA", .has_sys = 1, .sys_cutoff_ry = 60.0,
       .want_functional = "LDA", .want_cutoff_ry = 60.0, .want_mult = 1},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    failures += check_precedence_case(&cases[i]);
  assert_int_equal(failures, 0);
}

/* The overlay cutoff is in eV and the deck cutoff is in Ry. */
static void test_precedence_overlay_cutoff_unit_conversion(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDInputSection_list secs = new_sections(&arena, 1);
  CPMDDirective_list none;
  memset(&none, 0, sizeof(none));
  put_system(&arena, secs, 0, 70.0, 0, 0, none, 0);
  arm.inputSections = secs;
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.planewaveCutoffEv = 1088.45544983952;
  size_t sz = 0;
  unsigned char *buf = serialize_config(
      &arena, make_cpmd_plus_common(&arena, &arm, &com, 1), &sz);
  capn_free(&arena);
  char *deck = render_deck_from_config(buf, sz);
  assert_non_null(strstr(deck, " CUTOFF\n  80\n"));
  assert_null(strstr(deck, " CUTOFF\n  70\n"));
  assert_null(strstr(deck, " CUTOFF\n  1088"));
  free(deck);
  free(buf);
}

enum SingletCarrier {
  SINGLET_NO_CARRIER,
  SINGLET_INPUT_BLOCK,
  SINGLET_RAW_SYSTEM,
  SINGLET_GENERIC_SYSTEM,
  SINGLET_SET_MULTIPLICITY,
  SINGLET_TYPED_SYSTEM_DIRECTIVE,
};

/* spinMultiplicity 0 is unset, so 1 is an explicit singlet and goes through
 * the same opaque carrier checks as any other multiplicity. */
static unsigned char *build_singlet_config(enum SingletCarrier carrier,
                                           size_t *size) {
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  switch (carrier) {
  case SINGLET_NO_CARRIER:
    break;
  case SINGLET_INPUT_BLOCK: {
    capn_ptr blocks = capn_new_ptr_list(root.seg, 1);
    capn_set_text(blocks, 0,
                  capn_str("&SYSTEM\n MULTIPLICITY\n  3\n&END"));
    arm.inputBlocks = blocks;
    break;
  }
  case SINGLET_RAW_SYSTEM:
    sec.which = CPMDInputSection_raw;
    sec.raw = capn_str("&SYSTEM\n MULTIPLICITY\n  3\n&END");
    attach_section(&arena, &arm, &sec);
    break;
  case SINGLET_GENERIC_SYSTEM: {
    struct CPMDGenericSection gen;
    memset(&gen, 0, sizeof(gen));
    gen.name = capn_str("SYSTEM");
    CPMDGenericSection_ptr gen_ptr = new_CPMDGenericSection(root.seg);
    write_CPMDGenericSection(&gen, gen_ptr);
    sec.which = CPMDInputSection_generic;
    sec.generic = gen_ptr;
    attach_section(&arena, &arm, &sec);
    break;
  }
  case SINGLET_SET_MULTIPLICITY: {
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = capn_str("SYSTEM.MULTIPLICITY");
    setd.value = capn_str("3");
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    attach_section(&arena, &arm, &sec);
    break;
  }
  case SINGLET_TYPED_SYSTEM_DIRECTIVE: {
    CPMDDirective_list dirs = new_CPMDDirective_list(root.seg, 1);
    struct CPMDDirective d;
    memset(&d, 0, sizeof(d));
    d.keyword = capn_str("MULTIPLICITY");
    capn_ptr args = capn_new_ptr_list(root.seg, 1);
    capn_set_text(args, 0, capn_str("3"));
    d.args = args;
    set_CPMDDirective(&d, dirs, 0);
    CPMDInputSection_list secs = new_sections(&arena, 1);
    put_system(&arena, secs, 0, 0.0, 0, 0, dirs, 1);
    arm.inputSections = secs;
    break;
  }
  }
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  com.spinMultiplicity = 1;
  unsigned char *buf = serialize_config(
      &arena, make_cpmd_plus_common(&arena, &arm, &com, 1), size);
  capn_free(&arena);
  return buf;
}

static void test_merge_rejects_explicit_singlet_next_to_carriers(void **state) {
  (void)state;
  const enum SingletCarrier carriers[] = {
      SINGLET_INPUT_BLOCK,
      SINGLET_RAW_SYSTEM,
      SINGLET_GENERIC_SYSTEM,
      SINGLET_SET_MULTIPLICITY,
      SINGLET_TYPED_SYSTEM_DIRECTIVE,
  };
  for (size_t i = 0; i < sizeof(carriers) / sizeof(carriers[0]); ++i) {
    size_t sz = 0;
    unsigned char *buf = build_singlet_config(carriers[i], &sz);
    expect_config_rejected(buf, sz, "spinMultiplicity");
    free(buf);
  }
}

static void test_merge_accepts_explicit_singlet(void **state) {
  (void)state;
  size_t sz = 0;
  unsigned char *buf = build_singlet_config(SINGLET_NO_CARRIER, &sz);
  assert_int_equal(cpmdc_configure(buf, sz), 0);
  struct Recorded rec;
  read_recorded(&rec);
  assert_int_equal(rec.mult, 1);
  char *deck = render_deck_from_config(buf, sz);
  struct Recorded from_deck;
  read_deck_scalars(deck, &from_deck);
  assert_int_equal(from_deck.mult, 1);
  free(deck);
  free(buf);
}

enum OverlayScalar {
  SCALAR_FUNCTIONAL,
  SCALAR_CUTOFF,
  SCALAR_CHARGE,
  SCALAR_MULTIPLICITY,
};

enum OpaqueCarrier {
  CARRIER_INPUT_BLOCK,
  CARRIER_RAW_SECTION,
  CARRIER_GENERIC_SECTION,
  CARRIER_SET_DIRECTIVE,
  CARRIER_TYPED_DIRECTIVE,
};

struct ScalarCarrier {
  const char *section;
  const char *keyword;
  const char *value;
  const char *raw;
  const char *field;
};

static const struct ScalarCarrier k_scalar_carriers[] = {
    [SCALAR_FUNCTIONAL] = {"DFT", "FUNCTIONAL", "PBE0",
                           "&DFT\n FUNCTIONAL PBE0\n&END", "xcFunctionals"},
    [SCALAR_CUTOFF] = {"SYSTEM", "CUTOFF", "90",
                       "&SYSTEM\n CUTOFF\n  90\n&END", "planewaveCutoffEv"},
    [SCALAR_CHARGE] = {"SYSTEM", "CHARGE", "3", "&SYSTEM\n CHARGE\n  3\n&END",
                       "charge"},
    [SCALAR_MULTIPLICITY] = {"SYSTEM", "MULTIPLICITY", "3",
                             "&SYSTEM\n MULTIPLICITY\n  3\n&END",
                             "spinMultiplicity"},
};

static const enum OpaqueCarrier k_all_carriers[] = {
    CARRIER_INPUT_BLOCK,     CARRIER_RAW_SECTION,   CARRIER_GENERIC_SECTION,
    CARRIER_SET_DIRECTIVE,   CARRIER_TYPED_DIRECTIVE,
};

/* A cpmd arm whose typed &DFT or &SYSTEM section sets the scalar to a value
 * that beats the overlay (functional LDA, cutoff 50 Ry, charge 2,
 * multiplicity 2) when typed is non-zero, an opaque carrier for the same
 * setting, and an overlay that sets the scalar. */
static unsigned char *build_typed_carrier_config(enum OverlayScalar scalar,
                                                 enum OpaqueCarrier carrier,
                                                 int typed, int overlay_mult,
                                                 size_t *size) {
  const struct ScalarCarrier *sc = &k_scalar_carriers[scalar];
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);

  CPMDDirective_list dirs;
  memset(&dirs, 0, sizeof(dirs));
  int has_dirs = carrier == CARRIER_TYPED_DIRECTIVE;
  if (has_dirs) {
    dirs = new_CPMDDirective_list(root.seg, 1);
    struct CPMDDirective d;
    memset(&d, 0, sizeof(d));
    d.keyword = capn_str(sc->keyword);
    capn_ptr args = capn_new_ptr_list(root.seg, 1);
    capn_set_text(args, 0, capn_str(sc->value));
    d.args = args;
    set_CPMDDirective(&d, dirs, 0);
  }
  int extra = carrier == CARRIER_RAW_SECTION ||
              carrier == CARRIER_GENERIC_SECTION ||
              carrier == CARRIER_SET_DIRECTIVE;
  CPMDInputSection_list secs = new_sections(&arena, 1 + extra);
  if (scalar == SCALAR_FUNCTIONAL)
    put_dft(&arena, secs, 0, typed ? "LDA" : NULL, dirs, has_dirs);
  else
    put_system(&arena, secs, 0,
               typed && scalar == SCALAR_CUTOFF ? 50.0 : 0.0,
               typed && scalar == SCALAR_CHARGE ? 2 : 0,
               typed && scalar == SCALAR_MULTIPLICITY ? 2 : 0, dirs,
               has_dirs);

  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  switch (carrier) {
  case CARRIER_INPUT_BLOCK: {
    capn_ptr blocks = capn_new_ptr_list(root.seg, 1);
    capn_set_text(blocks, 0, capn_str(sc->raw));
    arm.inputBlocks = blocks;
    break;
  }
  case CARRIER_RAW_SECTION:
    sec.which = CPMDInputSection_raw;
    sec.raw = capn_str(sc->raw);
    set_CPMDInputSection(&sec, secs, 1);
    break;
  case CARRIER_GENERIC_SECTION: {
    struct CPMDGenericSection gen;
    memset(&gen, 0, sizeof(gen));
    gen.name = capn_str(sc->section);
    CPMDGenericSection_ptr gen_ptr = new_CPMDGenericSection(root.seg);
    write_CPMDGenericSection(&gen, gen_ptr);
    sec.which = CPMDInputSection_generic;
    sec.generic = gen_ptr;
    set_CPMDInputSection(&sec, secs, 1);
    break;
  }
  case CARRIER_SET_DIRECTIVE: {
    char key[64];
    snprintf(key, sizeof(key), "%s.%s", sc->section, sc->keyword);
    struct CPMDSetDirective setd;
    memset(&setd, 0, sizeof(setd));
    setd.key = capn_str(key);
    setd.value = capn_str(sc->value);
    CPMDSetDirective_ptr set_ptr = new_CPMDSetDirective(root.seg);
    write_CPMDSetDirective(&setd, set_ptr);
    sec.which = CPMDInputSection_set;
    sec.set = set_ptr;
    set_CPMDInputSection(&sec, secs, 1);
    break;
  }
  case CARRIER_TYPED_DIRECTIVE:
    break;
  }
  arm.inputSections = secs;

  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  switch (scalar) {
  case SCALAR_FUNCTIONAL:
    set_xc(&arena, &com, "PBE");
    break;
  case SCALAR_CUTOFF:
    com.planewaveCutoffEv = 80.0 * RY_EV;
    break;
  case SCALAR_CHARGE:
    com.charge = -1;
    break;
  case SCALAR_MULTIPLICITY:
    com.spinMultiplicity = overlay_mult;
    break;
  }
  unsigned char *buf = serialize_config(
      &arena, make_cpmd_plus_common(&arena, &arm, &com, 1), size);
  capn_free(&arena);
  return buf;
}

/* Both entry points fail with a cpmdc_last_error() naming the field, and
 * neither writes a deck to CPMDC_DECK_OUT. */
static int expect_rejected_without_deck(const unsigned char *config,
                                        size_t config_size, const char *want,
                                        const char *label) {
  int failures = 0;
  char deck_path[] = "overlay_carrier_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  int rc = cpmdc_configure(config, config_size);
  const char *err = cpmdc_last_error();
  if (rc != -1 || !strstr(err, want)) {
    fprintf(stderr, "%s: cpmdc_configure rc %d, want \"%s\", got \"%s\"\n",
            label, rc, want, err);
    failures++;
  }
  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  err = cpmdc_last_error();
  if (session || !strstr(err, want)) {
    fprintf(stderr,
            "%s: cpmdc_session_create_from_config %p, want \"%s\", got "
            "\"%s\"\n",
            label, (void *)session, want, err);
    failures++;
  }
  if (session)
    cpmdc_session_destroy(session);
  unsetenv("CPMDC_DECK_OUT");
  size_t deck_size = 0;
  unsigned char *deck = read_file(deck_path, &deck_size);
  remove(deck_path);
  if (deck) {
    fprintf(stderr, "%s: a deck of %zu bytes was written\n", label,
            deck_size);
    free(deck);
    failures++;
  }
  return failures;
}

static const char *carrier_name(enum OpaqueCarrier carrier) {
  switch (carrier) {
  case CARRIER_INPUT_BLOCK:
    return "inputBlocks";
  case CARRIER_RAW_SECTION:
    return "raw section";
  case CARRIER_GENERIC_SECTION:
    return "generic section";
  case CARRIER_SET_DIRECTIVE:
    return "set directive";
  case CARRIER_TYPED_DIRECTIVE:
    return "typed-section directive";
  }
  return "?";
}

/* A typed value that wins over the overlay does not exempt the overlay
 * field from the carrier check, and neither does the absence of one. */
static void check_overlay_scalar_carriers(enum OverlayScalar scalar,
                                          int overlay_mult) {
  int failures = 0;
  for (int typed = 1; typed >= 0; --typed) {
    for (size_t i = 0; i < sizeof(k_all_carriers) / sizeof(k_all_carriers[0]);
         ++i) {
      char label[128];
      snprintf(label, sizeof(label), "%s, %s, %s",
               k_scalar_carriers[scalar].field,
               typed ? "typed value" : "no typed value",
               carrier_name(k_all_carriers[i]));
      size_t sz = 0;
      unsigned char *buf = build_typed_carrier_config(
          scalar, k_all_carriers[i], typed, overlay_mult, &sz);
      failures += expect_rejected_without_deck(
          buf, sz, k_scalar_carriers[scalar].field, label);
      free(buf);
    }
  }
  assert_int_equal(failures, 0);
}

static void test_merge_rejects_typed_functional_with_carrier(void **state) {
  (void)state;
  check_overlay_scalar_carriers(SCALAR_FUNCTIONAL, 0);
}

static void test_merge_rejects_typed_cutoff_with_carrier(void **state) {
  (void)state;
  check_overlay_scalar_carriers(SCALAR_CUTOFF, 0);
}

static void test_merge_rejects_typed_charge_with_carrier(void **state) {
  (void)state;
  check_overlay_scalar_carriers(SCALAR_CHARGE, 0);
}

static void test_merge_rejects_typed_multiplicity_with_carrier(void **state) {
  (void)state;
  check_overlay_scalar_carriers(SCALAR_MULTIPLICITY, 3);
}

static void test_merge_rejects_typed_multiplicity_singlet_with_carrier(
    void **state) {
  (void)state;
  check_overlay_scalar_carriers(SCALAR_MULTIPLICITY, 1);
}

enum ArmOverlayField {
  ARM_OV_XC,
  ARM_OV_CUTOFF,
  ARM_OV_CHARGE,
  ARM_OV_MULT,
  ARM_OV_SMEARING,
  ARM_OV_VDW_S6,
};

enum ArmBadList {
  ARM_BAD_INPUT_BLOCKS,
  ARM_BAD_INPUT_SECTIONS,
  ARM_BAD_DFT_DIRECTIVES,
  ARM_BAD_SYSTEM_DIRECTIVES,
  ARM_BAD_CPMD_DIRECTIVES,
  ARM_BAD_VDW_DIRECTIVES,
};

struct ArmFailureCase {
  enum ArmOverlayField field;
  enum ArmBadList bad;
  const char *want_field;
  const char *want_list;
};

/* A list of the wrong kind where the arm checks expect inputBlocks text, an
 * inputSections struct list, or a directives list, next to an overlay that
 * sets one field. */
static unsigned char *build_arm_failure_config(const struct ArmFailureCase *fc,
                                               size_t *size) {
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  struct CPMDParams arm;
  init_default_arm(&arm);
  CPMDDirective_list bad_dirs;
  bad_dirs.p = capn_new_ptr_list(root.seg, 1);
  struct CPMDInputSection sec;
  memset(&sec, 0, sizeof(sec));
  switch (fc->bad) {
  case ARM_BAD_INPUT_BLOCKS:
    arm.inputBlocks = capn_new_list32(root.seg, 1).p;
    break;
  case ARM_BAD_INPUT_SECTIONS:
    arm.inputSections.p = capn_new_ptr_list(root.seg, 1);
    break;
  case ARM_BAD_DFT_DIRECTIVES: {
    CPMDInputSection_list secs = new_sections(&arena, 1);
    put_dft(&arena, secs, 0, NULL, bad_dirs, 1);
    arm.inputSections = secs;
    break;
  }
  case ARM_BAD_SYSTEM_DIRECTIVES: {
    CPMDInputSection_list secs = new_sections(&arena, 1);
    put_system(&arena, secs, 0, 0.0, 0, 0, bad_dirs, 1);
    arm.inputSections = secs;
    break;
  }
  case ARM_BAD_CPMD_DIRECTIVES: {
    struct CPMDCpmdSection cpmd;
    memset(&cpmd, 0, sizeof(cpmd));
    cpmd.directives = bad_dirs;
    CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
    write_CPMDCpmdSection(&cpmd, cpmd_ptr);
    sec.which = CPMDInputSection_cpmd;
    sec.cpmd = cpmd_ptr;
    attach_section(&arena, &arm, &sec);
    break;
  }
  case ARM_BAD_VDW_DIRECTIVES: {
    struct CPMDVdwSection vdw;
    memset(&vdw, 0, sizeof(vdw));
    vdw.empiricalCorrection = 1;
    vdw.grimme = capn_str("D3");
    vdw.directives = bad_dirs;
    CPMDVdwSection_ptr vdw_ptr = new_CPMDVdwSection(root.seg);
    write_CPMDVdwSection(&vdw, vdw_ptr);
    sec.which = CPMDInputSection_vdwParams;
    sec.vdwParams = vdw_ptr;
    attach_section(&arena, &arm, &sec);
    break;
  }
  }
  struct CommonMethodSpec com;
  memset(&com, 0, sizeof(com));
  switch (fc->field) {
  case ARM_OV_XC:
    set_xc(&arena, &com, "PBE");
    break;
  case ARM_OV_CUTOFF:
    com.planewaveCutoffEv = 80.0 * RY_EV;
    break;
  case ARM_OV_CHARGE:
    com.charge = -1;
    break;
  case ARM_OV_MULT:
    com.spinMultiplicity = 3;
    break;
  case ARM_OV_SMEARING: {
    CommonMethodSpec_Smearing_ptr sm_ptr =
        new_CommonMethodSpec_Smearing(root.seg);
    struct CommonMethodSpec_Smearing sm;
    memset(&sm, 0, sizeof(sm));
    sm.kind = CommonMethodSpec_Smearing_Kind_fermi;
    sm.widthEv = 0.025851999786;
    write_CommonMethodSpec_Smearing(&sm, sm_ptr);
    com.smearing = sm_ptr;
    break;
  }
  case ARM_OV_VDW_S6:
    com.vanDerWaalsS6 = 1.05;
    break;
  }
  unsigned char *buf = serialize_config(
      &arena, make_cpmd_plus_common(&arena, &arm, &com, 1), size);
  capn_free(&arena);
  return buf;
}

static void arm_failure_message(const struct ArmFailureCase *fc, char *dst,
                                size_t dst_size) {
  snprintf(dst, dst_size,
           "%s cannot merge: the cpmd arm has an invalid %s list",
           fc->want_field, fc->want_list);
}

/* Every arm check that rejects an overlay leaves a cpmdc_last_error() that
 * names the overlay field and the list it could not read. */
static void test_merge_arm_check_failures_name_cause(void **state) {
  (void)state;
  const struct ArmFailureCase cases[] = {
      {ARM_OV_SMEARING, ARM_BAD_INPUT_BLOCKS, "smearing", "inputBlocks"},
      {ARM_OV_SMEARING, ARM_BAD_INPUT_SECTIONS, "smearing", "inputSections"},
      {ARM_OV_SMEARING, ARM_BAD_CPMD_DIRECTIVES,
       "smearing", "&CPMD directives"},
      {ARM_OV_VDW_S6, ARM_BAD_INPUT_BLOCKS, "vanDerWaalsS6", "inputBlocks"},
      {ARM_OV_VDW_S6, ARM_BAD_INPUT_SECTIONS, "vanDerWaalsS6", "inputSections"},
      {ARM_OV_VDW_S6, ARM_BAD_VDW_DIRECTIVES,
       "vanDerWaalsS6", "&VDW directives"},
      {ARM_OV_XC, ARM_BAD_INPUT_BLOCKS, "xcFunctionals", "inputBlocks"},
      {ARM_OV_XC, ARM_BAD_INPUT_SECTIONS, "xcFunctionals", "inputSections"},
      {ARM_OV_XC, ARM_BAD_DFT_DIRECTIVES, "xcFunctionals", "&DFT directives"},
      {ARM_OV_CUTOFF, ARM_BAD_INPUT_BLOCKS, "planewaveCutoffEv", "inputBlocks"},
      {ARM_OV_CUTOFF, ARM_BAD_INPUT_SECTIONS,
       "planewaveCutoffEv", "inputSections"},
      {ARM_OV_CUTOFF, ARM_BAD_SYSTEM_DIRECTIVES,
       "planewaveCutoffEv", "&SYSTEM directives"},
      {ARM_OV_CHARGE, ARM_BAD_INPUT_BLOCKS, "charge", "inputBlocks"},
      {ARM_OV_CHARGE, ARM_BAD_INPUT_SECTIONS, "charge", "inputSections"},
      {ARM_OV_CHARGE, ARM_BAD_SYSTEM_DIRECTIVES,
       "charge", "&SYSTEM directives"},
      {ARM_OV_MULT, ARM_BAD_INPUT_BLOCKS, "spinMultiplicity", "inputBlocks"},
      {ARM_OV_MULT, ARM_BAD_INPUT_SECTIONS,
       "spinMultiplicity", "inputSections"},
      {ARM_OV_MULT, ARM_BAD_SYSTEM_DIRECTIVES,
       "spinMultiplicity", "&SYSTEM directives"},
  };
  int failures = 0;
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    const struct ArmFailureCase *fc = &cases[i];
    char want[256];
    arm_failure_message(fc, want, sizeof(want));
    size_t sz = 0;
    unsigned char *buf = build_arm_failure_config(fc, &sz);
    int rc = cpmdc_configure(buf, sz);
    const char *err = cpmdc_last_error();
    if (rc != -1 || !strstr(err, want)) {
      fprintf(stderr,
              "cpmdc_configure: want \"%s\", rc %d, last_error \"%s\"\n",
              want, rc, err);
      failures++;
    }
    CPMDCSession *session = cpmdc_session_create_from_config(buf, sz);
    err = cpmdc_last_error();
    if (session || !strstr(err, want)) {
      fprintf(stderr,
              "cpmdc_session_create_from_config: want \"%s\", session %p, "
              "last_error \"%s\"\n",
              want, (void *)session, err);
      failures++;
    }
    if (session)
      cpmdc_session_destroy(session);
    free(buf);
  }
  assert_int_equal(failures, 0);
}

int main(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr,
            "usage: %s OVERLAY_BIN REJECT_BIN STEP_BIN BAD_PARAMS_BIN "
            "ARM_OVERLAY_BIN\n",
            argv[0]);
    return 2;
  }
  g_overlay = argv[1];
  g_reject = argv[2];
  g_step = argv[3];
  g_bad_params = argv[4];
  g_arm_overlay = argv[5];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_overlay_session_evaluates),
      cmocka_unit_test(test_overlay_rejects_unlowerable_field),
      cmocka_unit_test(test_overlay_kmesh_reaches_deck),
      cmocka_unit_test(test_failed_set_params_keeps_overlay),
      cmocka_unit_test(test_merge_xc_functionals),
      cmocka_unit_test(test_merge_rejects_unmapped_xc_when_arm_wins),
      cmocka_unit_test(test_merge_unset_system_charge_keeps_root),
      cmocka_unit_test(test_merge_omitted_dft_functional_uses_overlay),
      cmocka_unit_test(test_merge_omitted_system_cutoff_uses_overlay),
      cmocka_unit_test(test_merge_omitted_system_multiplicity_keeps_arm),
      cmocka_unit_test(test_merge_rejects_typed_dft_directive_overlay),
      cmocka_unit_test(test_merge_rejects_raw_dft_overlay),
      cmocka_unit_test(test_merge_rejects_inputblocks_dft_overlay),
      cmocka_unit_test(test_merge_rejects_generic_dft_overlay),
      cmocka_unit_test(test_merge_rejects_set_functional_overlay),
      cmocka_unit_test(test_merge_rejects_raw_system_charge_overlay),
      cmocka_unit_test(test_merge_rejects_set_charge_overlay),
      cmocka_unit_test(test_merge_set_other_dft_keeps_overlay_functional),
      cmocka_unit_test(test_merge_wire_default_system_cutoff_takes_common),
      cmocka_unit_test(
          test_merge_wire_default_system_multiplicity_takes_common),
      cmocka_unit_test(test_merge_wire_default_system_charge_takes_common),
      cmocka_unit_test(test_merge_wire_default_dft_functional_takes_common),
      cmocka_unit_test(test_merge_unset_system_cutoff_keeps_root),
      cmocka_unit_test(test_merge_unset_system_multiplicity_keeps_root),
      cmocka_unit_test(test_merge_unset_dft_functional_keeps_root),
      cmocka_unit_test(test_merge_rejects_dft_directive_functional_overlay),
      cmocka_unit_test(test_merge_rejects_system_directive_scalar_overlays),
      cmocka_unit_test(test_merge_rejects_explicit_singlet_against_opaque_system),
      cmocka_unit_test(test_merge_rejects_split_temperature_electron_smearing),
      cmocka_unit_test(test_merge_basis_set),
      cmocka_unit_test(test_merge_cutoff_charge_multiplicity),
      cmocka_unit_test(test_merge_rejects_kmesh_scf_vdw_method_relativity),
      cmocka_unit_test(test_merge_smearing),
      cmocka_unit_test(test_merge_vdw_s6),
      cmocka_unit_test(test_merge_raw_cpmd_and_session_configure),
      cmocka_unit_test(test_precedence_overlay_fills_typed_defaults),
      cmocka_unit_test(test_precedence_typed_value_beats_overlay),
      cmocka_unit_test(test_precedence_root_value_with_default_sections),
      cmocka_unit_test(test_precedence_overlay_cutoff_unit_conversion),
      cmocka_unit_test(test_merge_rejects_explicit_singlet_next_to_carriers),
      cmocka_unit_test(test_merge_accepts_explicit_singlet),
      cmocka_unit_test(test_merge_rejects_typed_functional_with_carrier),
      cmocka_unit_test(test_merge_rejects_typed_cutoff_with_carrier),
      cmocka_unit_test(test_merge_rejects_typed_charge_with_carrier),
      cmocka_unit_test(test_merge_rejects_typed_multiplicity_with_carrier),
      cmocka_unit_test(
          test_merge_rejects_typed_multiplicity_singlet_with_carrier),
      cmocka_unit_test(test_merge_arm_check_failures_name_cause),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
