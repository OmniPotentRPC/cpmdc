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
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
