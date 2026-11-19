/**
 * PotentialConfig common-overlay lowering: session creation, evaluation, and
 * loud rejection of unsound combinations.
 */
#define _GNU_SOURCE
#include "cpmdc.h"
#include "Potentials.capnp.h"

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
      cmocka_unit_test(test_merge_basis_set),
      cmocka_unit_test(test_merge_cutoff_charge_multiplicity),
      cmocka_unit_test(test_merge_rejects_kmesh_scf_vdw_method_relativity),
      cmocka_unit_test(test_merge_smearing),
      cmocka_unit_test(test_merge_vdw_s6),
      cmocka_unit_test(test_merge_raw_cpmd_and_session_configure),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
