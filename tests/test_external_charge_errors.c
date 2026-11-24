/**
 * cpmdc_session_last_external_charge_forces writes a code and a message on
 * every failure and clears both after a successful copy. A deck that already
 * defines EXTERNAL POTENTIAL, with or without ADD, is its own error and is
 * not an evaluation failure.
 */
#include "cpmdc.h"
#include "schema/Potentials.capnp.h"
#include <capn.h>

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;
static const char *g_potential = NULL;
static const char *g_potential_add = NULL;
static const char *g_periodic = NULL;

static const double k_water[] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_numbers[] = {8, 1, 1};

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

static void expect_invalid(const char *fragment) {
  const char *err = cpmdc_last_error();
  assert_non_null(err);
  assert_non_null(strstr(err, fragment));
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
}

static CPMDCSession *open_params(const char *path) {
  size_t n = 0;
  unsigned char *bytes = read_file(path, &n);
  CPMDCSession *session;
  assert_non_null(bytes);
  session = cpmdc_session_create(bytes, n);
  free(bytes);
  assert_non_null(session);
  return session;
}

static void test_null_session(void **state) {
  double forces[3] = {0.0, 0.0, 0.0};
  (void)state;
  assert_int_equal(cpmdc_session_last_external_charge_forces(NULL, forces, 3),
                   -1);
  expect_invalid("null session");
}

static void test_missing_result_and_buffer_paths(void **state) {
  CPMDCSession *session = open_params(g_params);
  double positions[3] = {4.0, 0.0, 0.0};
  double charges[1] = {0.1};
  double widths[1] = {0.30};
  double qm_forces[9];
  double mm_forces[3];
  CPMDCResult step;
  (void)state;

  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 3), -1);
  expect_invalid("not available");

  assert_int_equal(cpmdc_session_set_external_charges(session, 1, positions,
                                                      charges, widths),
                   0);
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 3), -1);
  expect_invalid("not available");

  step = cpmdc_session_energy_forces(session, 3, k_water, k_numbers, qm_forces);
  assert_int_equal(step.ok, 1);

  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 0), -1);
  expect_invalid("wrong length");
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, NULL, 3), -1);
  expect_invalid("null");

  mm_forces[0] = mm_forces[1] = mm_forces[2] = NAN;
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 3), 0);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  assert_int_equal(cpmdc_last_error()[0], '\0');
  assert_true(isfinite(mm_forces[0]));
  assert_true(isfinite(mm_forces[1]));
  assert_true(isfinite(mm_forces[2]));

  cpmdc_session_destroy(session);
}

static capn_text capn_str(const char *s) {
  capn_text t;
  t.str = s;
  t.len = s ? (int)strlen(s) : 0;
  t.seg = NULL;
  return t;
}

static CPMDCSession *open_custom_water_params(const char *system_raw) {
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  CPMDParams_ptr ptr = new_CPMDParams(root.seg);
  struct CPMDParams view;
  memset(&view, 0, sizeof(view));
  view.functional = capn_str("BLYP");
  view.cutOffRy = 50.0;
  view.charge = 0;
  view.multiplicity = 1;
  view.task = capn_str("gradient");

  int nsec = system_raw ? 4 : 3;
  CPMDInputSection_list secs = new_CPMDInputSection_list(root.seg, nsec);
  int idx = 0;

  struct CPMDCpmdSection cpmd;
  memset(&cpmd, 0, sizeof(cpmd));
  cpmd.optimizeWavefunction = 1;
  cpmd.convergenceOrbitals = 1.0e-5;
  cpmd.maxStep = 40;
  CPMDCpmdSection_ptr cp = new_CPMDCpmdSection(root.seg);
  write_CPMDCpmdSection(&cpmd, cp);
  struct CPMDInputSection sc;
  memset(&sc, 0, sizeof(sc));
  sc.which = CPMDInputSection_cpmd;
  sc.cpmd = cp;
  set_CPMDInputSection(&sc, secs, idx++);

  struct CPMDDftSection dft;
  memset(&dft, 0, sizeof(dft));
  dft.functional = capn_str("BLYP");
  CPMDDftSection_ptr dp = new_CPMDDftSection(root.seg);
  write_CPMDDftSection(&dft, dp);
  struct CPMDInputSection sd;
  memset(&sd, 0, sizeof(sd));
  sd.which = CPMDInputSection_dft;
  sd.dft = dp;
  set_CPMDInputSection(&sd, secs, idx++);

  if (system_raw) {
    struct CPMDInputSection sr;
    memset(&sr, 0, sizeof(sr));
    sr.which = CPMDInputSection_raw;
    sr.raw = capn_str(system_raw);
    set_CPMDInputSection(&sr, secs, idx++);
  }

  struct CPMDAtomsSection atoms;
  memset(&atoms, 0, sizeof(atoms));
  CPMDAtomsPseudopotential_list pps =
      new_CPMDAtomsPseudopotential_list(root.seg, 2);
  struct CPMDAtomsPseudopotential pp0;
  memset(&pp0, 0, sizeof(pp0));
  pp0.element = capn_str("O");
  pp0.path = capn_str("O_MT_BLYP.psp");
  pp0.lmax = 1;
  set_CPMDAtomsPseudopotential(&pp0, pps, 0);
  struct CPMDAtomsPseudopotential pp1;
  memset(&pp1, 0, sizeof(pp1));
  pp1.element = capn_str("H");
  pp1.path = capn_str("H_CVB_BLYP.psp");
  pp1.lmax = 0;
  set_CPMDAtomsPseudopotential(&pp1, pps, 1);
  atoms.pseudopotentials = pps;
  CPMDAtomsSection_ptr ap = new_CPMDAtomsSection(root.seg);
  write_CPMDAtomsSection(&atoms, ap);
  struct CPMDInputSection sa;
  memset(&sa, 0, sizeof(sa));
  sa.which = CPMDInputSection_atoms;
  sa.atoms = ap;
  set_CPMDInputSection(&sa, secs, idx++);

  view.inputSections = secs;
  write_CPMDParams(&view, ptr);
  capn_setp(capn_root(&arena), 0, ptr.p);
  uint8_t buf[8192];
  int64_t len = capn_write_mem(&arena, buf, sizeof(buf), 0);
  capn_free(&arena);
  assert_true(len > 0);
  CPMDCSession *session = cpmdc_session_create(buf, (size_t)len);
  assert_non_null(session);
  return session;
}

static void test_typed_failure_replaces_previous_code(void **state) {
  CPMDCSession *session = open_params(g_periodic);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double water_bohr[9];
  double cell_bohr[9] = {20.0, 0.0, 0.0, 0.0, 20.0,
                         0.0,  0.0, 0.0, 20.0};
  double energy;
  double forces[9];
  (void)state;

  for (size_t i = 0; i < 9; ++i)
    water_bohr[i] = k_water[i] / 0.529177210903;
  assert_int_equal(cpmdc_session_set_external_charges(
                       session, 1, positions, charges, widths),
                   0);
  assert_int_equal(cpmdc_session_calculate_au(
                       session, 3, water_bohr, k_numbers, cell_bohr, &energy,
                       forces, NULL),
                   -1);
  assert_int_equal(cpmdc_last_error_code(),
                   CPMDC_ERROR_UNSUPPORTED_PERIODIC_EMBEDDING);

  assert_null(cpmdc_session_create(NULL, 0));
  assert_non_null(strstr(cpmdc_last_error(), "buffer is empty"));
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
  cpmdc_session_destroy(session);
}

static void test_isolated_cell_is_accepted_and_evaluates(void **state) {
  CPMDCSession *session = open_params(g_params);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double water_bohr[9];
  double cell_bohr[9] = {20.0, 0.0, 0.0, 0.0, 22.0,
                         0.0,  0.0, 0.0, 24.0};
  double zero_cell[9] = {0.0};
  double energy_box = 0.0;
  double energy_zero = 0.0;
  double qm_forces[9];
  double mm_forces[3];
  CPMDCStressTensor stress;
  (void)state;

  for (size_t i = 0; i < 9; ++i)
    water_bohr[i] = k_water[i] / 0.529177210903;
  assert_int_equal(cpmdc_session_set_external_charges(
                       session, 1, positions, charges, widths),
                   0);
  assert_int_equal(cpmdc_session_calculate_au(
                       session, 3, water_bohr, k_numbers, cell_bohr,
                       &energy_box, qm_forces, NULL),
                   0);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  assert_true(isfinite(energy_box));
  for (size_t i = 0; i < 9; ++i)
    assert_true(isfinite(qm_forces[i]));
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 3), 0);
  for (size_t i = 0; i < 3; ++i)
    assert_true(isfinite(mm_forces[i]));
  if (cpmdc_session_last_stress(session, &stress) == 0)
    assert_int_equal(stress.valid, 0);

  assert_int_equal(cpmdc_session_calculate_au(
                       session, 3, water_bohr, k_numbers, zero_cell,
                       &energy_zero, qm_forces, NULL),
                   0);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  assert_true(isfinite(energy_zero));
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(session, mm_forces, 3), 0);
  for (size_t i = 0; i < 3; ++i)
    assert_true(isfinite(mm_forces[i]));

  cpmdc_session_destroy(session);
}

static void test_cluster_cell_accepted_and_unset_cell_refused(void **state) {
  CPMDCSession *cluster = open_custom_water_params(
      "&SYSTEM\n CLUSTER\n POISSON SOLVER HOCKNEY\n ANGSTROM\n"
      " CELL\n  10.0 1.0 1.0 0.0 0.0 0.0\n CUTOFF\n  50.0\n&END\n");
  CPMDCSession *unset = open_custom_water_params(NULL);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double water_bohr[9];
  double cell_bohr[9] = {20.0, 0.0, 0.0, 0.0, 20.0,
                         0.0,  0.0, 0.0, 20.0};
  double zero_cell[9] = {0.0};
  double energy = 0.0;
  double qm_forces[9];
  double mm_forces[3];
  (void)state;

  for (size_t i = 0; i < 9; ++i)
    water_bohr[i] = k_water[i] / 0.529177210903;

  assert_int_equal(cpmdc_session_set_external_charges(
                       cluster, 1, positions, charges, widths),
                   0);
  assert_int_equal(cpmdc_session_calculate_au(
                       cluster, 3, water_bohr, k_numbers, cell_bohr, &energy,
                       qm_forces, NULL),
                   0);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  assert_true(isfinite(energy));
  assert_int_equal(
      cpmdc_session_last_external_charge_forces(cluster, mm_forces, 3), 0);
  for (size_t i = 0; i < 3; ++i)
    assert_true(isfinite(mm_forces[i]));

  assert_int_equal(cpmdc_session_set_external_charges(
                       unset, 1, positions, charges, widths),
                   0);
  assert_int_equal(cpmdc_session_calculate_au(
                       unset, 3, water_bohr, k_numbers, zero_cell, &energy,
                       qm_forces, NULL),
                   0);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  assert_true(isfinite(energy));

  assert_int_equal(cpmdc_session_calculate_au(
                       unset, 3, water_bohr, k_numbers, cell_bohr, &energy,
                       qm_forces, NULL),
                   -1);
  assert_non_null(strstr(cpmdc_last_error(), "periodic electrostatic embedding"));
  assert_int_equal(cpmdc_last_error_code(),
                   CPMDC_ERROR_UNSUPPORTED_PERIODIC_EMBEDDING);

  cpmdc_session_destroy(unset);
  cpmdc_session_destroy(cluster);
}

static void expect_deck_clash(const char *path) {
  CPMDCSession *session = open_params(path);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double qm_forces[9];
  CPMDCResult clash;

  assert_int_equal(cpmdc_session_set_external_charges(session, 1, positions,
                                                      charges, widths),
                   0);
  clash =
      cpmdc_session_energy_forces(session, 3, k_water, k_numbers, qm_forces);
  assert_int_equal(clash.ok, 0);
  assert_non_null(strstr(clash.message, "EXTERNAL POTENTIAL"));
  assert_string_equal(cpmdc_last_error(), clash.message);
  assert_int_equal(cpmdc_last_error_code(),
                   CPMDC_ERROR_DECK_EXTERNAL_POTENTIAL);
  assert_int_not_equal(cpmdc_last_error_code(), CPMDC_ERROR_EVALUATION_FAILED);

  cpmdc_session_destroy(session);
}

static void test_deck_external_potential(void **state) {
  (void)state;
  expect_deck_clash(g_potential);
}

static void test_deck_external_potential_add(void **state) {
  (void)state;
  expect_deck_clash(g_potential_add);
}

/* A periodic deck is periodic whether or not the call passes a cell, so a
 * session with charges refuses it on every evaluation path. The refusal is
 * a different variant from the deck clash and has to replace that code. */
static void expect_periodic_deck_refused(CPMDCResult r) {
  assert_int_equal(r.ok, 0);
  assert_non_null(strstr(r.message, "periodic electrostatic embedding"));
  assert_string_equal(cpmdc_last_error(), r.message);
  assert_int_equal(cpmdc_last_error_code(),
                   CPMDC_ERROR_UNSUPPORTED_PERIODIC_EMBEDDING);
}

static void test_periodic_deck_is_refused_without_a_cell(void **state) {
  CPMDCSession *session = open_params(g_periodic);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double qm_forces[9];
  double gradient[9];
  CPMDCResult bare;
  (void)state;

  bare = cpmdc_session_energy_forces(session, 3, k_water, k_numbers, qm_forces);
  assert_int_equal(bare.ok, 1);

  assert_int_equal(cpmdc_session_set_external_charges(session, 1, positions,
                                                      charges, widths),
                   0);
  expect_periodic_deck_refused(
      cpmdc_session_energy_forces(session, 3, k_water, k_numbers, qm_forces));
  expect_periodic_deck_refused(
      cpmdc_session_energy_gradient(session, 3, k_water, k_numbers, gradient));
  expect_periodic_deck_refused(
      cpmdc_session_energy(session, 3, k_water, k_numbers));

  assert_int_equal(
      cpmdc_session_set_external_charges(session, 0, NULL, NULL, NULL), 0);
  bare = cpmdc_session_energy_forces(session, 3, k_water, k_numbers, qm_forces);
  assert_int_equal(bare.ok, 1);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
  cpmdc_session_destroy(session);
}

static void test_periodic_refusal_replaces_the_deck_clash_code(void **state) {
  CPMDCSession *clash_session = open_params(g_potential);
  CPMDCSession *periodic_session = open_params(g_periodic);
  double positions[3] = {3.0, 0.0, 0.0};
  double charges[1] = {1.0};
  double widths[1] = {0.30};
  double qm_forces[9];
  (void)state;

  assert_int_equal(cpmdc_session_set_external_charges(
                       clash_session, 1, positions, charges, widths),
                   0);
  assert_int_equal(cpmdc_session_energy_forces(clash_session, 3, k_water,
                                               k_numbers, qm_forces)
                       .ok,
                   0);
  assert_int_equal(cpmdc_last_error_code(),
                   CPMDC_ERROR_DECK_EXTERNAL_POTENTIAL);

  assert_int_equal(cpmdc_session_set_external_charges(
                       periodic_session, 1, positions, charges, widths),
                   0);
  expect_periodic_deck_refused(cpmdc_session_energy_forces(
      periodic_session, 3, k_water, k_numbers, qm_forces));

  cpmdc_session_destroy(periodic_session);
  cpmdc_session_destroy(clash_session);
}

int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: %s PARAMS POTENTIAL POTENTIAL_ADD PERIODIC\n",
            argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_potential = argv[2];
  g_potential_add = argv[3];
  g_periodic = argv[4];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_null_session),
      cmocka_unit_test(test_missing_result_and_buffer_paths),
      cmocka_unit_test(test_typed_failure_replaces_previous_code),
      cmocka_unit_test(test_isolated_cell_is_accepted_and_evaluates),
      cmocka_unit_test(test_cluster_cell_accepted_and_unset_cell_refused),
      cmocka_unit_test(test_deck_external_potential),
      cmocka_unit_test(test_deck_external_potential_add),
      cmocka_unit_test(test_periodic_deck_is_refused_without_a_cell),
      cmocka_unit_test(test_periodic_refusal_replaces_the_deck_clash_code),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
