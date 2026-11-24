/**
 * cpmdc_session_last_external_charge_forces writes a code and a message on
 * every failure and clears both after a successful copy. A deck that already
 * defines EXTERNAL POTENTIAL, with or without ADD, is its own error and is
 * not an evaluation failure.
 */
#include "cpmdc.h"

#include <math.h>
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

static void test_typed_failure_replaces_previous_code(void **state) {
  CPMDCSession *session = open_params(g_params);
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
      cmocka_unit_test(test_deck_external_potential),
      cmocka_unit_test(test_deck_external_potential_add),
      cmocka_unit_test(test_periodic_deck_is_refused_without_a_cell),
      cmocka_unit_test(test_periodic_refusal_replaces_the_deck_clash_code),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
