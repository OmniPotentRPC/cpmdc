/**
 * A rejected array one-shot must not install parameters or disturb the
 * active session, the selected orbital key, or a result already stored.
 * The next identical session step is a cache hit: an empty message, the
 * same energy and forces, and no new evaluation.
 */
#include "cpmdc.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

int cpmdc_orbital_ledger_selected(long long *key);

static const char *g_params_path = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
static const long long k_key = 42;
static const char k_nan_msg[] = "positions contain a non-finite component";
static const char k_args_msg[] = "invalid arguments";
static const char k_z_msg[] = "atomic numbers must be in 1..118";
static const char k_limit_msg[] =
    "property gradient snapshot supports at most 1365 atoms";
static const char k_empty_msg[] = "CPMDParams buffer is empty";
static const int k_oversize = 1366;

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

static void expect_refusal(CPMDCResult r, const char *message) {
  assert_int_equal(r.ok, 0);
  assert_true(r.energy_h == 0.0);
  assert_string_equal(r.message, message);
  assert_string_equal(cpmdc_last_error(), message);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
}

static void expect_process_result(const double *gradient,
                                  const CPMDCEnergyComponents *parts) {
  double got[9];
  size_t count = 0;
  CPMDCEnergyComponents now;
  memset(&now, 0, sizeof(now));
  assert_int_equal(cpmdc_last_gradient(got, 9, &count), 0);
  assert_int_equal(count, 9);
  assert_memory_equal(got, gradient, sizeof(got));
  assert_int_equal(cpmdc_last_energy_components(&now), 0);
  assert_int_equal(now.valid, 1);
  assert_true(now.etot == parts->etot);
}

struct warm {
  unsigned char *params;
  size_t params_size;
  CPMDCSession *session;
  double energy;
  double forces[9];
  double gradient[9];
  CPMDCEnergyComponents parts;
};

static void warm_session(struct warm *w) {
  size_t count = 0;
  long long selected = 0;
  CPMDCResult first;
  CPMDCResult again;
  memset(w, 0, sizeof(*w));
  w->params = read_file(g_params_path, &w->params_size);
  assert_non_null(w->params);
  w->session = cpmdc_session_create(w->params, w->params_size);
  assert_non_null(w->session);
  assert_int_equal(cpmdc_session_select_orbitals(w->session, k_key), 0);
  first = cpmdc_session_energy_forces(w->session, 3, k_pos, k_z, w->forces);
  assert_int_equal(first.ok, 1);
  assert_string_equal(first.message, "ok");
  assert_true(isfinite(first.energy_h));
  w->energy = first.energy_h;
  assert_int_equal(cpmdc_last_gradient(w->gradient, 9, &count), 0);
  assert_int_equal(count, 9);
  assert_int_equal(cpmdc_last_energy_components(&w->parts), 0);
  assert_int_equal(w->parts.valid, 1);
  assert_true(w->parts.etot == w->energy);
  assert_int_equal(cpmdc_orbital_ledger_selected(&selected), 1);
  assert_true(selected == k_key);
  again = cpmdc_session_energy_forces(w->session, 3, k_pos, k_z, w->forces);
  assert_int_equal(again.ok, 1);
  assert_int_equal(again.message[0], '\0');
  assert_true(again.energy_h == w->energy);
}

static void expect_session_kept(const struct warm *w) {
  double forces[9];
  double gradient[9];
  size_t count = 0;
  long long selected = 0;
  CPMDCResult again;
  CPMDCEnergyComponents now;
  memset(&now, 0, sizeof(now));
  assert_int_equal(cpmdc_last_gradient(gradient, 9, &count), 0);
  assert_int_equal(count, 9);
  assert_memory_equal(gradient, w->gradient, sizeof(gradient));
  assert_int_equal(cpmdc_session_last_gradient(w->session, gradient, 9, &count),
                   0);
  assert_memory_equal(gradient, w->gradient, sizeof(gradient));
  assert_int_equal(cpmdc_last_energy_components(&now), 0);
  assert_true(now.etot == w->parts.etot);
  assert_int_equal(cpmdc_orbital_ledger_selected(&selected), 1);
  assert_true(selected == k_key);
  again = cpmdc_session_energy_forces(w->session, 3, k_pos, k_z, forces);
  assert_int_equal(again.ok, 1);
  /* An empty message is the cache hit. A recompute would say "ok". */
  assert_int_equal(again.message[0], '\0');
  assert_true(again.energy_h == w->energy);
  assert_memory_equal(forces, w->forces, sizeof(forces));
  assert_string_equal(cpmdc_last_error(), "");
}

static void release_warm(struct warm *w) {
  cpmdc_session_destroy(w->session);
  free(w->params);
}

static void refuse_all_wrappers(const struct warm *w, int n_atoms,
                                const double *pos, const int *z,
                                const unsigned char *params, size_t params_size,
                                double *scratch, const char *message) {
  expect_refusal(cpmdc_energy(n_atoms, pos, z, params, params_size), message);
  expect_session_kept(w);
  expect_refusal(cpmdc_energy_gradient(n_atoms, pos, z, params, params_size,
                                       scratch),
                 message);
  expect_session_kept(w);
  expect_refusal(
      cpmdc_energy_forces(n_atoms, pos, z, params, params_size, scratch),
      message);
  expect_session_kept(w);
}

static void test_session_survives_nan(void **state) {
  struct warm w;
  double pos[9];
  double scratch[9];
  (void)state;
  warm_session(&w);
  memcpy(pos, k_pos, sizeof(pos));
  pos[0] = nan("");
  refuse_all_wrappers(&w, 3, pos, k_z, w.params, w.params_size, scratch,
                      k_nan_msg);
  release_warm(&w);
}

static void test_session_survives_bad_length(void **state) {
  struct warm w;
  double scratch[9];
  (void)state;
  warm_session(&w);
  refuse_all_wrappers(&w, 0, k_pos, k_z, w.params, w.params_size, scratch,
                      k_args_msg);
  release_warm(&w);
}

static void test_session_survives_bad_atomic_number(void **state) {
  struct warm w;
  int z[3] = {8, 0, 1};
  double scratch[9];
  (void)state;
  warm_session(&w);
  refuse_all_wrappers(&w, 3, k_pos, z, w.params, w.params_size, scratch,
                      k_z_msg);
  release_warm(&w);
}

static void test_session_survives_null_positions(void **state) {
  struct warm w;
  double scratch[9];
  (void)state;
  warm_session(&w);
  refuse_all_wrappers(&w, 3, NULL, k_z, w.params, w.params_size, scratch,
                      k_args_msg);
  release_warm(&w);
}

static void test_session_survives_oversize(void **state) {
  struct warm w;
  double *pos;
  int *z;
  double scratch[3];
  int i;
  (void)state;
  warm_session(&w);
  pos = (double *)calloc((size_t)k_oversize * 3u, sizeof(*pos));
  z = (int *)malloc((size_t)k_oversize * sizeof(*z));
  assert_non_null(pos);
  assert_non_null(z);
  for (i = 0; i < k_oversize; ++i)
    z[i] = 1;
  /* The gradient buffer is not read when the atom count is refused. */
  refuse_all_wrappers(&w, k_oversize, pos, z, w.params, w.params_size, scratch,
                      k_limit_msg);
  free(pos);
  free(z);
  release_warm(&w);
}

static void warm_process(unsigned char **params, size_t *params_size,
                         double gradient[9], CPMDCEnergyComponents *parts) {
  double forces[9];
  size_t count = 0;
  CPMDCResult r;
  *params = read_file(g_params_path, params_size);
  assert_non_null(*params);
  r = cpmdc_energy_forces(3, k_pos, k_z, *params, *params_size, forces);
  assert_int_equal(r.ok, 1);
  assert_string_equal(r.message, "ok");
  assert_int_equal(cpmdc_last_gradient(gradient, 9, &count), 0);
  assert_int_equal(count, 9);
  assert_int_equal(cpmdc_last_energy_components(parts), 0);
  assert_true(parts->etot == r.energy_h);
}

static void test_process_result_survives_bad_length(void **state) {
  unsigned char *params = NULL;
  size_t params_size = 0;
  double gradient[9];
  double scratch[9];
  CPMDCEnergyComponents parts;
  (void)state;
  warm_process(&params, &params_size, gradient, &parts);
  expect_refusal(cpmdc_energy(0, k_pos, k_z, params, params_size), k_args_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(cpmdc_energy_gradient(0, k_pos, k_z, params, params_size,
                                       scratch),
                 k_args_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(
      cpmdc_energy_forces(0, k_pos, k_z, params, params_size, scratch),
      k_args_msg);
  expect_process_result(gradient, &parts);
  free(params);
}

static void test_process_result_survives_null_and_empty_params(void **state) {
  unsigned char *params = NULL;
  size_t params_size = 0;
  double gradient[9];
  double scratch[9];
  CPMDCEnergyComponents parts;
  (void)state;
  warm_process(&params, &params_size, gradient, &parts);
  expect_refusal(cpmdc_energy(3, NULL, k_z, params, params_size), k_args_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(
      cpmdc_energy_gradient(3, NULL, k_z, params, params_size, scratch),
      k_args_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(
      cpmdc_energy_forces(3, NULL, k_z, params, params_size, scratch),
      k_args_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(cpmdc_energy(3, k_pos, k_z, NULL, 0), k_empty_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(cpmdc_energy_gradient(3, k_pos, k_z, NULL, 0, scratch),
                 k_empty_msg);
  expect_process_result(gradient, &parts);
  expect_refusal(cpmdc_energy_forces(3, k_pos, k_z, NULL, 0, scratch),
                 k_empty_msg);
  expect_process_result(gradient, &parts);
  free(params);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_session_survives_nan),
      cmocka_unit_test(test_session_survives_bad_length),
      cmocka_unit_test(test_session_survives_bad_atomic_number),
      cmocka_unit_test(test_session_survives_null_positions),
      cmocka_unit_test(test_session_survives_oversize),
      cmocka_unit_test(test_process_result_survives_bad_length),
      cmocka_unit_test(test_process_result_survives_null_and_empty_params),
  };
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params_path = argv[1];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
