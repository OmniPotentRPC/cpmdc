/**
 * One-shot entry points reject a non-finite position or cell before
 * evaluation. ok is 0, the energy is 0, and the code is
 * CPMDC_ERROR_INVALID_ARGUMENT.
 */
#include "cpmdc.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

#include "Potentials.capnp.h"
#include <capn.h>

static unsigned char *g_params = NULL;
static size_t g_params_size = 0;
static unsigned char *g_config = NULL;
static size_t g_config_size = 0;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
static const double k_cell[9] = {
    10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0,
};
static const char k_pos_msg[] = "positions contain a non-finite component";
static const char k_cell_msg[] = "cell contains a non-finite component";

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

static void expect_rejected(CPMDCResult r, const char *message) {
  assert_int_equal(r.ok, 0);
  assert_true(r.energy_h == 0.0);
  assert_string_equal(r.message, message);
  assert_string_equal(cpmdc_last_error(), message);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
}

static void copy_pos(double out[9], double bad, int index, int use_bad) {
  memcpy(out, k_pos, 9 * sizeof(*out));
  if (use_bad)
    out[index] = bad;
}

static void reject_energy(double bad, int index) {
  double pos[9];
  copy_pos(pos, bad, index, 1);
  expect_rejected(cpmdc_energy(3, pos, k_z, g_params, g_params_size), k_pos_msg);
}

static void reject_gradient(double bad, int index) {
  double pos[9];
  double grad[9];
  copy_pos(pos, bad, index, 1);
  expect_rejected(cpmdc_energy_gradient(3, pos, k_z, g_params, g_params_size, grad),
                  k_pos_msg);
}

static void reject_forces(double bad, int index) {
  double pos[9];
  double forces[9];
  copy_pos(pos, bad, index, 1);
  expect_rejected(cpmdc_energy_forces(3, pos, k_z, g_params, g_params_size, forces),
                  k_pos_msg);
}

static int encode_force(const double *pos, const double *cell, unsigned char *out,
                        size_t cap) {
  struct capn c;
  capn_init_malloc(&c);
  capn_ptr root = capn_root(&c);
  ForceInput_ptr fp = new_ForceInput(root.seg);
  struct ForceInput v;
  memset(&v, 0, sizeof(v));
  v.pos = capn_new_list64(root.seg, 9);
  v.atmnrs = capn_new_list32(root.seg, 3);
  v.box = capn_new_list64(root.seg, cell ? 9 : 0);
  for (int i = 0; i < 9; ++i)
    capn_set64(v.pos, i, capn_from_f64(pos[i]));
  for (int i = 0; i < 3; ++i)
    capn_set32(v.atmnrs, i, (uint32_t)k_z[i]);
  for (int i = 0; cell && i < 9; ++i)
    capn_set64(v.box, i, capn_from_f64(cell[i]));
  v.lengthUnit = (capn_text){8, "angstrom", NULL};
  v.energyUnit = (capn_text){7, "hartree", NULL};
  write_ForceInput(&v, fp);
  int len =
      capn_setp(root, 0, fp.p) == 0 ? capn_write_mem(&c, out, cap, 0) : -1;
  capn_free(&c);
  return len;
}

static void reject_message(int from_config, const double *pos, const double *cell,
                           const char *message) {
  unsigned char raw[4096];
  int len = encode_force(pos, cell, raw, sizeof(raw));
  size_t need;
  unsigned char *out;
  size_t wrote = 0;
  CPMDCResult r;
  assert_true(len > 0);
  need = cpmdc_potential_result_size_for_force_input(raw, (size_t)len);
  assert_true(need > 0);
  out = (unsigned char *)calloc(1, need);
  assert_non_null(out);
  if (from_config)
    r = cpmdc_calculate_result_from_config(g_config, g_config_size, raw,
                                           (size_t)len, out, need, &wrote);
  else
    r = cpmdc_calculate_result(g_params, g_params_size, raw, (size_t)len, out,
                               need, &wrote);
  free(out);
  expect_rejected(r, message);
}

static void bad_position(double pos[9], double bad, int index) {
  copy_pos(pos, bad, index, 1);
}

static void bad_cell(double cell[9], double bad, int index) {
  memcpy(cell, k_cell, sizeof(k_cell));
  cell[index] = bad;
}

static void test_energy_nan(void **state) {
  (void)state;
  reject_energy(nan(""), 0);
}

static void test_energy_inf(void **state) {
  (void)state;
  reject_energy(INFINITY, 7);
}

static void test_gradient_nan(void **state) {
  (void)state;
  reject_gradient(nan(""), 2);
}

static void test_gradient_inf(void **state) {
  (void)state;
  reject_gradient(-INFINITY, 5);
}

static void test_forces_nan(void **state) {
  (void)state;
  reject_forces(nan(""), 8);
}

static void test_forces_inf(void **state) {
  (void)state;
  reject_forces(INFINITY, 1);
}

static void test_calculate_result_nan_position(void **state) {
  double pos[9];
  (void)state;
  bad_position(pos, nan(""), 3);
  reject_message(0, pos, NULL, k_pos_msg);
}

static void test_calculate_result_inf_position(void **state) {
  double pos[9];
  (void)state;
  bad_position(pos, INFINITY, 6);
  reject_message(0, pos, NULL, k_pos_msg);
}

static void test_calculate_result_nan_cell(void **state) {
  double cell[9];
  (void)state;
  bad_cell(cell, nan(""), 0);
  reject_message(0, k_pos, cell, k_cell_msg);
}

static void test_calculate_result_inf_cell(void **state) {
  double cell[9];
  (void)state;
  bad_cell(cell, INFINITY, 8);
  reject_message(0, k_pos, cell, k_cell_msg);
}

static void test_from_config_nan_position(void **state) {
  double pos[9];
  (void)state;
  bad_position(pos, nan(""), 4);
  reject_message(1, pos, NULL, k_pos_msg);
}

static void test_from_config_inf_position(void **state) {
  double pos[9];
  (void)state;
  bad_position(pos, -INFINITY, 0);
  reject_message(1, pos, NULL, k_pos_msg);
}

static void test_from_config_nan_cell(void **state) {
  double cell[9];
  (void)state;
  bad_cell(cell, nan(""), 4);
  reject_message(1, k_pos, cell, k_cell_msg);
}

static void test_from_config_inf_cell(void **state) {
  double cell[9];
  (void)state;
  bad_cell(cell, INFINITY, 1);
  reject_message(1, k_pos, cell, k_cell_msg);
}

static void test_finite_energy_still_runs(void **state) {
  CPMDCResult r;
  (void)state;
  r = cpmdc_energy(3, k_pos, k_z, g_params, g_params_size);
  assert_int_equal(r.ok, 1);
  assert_true(isfinite(r.energy_h));
  assert_string_equal(cpmdc_last_error(), "");
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_NONE);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s PARAMS_BIN CONFIG_BIN\n", argv[0]);
    return 2;
  }
  g_params = read_file(argv[1], &g_params_size);
  g_config = read_file(argv[2], &g_config_size);
  if (!g_params || !g_config)
    return 2;
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_energy_nan),
      cmocka_unit_test(test_energy_inf),
      cmocka_unit_test(test_gradient_nan),
      cmocka_unit_test(test_gradient_inf),
      cmocka_unit_test(test_forces_nan),
      cmocka_unit_test(test_forces_inf),
      cmocka_unit_test(test_calculate_result_nan_position),
      cmocka_unit_test(test_calculate_result_inf_position),
      cmocka_unit_test(test_calculate_result_nan_cell),
      cmocka_unit_test(test_calculate_result_inf_cell),
      cmocka_unit_test(test_from_config_nan_position),
      cmocka_unit_test(test_from_config_inf_position),
      cmocka_unit_test(test_from_config_nan_cell),
      cmocka_unit_test(test_from_config_inf_cell),
      cmocka_unit_test(test_finite_energy_still_runs),
  };
  int status = cmocka_run_group_tests(tests, NULL, NULL);
  free(g_params);
  free(g_config);
  return status;
}
