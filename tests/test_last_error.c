/**
 * Every public call that can fail leaves a non-empty cpmdc_last_error().
 * A successful configuration call clears it. Snapshot readers and the
 * capabilities size query do not write it.
 */
#include "cpmdc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;

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

static void expect_has(const char *fragment) {
  const char *err = cpmdc_last_error();
  assert_non_null(err);
  assert_non_null(strstr(err, fragment));
}

static void prime_config_error(void) {
  assert_int_equal(cpmdc_configure(NULL, 0), -1);
  expect_has("PotentialConfig");
}

static void test_set_params_failure(void **state) {
  (void)state;
  unsigned char garbage[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  prime_config_error();
  assert_int_equal(cpmdc_set_params(NULL, 0), -1);
  expect_has("CPMDParams buffer is empty");
  prime_config_error();
  assert_int_equal(cpmdc_set_params(garbage, sizeof(garbage)), -1);
  expect_has("CPMDParams parse failed");
}

static void test_session_create_failure(void **state) {
  (void)state;
  unsigned char garbage[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  prime_config_error();
  assert_null(cpmdc_session_create(NULL, 0));
  expect_has("CPMDParams buffer is empty");
  prime_config_error();
  assert_null(cpmdc_session_create(garbage, sizeof(garbage)));
  expect_has("CPMDParams parse failed");
}

static void test_session_set_params_failure(void **state) {
  (void)state;
  prime_config_error();
  assert_int_equal(cpmdc_session_set_params(NULL, NULL, 0), -1);
  expect_has("invalid session");
}

static void test_configure_failure(void **state) {
  (void)state;
  assert_int_equal(cpmdc_set_params(NULL, 0), -1);
  expect_has("CPMDParams");
  assert_int_equal(cpmdc_configure(NULL, 0), -1);
  expect_has("PotentialConfig buffer is empty");
}

static void test_session_create_from_config_failure(void **state) {
  (void)state;
  assert_int_equal(cpmdc_set_params(NULL, 0), -1);
  expect_has("CPMDParams");
  assert_null(cpmdc_session_create_from_config(NULL, 0));
  expect_has("PotentialConfig buffer is empty");
}

static void test_session_configure_failure(void **state) {
  (void)state;
  prime_config_error();
  assert_int_equal(cpmdc_session_configure(NULL, NULL, 0), -1);
  expect_has("invalid session");
}

static void test_bind_calculator_failure(void **state) {
  (void)state;
#if !defined(CPMDC_HAS_CPMD)
  prime_config_error();
  assert_int_equal(cpmdc_bind_calculator(1), -1);
  expect_has("calculator bind refused");
  prime_config_error();
  assert_int_equal(cpmdc_adopt_calculator_comm(NULL, 0, 1), -1);
  expect_has("calculator bind refused");
#else
  (void)state;
#endif
}

static void expect_result_error(CPMDCResult r, const char *fragment) {
  assert_int_equal(r.ok, 0);
  expect_has(fragment);
  assert_string_equal(cpmdc_last_error(), r.message);
}

static void test_energy_failures(void **state) {
  (void)state;
  prime_config_error();
  expect_result_error(cpmdc_energy(0, NULL, NULL, NULL, 0), "invalid arguments");
  prime_config_error();
  expect_result_error(
      cpmdc_energy_gradient(0, NULL, NULL, NULL, 0, NULL), "CPMDParams");
  prime_config_error();
  expect_result_error(cpmdc_energy_forces(0, NULL, NULL, NULL, 0, NULL),
                      "CPMDParams");
}

static void test_session_energy_failures(void **state) {
  (void)state;
  prime_config_error();
  expect_result_error(cpmdc_session_energy(NULL, 0, NULL, NULL),
                      "invalid arguments");
  prime_config_error();
  expect_result_error(
      cpmdc_session_energy_gradient(NULL, 0, NULL, NULL, NULL), "null session");
  prime_config_error();
  expect_result_error(cpmdc_session_energy_forces(NULL, 0, NULL, NULL, NULL),
                      "null session");
  prime_config_error();
  expect_result_error(
      cpmdc_session_calculate_forces(NULL, NULL, 0, NULL, 0),
      "invalid arguments");
  prime_config_error();
  expect_result_error(
      cpmdc_session_calculate_result(NULL, NULL, 0, NULL, 0, NULL),
      "invalid arguments");
}

static void test_one_shot_failures(void **state) {
  (void)state;
  prime_config_error();
  expect_result_error(cpmdc_calculate_result(NULL, 0, NULL, 0, NULL, 0, NULL),
                      "invalid arguments");
  prime_config_error();
  expect_result_error(
      cpmdc_calculate_result_from_config(NULL, 0, NULL, 0, NULL, 0, NULL),
      "invalid arguments");
}

static void test_force_input_size_failure(void **state) {
  (void)state;
  prime_config_error();
  assert_int_equal(cpmdc_potential_result_size_for_force_input(NULL, 0), 0);
  expect_has("ForceInput");
}

static void test_success_clears_last_error(void **state) {
  (void)state;
  size_t n = 0;
  unsigned char *params = read_file(g_params, &n);
  assert_non_null(params);
  prime_config_error();
  assert_int_equal(cpmdc_set_params(params, n), 0);
  assert_int_equal(cpmdc_last_error()[0], '\0');
  prime_config_error();
  CPMDCSession *session = cpmdc_session_create(params, n);
  assert_non_null(session);
  assert_int_equal(cpmdc_last_error()[0], '\0');
  cpmdc_session_destroy(session);
  free(params);
}

static void test_readers_do_not_write_last_error(void **state) {
  (void)state;
  assert_int_equal(cpmdc_set_params(NULL, 0), -1);
  expect_has("CPMDParams buffer is empty");
  CPMDCStressTensor stress;
  memset(&stress, 0x5a, sizeof(stress));
  assert_int_equal(cpmdc_last_stress(&stress), -1);
  expect_has("CPMDParams buffer is empty");
  assert_int_equal(cpmdc_adopted_comm(NULL, 0), -1);
  expect_has("CPMDParams buffer is empty");
  CPMDCEnergyComponents energy;
  assert_int_equal(cpmdc_last_energy_components(&energy), -1);
  expect_has("CPMDParams buffer is empty");
  size_t need = 0;
  assert_int_equal(cpmdc_capabilities_result(NULL, 0, &need), -1);
  assert_true(need > 0);
  expect_has("CPMDParams buffer is empty");
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_set_params_failure),
      cmocka_unit_test(test_session_create_failure),
      cmocka_unit_test(test_session_set_params_failure),
      cmocka_unit_test(test_configure_failure),
      cmocka_unit_test(test_session_create_from_config_failure),
      cmocka_unit_test(test_session_configure_failure),
      cmocka_unit_test(test_bind_calculator_failure),
      cmocka_unit_test(test_energy_failures),
      cmocka_unit_test(test_session_energy_failures),
      cmocka_unit_test(test_one_shot_failures),
      cmocka_unit_test(test_force_input_size_failure),
      cmocka_unit_test(test_success_clears_last_error),
      cmocka_unit_test(test_readers_do_not_write_last_error),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
