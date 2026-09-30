/* The one-shot entry points evaluate the same params more than once in
 * one process. Each call re-applies the params, so the second call must
 * find the pseudopotential directory from the host directory and leave
 * the process there afterwards. */
#define _DEFAULT_SOURCE
#include "cpmdc.h"

#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

static const char *g_params = NULL;

/* Si-N at 1.7 Angstrom: the two species the fixture names. */
static const double k_positions[6] = {0.0, 0.0, 0.0, 1.7, 0.0, 0.0};
static const int k_numbers[2] = {14, 7};

struct saved_env {
  int set;
  char value[PATH_MAX];
};

static void save_env(const char *name, struct saved_env *out) {
  const char *v = getenv(name);
  out->set = v != NULL;
  out->value[0] = '\0';
  if (v)
    snprintf(out->value, sizeof(out->value), "%s", v);
}

static void restore_env(const char *name, const struct saved_env *saved) {
  if (saved->set)
    setenv(name, saved->value, 1);
  else
    unsetenv(name);
}

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

static void report(const char *tag, const CPMDCResult *r) {
  char cwd[PATH_MAX];
  if (!getcwd(cwd, sizeof(cwd)))
    snprintf(cwd, sizeof(cwd), "?");
  print_message("[ %s ] ok=%d energy=%.8f message='%s' last_error='%s' "
                "cwd='%s' CPMDC_PSEUDO_DIR='%s' CPMD_PP_LIBRARY_PATH='%s'\n",
                tag, r->ok, r->energy_h, r->message, cpmdc_last_error(), cwd,
                getenv("CPMDC_PSEUDO_DIR") ? getenv("CPMDC_PSEUDO_DIR") : "",
                getenv("CPMD_PP_LIBRARY_PATH") ? getenv("CPMD_PP_LIBRARY_PATH")
                                               : "");
}

static void call_twice(const unsigned char *params, size_t params_size,
                       const char *tag, CPMDCResult *first,
                       CPMDCResult *second) {
  double grad_a[6];
  double grad_b[6];
  char before[PATH_MAX];
  char after[PATH_MAX];
  char label[64];
  assert_non_null(getcwd(before, sizeof(before)));
  *first = cpmdc_energy_gradient(2, k_positions, k_numbers, params,
                                 params_size, grad_a);
  snprintf(label, sizeof(label), "%s call 1", tag);
  report(label, first);
  assert_non_null(getcwd(after, sizeof(after)));
  assert_string_equal(after, before);
  *second = cpmdc_energy_gradient(2, k_positions, k_numbers, params,
                                  params_size, grad_b);
  snprintf(label, sizeof(label), "%s call 2", tag);
  report(label, second);
  assert_non_null(getcwd(after, sizeof(after)));
  assert_string_equal(after, before);
  assert_int_equal(first->ok, 1);
  assert_int_equal(second->ok, 1);
  assert_true(first->energy_h < -5.0);
  assert_true(fabs(first->energy_h - second->energy_h) < 1e-5);
  for (int i = 0; i < 6; ++i)
    assert_true(fabs(grad_a[i] - grad_b[i]) < 1e-4);
}

static void test_one_shot_twice_absolute_dir(void **state) {
  unsigned char *params = NULL;
  size_t params_size = 0;
  CPMDCResult first;
  CPMDCResult second;
  (void)state;
  if (!g_params || !getenv("CPMDC_PSEUDO_DIR"))
    skip();
  params = read_file(g_params, &params_size);
  assert_non_null(params);
  call_twice(params, params_size, "absolute", &first, &second);
  free(params);
}

/* Only CPMD_PP_LIBRARY_PATH is set, and it is relative to the host
 * directory. The library rewrite of that variable must still resolve
 * from the host directory on the next call. */
static void test_one_shot_twice_relative_library_path(void **state) {
  struct saved_env pseudo_env;
  struct saved_env library_env;
  struct saved_env pp_env;
  unsigned char *params = NULL;
  size_t params_size = 0;
  char original[PATH_MAX];
  char parent[PATH_MAX];
  char relative[PATH_MAX];
  const char *pseudo;
  const char *slash;
  CPMDCResult first;
  CPMDCResult second;
  (void)state;
  pseudo = getenv("CPMDC_PSEUDO_DIR");
  if (!g_params || !pseudo || pseudo[0] != '/')
    skip();
  params = read_file(g_params, &params_size);
  assert_non_null(params);
  assert_non_null(getcwd(original, sizeof(original)));
  save_env("CPMDC_PSEUDO_DIR", &pseudo_env);
  save_env("CPMD_PP_LIBRARY_PATH", &library_env);
  save_env("PP_LIBRARY_PATH", &pp_env);
  snprintf(parent, sizeof(parent), "%s", pseudo);
  while (strlen(parent) > 1 && parent[strlen(parent) - 1] == '/')
    parent[strlen(parent) - 1] = '\0';
  slash = strrchr(parent, '/');
  assert_non_null(slash);
  assert_true(slash > parent);
  snprintf(relative, sizeof(relative), "%s", slash + 1);
  parent[slash - parent] = '\0';
  assert_int_equal(chdir(parent), 0);
  unsetenv("CPMDC_PSEUDO_DIR");
  setenv("CPMD_PP_LIBRARY_PATH", relative, 1);
  memset(&first, 0, sizeof(first));
  memset(&second, 0, sizeof(second));
  call_twice(params, params_size, "relative", &first, &second);
  assert_int_equal(chdir(original), 0);
  restore_env("CPMDC_PSEUDO_DIR", &pseudo_env);
  restore_env("CPMD_PP_LIBRARY_PATH", &library_env);
  restore_env("PP_LIBRARY_PATH", &pp_env);
  free(params);
}

int main(int argc, char **argv) {
  if (argc >= 2)
    g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_one_shot_twice_absolute_dir),
      cmocka_unit_test(test_one_shot_twice_relative_library_path),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
