/* Pseudopotential directory errors, species-block force order, and the
 * directory CPMD writes into. */
#define _DEFAULT_SOURCE
#include "cpmdc.h"
#include "cpmdc_embed_host.h"

#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cmocka.h>

static const char *g_params = NULL;
static const char *g_encode_py = NULL;
static const char *g_capnp = NULL;
static const char *g_schema = NULL;
static const char *g_params_text = NULL;

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

static char *make_temp_dir(void) {
  char tmpl[] = "/tmp/cpmdc-pp-XXXXXX";
  char *dir = mkdtemp(tmpl);
  if (!dir)
    return NULL;
  return strdup(dir);
}

static void remove_tree(const char *dir) {
  DIR *handle;
  struct dirent *ent;
  if (!dir)
    return;
  handle = opendir(dir);
  if (!handle)
    return;
  while ((ent = readdir(handle)) != NULL) {
    char path[PATH_MAX];
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
      continue;
    snprintf(path, sizeof(path), "%s/%s", dir, ent->d_name);
    unlink(path);
  }
  closedir(handle);
  rmdir(dir);
}

static int file_exists(const char *dir, const char *name) {
  char path[PATH_MAX];
  struct stat st;
  snprintf(path, sizeof(path), "%s/%s", dir, name);
  return stat(path, &st) == 0;
}

static void test_missing_pseudopotential_directory_message(void **state) {
  struct saved_env pseudo;
  struct saved_env library;
  char buf[512];
  char *dir;
  (void)state;
  save_env("CPMDC_PSEUDO_DIR", &pseudo);
  save_env("CPMD_PP_LIBRARY_PATH", &library);

  unsetenv("CPMDC_PSEUDO_DIR");
  unsetenv("CPMD_PP_LIBRARY_PATH");
  assert_int_equal(cpmdc_pseudopotential_directory(buf, sizeof(buf)), -1);
  assert_non_null(strstr(buf, "CPMDC_PSEUDO_DIR"));
  assert_non_null(strstr(buf, "CPMD_PP_LIBRARY_PATH"));
  assert_non_null(strstr(buf, "not set"));
  assert_null(strstr(buf, "orbitals not converged"));

  setenv("CPMDC_PSEUDO_DIR", "/no/such/cpmdc-pseudopotential-dir", 1);
  assert_int_equal(cpmdc_pseudopotential_directory(buf, sizeof(buf)), -1);
  assert_non_null(strstr(buf, "not a directory"));
  assert_null(strstr(buf, "orbitals not converged"));

  dir = make_temp_dir();
  assert_non_null(dir);
  unsetenv("CPMDC_PSEUDO_DIR");
  setenv("CPMD_PP_LIBRARY_PATH", dir, 1);
  assert_int_equal(cpmdc_pseudopotential_directory(buf, sizeof(buf)), 0);
  assert_string_equal(buf, dir);

  setenv("CPMDC_PSEUDO_DIR", dir, 1);
  setenv("CPMD_PP_LIBRARY_PATH", "/no/such/cpmdc-pseudopotential-dir", 1);
  assert_int_equal(cpmdc_pseudopotential_directory(buf, sizeof(buf)), 0);
  assert_string_equal(buf, dir);

  remove_tree(dir);
  free(dir);
  restore_env("CPMDC_PSEUDO_DIR", &pseudo);
  restore_env("CPMD_PP_LIBRARY_PATH", &library);
}

static void test_species_order_scatters_forces(void **state) {
  int atomic_numbers[] = {1, 8, 1};
  int species_z[] = {8, 1};
  int species_n[] = {1, 2};
  int map[3];
  double species_grad[] = {0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 2.0, 2.1, 2.2};
  double grad[9];
  int grouped_z[] = {8, 1, 1};
  int identity[3];
  (void)state;

  assert_int_equal(
      cpmdc_species_order_map(3, atomic_numbers, 2, species_z, species_n, map),
      0);
  assert_int_equal(map[0], 1);
  assert_int_equal(map[1], 0);
  assert_int_equal(map[2], 2);

  memset(grad, 0, sizeof(grad));
  cpmdc_scatter_species_gradient(3, map, species_grad, grad);
  assert_true(fabs(grad[0] - 1.0) < 1e-15);
  assert_true(fabs(grad[1] - 1.1) < 1e-15);
  assert_true(fabs(grad[2] - 1.2) < 1e-15);
  assert_true(fabs(grad[3] - 0.1) < 1e-15);
  assert_true(fabs(grad[4] - 0.2) < 1e-15);
  assert_true(fabs(grad[5] - 0.3) < 1e-15);
  assert_true(fabs(grad[6] - 2.0) < 1e-15);
  assert_true(fabs(grad[7] - 2.1) < 1e-15);
  assert_true(fabs(grad[8] - 2.2) < 1e-15);

  assert_int_equal(
      cpmdc_species_order_map(3, grouped_z, 2, species_z, species_n, identity),
      0);
  assert_int_equal(identity[0], 0);
  assert_int_equal(identity[1], 1);
  assert_int_equal(identity[2], 2);

  species_n[0] = 2;
  assert_int_equal(
      cpmdc_species_order_map(3, atomic_numbers, 2, species_z, species_n, map),
      -1);
}

static void test_output_cwd_leaves_pseudopotential_directory(void **state) {
  struct saved_env pseudo_env;
  struct saved_env library_env;
  struct saved_env pp_env;
  char original[PATH_MAX];
  char *host = NULL;
  char *pseudo = NULL;
  char *scratch = NULL;
  char marker[PATH_MAX];
  FILE *fp;
  (void)state;

  assert_non_null(getcwd(original, sizeof(original)));
  save_env("CPMD_PP_LIBRARY_PATH", &pseudo_env);
  save_env("PP_LIBRARY_PATH", &pp_env);
  save_env("CPMDC_PSEUDO_DIR", &library_env);

  host = make_temp_dir();
  pseudo = make_temp_dir();
  scratch = make_temp_dir();
  assert_non_null(host);
  assert_non_null(pseudo);
  assert_non_null(scratch);
  snprintf(marker, sizeof(marker), "%s/KEEP", pseudo);
  fp = fopen(marker, "w");
  assert_non_null(fp);
  fputs("keep\n", fp);
  fclose(fp);

  assert_int_equal(chdir(host), 0);
  assert_int_equal(cpmdc_prepare_pp_cwd(pseudo), 0);
  assert_int_equal(cpmdc_enter_output_cwd(scratch), 0);
  fp = fopen("RESTART.1", "w");
  assert_non_null(fp);
  fputs("restart\n", fp);
  fclose(fp);
  fp = fopen("LATEST", "w");
  assert_non_null(fp);
  fputs("latest\n", fp);
  fclose(fp);
  fp = fopen("GEOMETRY", "w");
  assert_non_null(fp);
  fputs("geometry\n", fp);
  fclose(fp);
  fp = fopen("GEOMETRY.xyz", "w");
  assert_non_null(fp);
  fputs("xyz\n", fp);
  fclose(fp);
  assert_int_equal(cpmdc_restore_host_cwd(), 0);
  assert_int_equal(chdir(original), 0);

  assert_true(file_exists(pseudo, "KEEP"));
  assert_false(file_exists(pseudo, "RESTART.1"));
  assert_false(file_exists(pseudo, "LATEST"));
  assert_false(file_exists(pseudo, "GEOMETRY"));
  assert_false(file_exists(pseudo, "GEOMETRY.xyz"));
  assert_true(file_exists(scratch, "RESTART.1"));
  assert_true(file_exists(scratch, "LATEST"));
  assert_true(file_exists(scratch, "GEOMETRY"));
  assert_true(file_exists(scratch, "GEOMETRY.xyz"));

  remove_tree(host);
  remove_tree(pseudo);
  remove_tree(scratch);
  free(host);
  free(pseudo);
  free(scratch);
  restore_env("CPMD_PP_LIBRARY_PATH", &pseudo_env);
  restore_env("PP_LIBRARY_PATH", &pp_env);
  restore_env("CPMDC_PSEUDO_DIR", &library_env);
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

#if defined(CPMDC_HAS_CPMD)
static void test_energy_missing_pp_dir_sets_last_error(void **state) {
  struct saved_env pseudo;
  struct saved_env library;
  unsigned char *params = NULL;
  size_t params_size = 0;
  CPMDCSession *session = NULL;
  double positions[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692,
                         0.0, -0.7572, -0.4692};
  int atomic_numbers[3] = {8, 1, 1};
  double forces[9];
  CPMDCResult result;
  char message[512];
  char last[512];
  memset(&result, 0, sizeof(result));
  (void)state;
  if (!g_params)
    skip();
  save_env("CPMDC_PSEUDO_DIR", &pseudo);
  save_env("CPMD_PP_LIBRARY_PATH", &library);
  unsetenv("CPMDC_PSEUDO_DIR");
  unsetenv("CPMD_PP_LIBRARY_PATH");
  params = read_file(g_params, &params_size);
  if (params)
    session = cpmdc_session_create(params, params_size);
  if (session)
    result = cpmdc_session_energy_forces(session, 3, positions, atomic_numbers,
                                         forces);
  snprintf(message, sizeof(message), "%s", result.message);
  snprintf(last, sizeof(last), "%s", cpmdc_last_error());
  cpmdc_session_destroy(session);
  free(params);
  restore_env("CPMDC_PSEUDO_DIR", &pseudo);
  restore_env("CPMD_PP_LIBRARY_PATH", &library);
  assert_non_null(params);
  assert_non_null(session);

  assert_int_equal(result.ok, 0);
  assert_string_equal(message, last);
  assert_non_null(strstr(message, "CPMDC_PSEUDO_DIR"));
  assert_non_null(strstr(message, "CPMD_PP_LIBRARY_PATH"));
  assert_non_null(strstr(message, "not set"));
  assert_null(strstr(message, "orbitals not converged"));
}

struct dir_ent {
  char name[256];
  off_t size;
  time_t mtime;
};

static int cmp_ent(const void *a, const void *b) {
  return strcmp(((const struct dir_ent *)a)->name,
                ((const struct dir_ent *)b)->name);
}

static struct dir_ent *snapshot_dir(const char *dir, int *count) {
  DIR *handle = opendir(dir);
  struct dirent *ent;
  struct dir_ent *rows = NULL;
  int n = 0;
  int cap = 0;
  *count = -1;
  if (!handle)
    return NULL;
  while ((ent = readdir(handle)) != NULL) {
    char path[PATH_MAX];
    struct stat st;
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
      continue;
    if ((int)strlen(ent->d_name) >= 256)
      continue;
    snprintf(path, sizeof(path), "%s/%s", dir, ent->d_name);
    if (stat(path, &st) != 0)
      continue;
    if (n == cap) {
      int next = cap == 0 ? 32 : cap * 2;
      struct dir_ent *grown =
          (struct dir_ent *)realloc(rows, (size_t)next * sizeof(*rows));
      if (!grown) {
        closedir(handle);
        free(rows);
        return NULL;
      }
      rows = grown;
      cap = next;
    }
    snprintf(rows[n].name, sizeof(rows[n].name), "%s", ent->d_name);
    rows[n].size = st.st_size;
    rows[n].mtime = st.st_mtime;
    n++;
  }
  closedir(handle);
  if (n > 1)
    qsort(rows, (size_t)n, sizeof(*rows), cmp_ent);
  *count = n;
  return rows;
}

static int snapshots_equal(const struct dir_ent *a, int na,
                           const struct dir_ent *b, int nb) {
  int i;
  if (na != nb)
    return 0;
  for (i = 0; i < na; ++i) {
    if (strcmp(a[i].name, b[i].name) != 0 || a[i].size != b[i].size ||
        a[i].mtime != b[i].mtime)
      return 0;
  }
  return 1;
}

static int encode_with_scratch(const char *scratch, char *bin_path,
                               size_t bin_cap) {
  char text_path[] = "/tmp/cpmdc-params-XXXXXX";
  int fd;
  FILE *in;
  FILE *out;
  char *body;
  long n;
  const char *nl;
  char cmd[8192];
  int rc;
  fd = mkstemp(text_path);
  if (fd < 0)
    return -1;
  close(fd);
  in = fopen(g_params_text, "rb");
  if (!in)
    return -1;
  if (fseek(in, 0, SEEK_END) != 0) {
    fclose(in);
    return -1;
  }
  n = ftell(in);
  if (n <= 0) {
    fclose(in);
    return -1;
  }
  rewind(in);
  body = (char *)malloc((size_t)n + 1);
  if (!body) {
    fclose(in);
    return -1;
  }
  if (fread(body, 1, (size_t)n, in) != (size_t)n) {
    free(body);
    fclose(in);
    return -1;
  }
  fclose(in);
  body[n] = '\0';
  nl = strchr(body, '\n');
  if (!nl) {
    free(body);
    return -1;
  }
  out = fopen(text_path, "w");
  if (!out) {
    free(body);
    return -1;
  }
  fprintf(out, "(\n  scratchDir = \"%s\",\n%s", scratch, nl + 1);
  fclose(out);
  free(body);
  snprintf(bin_path, bin_cap, "/tmp/cpmdc-params-bin-XXXXXX");
  fd = mkstemp(bin_path);
  if (fd < 0) {
    unlink(text_path);
    return -1;
  }
  close(fd);
  snprintf(cmd, sizeof(cmd), "python3 '%s' '%s' '%s' CPMDParams '%s' '%s'",
           g_encode_py, g_capnp, g_schema, text_path, bin_path);
  rc = system(cmd);
  unlink(text_path);
  return rc == 0 ? 0 : -1;
}

static void test_split_species_and_output_directory(void **state) {
  const char *pseudo;
  char *scratch = NULL;
  char bin_path[PATH_MAX];
  unsigned char *params = NULL;
  size_t params_size = 0;
  struct dir_ent *before = NULL;
  struct dir_ent *after = NULL;
  int before_n = 0;
  int after_n = 0;
  char original[PATH_MAX];
  double grouped[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692,
                       0.0, -0.7572, -0.4692};
  int grouped_z[3] = {8, 1, 1};
  double split[9] = {0.0, 0.7572, -0.4692, 0.0, 0.0, 0.1173,
                     0.0, -0.7572, -0.4692};
  int split_z[3] = {1, 8, 1};
  double forces_grouped[9];
  double forces_split[9];
  CPMDCSession *session_grouped = NULL;
  CPMDCSession *session_split = NULL;
  CPMDCResult grouped_result;
  CPMDCResult split_result;
  double force_scale;
  int i;
  (void)state;

  pseudo = getenv("CPMDC_PSEUDO_DIR");
  if (!pseudo || !pseudo[0] || !g_params_text || !g_encode_py || !g_capnp ||
      !g_schema)
    skip();
  assert_non_null(getcwd(original, sizeof(original)));
  scratch = make_temp_dir();
  assert_non_null(scratch);
  assert_int_equal(encode_with_scratch(scratch, bin_path, sizeof(bin_path)), 0);
  params = read_file(bin_path, &params_size);
  unlink(bin_path);
  assert_non_null(params);
  before = snapshot_dir(pseudo, &before_n);
  assert_non_null(before);
  assert_true(before_n >= 0);

  session_grouped = cpmdc_session_create(params, params_size);
  assert_non_null(session_grouped);
  grouped_result = cpmdc_session_energy_forces(
      session_grouped, 3, grouped, grouped_z, forces_grouped);
  cpmdc_session_destroy(session_grouped);
  assert_int_equal(chdir(original), 0);
  assert_int_equal(grouped_result.ok, 1);

  session_split = cpmdc_session_create(params, params_size);
  assert_non_null(session_split);
  split_result = cpmdc_session_energy_forces(session_split, 3, split, split_z,
                                             forces_split);
  cpmdc_session_destroy(session_split);
  assert_int_equal(chdir(original), 0);

  after = snapshot_dir(pseudo, &after_n);
  assert_int_equal(split_result.ok, 1);
  assert_true(grouped_result.energy_h < -10.0);
  assert_true(fabs(grouped_result.energy_h - split_result.energy_h) < 1e-6);
  force_scale = 0.0;
  for (i = 0; i < 9; ++i)
    force_scale += fabs(forces_grouped[i]);
  assert_true(force_scale > 1e-4);
  for (i = 0; i < 3; ++i) {
    assert_true(fabs(forces_split[i] - forces_grouped[3 + i]) < 1e-5);
    assert_true(fabs(forces_split[3 + i] - forces_grouped[i]) < 1e-5);
    assert_true(fabs(forces_split[6 + i] - forces_grouped[6 + i]) < 1e-5);
  }
  assert_non_null(after);
  assert_true(snapshots_equal(before, before_n, after, after_n));
  assert_false(file_exists(pseudo, "RESTART.1"));
  assert_false(file_exists(pseudo, "LATEST"));
  assert_false(file_exists(pseudo, "GEOMETRY"));
  assert_false(file_exists(pseudo, "GEOMETRY.xyz"));
  assert_true(file_exists(scratch, "GEOMETRY.xyz") ||
              file_exists(scratch, "GEOMETRY") ||
              file_exists(scratch, "LATEST") ||
              file_exists(scratch, "RESTART.1"));

  free(before);
  free(after);
  free(params);
  remove_tree(scratch);
  free(scratch);
}
#endif

int main(int argc, char **argv) {
  if (argc >= 2)
    g_params = argv[1];
  if (argc >= 5) {
    g_encode_py = argv[2];
    g_capnp = argv[3];
    g_schema = argv[4];
  }
  if (argc >= 6)
    g_params_text = argv[5];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_missing_pseudopotential_directory_message),
      cmocka_unit_test(test_species_order_scatters_forces),
      cmocka_unit_test(test_output_cwd_leaves_pseudopotential_directory),
#if defined(CPMDC_HAS_CPMD)
      cmocka_unit_test(test_energy_missing_pp_dir_sets_last_error),
      cmocka_unit_test(test_split_species_and_output_directory),
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
