/**
 * A cold deck memfd is closed once the deck has been read. Thirty-two
 * session create and destroy cycles must not grow the process fd count.
 */
#include "cpmdc.h"

#include <dirent.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_params = NULL;
static const char *g_root = NULL;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n < 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n + 1u);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (n > 0 && fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  buf[n] = '\0';
  if (size)
    *size = (size_t)n;
  return buf;
}

static int count_fds(void) {
  DIR *dir = opendir("/proc/self/fd");
  struct dirent *ent;
  int n = 0;
  if (!dir)
    return -1;
  while ((ent = readdir(dir)) != NULL) {
    if (ent->d_name[0] == '.')
      continue;
    n++;
  }
  closedir(dir);
  return n;
}

static void test_deck_is_closed_and_fds_stay_flat(void **state) {
  char path[512];
  char *source;
  const char *routine;
  const char *vdw;
  const char *eval;
  const char *closed;
  unsigned char *params;
  size_t n = 0;
  int before;
  int after;
  int cycle;
  (void)state;

  snprintf(path, sizeof(path), "%s/src/cpmd_embed_c_api.F90", g_root);
  source = (char *)read_file(path, NULL);
  assert_non_null(source);
  routine = strstr(source, "SUBROUTINE run_embed_scf");
  assert_non_null(routine);
  vdw = strstr(routine, "CALL vdw_wf_alloc");
  assert_non_null(vdw);
  eval = strstr(vdw, "CALL embed_eval_energy_grad");
  assert_non_null(eval);
  closed = strstr(vdw, "close_cold_deck");
  assert_non_null(closed);
  assert_true(closed < eval);
  free(source);

  params = read_file(g_params, &n);
  assert_non_null(params);
  before = count_fds();
  assert_true(before > 0);
  for (cycle = 0; cycle < 32; ++cycle) {
    CPMDCSession *session = cpmdc_session_create(params, n);
    double forces[9];
    CPMDCResult result;
    assert_non_null(session);
    result = cpmdc_session_energy_forces(session, 3, k_pos, k_z, forces);
    assert_int_equal(result.ok, 1);
    cpmdc_session_destroy(session);
  }
  after = count_fds();
  assert_int_equal(after, before);
  free(params);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s PARAMS_BIN SOURCE_ROOT\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_root = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_deck_is_closed_and_fds_stay_flat),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
