/**
 * Two ranks, one calculator, a cold setup failure on rank 1.
 *
 * The rank that failed used to free the group while the other rank
 * entered the next setup collective, and the job waited. Every rank now
 * agrees first, tears the stand-in group down, and reports the same
 * text. The launcher timeout is what fails the run when that agreement
 * is missing. A stopgm, a cold deck close, and the output directory are
 * the three failures.
 */
#define _POSIX_C_SOURCE 200809L

#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_rank = 0;
static int g_fail = 0;

static const double k_pos[9] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692,
};
static const int k_z[3] = {8, 1, 1};
static const char k_stop[] = "CPMD stopgm during embed SCF";
static const char k_close[] = "cannot close the cold deck";
static const char k_output[] = "CPMD output directory is not a directory";
static const char *g_only = NULL;

static int selected(const char *name) { return !g_only || strcmp(g_only, name) == 0; }

static void check(int cond, const char *what);

static void expect_fault(const char *fault, const char *message,
                         unsigned char *params, size_t params_size) {
  double pos[9];
  double forces[9];
  CPMDCResult stopped;
  int oks[2];
  char messages[2][512];
  if (!selected(fault))
    return;
  setenv("CPMDC_TEST_FAULT", fault, 1);
  memcpy(pos, k_pos, sizeof(pos));
  stopped = cpmdc_energy_forces(3, pos, k_z, params, params_size, forces);
  check(stopped.ok == 0, "a one-rank failure fails the call");
  check(strcmp(stopped.message, message) == 0, "the failure text is reported");
  check(strcmp(cpmdc_last_error(), stopped.message) == 0, "last error matches");
  check(stopped.energy_h == 0.0, "a failure does not publish an energy");
  MPI_Allgather(&stopped.ok, 1, MPI_INT, oks, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(stopped.message, 512, MPI_CHAR, messages, 512, MPI_CHAR,
                MPI_COMM_WORLD);
  check(oks[0] == 0 && oks[1] == 0, "both ranks fail");
  check(strcmp(messages[0], message) == 0 && strcmp(messages[1], message) == 0,
        "both ranks report the failure");
}

static void check(int cond, const char *what) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_rank, what);
  g_fail = 1;
}

static unsigned char *read_params(const char *path, size_t *size) {
  FILE *file = path ? fopen(path, "rb") : NULL;
  long length;
  unsigned char *buf;
  if (!file)
    return NULL;
  if (fseek(file, 0, SEEK_END) != 0 || (length = ftell(file)) <= 0 ||
      fseek(file, 0, SEEK_SET) != 0) {
    fclose(file);
    return NULL;
  }
  buf = (unsigned char *)malloc((size_t)length);
  if (!buf || fread(buf, 1, (size_t)length, file) != (size_t)length) {
    free(buf);
    fclose(file);
    return NULL;
  }
  fclose(file);
  *size = (size_t)length;
  return buf;
}

int main(int argc, char **argv) {
  unsigned char *params;
  size_t params_size = 0;
  double pos[9];
  double forces[9];
  CPMDCResult again;
  int size = 0;
  int any = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    if (g_rank == 0)
      fprintf(stderr, "needs 2 ranks and a CPMDParams file\n");
    MPI_Finalize();
    return 1;
  }
  params = read_params(argc > 1 ? argv[1] : NULL, &params_size);
  if (!params) {
    if (g_rank == 0)
      fprintf(stderr, "needs 2 ranks and a CPMDParams file\n");
    MPI_Finalize();
    return 1;
  }
  if (argc > 2)
    g_only = argv[2];

  expect_fault("cold-setup-stop", k_stop, params, params_size);
  expect_fault("cold-deck-close", k_close, params, params_size);
  expect_fault("cold-output-dir", k_output, params, params_size);

  unsetenv("CPMDC_TEST_FAULT");
  memcpy(pos, k_pos, sizeof(pos));
  pos[0] += 0.02;
  again = cpmdc_energy_forces(3, pos, k_z, params, params_size, forces);
  check(again.ok == 1, "the next finite step still evaluates");
  check(strcmp(again.message, "ok") == 0, "the next step says ok");

  free(params);
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_rank == 0)
    printf("%s cold setup stop on one rank\n", any ? "FAIL" : "ok");
  MPI_Finalize();
  return any ? 1 : 0;
}
