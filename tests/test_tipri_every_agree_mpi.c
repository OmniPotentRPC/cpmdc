/**
 * CPMDC_TIPRI_EVERY is read on rank 0 and broadcast before tipri.
 *
 * The ranks below see different values. Every rank leaves the collective
 * with rank 0's interval, and a value longer than 16 characters fails
 * every rank. A hung rank ends the test through the alarm.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int cpmdc_agree_tipri_every(int *every, const char **err);

static int g_rank = 0;
static int g_fail = 0;

static void expect(int cond, const char *what) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_rank, what);
  g_fail = 1;
}

static void agree(const char *value, int want_rc, int want_every, const char *what) {
  int every = -1;
  const char *err = "unset";
  int rc;
  if (value)
    setenv("CPMDC_TIPRI_EVERY", value, 1);
  else
    unsetenv("CPMDC_TIPRI_EVERY");
  MPI_Barrier(MPI_COMM_WORLD);
  rc = cpmdc_agree_tipri_every(&every, &err);
  expect(rc == want_rc, what);
  expect(every == want_every, what);
  if (want_rc == 0)
    expect(err == NULL, what);
  else
    expect(err && strstr(err, "16") != NULL, what);
}

int main(int argc, char **argv) {
  int total = 0;
  (void)argc;
  (void)argv;
  alarm(60);
  if (MPI_Init(NULL, NULL) != MPI_SUCCESS)
    return 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  agree(g_rank == 0 ? "4" : "99", 0, 4, "ranks take rank 0's interval");
  agree(g_rank == 0 ? "1" : "12345678901234567", 0, 1,
        "a long value on another rank is ignored");
  agree(g_rank == 0 ? NULL : "7", 0, 0, "unset on rank 0 prints none");
  agree(g_rank == 0 ? "12345678901234567" : "2", -1, 0,
        "a value longer than 16 characters fails every rank");
  MPI_Allreduce(&g_fail, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
  return total ? 1 : 0;
}
