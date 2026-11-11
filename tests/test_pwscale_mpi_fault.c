/* SPDX-License-Identifier: MIT */
/* Bounded check that a rank-local setup failure still returns on every rank. */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc_pwscale.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int call_fault(const char *fault, int natoms, int npw) {
  cpmdc_pwscale_problem prob;
  double *pos, *forces;
  int *z;
  cpmdc_pwscale_stats stats;
  int rc, rank = 0;
  alarm(8);
  setenv("CPMDC_PWSCALE_MPI_FAULT", fault, 1);
  memset(&prob, 0, sizeof(prob));
  prob.natoms = natoms;
  prob.npw = npw;
  prob.cell_bohr = 10.0;
  prob.alpha = 0.5;
  prob.sigma = 1.0;
  pos = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  forces = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  z = (int *)calloc((size_t)natoms, sizeof(int));
  if (!pos || !forces || !z) {
    free(pos);
    free(forces);
    free(z);
    alarm(0);
    return 2;
  }
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, forces, 1, CPMDC_PWSCALE_DIST_PW,
                              &stats, NULL, 0);
  alarm(0);
  free(pos);
  free(forces);
  free(z);
  if (rc != CPMDC_PWSCALE_ERR_ARG) {
    fprintf(stderr, "rank %d fault %s rc %d\n", rank, fault, rc);
    return 1;
  }
  return 0;
}

int main(void) {
  int rank = 0, size = 1, local = 0, any = 0;
  cpmdc_pwscale_problem prob;
  double pos[6], forces[6];
  int z[2];
  cpmdc_pwscale_stats stats;
  int rc;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    if (rank == 0)
      fprintf(stderr, "mpi fault test needs at least 2 ranks\n");
    MPI_Finalize();
    return 2;
  }
  local = call_fault("gather", 2, 8);
  MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any) {
    MPI_Finalize();
    return 1;
  }
  local = call_fault("team", 2, 8);
  MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any) {
    MPI_Finalize();
    return 1;
  }
  unsetenv("CPMDC_PWSCALE_MPI_FAULT");
  memset(&prob, 0, sizeof(prob));
  prob.natoms = 2;
  prob.npw = 8;
  prob.cell_bohr = 10.0;
  prob.alpha = 0.5;
  prob.sigma = 1.0;
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  {
    cpmdc_pwscale_mpi_bind *bind = NULL;
    alarm(8);
    rc = cpmdc_pwscale_mpi_bind_create(&prob, CPMDC_PWSCALE_DIST_PW, NULL, 0,
                                       &bind);
    if (rc == CPMDC_PWSCALE_OK)
      rc = cpmdc_pwscale_eval_mpi_bound(bind, rank == size - 1 ? NULL : pos, z,
                                        forces, 1, &stats);
    alarm(0);
    cpmdc_pwscale_mpi_bind_free(bind);
    local = rc != CPMDC_PWSCALE_ERR_ARG;
    if (local)
      fprintf(stderr, "rank %d missing positions rc %d\n", rank, rc);
    MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (any) {
      MPI_Finalize();
      return 1;
    }
  }
  memset(&prob, 0, sizeof(prob));
  prob.natoms = 2;
  prob.npw = 8;
  prob.cell_bohr = 10.0;
  prob.alpha = 0.5;
  prob.sigma = 1.0;
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, forces, 1, CPMDC_PWSCALE_DIST_PW,
                              &stats, NULL, 0);
  if (rc != CPMDC_PWSCALE_OK || stats.mpi_ranks != size) {
    if (rank == 0)
      fprintf(stderr, "follow-up eval rc %d ranks %d\n", rc, stats.mpi_ranks);
    MPI_Finalize();
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf("pwscale mpi fault ok\n");
  MPI_Finalize();
  return 0;
}
