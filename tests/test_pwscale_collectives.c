/* SPDX-License-Identifier: MIT */
/* Copyright (c) cpmdc developers */
/* Counts the collectives of one rank force call through the PMPI shim, and
 * prints a digest of the result bytes so separate runs can be compared. */
#include "cpmdc_pwscale.h"
#include "pmpi_count.h"

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REPEATS 3

static int g_fail = 0;

static uint64_t digest(double energy, const double *forces, int n3) {
  uint64_t h = 14695981039346656037ULL;
  const unsigned char *p = (const unsigned char *)&energy;
  size_t i;
  for (i = 0; i < sizeof(energy); ++i)
    h = (h ^ p[i]) * 1099511628211ULL;
  p = (const unsigned char *)forces;
  for (i = 0; i < (size_t)n3 * sizeof(double); ++i)
    h = (h ^ p[i]) * 1099511628211ULL;
  return h;
}

static long long total_of(const long long *c) {
  long long t = 0;
  int k;
  for (k = 0; k < PMPI_COUNT_KINDS; ++k)
    t += c[k];
  return t;
}

/* Every rank must see the same per-call counts and the same result bytes. */
static void same_on_every_rank(const char *what, long long total, uint64_t h) {
  long long lo = 0, hi = 0;
  unsigned long long hlo = 0, hhi = 0, hv = (unsigned long long)h;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  PMPI_Allreduce(&total, &lo, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  PMPI_Allreduce(&total, &hi, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  PMPI_Allreduce(&hv, &hlo, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  PMPI_Allreduce(&hv, &hhi, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  if (lo != hi || hlo != hhi) {
    if (rank == 0)
      fprintf(stderr, "%s differs across ranks: calls %lld..%lld\n", what, lo,
              hi);
    g_fail = 1;
  }
}

/* Returns this path's collectives per call. */
static long long report(const char *path, int dist, int natoms, int npw,
                        const long long *before, const long long *after,
                        int calls, uint64_t h) {
  int rank = 0, size = 1, k;
  long long d[PMPI_COUNT_KINDS];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (k = 0; k < PMPI_COUNT_KINDS; ++k) {
    d[k] = after[k] - before[k];
    if (d[k] % calls != 0) {
      if (rank == 0)
        fprintf(stderr, "%s %s count %lld over %d calls is uneven\n", path,
                pmpi_count_name(k), d[k], calls);
      g_fail = 1;
    }
    d[k] /= calls;
  }
  same_on_every_rank(path, total_of(d), h);
  if (rank != 0)
    return total_of(d);
  printf("path %s np %d dist %s natoms %d npw %d per_call_total %lld", path,
         size, dist == CPMDC_PWSCALE_DIST_PW ? "pw" : "atoms", natoms, npw,
         total_of(d));
  for (k = 0; k < PMPI_COUNT_KINDS; ++k)
    if (d[k] != 0)
      printf(" %s=%lld", pmpi_count_name(k), d[k]);
  printf(" digest %016llx\n", (unsigned long long)h);
  return total_of(d);
}

static void one_shot(int natoms, int npw, int dist) {
  cpmdc_pwscale_problem prob;
  cpmdc_pwscale_stats stats;
  cpmdc_pwscale_mpi_bind *bind = NULL;
  long long before[PMPI_COUNT_KINDS], after[PMPI_COUNT_KINDS];
  long long per_one_shot, per_bind, per_bound;
  double *pos = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  double *f = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  int *z = (int *)calloc((size_t)natoms, sizeof(int));
  uint64_t h = 0, first = 0;
  int i, rc;
  if (!pos || !f || !z) {
    g_fail = 1;
    goto done;
  }
  prob.natoms = natoms;
  prob.npw = npw;
  prob.cell_bohr = 12.0;
  prob.alpha = 0.4;
  prob.sigma = 0.9;
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  pmpi_count_snapshot(before);
  for (i = 0; i < REPEATS; ++i) {
    rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, f, 1, dist, &stats, NULL, 0);
    if (rc != CPMDC_PWSCALE_OK) {
      g_fail = 1;
      goto done;
    }
    h = digest(stats.energy_h, f, natoms * 3);
    if (i == 0)
      first = h;
    else if (h != first)
      g_fail = 1;
  }
  pmpi_count_snapshot(after);
  per_one_shot = report("eval_mpi", dist, natoms, npw, before, after, REPEATS, h);

  pmpi_count_snapshot(before);
  rc = cpmdc_pwscale_mpi_bind_create(&prob, dist, NULL, 0, &bind);
  pmpi_count_snapshot(after);
  if (rc != CPMDC_PWSCALE_OK) {
    g_fail = 1;
    goto done;
  }
  per_bind = report("bind_create", dist, natoms, npw, before, after, 1, 0);
  pmpi_count_snapshot(before);
  for (i = 0; i < REPEATS; ++i) {
    rc = cpmdc_pwscale_eval_mpi_bound(bind, pos, z, f, 1, &stats);
    if (rc != CPMDC_PWSCALE_OK) {
      g_fail = 1;
      goto done;
    }
    if (digest(stats.energy_h, f, natoms * 3) != first)
      g_fail = 1;
  }
  pmpi_count_snapshot(after);
  per_bound = report("eval_mpi_bound", dist, natoms, npw, before, after,
                     REPEATS, digest(stats.energy_h, f, natoms * 3));
  if (per_bind != 1 || per_bound != 1 || per_bound + per_bind != per_one_shot) {
    fprintf(stderr, "bind %lld bound %lld one-shot %lld\n", per_bind,
            per_bound, per_one_shot);
    g_fail = 1;
  }
done:
  cpmdc_pwscale_mpi_bind_free(bind);
  free(pos);
  free(f);
  free(z);
}

int main(void) {
  int rank = 0, any = 0;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  one_shot(96, 32768, CPMDC_PWSCALE_DIST_PW);
  one_shot(96, 32768, CPMDC_PWSCALE_DIST_ATOMS);
  one_shot(17, 100, CPMDC_PWSCALE_DIST_PW);
  one_shot(17, 100, CPMDC_PWSCALE_DIST_ATOMS);
  PMPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0)
    printf(any ? "pwscale collectives failed\n" : "pwscale collectives ok\n");
  MPI_Finalize();
  return any ? 1 : 0;
}
