/* SPDX-License-Identifier: MIT */
#include "cpmdc_pwscale.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(CPMDC_PWSCALE_HAVE_MPI)
#include <mpi.h>
#endif

static int g_fail = 0;

static void expect(int cond, const char *msg) {
  if (!cond) {
    fprintf(stderr, "fail: %s\n", msg);
    g_fail = 1;
  }
}

static int forces_close(const double *got, const double *ref, int n3,
                        double *max_abs) {
  int i, ok = 1;
  double peak = 0.0;
  for (i = 0; i < n3; ++i) {
    double gap = fabs(got[i] - ref[i]);
    if (gap > peak)
      peak = gap;
    if (!cpmdc_pwscale_within(got[i], ref[i]))
      ok = 0;
  }
  if (max_abs)
    *max_abs = peak;
  return ok;
}

static void oracle(void) {
  cpmdc_pwscale_problem prob;
  double pos[3];
  int z[1];
  double forces[3];
  cpmdc_pwscale_stats stats;
  int rc;
  /* Independent sum in tests/pwscale_oracle.py for this geometry. */
  const double energy_ref = 1.7274791563162685e+02;
  const double force_ref[3] = {8.1165490211079074e+00, -7.4606561327671184e+00,
                               2.7471849298200780e+01};
  memset(&prob, 0, sizeof(prob));
  prob.natoms = 1;
  prob.npw = 8;
  prob.cell_bohr = 10.0;
  prob.alpha = 0.5;
  prob.sigma = 1.0;
  pos[0] = 0.2;
  pos[1] = -0.3;
  pos[2] = 0.4;
  z[0] = 3;
  rc = cpmdc_pwscale_eval(&prob, pos, z, forces, 1, CPMDC_PWSCALE_DIST_PW, &stats);
  expect(rc == CPMDC_PWSCALE_OK, "oracle eval");
  if (rc == CPMDC_PWSCALE_OK) {
    expect(fabs(stats.energy_h - energy_ref) < 1e-8, "oracle energy");
    expect(fabs(forces[0] - force_ref[0]) < 1e-8, "oracle fx");
    expect(fabs(forces[1] - force_ref[1]) < 1e-8, "oracle fy");
    expect(fabs(forces[2] - force_ref[2]) < 1e-8, "oracle fz");
    expect(stats.g_visited == 8, "oracle g count");
    expect(stats.atom_visits == 8, "oracle atom visits");
    if (fabs(stats.energy_h - energy_ref) >= 1e-8)
      fprintf(stderr, "energy got %.16e ref %.16e\n", stats.energy_h, energy_ref);
  }
}

static void compare_threads(int natoms, int npw, int threads, int dist) {
  cpmdc_pwscale_problem prob;
  double *pos, *f_serial, *f_par;
  int *z;
  cpmdc_pwscale_stats serial, parallel;
  double peak = 0.0;
  int rc;
  char msg[128];

  prob.natoms = natoms;
  prob.npw = npw;
  prob.cell_bohr = 12.0;
  prob.alpha = 0.4;
  prob.sigma = 0.9;
  pos = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  f_serial = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  f_par = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  z = (int *)calloc((size_t)natoms, sizeof(int));
  expect(pos && f_serial && f_par && z, "alloc");
  if (!pos || !f_serial || !f_par || !z)
    goto done;
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval(&prob, pos, z, f_serial, 1, dist, &serial);
  expect(rc == CPMDC_PWSCALE_OK, "serial eval");
  rc = cpmdc_pwscale_eval(&prob, pos, z, f_par, threads, dist, &parallel);
  expect(rc == CPMDC_PWSCALE_OK, "parallel eval");
  if (rc != CPMDC_PWSCALE_OK)
    goto done;
  snprintf(msg, sizeof(msg), "energy dist %d threads %d", dist, threads);
  expect(cpmdc_pwscale_within(parallel.energy_h, serial.energy_h), msg);
  expect(forces_close(f_par, f_serial, natoms * 3, &peak), "forces");
  expect(parallel.g_visited == serial.g_visited, "g count");
  expect(parallel.atom_visits == serial.atom_visits, "atom visits");
  expect(serial.g_visited == npw, "g visited is npw");
  expect(parallel.workers == threads, "actual OpenMP team");
  expect(parallel.mpi_ranks == 1, "threaded eval is one rank");
  if (parallel.workers != threads)
    fprintf(stderr, "dist %d requested %d actual team %d\n", dist, threads,
            parallel.workers);
  if (!cpmdc_pwscale_within(parallel.energy_h, serial.energy_h) || peak > 1e-6)
    fprintf(stderr, "dist %d threads %d dE %.3e peak |dF| %.3e\n", dist, threads,
            parallel.energy_h - serial.energy_h, peak);
done:
  free(pos);
  free(f_serial);
  free(f_par);
  free(z);
}

static void zero_mode_and_bound(void) {
  cpmdc_pwscale_problem prob;
  double *pos, *forces;
  int *z;
  cpmdc_pwscale_stats stats;
  int rc;
  prob.natoms = 4;
  prob.npw = 27;
  prob.cell_bohr = 12.0;
  prob.alpha = 0.4;
  prob.sigma = 0.9;
  pos = (double *)calloc(12, sizeof(double));
  forces = (double *)calloc(12, sizeof(double));
  z = (int *)calloc(4, sizeof(int));
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval(&prob, pos, z, forces, 4, CPMDC_PWSCALE_DIST_PW, &stats);
  expect(rc == CPMDC_PWSCALE_OK, "g0 eval");
  expect(stats.g_visited == 27, "g0 plane waves");
  expect(stats.atom_visits == 26L * 4L, "g0 skips the origin");
  prob.natoms = 16;
  prob.npw = 128;
  free(pos);
  free(forces);
  free(z);
  pos = (double *)calloc(48, sizeof(double));
  forces = (double *)calloc(48, sizeof(double));
  z = (int *)calloc(16, sizeof(int));
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval(&prob, pos, z, forces, 4, CPMDC_PWSCALE_DIST_PW, &stats);
  expect(rc == CPMDC_PWSCALE_OK, "bound eval");
  expect(stats.alloc_bytes < (long long)prob.npw * prob.natoms * (long long)sizeof(double),
         "scratch stays under a plane-wave by atom matrix");
  prob.npw = 0;
  rc = cpmdc_pwscale_eval(&prob, pos, z, forces, 1, CPMDC_PWSCALE_DIST_PW, &stats);
  expect(rc == CPMDC_PWSCALE_ERR_ARG, "reject empty grid");
  free(pos);
  free(forces);
  free(z);
}

#if defined(CPMDC_PWSCALE_HAVE_MPI)
static void mpi_case(int natoms, int npw, int dist) {
  cpmdc_pwscale_problem prob;
  double *pos, *f_serial, *f_mpi;
  int *z;
  cpmdc_pwscale_stats serial, parallel;
  int rc, rank = 0, size = 1;
  double peak = 0.0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  prob.natoms = natoms;
  prob.npw = npw;
  prob.cell_bohr = 12.0;
  prob.alpha = 0.4;
  prob.sigma = 0.9;
  pos = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  f_serial = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  f_mpi = (double *)calloc((size_t)natoms * 3u, sizeof(double));
  z = (int *)calloc((size_t)natoms, sizeof(int));
  if (!pos || !f_serial || !f_mpi || !z) {
    g_fail = 1;
    goto done;
  }
  cpmdc_pwscale_fill_geometry(&prob, pos, z);
  rc = cpmdc_pwscale_eval(&prob, pos, z, f_serial, 1, dist, &serial);
  if (rc != CPMDC_PWSCALE_OK)
    g_fail = 1;
  rc = cpmdc_pwscale_eval_mpi(&prob, pos, z, f_mpi, 1, dist, &parallel, NULL, 0);
  if (rc != CPMDC_PWSCALE_OK)
    g_fail = 1;
  if (rc == CPMDC_PWSCALE_OK && serial.energy_h == serial.energy_h) {
    if (!cpmdc_pwscale_within(parallel.energy_h, serial.energy_h) ||
        !forces_close(f_mpi, f_serial, natoms * 3, &peak) ||
        parallel.g_visited != serial.g_visited ||
        parallel.atom_visits != serial.atom_visits ||
        parallel.mpi_ranks != size || parallel.workers != 1) {
      if (rank == 0)
        fprintf(stderr, "mpi dist %d dE %.3e peak |dF| %.3e visits %lld %lld\n",
                dist, parallel.energy_h - serial.energy_h, peak,
                parallel.g_visited, parallel.atom_visits);
      g_fail = 1;
    }
  }
done:
  free(pos);
  free(f_serial);
  free(f_mpi);
  free(z);
}
#endif

static void tolerance_and_sample(void) {
  double de = -1.0, df = -1.0;
  int counts = -1;
  double got[1] = {1.0};
  double ref[1] = {1.0};
  int ok;
  expect(cpmdc_pwscale_within(1.0, 1.0), "equal values");
  expect(cpmdc_pwscale_within(1.0 + 5e-9, 1.0), "inside 1e-8 at unit scale");
  expect(!cpmdc_pwscale_within(1.0 + 2e-8, 1.0), "outside 1e-8 at unit scale");
  expect(cpmdc_pwscale_within(100.0 + 5e-7, 100.0), "scale uses the larger magnitude");
  expect(!cpmdc_pwscale_within(100.0 + 2e-6, 100.0), "larger magnitude stays bounded");
  expect(!cpmdc_pwscale_within(2.0, 0.0), "a large gap against zero fails");
  ok = cpmdc_pwscale_sample_check(1.0, got, 1.0, ref, 1, 8, 8, 8, 8, &de, &df,
                                  &counts);
  expect(ok && counts == 1 && de == 0.0 && df == 0.0, "matching sample");
  ok = cpmdc_pwscale_sample_check(2.0, got, 1.0, ref, 1, 8, 8, 8, 8, &de, &df,
                                  &counts);
  expect(!ok && counts == 1 && de == 1.0, "poisoned energy is this sample");
  ok = cpmdc_pwscale_sample_check(1.0, got, 1.0, ref, 1, 8, 9, 8, 8, &de, &df,
                                  &counts);
  expect(!ok && counts == 0, "a visit mismatch is not a hardcoded pass");
}

int main(void) {
  int threads;
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  int rank = 0, size = 1, any = 0;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size == 1) {
#endif
    oracle();
    tolerance_and_sample();
    zero_mode_and_bound();
    for (threads = 2; threads <= 8; threads *= 2) {
      compare_threads(24, 96, threads, CPMDC_PWSCALE_DIST_PW);
      compare_threads(24, 96, threads, CPMDC_PWSCALE_DIST_ATOMS);
    }
    compare_threads(17, 100, 3, CPMDC_PWSCALE_DIST_PW);
    compare_threads(17, 100, 3, CPMDC_PWSCALE_DIST_ATOMS);
#if defined(CPMDC_PWSCALE_HAVE_MPI)
  } else {
    mpi_case(16, 64, CPMDC_PWSCALE_DIST_PW);
    mpi_case(16, 64, CPMDC_PWSCALE_DIST_ATOMS);
    mpi_case(32, 256, CPMDC_PWSCALE_DIST_PW);
    mpi_case(32, 128, CPMDC_PWSCALE_DIST_ATOMS);
  }
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  g_fail = any;
  MPI_Finalize();
  if (rank != 0)
    return g_fail ? 1 : 0;
#endif
  if (g_fail) {
    fprintf(stderr, "pwscale agreement failed\n");
    return 1;
  }
  printf("pwscale agreement ok\n");
  return 0;
}
