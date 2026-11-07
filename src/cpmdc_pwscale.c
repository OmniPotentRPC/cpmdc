/* SPDX-License-Identifier: MIT */
#define _POSIX_C_SOURCE 200809L

#include "cpmdc_pwscale.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(CPMDC_PWSCALE_HAVE_MPI)
#include <mpi.h>
#endif

static const double PW_PI = 3.14159265358979323846;

static double pw_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int grid_n(int npw) {
  int n = 1;
  while ((long long)n * n * n < (long long)npw)
    n += 2;
  return n;
}

static void miller_at(int k, int n, int *ix, int *iy, int *iz) {
  int half = n / 2;
  int x = k % n;
  int y = (k / n) % n;
  int z = k / (n * n);
  *ix = x - half;
  *iy = y - half;
  *iz = z - half;
}

static void split_range(int lo, int hi, int index, int pieces, int *begin,
                        int *end) {
  int n = hi - lo;
  if (pieces < 1)
    pieces = 1;
  *begin = lo + (int)(((long long)index * n) / pieces);
  *end = lo + (int)(((long long)(index + 1) * n) / pieces);
}

static int problem_ok(const cpmdc_pwscale_problem *prob, const double *positions,
                      const int *atomic_numbers, const double *forces) {
  if (!prob || !positions || !atomic_numbers || !forces)
    return 0;
  if (prob->natoms < 1 || prob->npw < 1)
    return 0;
  if (prob->natoms > 1000000 || prob->npw > 100000000)
    return 0;
  if (!(prob->cell_bohr > 0.0) || !(prob->alpha > 0.0) || !(prob->sigma > 0.0))
    return 0;
  if (prob->natoms > 2000000000 / 3)
    return 0;
  return 1;
}

/* Adds into energy and forces. G = 0 is visited and adds nothing. */
static void accum_range(const cpmdc_pwscale_problem *prob, const double *pos,
                        const int *atomic_numbers, int k0, int k1, int a0,
                        int a1, double *energy, double *forces,
                        long long *g_visited, long long *atom_visits) {
  const double scale = (2.0 * PW_PI) / prob->cell_bohr;
  const int n = grid_n(prob->npw);
  const double alpha = prob->alpha;
  const double sigma = prob->sigma;
  double e = 0.0;
  long long gv = 0;
  long long av = 0;
  int k;

  for (k = k0; k < k1; ++k) {
    int ix, iy, iz, a;
    double gx, gy, gz, g2, weight;
    miller_at(k, n, &ix, &iy, &iz);
    gx = scale * (double)ix;
    gy = scale * (double)iy;
    gz = scale * (double)iz;
    g2 = (gx * gx + gy * gy) + gz * gz;
    gv += 1;
    if (g2 == 0.0)
      continue;
    weight = ((4.0 * PW_PI) / g2) * exp(-g2 / (4.0 * alpha)) *
              exp(-g2 / (2.0 * sigma * sigma));
    for (a = a0; a < a1; ++a) {
      double phase =
          (gx * pos[3 * a] + gy * pos[3 * a + 1]) + gz * pos[3 * a + 2];
      double c = cos(phase);
      double s = sin(phase);
      double wz = weight * (double)atomic_numbers[a];
      e += wz * c;
      forces[3 * a] += wz * s * gx;
      forces[3 * a + 1] += wz * s * gy;
      forces[3 * a + 2] += wz * s * gz;
      av += 1;
    }
  }
  *energy += e;
  *g_visited += gv;
  *atom_visits += av;
}

/* Threads split either the plane-wave box or the atom box.
 * Plane-wave splits keep a private force vector per thread and sum in
 * thread order. Atom splits write disjoint force slots.
 * Scratch is the private force vectors, not a plane-wave by atom matrix.
 */
static int team_box(const cpmdc_pwscale_problem *prob, const double *pos,
                    const int *atomic_numbers, int k0, int k1, int a0, int a1,
                    int nthreads, int split_pw, double *energy, double *forces,
                    long long *g_visited, long long *atom_visits,
                    long long *alloc_bytes, int *workers) {
  int nt, t, n3, used;
  double *partial_e;
  long long *partial_g;
  long long *partial_a;
  double **partial_f = NULL;

  if (nthreads < 1)
    nthreads = 1;
  if (nthreads > 4096)
    return CPMDC_PWSCALE_ERR_ARG;
  nt = nthreads;
  n3 = prob->natoms * 3;
  *energy = 0.0;
  *g_visited = 0;
  *atom_visits = 0;
  *workers = 1;

  partial_e = (double *)calloc((size_t)nt, sizeof(double));
  partial_g = (long long *)calloc((size_t)nt, sizeof(long long));
  partial_a = (long long *)calloc((size_t)nt, sizeof(long long));
  if (!partial_e || !partial_g || !partial_a) {
    free(partial_e);
    free(partial_g);
    free(partial_a);
    return CPMDC_PWSCALE_ERR_ARG;
  }
  *alloc_bytes += (long long)nt * (long long)sizeof(double);
  *alloc_bytes += 2 * (long long)nt * (long long)sizeof(long long);

  if (split_pw) {
    partial_f = (double **)calloc((size_t)nt, sizeof(double *));
    if (!partial_f) {
      free(partial_e);
      free(partial_g);
      free(partial_a);
      return CPMDC_PWSCALE_ERR_ARG;
    }
    *alloc_bytes += (long long)nt * (long long)sizeof(double *);
    for (t = 0; t < nt; ++t) {
      partial_f[t] = (double *)calloc((size_t)n3, sizeof(double));
      if (!partial_f[t]) {
        int j;
        for (j = 0; j < t; ++j)
          free(partial_f[j]);
        free(partial_f);
        free(partial_e);
        free(partial_g);
        free(partial_a);
        return CPMDC_PWSCALE_ERR_ARG;
      }
      *alloc_bytes += (long long)n3 * (long long)sizeof(double);
    }
  }

  used = 1;
#if defined(_OPENMP)
  omp_set_dynamic(0);
#pragma omp parallel num_threads(nt)
  {
    int tid = omp_get_thread_num();
    int nuse = omp_get_num_threads();
    int b, e;
    if (tid == 0)
      used = nuse;
    if (split_pw) {
      split_range(k0, k1, tid, nuse, &b, &e);
      accum_range(prob, pos, atomic_numbers, b, e, a0, a1, &partial_e[tid],
                  partial_f[tid], &partial_g[tid], &partial_a[tid]);
    } else {
      split_range(a0, a1, tid, nuse, &b, &e);
      accum_range(prob, pos, atomic_numbers, k0, k1, b, e, &partial_e[tid],
                  forces, &partial_g[tid], &partial_a[tid]);
    }
  }
#else
  if (split_pw) {
    accum_range(prob, pos, atomic_numbers, k0, k1, a0, a1, &partial_e[0],
                partial_f[0], &partial_g[0], &partial_a[0]);
  } else {
    accum_range(prob, pos, atomic_numbers, k0, k1, a0, a1, &partial_e[0],
                forces, &partial_g[0], &partial_a[0]);
  }
  used = 1;
#endif

  for (t = 0; t < nt; ++t)
    *energy += partial_e[t];
  if (split_pw) {
    for (t = 0; t < nt; ++t) {
      int i;
      for (i = 0; i < n3; ++i)
        forces[i] += partial_f[t][i];
      *g_visited += partial_g[t];
      *atom_visits += partial_a[t];
      free(partial_f[t]);
    }
    free(partial_f);
  } else {
    /* Every thread walks the same plane waves, so the count is one copy. */
    *g_visited = partial_g[0];
    for (t = 0; t < nt; ++t)
      *atom_visits += partial_a[t];
  }
  free(partial_e);
  free(partial_g);
  free(partial_a);
  *workers = used;
  return CPMDC_PWSCALE_OK;
}

void cpmdc_pwscale_fill_geometry(const cpmdc_pwscale_problem *prob,
                                 double *positions_bohr, int *atomic_numbers) {
  int i;
  if (!prob || !positions_bohr || !atomic_numbers)
    return;
  for (i = 0; i < prob->natoms; ++i) {
    double s = (double)(i + 1);
    positions_bohr[3 * i] = fmod(0.11 * s, 0.73) * prob->cell_bohr;
    positions_bohr[3 * i + 1] = fmod(0.17 * s, 0.71) * prob->cell_bohr;
    positions_bohr[3 * i + 2] = fmod(0.23 * s, 0.69) * prob->cell_bohr;
    atomic_numbers[i] = 1 + (i % 8);
  }
}

int cpmdc_pwscale_within(double got, double ref) {
  double gap = fabs(got - ref);
  double scale = fabs(ref);
  if (fabs(got) > scale)
    scale = fabs(got);
  if (scale < 1.0)
    scale = 1.0;
  return gap <= CPMDC_PWSCALE_RTOL * scale;
}

int cpmdc_pwscale_sample_check(double got_energy, const double *got_forces,
                               double ref_energy, const double *ref_forces,
                               int n3, long long g_got, long long g_ref,
                               long long a_got, long long a_ref,
                               double *max_abs_dE, double *max_abs_dF,
                               int *counts_match) {
  int i, forces_ok = 1, counts, energy_ok;
  double de, df = 0.0;
  if (!got_forces || !ref_forces || n3 < 1 || !max_abs_dE || !max_abs_dF ||
      !counts_match)
    return 0;
  de = fabs(got_energy - ref_energy);
  for (i = 0; i < n3; ++i) {
    double gap = fabs(got_forces[i] - ref_forces[i]);
    if (gap > df)
      df = gap;
    if (!cpmdc_pwscale_within(got_forces[i], ref_forces[i]))
      forces_ok = 0;
  }
  counts = g_got == g_ref && a_got == a_ref;
  energy_ok = cpmdc_pwscale_within(got_energy, ref_energy);
  *max_abs_dE = de;
  *max_abs_dF = df;
  *counts_match = counts;
  return energy_ok && forces_ok && counts;
}

int cpmdc_pwscale_eval(const cpmdc_pwscale_problem *prob,
                       const double *positions_bohr, const int *atomic_numbers,
                       double *forces_h_bohr, int nthreads, int dist,
                       cpmdc_pwscale_stats *stats) {
  int rc, workers = 1;
  double t0, t1, energy = 0.0;
  long long g_visited = 0, atom_visits = 0, alloc_bytes = 0;

  if (stats)
    memset(stats, 0, sizeof(*stats));
  if (!problem_ok(prob, positions_bohr, atomic_numbers, forces_h_bohr))
    return CPMDC_PWSCALE_ERR_ARG;
  if (dist != CPMDC_PWSCALE_DIST_PW && dist != CPMDC_PWSCALE_DIST_ATOMS)
    return CPMDC_PWSCALE_ERR_ARG;

  memset(forces_h_bohr, 0, (size_t)prob->natoms * 3u * sizeof(double));
  t0 = pw_now();
  rc = team_box(prob, positions_bohr, atomic_numbers, 0, prob->npw, 0,
                prob->natoms, nthreads, dist == CPMDC_PWSCALE_DIST_PW, &energy,
                forces_h_bohr, &g_visited, &atom_visits, &alloc_bytes, &workers);
  t1 = pw_now();
  if (rc != CPMDC_PWSCALE_OK)
    return rc;
  if (stats) {
    stats->energy_h = energy;
    stats->seconds = t1 - t0;
    stats->g_visited = g_visited;
    stats->atom_visits = atom_visits;
    stats->alloc_bytes = alloc_bytes;
    stats->workers = workers;
    stats->mpi_ranks = 1;
  }
  return CPMDC_PWSCALE_OK;
}

#if defined(CPMDC_PWSCALE_HAVE_MPI)

/* Test hook. Unset in ordinary runs. "gather" drops rank 0's buffers.
 * "team" fails the last rank's team. Every rank still joins the allreduce. */
static int mpi_fault(const char *name) {
  const char *env = getenv("CPMDC_PWSCALE_MPI_FAULT");
  return env && strcmp(env, name) == 0;
}

int cpmdc_pwscale_eval_mpi(const cpmdc_pwscale_problem *prob,
                           const double *positions_bohr,
                           const int *atomic_numbers, double *forces_h_bohr,
                           int nthreads, int dist, cpmdc_pwscale_stats *stats,
                           const void *comm, size_t comm_bytes) {
  MPI_Comm use;
  int rank = 0, size = 1, rc, b, e, i, workers = 1;
  int key[2], minv[2], maxv[2];
  int natoms, npw, n3;
  double t0, t1, local_e = 0.0, energy = 0.0;
  double *local_f = NULL;
  double *gathered_e = NULL;
  double *gathered_f = NULL;
  long long local_g = 0, local_a = 0, g_visited = 0, atom_visits = 0;
  long long alloc_bytes = 0;

  if (stats)
    memset(stats, 0, sizeof(*stats));
  if (!problem_ok(prob, positions_bohr, atomic_numbers, forces_h_bohr))
    return CPMDC_PWSCALE_ERR_ARG;
  if (dist != CPMDC_PWSCALE_DIST_PW && dist != CPMDC_PWSCALE_DIST_ATOMS)
    return CPMDC_PWSCALE_ERR_ARG;

  if (comm) {
    if (comm_bytes != sizeof(MPI_Comm))
      return CPMDC_PWSCALE_ERR_ARG;
    memcpy(&use, comm, sizeof(use));
  } else {
    use = MPI_COMM_WORLD;
  }
  if (MPI_Comm_rank(use, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(use, &size) != MPI_SUCCESS || size < 1)
    return CPMDC_PWSCALE_ERR_ARG;

  natoms = prob->natoms;
  npw = prob->npw;
  n3 = natoms * 3;
  key[0] = npw;
  key[1] = natoms;
  if (MPI_Allreduce(key, minv, 2, MPI_INT, MPI_MIN, use) != MPI_SUCCESS ||
      MPI_Allreduce(key, maxv, 2, MPI_INT, MPI_MAX, use) != MPI_SUCCESS)
    return CPMDC_PWSCALE_ERR_ARG;
  if (minv[0] != maxv[0] || minv[1] != maxv[1])
    return CPMDC_PWSCALE_ERR_ARG;

  local_f = (double *)calloc((size_t)n3, sizeof(double));
  if (local_f)
    alloc_bytes += (long long)n3 * (long long)sizeof(double);
  if (rank == 0 && !mpi_fault("gather")) {
    gathered_e = (double *)calloc((size_t)size, sizeof(double));
    gathered_f = (double *)calloc((size_t)size * (size_t)n3, sizeof(double));
    if (gathered_e && gathered_f) {
      alloc_bytes += (long long)size * (long long)sizeof(double);
      alloc_bytes += (long long)size * (long long)n3 * (long long)sizeof(double);
    }
  }
  {
    int local_ok = local_f != NULL && (rank != 0 || (gathered_e && gathered_f));
    int all_ok = 0;
    if (MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, use) != MPI_SUCCESS ||
        !all_ok) {
      free(local_f);
      free(gathered_e);
      free(gathered_f);
      return CPMDC_PWSCALE_ERR_ARG;
    }
  }

  if (MPI_Barrier(use) != MPI_SUCCESS) {
    free(local_f);
    free(gathered_e);
    free(gathered_f);
    return CPMDC_PWSCALE_ERR_ARG;
  }
  t0 = pw_now();

  if (dist == CPMDC_PWSCALE_DIST_PW) {
    split_range(0, npw, rank, size, &b, &e);
    rc = team_box(prob, positions_bohr, atomic_numbers, b, e, 0, natoms,
                  nthreads, 1, &local_e, local_f, &local_g, &local_a,
                  &alloc_bytes, &workers);
  } else {
    split_range(0, natoms, rank, size, &b, &e);
    rc = team_box(prob, positions_bohr, atomic_numbers, 0, npw, b, e, nthreads,
                  0, &local_e, local_f, &local_g, &local_a, &alloc_bytes,
                  &workers);
  }
  if (mpi_fault("team") && size > 1 && rank == size - 1)
    rc = CPMDC_PWSCALE_ERR_ARG;
  {
    int rc_all = rc;
    if (MPI_Allreduce(&rc, &rc_all, 1, MPI_INT, MPI_MIN, use) != MPI_SUCCESS) {
      free(local_f);
      free(gathered_e);
      free(gathered_f);
      return CPMDC_PWSCALE_ERR_ARG;
    }
    if (rc_all != CPMDC_PWSCALE_OK) {
      free(local_f);
      free(gathered_e);
      free(gathered_f);
      return rc_all;
    }
  }

  if (MPI_Gather(&local_e, 1, MPI_DOUBLE, gathered_e, 1, MPI_DOUBLE, 0, use) !=
          MPI_SUCCESS ||
      MPI_Gather(local_f, n3, MPI_DOUBLE, gathered_f, n3, MPI_DOUBLE, 0, use) !=
          MPI_SUCCESS) {
    free(local_f);
    free(gathered_e);
    free(gathered_f);
    return CPMDC_PWSCALE_ERR_ARG;
  }

  memset(forces_h_bohr, 0, (size_t)n3 * sizeof(double));
  if (rank == 0) {
    int r;
    for (r = 0; r < size; ++r) {
      energy += gathered_e[r];
      for (i = 0; i < n3; ++i)
        forces_h_bohr[i] += gathered_f[(size_t)r * (size_t)n3 + (size_t)i];
    }
  }
  if (MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, use) != MPI_SUCCESS ||
      MPI_Bcast(forces_h_bohr, n3, MPI_DOUBLE, 0, use) != MPI_SUCCESS) {
    free(local_f);
    free(gathered_e);
    free(gathered_f);
    return CPMDC_PWSCALE_ERR_ARG;
  }

  if (dist == CPMDC_PWSCALE_DIST_PW) {
    if (MPI_Allreduce(&local_g, &g_visited, 1, MPI_LONG_LONG, MPI_SUM, use) !=
            MPI_SUCCESS ||
        MPI_Allreduce(&local_a, &atom_visits, 1, MPI_LONG_LONG, MPI_SUM, use) !=
            MPI_SUCCESS) {
      free(local_f);
      free(gathered_e);
      free(gathered_f);
      return CPMDC_PWSCALE_ERR_ARG;
    }
  } else {
    g_visited = npw;
    if (MPI_Allreduce(&local_a, &atom_visits, 1, MPI_LONG_LONG, MPI_SUM, use) !=
        MPI_SUCCESS) {
      free(local_f);
      free(gathered_e);
      free(gathered_f);
      return CPMDC_PWSCALE_ERR_ARG;
    }
  }

  if (MPI_Barrier(use) != MPI_SUCCESS) {
    free(local_f);
    free(gathered_e);
    free(gathered_f);
    return CPMDC_PWSCALE_ERR_ARG;
  }
  t1 = pw_now();
  free(local_f);
  free(gathered_e);
  free(gathered_f);

  if (stats) {
    stats->energy_h = energy;
    stats->seconds = t1 - t0;
    stats->g_visited = g_visited;
    stats->atom_visits = atom_visits;
    stats->alloc_bytes = alloc_bytes;
    stats->workers = workers;
    stats->mpi_ranks = size;
  }
  return CPMDC_PWSCALE_OK;
}

#else

int cpmdc_pwscale_eval_mpi(const cpmdc_pwscale_problem *prob,
                           const double *positions_bohr,
                           const int *atomic_numbers, double *forces_h_bohr,
                           int nthreads, int dist, cpmdc_pwscale_stats *stats,
                           const void *comm, size_t comm_bytes) {
  (void)prob;
  (void)positions_bohr;
  (void)atomic_numbers;
  (void)forces_h_bohr;
  (void)nthreads;
  (void)dist;
  (void)comm;
  (void)comm_bytes;
  if (stats)
    memset(stats, 0, sizeof(*stats));
  return CPMDC_PWSCALE_ERR_NOMPI;
}

#endif
