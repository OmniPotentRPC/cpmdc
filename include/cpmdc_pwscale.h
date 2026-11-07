/* SPDX-License-Identifier: MIT */
#ifndef CPMDC_PWSCALE_H
#define CPMDC_PWSCALE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Plane-wave structure-factor scaffold.
 * This is not the embed force call and not the reference potential.
 * The ordinary build leaves this scaffold out. Turning it on uses the CPU
 * path unless a separate flag compiles the CUDA backend.
 *
 * Parallel and serial energies and forces agree when every absolute gap is
 * at most CPMDC_PWSCALE_RTOL * max(1, |got|, |reference|).
 * The scale follows the larger magnitude, then the floor of 1.
 * Visit counts match exactly.
 * The ordered reduction is what the tolerance covers. It is not a bitwise claim.
 */

#define CPMDC_PWSCALE_OK 0
#define CPMDC_PWSCALE_ERR_ARG -1
#define CPMDC_PWSCALE_ERR_NOMPI -2
#define CPMDC_PWSCALE_ERR_NOCUDA -3
#define CPMDC_PWSCALE_ERR_CUDA -4
#define CPMDC_PWSCALE_ERR_NOFFTW -5

#define CPMDC_PWSCALE_DIST_PW 0
#define CPMDC_PWSCALE_DIST_ATOMS 1

#define CPMDC_PWSCALE_RTOL 1e-8

typedef struct cpmdc_pwscale_problem {
  int natoms;
  int npw;
  double cell_bohr;
  double alpha;
  double sigma;
} cpmdc_pwscale_problem;

typedef struct cpmdc_pwscale_stats {
  double energy_h;
  double seconds;
  double h2d_seconds;
  double d2h_seconds;
  long long g_visited;
  long long atom_visits;
  long long alloc_bytes;
  int workers;   /* OpenMP team that ran, from omp_get_num_threads, else 1 */
  int mpi_ranks; /* communicator size, or 1 when the call is not MPI */
} cpmdc_pwscale_stats;

void cpmdc_pwscale_fill_geometry(const cpmdc_pwscale_problem *prob,
                                 double *positions_bohr,
                                 int *atomic_numbers);

int cpmdc_pwscale_within(double got, double ref);

/* Compare one finished sample with the serial reference.
 * Writes this sample's gaps and whether the visit counts match.
 * Returns 1 only when the energy, every force, and both counts agree.
 */
int cpmdc_pwscale_sample_check(double got_energy, const double *got_forces,
                               double ref_energy, const double *ref_forces,
                               int n3, long long g_got, long long g_ref,
                               long long a_got, long long a_ref,
                               double *max_abs_dE, double *max_abs_dF,
                               int *counts_match);

/* dist selects the split: plane waves, or atoms.
 * nthreads is the OpenMP team size. One thread is the serial loop.
 * Wall time covers scratch allocation, the arithmetic, and the ordered sum.
 */
int cpmdc_pwscale_eval(const cpmdc_pwscale_problem *prob,
                       const double *positions_bohr,
                       const int *atomic_numbers, double *forces_h_bohr,
                       int nthreads, int dist, cpmdc_pwscale_stats *stats);

/* Every rank passes the same problem. The library sums partials in rank
 * order and broadcasts the parent result. nthreads is the team inside a rank.
 * comm NULL selects MPI_COMM_WORLD. A build without MPI returns ERR_NOMPI.
 */
int cpmdc_pwscale_eval_mpi(const cpmdc_pwscale_problem *prob,
                           const double *positions_bohr,
                           const int *atomic_numbers, double *forces_h_bohr,
                           int nthreads, int dist, cpmdc_pwscale_stats *stats,
                           const void *comm, size_t comm_bytes);

int cpmdc_pwscale_cuda_compiled(void);

/* Writes a device name when one is present. Returns 1 when a device exists. */
int cpmdc_pwscale_cuda_device(char *name, int namelen);

/* Total seconds include host-to-device, device-to-host, and synchronization.
 * A missing build or a missing device is an error and writes no timing.
 */
int cpmdc_pwscale_eval_cuda(const cpmdc_pwscale_problem *prob,
                            const double *positions_bohr,
                            const int *atomic_numbers, double *forces_h_bohr,
                            cpmdc_pwscale_stats *stats);

#ifdef __cplusplus
}
#endif

#endif
