/* SPDX-License-Identifier: MIT */
#ifndef CPMDC_PWSCALE_H
#define CPMDC_PWSCALE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc_pwscale.h
 * @brief Plane-wave structure-factor scaffold.
 *
 * This is not the embed force call and not the reference potential.
 * The ordinary build leaves this scaffold out. Turning it on uses the CPU
 * path unless a separate flag compiles the CUDA backend.
 *
 * Parallel and serial energies and forces agree when every absolute gap is
 * at most CPMDC_PWSCALE_RTOL * max(1, |got|, |reference|).
 * The scale follows the larger magnitude, then the floor of 1.
 * Visit counts match exactly.
 * The ordered reduction is what the tolerance covers. It is not a bitwise claim.
 * A rank call is bitwise the same on every rank and from run to run for one
 * rank count, because every rank sums the same partials in rank order.
 */

/** @brief Status and distribution codes for the scaffold. */
#define CPMDC_PWSCALE_OK 0
/** Argument or collective failure. */
#define CPMDC_PWSCALE_ERR_ARG -1
/** This build has no MPI. */
#define CPMDC_PWSCALE_ERR_NOMPI -2
/** This build has no CUDA. */
#define CPMDC_PWSCALE_ERR_NOCUDA -3
/** A CUDA call failed. */
#define CPMDC_PWSCALE_ERR_CUDA -4
/** This build has no FFTW. */
#define CPMDC_PWSCALE_ERR_NOFFTW -5

/** Split the plane-wave loop across workers. */
#define CPMDC_PWSCALE_DIST_PW 0
/** Split the atom loop across workers. */
#define CPMDC_PWSCALE_DIST_ATOMS 1

/**
 * Relative tolerance for energy and force agreement.
 * A gap is accepted when it is at most this times max(1, |got|, |reference|).
 */
#define CPMDC_PWSCALE_RTOL 1e-8

/** @brief One scaffold problem: atom count, plane waves, and the Gaussian width. */
typedef struct cpmdc_pwscale_problem {
  /** Number of atoms. */
  int natoms;
  /** Number of plane waves. */
  int npw;
  /** Cubic cell edge, Bohr. */
  double cell_bohr;
  /** Gaussian exponent. */
  double alpha;
  /** Gaussian width. */
  double sigma;
} cpmdc_pwscale_problem;

/** @brief Timers and visit counts from one scaffold evaluation. */
typedef struct cpmdc_pwscale_stats {
  /** Energy, Hartree. */
  double energy_h;
  /** Wall time of the evaluated region, seconds. */
  double seconds;
  /** Host-to-device time, seconds. Zero on the CPU path. */
  double h2d_seconds;
  /** Device-to-host time, seconds. Zero on the CPU path. */
  double d2h_seconds;
  /** Plane waves visited. */
  long long g_visited;
  /** Atom visits. */
  long long atom_visits;
  /** Bytes allocated for the evaluation. */
  long long alloc_bytes;
  /** OpenMP team that ran, from omp_get_num_threads, else 1. */
  int workers;
  /** Communicator size, or 1 when the call is not MPI. */
  int mpi_ranks;
} cpmdc_pwscale_stats;

/**
 * @brief Fill a deterministic geometry for @p prob.
 *
 * A null pointer returns without writing. Positions are Bohr.
 */
void cpmdc_pwscale_fill_geometry(const cpmdc_pwscale_problem *prob,
                                 double *positions_bohr,
                                 int *atomic_numbers);

/**
 * @brief Return 1 when @p got agrees with @p ref under CPMDC_PWSCALE_RTOL.
 */
int cpmdc_pwscale_within(double got, double ref);

/**
 * @brief Compare one finished sample with the serial reference.
 *
 * Writes this sample's gaps and whether the visit counts match.
 * Returns 1 only when the energy, every force, and both counts agree.
 * A null buffer or a non-positive @p n3 returns 0 and writes nothing.
 */
int cpmdc_pwscale_sample_check(double got_energy, const double *got_forces,
                               double ref_energy, const double *ref_forces,
                               int n3, long long g_got, long long g_ref,
                               long long a_got, long long a_ref,
                               double *max_abs_dE, double *max_abs_dF,
                               int *counts_match);

/**
 * @brief Evaluate the scaffold on one rank.
 *
 * @p dist selects the split: plane waves, or atoms.
 * @p nthreads is the OpenMP team size. One thread is the serial loop.
 * Wall time covers scratch allocation, the arithmetic, and the ordered sum.
 * A null buffer or an unknown @p dist returns CPMDC_PWSCALE_ERR_ARG.
 */
int cpmdc_pwscale_eval(const cpmdc_pwscale_problem *prob,
                       const double *positions_bohr,
                       const int *atomic_numbers, double *forces_h_bohr,
                       int nthreads, int dist, cpmdc_pwscale_stats *stats);

/**
 * @brief Evaluate the scaffold on every rank of one communicator.
 *
 * Every rank passes the same problem. Every rank receives every rank's
 * partial and sums them in rank order. @p nthreads is the team inside a rank.
 * @p comm NULL selects MPI_COMM_WORLD. A build without MPI returns
 * CPMDC_PWSCALE_ERR_NOMPI.
 * seconds runs from this rank's entry to the arrival of every partial.
 * Each rank holds (ranks + 1) * (3 * natoms + 4) doubles for the sum.
 */
int cpmdc_pwscale_eval_mpi(const cpmdc_pwscale_problem *prob,
                           const double *positions_bohr,
                           const int *atomic_numbers, double *forces_h_bohr,
                           int nthreads, int dist, cpmdc_pwscale_stats *stats,
                           const void *comm, size_t comm_bytes);

/**
 * @brief A problem bound to a communicator for repeated rank force calls.
 *
 * cpmdc_pwscale_eval_mpi binds, evaluates once, and frees on every call.
 */
typedef struct cpmdc_pwscale_mpi_bind cpmdc_pwscale_mpi_bind;

/**
 * @brief Bind a problem to a communicator.
 *
 * Collective over @p comm. Every rank passes the same problem and @p dist.
 * One collective checks that agreement and every rank's buffers, so the
 * bound calls repeat neither check. @p comm NULL selects MPI_COMM_WORLD, and
 * the communicator must outlive the binding. On failure every rank
 * returns the same code and @p out is NULL. A build without MPI returns
 * CPMDC_PWSCALE_ERR_NOMPI.
 */
int cpmdc_pwscale_mpi_bind_create(const cpmdc_pwscale_problem *prob, int dist,
                                  const void *comm, size_t comm_bytes,
                                  cpmdc_pwscale_mpi_bind **out);

/**
 * @brief One bound force call.
 *
 * Collective over the bound communicator, with the result of
 * cpmdc_pwscale_eval_mpi for the same arguments. The call runs one
 * collective, an Allgather of each rank's energy, forces, visit counts, and
 * return code. A failure on one rank, including a missing array, returns
 * an error on every rank.
 */
int cpmdc_pwscale_eval_mpi_bound(cpmdc_pwscale_mpi_bind *bind,
                                 const double *positions_bohr,
                                 const int *atomic_numbers,
                                 double *forces_h_bohr, int nthreads,
                                 cpmdc_pwscale_stats *stats);

/** @brief Release a binding. Local. Leaves the communicator alone. */
void cpmdc_pwscale_mpi_bind_free(cpmdc_pwscale_mpi_bind *bind);

/** @brief Return 1 when this build compiled the CUDA scaffold. */
int cpmdc_pwscale_cuda_compiled(void);

/**
 * @brief Write a device name when one is present.
 *
 * Returns 1 when a device exists.
 */
int cpmdc_pwscale_cuda_device(char *name, int namelen);

/**
 * @brief Evaluate the scaffold on the CUDA device.
 *
 * Total seconds include host-to-device, device-to-host, and synchronization.
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
