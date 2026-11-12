/* SPDX-License-Identifier: MIT */
#ifndef CPMDC_CP_H
#define CPMDC_CP_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Opt-in Car-Parrinello velocity Verlet.
 *
 * A force evaluation still optimises the wavefunction. Nothing here runs
 * unless the caller asks for a step. Coefficient matrices are column-major,
 * one state per column, in atomic units. Positions are Bohr, energies are
 * Hartree, and ionic masses are electron masses.
 *
 * CPMDC_CP_ORTHO_LOWDIN diagonalises the overlap, which is the Löwdin step
 * OpenCPMD takes when LOWDIN is set. CPMDC_CP_ORTHO_RATTLE solves the same
 * orthonormality constraints iteratively: SHAKE on the coefficients and
 * RATTLE on their velocities.
 */
#define CPMDC_CP_ORTHO_LOWDIN 0
#define CPMDC_CP_ORTHO_RATTLE 1

#define CPMDC_CP_WATER_ATOMS 3
#define CPMDC_CP_WATER_BASIS 3
#define CPMDC_CP_WATER_STATES 2

typedef struct cpmdc_cp_system {
    int n_basis;
    int n_state;
    int n_atom;
    /** Fictitious electronic mass, atomic units. */
    double mu;
    /** Ionic masses, electron masses, length n_atom. */
    const double *mass;
} cpmdc_cp_system;

/**
 * Partial derivatives at fixed coefficients.
 * f_elec has the layout of C and is minus the derivative of the energy.
 * F_ion is minus the derivative with respect to R, length 3 * n_atom.
 * Returns 0 on success.
 */
typedef int (*cpmdc_cp_force_fn)(void *ctx, const double *R, const double *C,
                                 double *energy, double *F_ion, double *f_elec);

typedef struct cpmdc_cp_step_opts {
    /** Time step in atomic units. */
    double dt;
    /** CPMDC_CP_ORTHO_LOWDIN or CPMDC_CP_ORTHO_RATTLE. */
    int ortho;
    /** Largest allowed absolute overlap error after SHAKE. */
    double shake_tol;
    /** Iteration cap for SHAKE. Ignored by the Löwdin path. */
    int shake_max;
} cpmdc_cp_step_opts;

typedef struct cpmdc_cp_step_info {
    /**
     * SHAKE iterations on the RATTLE path, or Jacobi sweeps of the overlap
     * on the Löwdin path.
     */
    int iterations;
    /** A Car-Parrinello step does not run an SCF. This is 0. */
    int scf_iterations;
    double energy;
    double ekin_ion;
    /** Fictitious electronic kinetic energy. */
    double ekinc;
    /** energy + ekin_ion + ekinc. */
    double econs;
    /** Largest absolute entry of C^T C - I after the step. */
    double ortho_err;
} cpmdc_cp_step_info;

typedef struct cpmdc_cp_work cpmdc_cp_work;

cpmdc_cp_work *cpmdc_cp_work_new(int n_basis, int n_state, int n_atom);
void cpmdc_cp_work_free(cpmdc_cp_work *work);

/** Symmetric Löwdin orthonormalisation. C_out may alias C. Returns 0 on success. */
int cpmdc_cp_lowdin(int n_basis, int n_state, const double *C, double *C_out,
                    int *sweeps, double *err);

/**
 * Iterative SHAKE. C_new = C_pred + C_old * lambda, with lambda symmetric,
 * and C_new orthonormal. C_new must not alias C_old. Returns 0 on success.
 */
int cpmdc_cp_shake(int n_basis, int n_state, const double *C_old,
                   const double *C_pred, double *C_new, double tol, int max_iter,
                   int *iterations, double *err);

/** Remove the part of U that changes the orthonormality of C. Returns 0 on success. */
int cpmdc_cp_rattle_velocities(int n_basis, int n_state, const double *C, double *U);

int cpmdc_cp_conserved(const cpmdc_cp_system *sys, cpmdc_cp_force_fn force, void *ctx,
                       const double *R, const double *V, const double *C, const double *U,
                       cpmdc_cp_step_info *info);

int cpmdc_cp_step(const cpmdc_cp_system *sys, const cpmdc_cp_step_opts *opt,
                  cpmdc_cp_force_fn force, void *ctx, cpmdc_cp_work *work, double *R,
                  double *V, double *C, double *U, cpmdc_cp_step_info *info);

/** Reference water geometry from tests/force_input_water.capnp.txt, in Bohr. */
int cpmdc_cp_water_positions(double R[9]);
int cpmdc_cp_water_masses(double M[3]);
int cpmdc_cp_water_hamiltonian(const double *R, double H[9]);
int cpmdc_cp_water_ground_state(const double *R, double *C, double *energy);
int cpmdc_cp_water_band_energy(const double *R, double *energy);
int cpmdc_cp_water_force(void *ctx, const double *R, const double *C, double *energy,
                         double *F_ion, double *f_elec);

#ifdef __cplusplus
}
#endif

#endif
