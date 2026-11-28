/* SPDX-License-Identifier: MIT */
#ifndef CPMDC_DAMP_H
#define CPMDC_DAMP_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc_damp.h
 * @brief Opt-in damped electronic dynamics.
 *
 * Nothing in a session evaluation calls these functions, and libcpmdc does
 * not export them. A host that already has fictitious orbital coefficients
 * and the force on those coefficients calls them between force evaluations.
 * Omitting the call leaves the coefficients unchanged.
 *
 * The step is the damped Verlet update in Quantum ESPRESSO CPV: wave_verlet
 * in CPV/src/wave_base.f90, with the coefficients set in CPV/src/runcp.f90
 * from fccc in CPV/src/cpr.f90. input.f90 copies electron_damping into frice
 * when electron_dynamics is damp, and sets frice to 0 for verlet. Friction 0
 * is that Verlet update. This step does not compute the force, does not
 * orthonormalize, and does not clear a gamma-point imaginary part.
 *
 * Units are the atomic units of that update. Friction is the dimensionless
 * electron_damping value. dt2_over_emass is dt squared over the fictitious
 * electronic mass. ema0bg is the per-coefficient mass preconditioner, or
 * all ones when the pointer is NULL. The force is the grad argument of
 * wave_verlet. A return of 0 writes every output. A return of -1 writes
 * nothing. cpmdc_damp_grease does that when the pointer is null, a factor
 * is not finite, the input friction is not greater than -1, or the
 * product of that friction and grease is not finite or is not greater
 * than -1. The functions do not call cpmdc_last_error. They keep no
 * process state.
 */

/**
 * Verlet weights for one friction.
 * fccc and verl3 are 1 / (1 + friction). verl1 is 2 * fccc.
 * verl2 is 1 - verl1.
 */
int cpmdc_damp_coefficients(double friction, double *fccc, double *verl1,
                            double *verl2, double *verl3);

/**
 * Write c(t+dt) into the cm arrays.
 * c0 is c(t) and is not modified. cm enters as c(t-dt).
 * n is the number of complex coefficients. ema0bg has n entries, or is
 * NULL for ones. A friction of 0 is undamped Verlet.
 */
int cpmdc_damp_step(int n, const double *c0_re, const double *c0_im,
                    double *cm_re, double *cm_im, const double *force_re,
                    const double *force_im, const double *ema0bg,
                    double friction, double dt2_over_emass);

/**
 * Exchange c0 and cm.
 * After cpmdc_damp_step, this is the dswap that leaves c0 at c(t+dt)
 * and cm at c(t).
 */
int cpmdc_damp_swap(int n, double *c0_re, double *c0_im, double *cm_re,
                    double *cm_im);

/**
 * Multiply friction by grease.
 * cpr.f90 does this after the electronic step. grease 1 leaves a finite
 * friction that is greater than -1 unchanged.
 * Returns -1 and leaves *friction unchanged when friction is null, when
 * *friction or grease is not finite, when *friction is not greater than
 * -1, or when the product *friction * grease is not finite or is not
 * greater than -1. A return of 0 stores that product.
 */
int cpmdc_damp_grease(double *friction, double grease);

/**
 * One electronic step, then the swap.
 * On success c0 is c(t+dt) and cm is c(t). On failure both are unchanged.
 */
int cpmdc_damp_advance(int n, double *c0_re, double *c0_im, double *cm_re,
                       double *cm_im, const double *force_re,
                       const double *force_im, const double *ema0bg,
                       double friction, double dt2_over_emass);

#ifdef __cplusplus
}
#endif

#endif
