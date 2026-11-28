#ifndef CPMDC_CSVR_H
#define CPMDC_CSVR_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Caller-owned generator for one CSVR trajectory.
 *
 * The four state words are a xoshiro256** generator. spare_normal and
 * spare_ready hold one unused normal from the Box-Muller pair. Seed the
 * object before a scale call. An all-zero object is rejected. The host
 * owns the memory and may copy the bytes to replay a trajectory.
 */
typedef struct cpmdc_csvr_rng {
    uint64_t state[4];
    double spare_normal;
    int spare_ready;
} cpmdc_csvr_rng;

/**
 * Fill rng from seed and clear the stored normal.
 *
 * Any uint64_t seed is accepted, including zero. Returns 0, or -1 when
 * rng is NULL. A failure writes nothing.
 */
int cpmdc_csvr_rng_seed(cpmdc_csvr_rng *rng, uint64_t seed);

/**
 * Positive velocity scale for one Bussi CSVR step, and the thermostat work.
 *
 * kinetic is the kinetic energy of the degrees being thermostatted, in
 * any energy unit. ndof is how many of those quadratic degrees the host
 * wants in the ensemble. temperature is k_B T in that same energy unit,
 * not a Kelvin count. time_step and time_constant share one time unit.
 * A time_constant of positive infinity, or a zero time_step, is NVE:
 * the scale is 1 and the work is 0. The host multiplies each velocity
 * component by the scale. This library does not see the velocities.
 *
 * On success the scale is non-negative and the work is
 * kinetic * (1 - scale * scale), the energy removed from the ions.
 * Adding the work onto the Hamiltonian is the effective energy of
 * Bussi, Donadio, and Parrinello. Returns 0 and writes scale, work, and,
 * when a draw was required, rng. Returns -1 for a NULL pointer, an
 * unseeded rng, a non-positive kinetic energy, a negative temperature,
 * a negative time step, a time constant that is not positive, a degree
 * count below 1, or a non-finite result. A failure writes nothing and
 * does not advance rng.
 */
int cpmdc_csvr_scale(double kinetic, int ndof, double temperature,
                     double time_step, double time_constant,
                     cpmdc_csvr_rng *rng, double *scale, double *work);

#ifdef __cplusplus
}
#endif

#endif
