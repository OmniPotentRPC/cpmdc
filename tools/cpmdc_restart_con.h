#ifndef CPMDC_RESTART_CON_H
#define CPMDC_RESTART_CON_H

#include "cpmdc_restart.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc_restart_con.h
 * @brief Optional CON sidecar for a CPMD RESTART image.
 *
 * The binary RESTART stays the wavefunction. The CON file is the geometry
 * a later step resumes. This bridge does not link readcon-db and does not
 * store coefficients, cutoffs, or electronic states in the CON file.
 *
 * Bohr values become angstroms by multiplying CPMDC_RESTART_ANGSTROM_PER_BOHR
 * in cpmdc_restart_bohr_to_angstrom. The inverse divides by that same macro.
 * The CON writer is create_writer_from_path_round_trip_c, canonical left off,
 * with no decimal precision argument.
 */

/** Multiply `bohr` by CPMDC_RESTART_ANGSTROM_PER_BOHR. */
double cpmdc_restart_bohr_to_angstrom(double bohr);

/** Divide `angstrom` by CPMDC_RESTART_ANGSTROM_PER_BOHR. */
double cpmdc_restart_angstrom_to_bohr(double angstrom);

/**
 * @brief Map CPMD celldm onto CON lengths (angstrom) and angles (degrees).
 *
 * celldm[0] is a in Bohr. celldm[1] and celldm[2] are b/a and c/a, and a
 * stored 0 means 1. celldm[3], celldm[4], and celldm[5] are cos(bc),
 * cos(ac), and cos(ab).
 *
 * @return 0 on success.
 */
int cpmdc_restart_celldm_to_con(const double celldm[6], double lengths[3],
                                double angles[3]);

/**
 * @brief Map CON lengths and angles back onto CPMD celldm.
 *
 * Lengths are divided by CPMDC_RESTART_ANGSTROM_PER_BOHR. The ratios are
 * b/a and c/a in Bohr.
 *
 * @return 0 on success.
 */
int cpmdc_restart_con_to_celldm(const double lengths[3], const double angles[3],
                                double celldm[6]);

/**
 * @brief Write one RoundTrip CON for this image.
 *
 * `atom_symbols` has one symbol per coordinate, in CPMD order. `atom_id`
 * is that index. Velocities are written when section 5 has records, with
 * units.time set to (amu*angstrom^2/eV)^0.5.
 *
 * @return 0 on success.
 */
int cpmdc_restart_con_export(const cpmdc_restart *file, const char *con_path,
                             const char *const *atom_symbols, int n_atoms,
                             char *err, size_t err_cap);

/**
 * @brief Patch coordinates, velocities, and the cell from a CON file.
 *
 * Positions and velocities are scattered by `atom_id` because `build`
 * groups atoms by symbol. The caller writes the image with
 * cpmdc_restart_write_path.
 *
 * @return 0 on success.
 */
int cpmdc_restart_con_import(cpmdc_restart *file, const char *con_path, char *err,
                             size_t err_cap);

#ifdef __cplusplus
}
#endif

#endif
