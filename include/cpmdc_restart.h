#ifndef CPMDC_RESTART_H
#define CPMDC_RESTART_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc_restart.h
 * @brief Read and write an OpenCPMD / CPMD `RESTART.n` file.
 *
 * The file is a gfortran sequential unformatted stream (4-byte record
 * markers, negative marker means the record continues). Section 1 is the
 * 80-character "STUTTGART VERSION 3.0" header. Each later section starts
 * with a record whose first integer is the number of records that follow
 * (negative for the wavefunction writers). A current OpenCPMD writer also
 * stores an 8-byte file position in that header record; a classic header
 * that is only the integer is accepted too.
 *
 * Coordinates, velocities, and the initial geometry are Bohr, in the order
 * wv30 writes tau0: species, then atom, then x,y,z. The angstrom conversion
 * matches OpenCPMD cnst%fbohr (1/0.529177210859).
 *
 * Replacing coordinates rewrites those 24-byte records and leaves every
 * other record, including plane-wave coefficients, byte for byte. The host
 * endianness must match the Fortran file (little endian on x86_64 gfortran).
 */

/** Angstroms per Bohr, matching OpenCPMD `cnst%fbohr`. */
#define CPMDC_RESTART_ANGSTROM_PER_BOHR 0.529177210859

/** Opaque RESTART image. */
typedef struct cpmdc_restart cpmdc_restart;

/**
 * @brief Read a RESTART image from a memory buffer.
 *
 * `err` receives a diagnostic when the read fails and may be NULL.
 * @return The image, or NULL on failure.
 */
cpmdc_restart *cpmdc_restart_read_mem(const void *bytes, size_t nbytes,
                                      char *err, size_t err_cap);

/**
 * @brief Read a RESTART image from a filesystem path.
 *
 * `err` receives a diagnostic when the read fails and may be NULL.
 * @return The image, or NULL on failure.
 */
cpmdc_restart *cpmdc_restart_read_path(const char *path, char *err,
                                       size_t err_cap);

/** @brief Release an image returned by the readers. */
void cpmdc_restart_free(cpmdc_restart *file);

/** @brief 1 when @p file uses the continuing-record stream layout. */
int cpmdc_restart_is_stream(const cpmdc_restart *file);

/** @brief Section 1 header text, or NULL when the image has none. */
const char *cpmdc_restart_header(const cpmdc_restart *file);

/**
 * @brief Report one section's record count.
 *
 * `section` is 1-based. `count` is the integer in the section header,
 * or 0 when the section is empty.
 * @return 0 when the section is present.
 */
int cpmdc_restart_section(const cpmdc_restart *file, int section, int *count);

/**
 * @brief Read the cell description.
 *
 * @return 0 when the cell section is present.
 */
int cpmdc_restart_cell(const cpmdc_restart *file, int *ibrav, int *indpg,
                       double celldm[6]);

/**
 * @brief Read the species counts.
 *
 * `na_per_species` points at the image's own array until the image is
 * freed or rewritten.
 * @return 0 when the species section is present.
 */
int cpmdc_restart_species(const cpmdc_restart *file, int *nsp,
                          const int **na_per_species);

/** @brief Number of coordinate triples in section 4, or 0 when that section is empty. */
int cpmdc_restart_ncoords(const cpmdc_restart *file);

/**
 * @brief Copy coordinates in Bohr.
 *
 * `xyz` receives `ncoords * 3` doubles.
 * @return 0 on success.
 */
int cpmdc_restart_coordinates(const cpmdc_restart *file, double *xyz, int n3);

/**
 * @brief Copy velocities in Bohr.
 *
 * `xyz` receives `ncoords * 3` doubles.
 * @return 0 on success.
 */
int cpmdc_restart_velocities(const cpmdc_restart *file, double *xyz, int n3);

/**
 * @brief Copy the initial geometry in Bohr.
 *
 * `xyz` receives `ncoords * 3` doubles.
 * @return 0 on success.
 */
int cpmdc_restart_initial_coordinates(const cpmdc_restart *file, double *xyz,
                                      int n3);

/**
 * @brief Replace coordinates, leaving every other record unchanged.
 *
 * `xyz` holds `n3` doubles, `ncoords * 3`.
 * @return 0 on success.
 */
int cpmdc_restart_set_coordinates(cpmdc_restart *file, const double *xyz,
                                  int n3);

/**
 * @brief Replace velocities, leaving every other record unchanged.
 *
 * @return 0 on success.
 */
int cpmdc_restart_set_velocities(cpmdc_restart *file, const double *xyz,
                                 int n3);

/**
 * @brief Replace the cell lattice parameters.
 *
 * @return 0 on success.
 */
int cpmdc_restart_set_cell(cpmdc_restart *file, const double celldm[6]);

/**
 * @brief Read the cutoff section.
 *
 * `dual_flag` is the logical `dual` stored in section 7 (0 or 1), not the
 * real dual factor. That factor is `cdual`.
 *
 * @return 0 when the section is present.
 */
int cpmdc_restart_cutoff(const cpmdc_restart *file, double *ecut, double *cdual,
                         int *dual_flag, int *nel, int *nr1s, int *nr2s,
                         int *nr3s);

/**
 * @brief Read the electronic-state counts.
 *
 * @return 0 when the section is present.
 */
int cpmdc_restart_states(const cpmdc_restart *file, int *n, int *nkpts,
                         int *ngw, int *ngwl, int *nhg, int *nhgl);

/**
 * @brief Write the image to a newly allocated buffer.
 *
 * `bytes` receives a buffer the caller frees. `err` receives a diagnostic
 * when the write fails and may be NULL.
 * @return 0 on success.
 */
int cpmdc_restart_write_mem(const cpmdc_restart *file, void **bytes,
                            size_t *nbytes, char *err, size_t err_cap);

/**
 * @brief Write the image to a filesystem path.
 *
 * `err` receives a diagnostic when the write fails and may be NULL.
 * @return 0 on success.
 */
int cpmdc_restart_write_path(const cpmdc_restart *file, const char *path,
                             char *err, size_t err_cap);

#ifdef __cplusplus
}
#endif

#endif
