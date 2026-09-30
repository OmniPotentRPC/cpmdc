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
 * @param bytes RESTART image bytes.
 * @param nbytes Number of bytes in `bytes`.
 * @param err Diagnostic buffer. May be NULL.
 * @param err_cap Capacity of `err`, including the terminating NUL.
 * @return The image, or NULL on failure.
 */
cpmdc_restart *cpmdc_restart_read_mem(const void *bytes, size_t nbytes,
                                      char *err, size_t err_cap);

/**
 * @brief Read a RESTART image from a filesystem path.
 *
 * @param path Filesystem path of the RESTART file.
 * @param err Diagnostic buffer. May be NULL.
 * @param err_cap Capacity of `err`, including the terminating NUL.
 * @return The image, or NULL on failure.
 */
cpmdc_restart *cpmdc_restart_read_path(const char *path, char *err,
                                       size_t err_cap);

/**
 * @brief Release an image returned by the readers.
 *
 * @param file Image to free. May be NULL.
 */
void cpmdc_restart_free(cpmdc_restart *file);

/**
 * @brief Report whether the image uses the continuing-record stream layout.
 *
 * @param file Image to query.
 * @return 1 when `file` uses that layout, 0 otherwise.
 */
int cpmdc_restart_is_stream(const cpmdc_restart *file);

/**
 * @brief Section 1 header text.
 *
 * @param file Image to query.
 * @return The header, or NULL when the image has none.
 */
const char *cpmdc_restart_header(const cpmdc_restart *file);

/**
 * @brief Report one section's record count.
 *
 * @param file Image to query.
 * @param section 1-based section index.
 * @param count Receives the integer in the section header, or 0 when the
 *        section is empty.
 * @return 0 when the section is present.
 */
int cpmdc_restart_section(const cpmdc_restart *file, int section, int *count);

/**
 * @brief Read the cell description.
 *
 * @param file Image to query.
 * @param ibrav Receives the stored `ibrav`. May be NULL.
 * @param indpg Receives the stored `indpg`. May be NULL.
 * @param celldm Receives the six lattice parameters. Required.
 * @return 0 when the cell section is present.
 */
int cpmdc_restart_cell(const cpmdc_restart *file, int *ibrav, int *indpg,
                       double celldm[6]);

/**
 * @brief Read the species counts.
 *
 * @param file Image to query.
 * @param nsp Receives the stored species count. May be NULL.
 * @param na_per_species Points at the image's own array. May be NULL.
 *        Valid until the image is freed or rewritten.
 * @return 0 when the species section is present.
 */
int cpmdc_restart_species(const cpmdc_restart *file, int *nsp,
                          const int **na_per_species);

/**
 * @brief Number of coordinate triples in section 4.
 *
 * @param file Image to query.
 * @return The triple count, or 0 when that section is empty.
 */
int cpmdc_restart_ncoords(const cpmdc_restart *file);

/**
 * @brief Copy coordinates in Bohr.
 *
 * @param file Image to query.
 * @param xyz Receives `ncoords * 3` doubles.
 * @param n3 Capacity of `xyz` in doubles.
 * @return 0 on success.
 */
int cpmdc_restart_coordinates(const cpmdc_restart *file, double *xyz, int n3);

/**
 * @brief Copy velocities in Bohr.
 *
 * @param file Image to query.
 * @param xyz Receives `ncoords * 3` doubles.
 * @param n3 Capacity of `xyz` in doubles.
 * @return 0 on success.
 */
int cpmdc_restart_velocities(const cpmdc_restart *file, double *xyz, int n3);

/**
 * @brief Copy the initial geometry in Bohr.
 *
 * @param file Image to query.
 * @param xyz Receives `ncoords * 3` doubles.
 * @param n3 Capacity of `xyz` in doubles.
 * @return 0 on success.
 */
int cpmdc_restart_initial_coordinates(const cpmdc_restart *file, double *xyz,
                                      int n3);

/**
 * @brief Replace coordinates, leaving every other record unchanged.
 *
 * @param file Image to modify.
 * @param xyz `n3` doubles, `ncoords * 3`, in Bohr.
 * @param n3 Number of doubles in `xyz`.
 * @return 0 on success.
 */
int cpmdc_restart_set_coordinates(cpmdc_restart *file, const double *xyz,
                                  int n3);

/**
 * @brief Replace velocities, leaving every other record unchanged.
 *
 * @param file Image to modify.
 * @param xyz `n3` doubles, in Bohr.
 * @param n3 Number of doubles in `xyz`.
 * @return 0 on success.
 */
int cpmdc_restart_set_velocities(cpmdc_restart *file, const double *xyz,
                                 int n3);

/**
 * @brief Replace the cell lattice parameters.
 *
 * @param file Image to modify.
 * @param celldm Six lattice parameters.
 * @return 0 on success.
 */
int cpmdc_restart_set_cell(cpmdc_restart *file, const double celldm[6]);

/**
 * @brief Read the cutoff section.
 *
 * `dual_flag` is the logical `dual` stored in section 7 (0 or 1), not the
 * real dual factor. That factor is `cdual`.
 *
 * @param file Image to query.
 * @param ecut Receives the stored cutoff. May be NULL.
 * @param cdual Receives the stored real dual factor. May be NULL.
 * @param dual_flag Receives the stored logical dual, 0 or 1. May be NULL.
 * @param nel Receives the stored `nel`. May be NULL.
 * @param nr1s Receives the stored `nr1s`. May be NULL.
 * @param nr2s Receives the stored `nr2s`. May be NULL.
 * @param nr3s Receives the stored `nr3s`. May be NULL.
 * @return 0 when the section is present.
 */
int cpmdc_restart_cutoff(const cpmdc_restart *file, double *ecut, double *cdual,
                         int *dual_flag, int *nel, int *nr1s, int *nr2s,
                         int *nr3s);

/**
 * @brief Read the electronic-state counts.
 *
 * @param file Image to query.
 * @param n Receives the stored state count. May be NULL.
 * @param nkpts Receives the stored k-point count. May be NULL.
 * @param ngw Receives the stored `ngw`. May be NULL.
 * @param ngwl Receives the stored `ngwl`. May be NULL.
 * @param nhg Receives the stored `nhg`. May be NULL.
 * @param nhgl Receives the stored `nhgl`. May be NULL.
 * @return 0 when the section is present.
 */
int cpmdc_restart_states(const cpmdc_restart *file, int *n, int *nkpts,
                         int *ngw, int *ngwl, int *nhg, int *nhgl);

/**
 * @brief Write the image to a newly allocated buffer.
 *
 * @param file Image to write.
 * @param bytes Receives a buffer the caller frees.
 * @param nbytes Receives the number of bytes written.
 * @param err Diagnostic buffer. May be NULL.
 * @param err_cap Capacity of `err`, including the terminating NUL.
 * @return 0 on success.
 */
int cpmdc_restart_write_mem(const cpmdc_restart *file, void **bytes,
                            size_t *nbytes, char *err, size_t err_cap);

/**
 * @brief Write the image to a filesystem path.
 *
 * @param file Image to write.
 * @param path Destination path.
 * @param err Diagnostic buffer. May be NULL.
 * @param err_cap Capacity of `err`, including the terminating NUL.
 * @return 0 on success.
 */
int cpmdc_restart_write_path(const cpmdc_restart *file, const char *path,
                             char *err, size_t err_cap);

#ifdef __cplusplus
}
#endif

#endif
