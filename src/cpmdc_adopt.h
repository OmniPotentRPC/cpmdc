#ifndef CPMDC_ADOPT_H
#define CPMDC_ADOPT_H

#include <stddef.h>

/* Query the adopted communicator and choose the calculator index.
 * Rank and size come from comm. When ranks_per_calc > 0 that size must
 * equal it. The group must be one aligned block of that many world ranks.
 * The index is the lowest world rank in the block divided by the size.
 * live is nonzero when this rank has a CPMDCSession. The call reduces that
 * flag with the rank bounds, and any live session fails the whole group.
 * Returns 0 and writes *index_out on success. Returns -1 and writes err
 * on failure. Does not store the handle.
 */
int cpmdc_adopt_resolve(const void *comm, size_t comm_bytes, int ranks_per_calc,
                        int live, int *index_out, char *err, size_t errcap);

/* Remember the handle of the last successful adopt. */
int cpmdc_adopt_commit(const void *comm, size_t comm_bytes);

/* Copy the remembered handle into out. Returns 0, or -1 when none is
 * stored or the buffer is the wrong size. Does not write an error string.
 */
int cpmdc_adopt_copy(void *out, size_t nbytes);

#endif
