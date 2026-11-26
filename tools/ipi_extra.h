#ifndef CPMDC_IPI_EXTRA_H
#define CPMDC_IPI_EXTRA_H

#include <stddef.h>
#include <stdint.h>

/* `{"dipole":[%.16e,%.16e,%.16e]}` for three finite binary64 values, plus NUL.
 * %.16e of DBL_MAX is 23 bytes, so three components and the JSON wrapper fit. */
#define CPMDC_IPI_DIPOLE_EXTRA_CAP 128

/* Write the i-PI extra string for a molecular dipole in electron Bohr.
 *
 * valid == 0 leaves *out_len at 0 and does not read values: the library did
 * not recompute a dipole, and the extra string stays empty. A finite dipole
 * is UTF-8 JSON {"dipole":[x,y,z]} with each component printed as %.16e in
 * the C numeric locale, the same text Python's format(value, ".16e")
 * produces. That width is the number of digits in the serialization, not an
 * acceptance tolerance. Returns 0 on success and -1 when a component is not
 * finite, the buffer is too small, or an argument is missing. On failure
 * *out_len is 0 when out_len itself is non-NULL. */
int cpmdc_ipi_format_dipole(int valid, const double *values, char *buf,
                            size_t cap, int32_t *out_len);

#endif
