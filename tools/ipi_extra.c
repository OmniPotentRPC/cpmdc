#define _POSIX_C_SOURCE 200809L

#include "ipi_extra.h"

#include <locale.h>
#include <math.h>
#include <stdio.h>

int cpmdc_ipi_format_dipole(int valid, const double *values, char *buf,
                            size_t cap, int32_t *out_len) {
  locale_t loc;
  locale_t prev;
  int wrote;
  int i;

  if (out_len == NULL)
    return -1;
  *out_len = 0;
  if (!valid)
    return 0;
  if (values == NULL || buf == NULL || cap == 0)
    return -1;
  for (i = 0; i < 3; ++i)
    if (!isfinite(values[i]))
      return -1;

  loc = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
  if (loc == (locale_t)0)
    return -1;
  prev = uselocale(loc);
  if (prev == (locale_t)0) {
    freelocale(loc);
    return -1;
  }
  wrote = snprintf(buf, cap, "{\"dipole\":[%.16e,%.16e,%.16e]}", values[0],
                   values[1], values[2]);
  uselocale(prev);
  freelocale(loc);
  if (wrote < 0 || (size_t)wrote >= cap)
    return -1;
  *out_len = (int32_t)wrote;
  return 0;
}
