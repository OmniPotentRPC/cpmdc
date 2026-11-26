#ifndef CPMDC_XLBOMD_INTERNAL_H
#define CPMDC_XLBOMD_INTERNAL_H

#include "cpmdc.h"

#include <string.h>

/* Auxiliary extended-Lagrangian step. Public wrappers live in cpmdc.c.
 * Order 0 is the reversible Verlet map. Orders 3 through 9 are the
 * dissipative rows.
 * Classification is inline so a deck renderer that compiles cpmdc_params.c
 * alone does not need the history store. */

__attribute__((unused)) static inline int
xlbomd_classify(const char *text, int len, int order, int *order_out) {
  if (!text || len != 6 || memcmp(text, "XLBOMD", 6) != 0)
    return CPMDC_WFN_EXTRAP_DECK;
  if (order != 0 && (order < 3 || order > 9))
    return -1;
  if (order_out)
    *order_out = order;
  return CPMDC_WFN_EXTRAP_XLBOMD;
}

int xlbomd_coefficients(int order, double *kappa, double *alpha, double *c,
                        int c_len);

int xlbomd_propagate(int order, const double *density, int n,
                     const double *history, int nhist, double *p_next);

int xlbomd_observe(int order, long long key, int has_key, const double *density,
                   int n, double *p_next);

void xlbomd_clear_key(long long key, int has_key);

/* have_key selects that key's samples. Otherwise drop_unkeyed clears the
 * unkeyed trajectory, matching embed_clear_orbital_key. */
void xlbomd_clear_current(int have_key, long long key, int drop_unkeyed);

/* A copy onto the same key keeps the auxiliary samples. Any other
 * destination loses them. The source is left as it was. */
void xlbomd_note_copy(long long from_key, long long to_key);

void xlbomd_reset(void);

void xlbomd_set_configured(int order);

int xlbomd_configured(void);

#endif
