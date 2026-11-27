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

/* One session owns one store. The slot table is not process global. */
typedef struct XlStore XlStore;

XlStore *xlbomd_store_new(void);
void xlbomd_store_free(XlStore *store);
void xlbomd_store_reset(XlStore *store);

int xlbomd_store_observe(XlStore *store, int order, long long key, int has_key,
                         const double *density, int n, double *p_next);

/* 1 when this key holds a sample of that order and length. The newest
 * sample, which is the next SCF guess, is copied to out. 0 when that
 * sample is absent. -1 when the store or out is missing. */
int xlbomd_store_newest(const XlStore *store, int order, long long key,
                        int has_key, int n, double *out);

void xlbomd_store_clear_key(XlStore *store, long long key, int has_key);

/* have_key selects that key's samples. Otherwise drop_unkeyed clears the
 * unkeyed trajectory, matching embed_clear_orbital_key. */
void xlbomd_store_clear_current(XlStore *store, int have_key, long long key,
                                int drop_unkeyed);

/* A copy onto the same key keeps the auxiliary samples. Any other
 * destination loses them. The source is left as it was. */
void xlbomd_store_note_copy(XlStore *store, long long from_key,
                            long long to_key);

/* Host-endian blob. Empty history writes *bytes NULL and *nbytes 0.
 * The caller frees a non-NULL buffer. Import replaces the store, and a
 * NULL buffer of length 0 clears it. A bad blob leaves the store as it was. */
int xlbomd_store_export(const XlStore *store, void **bytes, size_t *nbytes);
int xlbomd_store_import(XlStore *store, const void *bytes, size_t nbytes);

/* The orbital ledger is still process global. It acts on the store bound
 * by the session that is evaluating, and on no store when none is bound. */
void xlbomd_bind(XlStore *store);
void xlbomd_reset_bound(void);
void xlbomd_clear_bound(int have_key, long long key, int drop_unkeyed);

void xlbomd_set_configured(int order);

int xlbomd_configured(void);

#endif
