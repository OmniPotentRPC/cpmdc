/* Extended-Lagrangian Born-Oppenheimer propagation of one auxiliary vector.
 *
 * The reversible step is the Verlet map of Niklasson, Phys. Rev. Lett. 100,
 * 123004 (2008), with kappa = delta_t^2 * omega^2 = 2. The dissipative step
 * is equation (22) of Niklasson, Steneteg, Odell, Bock, Challacombe, Tymczak,
 * Holmstrom, Zheng, and Weber, J. Chem. Phys. 130, 214109 (2009):
 *
 *   P_{n+1} = 2 P_n - P_{n-1} + kappa (D_n - P_n)
 *             + alpha * sum_{k=0}^{K} c_k P_{n-k}
 *
 * kappa, alpha, and c_k for K = 3..9 are that paper's Table I. The printed
 * alpha is alpha times 10^3, and this file stores the unscaled alpha.
 * Samples are most recent first. Each orbital key keeps its own samples.
 */
#include "cpmdc_xlbomd.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double k_kappa[] = {1.69, 1.75, 1.82, 1.84, 1.86, 1.88, 1.89};

/* Table I alpha times 10^3. */
static const double k_alpha_milli[] = {150.0, 57.0, 18.0, 5.5, 1.6, 0.44, 0.12};

/* c_0 .. c_K, then unused zeros. Row 0 is K = 3. */
static const int k_c[7][10] = {
    {-2, 3, 0, -1, 0, 0, 0, 0, 0, 0},
    {-3, 6, -2, -2, 1, 0, 0, 0, 0, 0},
    {-6, 14, -8, -3, 4, -1, 0, 0, 0, 0},
    {-14, 36, -27, -2, 12, -6, 1, 0, 0, 0},
    {-36, 99, -88, 11, 32, -25, 8, -1, 0, 0},
    {-99, 286, -286, 78, 78, -90, 42, -10, 1, 0},
    {-286, 858, -936, 364, 168, -300, 184, -63, 12, -1},
};

static int order_row(int order) {
  if (order < 3 || order > 9)
    return -1;
  return order - 3;
}

static int samples_needed(int order) {
  if (order == 0)
    return 2;
  if (order_row(order) < 0)
    return -1;
  return order + 1;
}

int xlbomd_coefficients(int order, double *kappa, double *alpha, double *c,
                        int c_len) {
  int row;
  int k;
  if (!kappa)
    return -1;
  if (order == 0) {
    *kappa = 2.0;
    if (alpha)
      *alpha = 0.0;
    return 0;
  }
  row = order_row(order);
  if (row < 0)
    return -1;
  if (c && c_len < order + 1)
    return -1;
  *kappa = k_kappa[row];
  if (alpha)
    *alpha = k_alpha_milli[row] * 1.0e-3;
  if (c) {
    for (k = 0; k <= order; ++k)
      c[k] = (double)k_c[row][k];
  }
  return 0;
}

int xlbomd_propagate(int order, const double *density, int n,
                     const double *history, int nhist, double *p_next) {
  double kappa = 0.0;
  double alpha = 0.0;
  double coeff[10];
  int need;
  int i;
  int k;
  if (n < 0)
    return -1;
  if (n == 0)
    return 0;
  if (!density || !history || !p_next)
    return -1;
  need = samples_needed(order);
  if (need < 0 || nhist < need)
    return -1;
  if (xlbomd_coefficients(order, &kappa, &alpha, order == 0 ? NULL : coeff,
                          order == 0 ? 0 : 10) != 0)
    return -1;
  for (i = 0; i < n; ++i) {
    double pn = history[i];
    double pnm1 = history[n + i];
    double acc = (2.0 - kappa) * pn - pnm1 + kappa * density[i];
    if (order != 0) {
      double diss = 0.0;
      for (k = 0; k <= order; ++k)
        diss += coeff[k] * history[(size_t)k * (size_t)n + (size_t)i];
      acc += alpha * diss;
    }
    p_next[i] = acc;
  }
  return 0;
}

typedef struct XlSlot {
  int used;
  int has_key;
  long long key;
  int order;
  int n;
  int count;
  int cap;
  double *hist;
} XlSlot;

static XlSlot *g_slots = NULL;
static int g_nslots = 0;
static int g_configured = -1;

static void slot_release(XlSlot *slot) {
  free(slot->hist);
  memset(slot, 0, sizeof(*slot));
}

static int slot_bytes(int samples, int n, size_t *out) {
  if (samples < 0 || n < 0)
    return -1;
  if (samples == 0 || n == 0) {
    *out = 0;
    return 0;
  }
  if ((size_t)samples > SIZE_MAX / (size_t)n)
    return -1;
  if ((size_t)samples * (size_t)n > SIZE_MAX / sizeof(double))
    return -1;
  *out = (size_t)samples * (size_t)n * sizeof(double);
  return 0;
}

static XlSlot *slot_find(long long key, int has_key, int create) {
  int i;
  XlSlot *grown;
  for (i = 0; i < g_nslots; ++i) {
    if (g_slots[i].used && g_slots[i].has_key == has_key &&
        (has_key == 0 || g_slots[i].key == key))
      return &g_slots[i];
  }
  if (!create)
    return NULL;
  for (i = 0; i < g_nslots; ++i) {
    if (!g_slots[i].used) {
      memset(&g_slots[i], 0, sizeof(g_slots[i]));
      g_slots[i].used = 1;
      g_slots[i].has_key = has_key;
      g_slots[i].key = key;
      g_slots[i].order = -1;
      return &g_slots[i];
    }
  }
  grown = (XlSlot *)realloc(g_slots, (size_t)(g_nslots + 1) * sizeof(*grown));
  if (!grown)
    return NULL;
  g_slots = grown;
  memset(&g_slots[g_nslots], 0, sizeof(g_slots[g_nslots]));
  g_slots[g_nslots].used = 1;
  g_slots[g_nslots].has_key = has_key;
  g_slots[g_nslots].key = key;
  g_slots[g_nslots].order = -1;
  g_nslots += 1;
  return &g_slots[g_nslots - 1];
}

static int slot_reserve(XlSlot *slot, int samples, int n) {
  size_t bytes = 0;
  double *hist;
  if (slot->cap >= samples && slot->n == n && (n == 0 || slot->hist))
    return 0;
  if (slot_bytes(samples, n, &bytes) != 0)
    return -1;
  hist = (double *)realloc(slot->hist, bytes == 0 ? 1 : bytes);
  if (!hist && bytes != 0)
    return -1;
  slot->hist = hist;
  slot->cap = samples;
  slot->n = n;
  if (slot->count > samples)
    slot->count = samples;
  return 0;
}

void xlbomd_clear_key(long long key, int has_key) {
  XlSlot *slot = slot_find(key, has_key, 0);
  if (slot)
    slot_release(slot);
}

void xlbomd_clear_current(int have_key, long long key, int drop_unkeyed) {
  if (have_key)
    xlbomd_clear_key(key, 1);
  else if (drop_unkeyed)
    xlbomd_clear_key(0, 0);
}

void xlbomd_note_copy(long long from_key, long long to_key) {
  if (from_key == to_key)
    return;
  xlbomd_clear_key(to_key, 1);
}

void xlbomd_reset(void) {
  int i;
  for (i = 0; i < g_nslots; ++i)
    free(g_slots[i].hist);
  free(g_slots);
  g_slots = NULL;
  g_nslots = 0;
}

void xlbomd_set_configured(int order) { g_configured = order; }

int xlbomd_configured(void) { return g_configured; }

int xlbomd_observe(int order, long long key, int has_key, const double *density,
                   int n, double *p_next) {
  XlSlot *slot;
  int need;
  size_t step;
  if (order != 0 && order_row(order) < 0)
    return -1;
  if (n < 0)
    return -1;
  if (n == 0)
    return 0;
  if (!density || !p_next)
    return -1;
  need = samples_needed(order);
  slot = slot_find(key, has_key, 1);
  if (!slot)
    return -1;
  if (slot->order != order || slot->n != n) {
    free(slot->hist);
    slot->hist = NULL;
    slot->count = 0;
    slot->cap = 0;
    slot->order = order;
    slot->n = n;
  }
  if (slot_reserve(slot, need, n) != 0)
    return -1;
  step = (size_t)n * sizeof(double);
  if (slot->count == 0) {
    memcpy(slot->hist, density, step);
    slot->count = 1;
    memcpy(p_next, density, step);
    return 0;
  }
  if (order == 0 && slot->count == 1) {
    memcpy(slot->hist + n, slot->hist, step);
    if (xlbomd_propagate(0, density, n, slot->hist, 2, p_next) != 0)
      return -1;
    memmove(slot->hist + n, slot->hist, step);
    memcpy(slot->hist, p_next, step);
    slot->count = 2;
    return 1;
  }
  if (slot->count < need) {
    memmove(slot->hist + n, slot->hist, (size_t)slot->count * step);
    memcpy(slot->hist, density, step);
    slot->count += 1;
    memcpy(p_next, density, step);
    return 0;
  }
  if (xlbomd_propagate(order, density, n, slot->hist, slot->count, p_next) != 0)
    return -1;
  memmove(slot->hist + n, slot->hist, (size_t)(need - 1) * step);
  memcpy(slot->hist, p_next, step);
  slot->count = need;
  return 1;
}
