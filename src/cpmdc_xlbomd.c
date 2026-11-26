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
 * Each session owns its own slot table.
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

/* Host-endian history blob. Magic bytes are 'XLB1'. */
static const uint32_t k_blob_magic = 0x31424C58u;
static const uint32_t k_blob_version = 1u;
static const int k_max_slots = 1024;
static const int k_max_n = 4096;
static const int k_max_count = 16;

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

struct XlStore {
  XlSlot *slots;
  int nslots;
};

static XlStore *g_bound = NULL;
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

static XlSlot *slot_find(XlStore *store, long long key, int has_key,
                         int create) {
  int i;
  XlSlot *grown;
  if (!store)
    return NULL;
  for (i = 0; i < store->nslots; ++i) {
    if (store->slots[i].used && store->slots[i].has_key == has_key &&
        (has_key == 0 || store->slots[i].key == key))
      return &store->slots[i];
  }
  if (!create)
    return NULL;
  for (i = 0; i < store->nslots; ++i) {
    if (!store->slots[i].used) {
      memset(&store->slots[i], 0, sizeof(store->slots[i]));
      store->slots[i].used = 1;
      store->slots[i].has_key = has_key;
      store->slots[i].key = key;
      store->slots[i].order = -1;
      return &store->slots[i];
    }
  }
  grown = (XlSlot *)realloc(store->slots,
                            (size_t)(store->nslots + 1) * sizeof(*grown));
  if (!grown)
    return NULL;
  store->slots = grown;
  memset(&store->slots[store->nslots], 0, sizeof(store->slots[store->nslots]));
  store->slots[store->nslots].used = 1;
  store->slots[store->nslots].has_key = has_key;
  store->slots[store->nslots].key = key;
  store->slots[store->nslots].order = -1;
  store->nslots += 1;
  return &store->slots[store->nslots - 1];
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

XlStore *xlbomd_store_new(void) {
  XlStore *store = (XlStore *)calloc(1, sizeof(*store));
  return store;
}

void xlbomd_store_reset(XlStore *store) {
  int i;
  if (!store)
    return;
  for (i = 0; i < store->nslots; ++i)
    free(store->slots[i].hist);
  free(store->slots);
  store->slots = NULL;
  store->nslots = 0;
}

void xlbomd_store_free(XlStore *store) {
  if (!store)
    return;
  if (g_bound == store)
    g_bound = NULL;
  xlbomd_store_reset(store);
  free(store);
}

void xlbomd_store_clear_key(XlStore *store, long long key, int has_key) {
  XlSlot *slot = slot_find(store, key, has_key, 0);
  if (slot)
    slot_release(slot);
}

void xlbomd_store_clear_current(XlStore *store, int have_key, long long key,
                                int drop_unkeyed) {
  if (have_key)
    xlbomd_store_clear_key(store, key, 1);
  else if (drop_unkeyed)
    xlbomd_store_clear_key(store, 0, 0);
}

void xlbomd_store_note_copy(XlStore *store, long long from_key,
                            long long to_key) {
  if (from_key == to_key)
    return;
  xlbomd_store_clear_key(store, to_key, 1);
}

int xlbomd_store_observe(XlStore *store, int order, long long key, int has_key,
                         const double *density, int n, double *p_next) {
  XlSlot *slot;
  int need;
  size_t step;
  if (!store)
    return -1;
  if (order != 0 && order_row(order) < 0)
    return -1;
  if (n < 0)
    return -1;
  if (n == 0)
    return 0;
  if (!density || !p_next)
    return -1;
  need = samples_needed(order);
  slot = slot_find(store, key, has_key, 1);
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

int xlbomd_store_newest(const XlStore *store, int order, long long key,
                        int has_key, int n, double *out) {
  int i;
  if (!store || !out || n < 0)
    return -1;
  for (i = 0; i < store->nslots; ++i) {
    const XlSlot *slot = &store->slots[i];
    if (!slot->used || slot->has_key != has_key)
      continue;
    if (has_key && slot->key != key)
      continue;
    if (slot->order != order || slot->n != n || slot->count < 1 || !slot->hist)
      return 0;
    memcpy(out, slot->hist, (size_t)n * sizeof(double));
    return 1;
  }
  return 0;
}

static int put_bytes(unsigned char **buf, size_t *len, size_t *cap,
                     const void *src, size_t n) {
  unsigned char *grown;
  size_t next;
  if (n > SIZE_MAX - *len)
    return -1;
  next = *len + n;
  if (next > *cap) {
    size_t want = *cap ? *cap : 64;
    while (want < next) {
      if (want > SIZE_MAX / 2)
        return -1;
      want *= 2;
    }
    grown = (unsigned char *)realloc(*buf, want);
    if (!grown)
      return -1;
    *buf = grown;
    *cap = want;
  }
  if (n)
    memcpy(*buf + *len, src, n);
  *len = next;
  return 0;
}

static int put_u32(unsigned char **buf, size_t *len, size_t *cap, uint32_t v) {
  return put_bytes(buf, len, cap, &v, sizeof(v));
}

static int put_i32(unsigned char **buf, size_t *len, size_t *cap, int32_t v) {
  return put_bytes(buf, len, cap, &v, sizeof(v));
}

int xlbomd_store_export(const XlStore *store, void **bytes, size_t *nbytes) {
  unsigned char *buf = NULL;
  size_t len = 0;
  size_t cap = 0;
  int i;
  int used = 0;
  if (!bytes || !nbytes)
    return -1;
  *bytes = NULL;
  *nbytes = 0;
  if (!store)
    return -1;
  for (i = 0; i < store->nslots; ++i) {
    if (store->slots[i].used && store->slots[i].count > 0)
      used += 1;
  }
  if (used == 0)
    return 0;
  if (put_u32(&buf, &len, &cap, k_blob_magic) != 0 ||
      put_u32(&buf, &len, &cap, k_blob_version) != 0 ||
      put_u32(&buf, &len, &cap, (uint32_t)used) != 0) {
    free(buf);
    return -1;
  }
  for (i = 0; i < store->nslots; ++i) {
    const XlSlot *slot = &store->slots[i];
    size_t samples;
    int64_t key;
    if (!slot->used || slot->count <= 0)
      continue;
    key = (int64_t)slot->key;
    samples = (size_t)slot->count * (size_t)slot->n * sizeof(double);
    if (put_i32(&buf, &len, &cap, slot->has_key ? 1 : 0) != 0 ||
        put_i32(&buf, &len, &cap, slot->order) != 0 ||
        put_i32(&buf, &len, &cap, slot->n) != 0 ||
        put_i32(&buf, &len, &cap, slot->count) != 0 ||
        put_bytes(&buf, &len, &cap, &key, sizeof(key)) != 0 ||
        put_bytes(&buf, &len, &cap, slot->hist, samples) != 0) {
      free(buf);
      return -1;
    }
  }
  *bytes = buf;
  *nbytes = len;
  return 0;
}

static int take_bytes(const unsigned char **p, size_t *left, void *dst,
                      size_t n) {
  if (*left < n)
    return -1;
  if (n)
    memcpy(dst, *p, n);
  *p += n;
  *left -= n;
  return 0;
}

int xlbomd_store_import(XlStore *store, const void *bytes, size_t nbytes) {
  const unsigned char *p;
  size_t left;
  XlStore *tmp;
  uint32_t magic = 0;
  uint32_t version = 0;
  uint32_t nslots = 0;
  uint32_t s;
  if (!store)
    return -1;
  if (!bytes || nbytes == 0) {
    xlbomd_store_reset(store);
    return 0;
  }
  p = (const unsigned char *)bytes;
  left = nbytes;
  if (take_bytes(&p, &left, &magic, sizeof(magic)) != 0 ||
      take_bytes(&p, &left, &version, sizeof(version)) != 0 ||
      take_bytes(&p, &left, &nslots, sizeof(nslots)) != 0)
    return -1;
  if (magic != k_blob_magic || version != k_blob_version || nslots > (uint32_t)k_max_slots)
    return -1;
  tmp = xlbomd_store_new();
  if (!tmp)
    return -1;
  for (s = 0; s < nslots; ++s) {
    int32_t has_key = 0;
    int32_t order = 0;
    int32_t n = 0;
    int32_t count = 0;
    int64_t key = 0;
    XlSlot *slot;
    size_t samples;
    if (take_bytes(&p, &left, &has_key, sizeof(has_key)) != 0 ||
        take_bytes(&p, &left, &order, sizeof(order)) != 0 ||
        take_bytes(&p, &left, &n, sizeof(n)) != 0 ||
        take_bytes(&p, &left, &count, sizeof(count)) != 0 ||
        take_bytes(&p, &left, &key, sizeof(key)) != 0)
      goto fail;
    if ((has_key != 0 && has_key != 1) || n < 1 || n > k_max_n || count < 1 ||
        count > k_max_count)
      goto fail;
    if (order != 0 && order_row(order) < 0)
      goto fail;
    if (count > samples_needed(order))
      goto fail;
    if ((size_t)count > SIZE_MAX / (size_t)n ||
        (size_t)count * (size_t)n > SIZE_MAX / sizeof(double))
      goto fail;
    samples = (size_t)count * (size_t)n * sizeof(double);
    if (left < samples)
      goto fail;
    slot = slot_find(tmp, (long long)key, has_key, 1);
    if (!slot || slot->count != 0)
      goto fail;
    slot->hist = (double *)malloc(samples);
    if (!slot->hist)
      goto fail;
    memcpy(slot->hist, p, samples);
    p += samples;
    left -= samples;
    slot->order = order;
    slot->n = n;
    slot->count = count;
    slot->cap = count;
  }
  if (left != 0)
    goto fail;
  xlbomd_store_reset(store);
  store->slots = tmp->slots;
  store->nslots = tmp->nslots;
  tmp->slots = NULL;
  tmp->nslots = 0;
  xlbomd_store_free(tmp);
  return 0;
fail:
  xlbomd_store_free(tmp);
  return -1;
}

void xlbomd_bind(XlStore *store) { g_bound = store; }

void xlbomd_reset_bound(void) { xlbomd_store_reset(g_bound); }

void xlbomd_clear_bound(int have_key, long long key, int drop_unkeyed) {
  xlbomd_store_clear_current(g_bound, have_key, key, drop_unkeyed);
}

void xlbomd_set_configured(int order) { g_configured = order; }

int xlbomd_configured(void) { return g_configured; }
