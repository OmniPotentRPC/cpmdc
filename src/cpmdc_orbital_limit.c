/* Keyed-orbital cap for CPMDC_ORBITAL_SLOT_LIMIT.
 *
 * The OpenCPMD table in tools/opencpmd_embed_rwfopt.patch is the store when
 * this build links CPMD; that table counts parked c0 plus populated
 * extrapolation histories. This file agrees the cap across calculator ranks
 * and, on the reference evaluator, applies the same rule to keyed evaluations.
 * The reference evaluator has no orbital or extrapolation storage, so its
 * ledger counts parked keys only. Exceeding the cap reports an error and
 * leaves the parked keys in place. The codes ok, invalid and exceeded match
 * the Fortran parameters embed_slot_limit_ok, embed_slot_limit_invalid and
 * embed_slot_limit_exceeded. The no memory code has no Fortran twin. */
#include "cpmdc.h"

#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(CPMDC_HAVE_MPI)
#include <mpi.h>
#endif

enum {
  CPMDC_ORBITAL_SLOT_READ_ENV = -2,
  CPMDC_ORBITAL_SLOT_UNLIMITED = -1,
  CPMDC_ORBITAL_SLOT_OK = 0,
  CPMDC_ORBITAL_SLOT_INVALID = 1,
  CPMDC_ORBITAL_SLOT_EXCEEDED = 2,
  CPMDC_ORBITAL_SLOT_NO_MEMORY = 3,
  CPMDC_ORBITAL_SLOT_DECIMAL_BASE = 10,
  CPMDC_ORBITAL_SLOT_FIRST_CAP = 8,
  CPMDC_ORBITAL_SLOT_KIND_UNSET = 0,
  CPMDC_ORBITAL_SLOT_KIND_SET = 1,
  CPMDC_ORBITAL_SLOT_KIND_INVALID = 2
};

static const char k_orbital_limit_invalid[] =
    "CPMDC_ORBITAL_SLOT_LIMIT must be a positive integer";
static const char k_orbital_limit_exceeded[] =
    "CPMDC_ORBITAL_SLOT_LIMIT exceeded; parked orbital and extrapolation-history storage was left unchanged";
static const char k_orbital_limit_differs[] =
    "CPMDC_ORBITAL_SLOT_LIMIT differs across calculator ranks";
static const char k_orbital_limit_mpi[] =
    "CPMDC_ORBITAL_SLOT_LIMIT: MPI collective failed";
static const char k_orbital_limit_memory[] = "out of memory";
static const char k_orbital_limit_peer[] =
    "another rank of the calculator rejected this call";
#if defined(CPMDC_HAVE_MPI)
static const char k_orbital_key_differs[] =
    "orbital key selection differs across calculator ranks";
#endif

/* Reference-evaluator ledger. Parked keys exclude the selected key.
 * g_selected_has_store is the live copy of the selected key. */
static long long *g_parked = NULL;
static size_t g_parked_n = 0;
static size_t g_parked_cap = 0;
static int g_have_selected = 0;
static long long g_selected_key = 0;
static int g_selected_has_store = 0;
static int g_live_store = 0;
static int g_agreed_limit = CPMDC_ORBITAL_SLOT_READ_ENV;
static int g_limit_ready = 0;

void cpmdc_embed_set_orbital_slot_limit(int limit);
void cpmdc_embed_set_cell_in_place(double limit);
int cpmdc_embed_slot_limit_code(void);
int cpmdc_embed_prepare_orbitals(long long key, int fail_allocation);
void cpmdc_embed_abort_orbitals(void);
void cpmdc_embed_commit_orbitals(void);

const char *cpmdc_orbital_limit_message(int code) {
  if (code == CPMDC_ORBITAL_SLOT_INVALID)
    return k_orbital_limit_invalid;
  if (code == CPMDC_ORBITAL_SLOT_EXCEEDED)
    return k_orbital_limit_exceeded;
  return NULL;
}

void cpmdc_orbital_ledger_reset(void) {
  g_parked_n = 0;
  g_selected_has_store = 0;
  g_live_store = 0;
}

void cpmdc_orbital_ledger_finalize(void) {
  free(g_parked);
  g_parked = NULL;
  g_parked_n = 0;
  g_parked_cap = 0;
  g_have_selected = 0;
  g_selected_key = 0;
  g_selected_has_store = 0;
  g_live_store = 0;
  g_limit_ready = 0;
  g_agreed_limit = CPMDC_ORBITAL_SLOT_READ_ENV;
}

/* Slots allocated for parked keys. For tests that bound the table. */
size_t cpmdc_orbital_ledger_capacity(void) { return g_parked_cap; }

#if !defined(CPMDC_HAS_CPMD)
static const size_t k_parked_not_found = SIZE_MAX;

static size_t parked_index(long long key) {
  for (size_t i = 0; i < g_parked_n; ++i) {
    if (g_parked[i] == key)
      return i;
  }
  return k_parked_not_found;
}

/* Makes room for one more parked key. Grows the table only when it is
 * full, and leaves it as it was when the growth cannot be sized or
 * allocated. */
static int parked_reserve(void) {
  long long *grown;
  size_t cap;
  if (g_parked_n < g_parked_cap)
    return 0;
  if (g_parked_cap == 0)
    cap = CPMDC_ORBITAL_SLOT_FIRST_CAP;
  else if (g_parked_cap > SIZE_MAX / (2u * sizeof(*grown)))
    return -1;
  else
    cap = g_parked_cap * 2u;
  grown = (long long *)realloc(g_parked, cap * sizeof(*grown));
  if (!grown)
    return -1;
  g_parked = grown;
  g_parked_cap = cap;
  return 0;
}

static void parked_remove_at(size_t at) {
  if (at >= g_parked_n)
    return;
  if (at + 1u < g_parked_n) {
    memmove(&g_parked[at], &g_parked[at + 1u],
            (g_parked_n - at - 1u) * sizeof(*g_parked));
  }
  g_parked_n--;
}

/* How many parked keys a select of `key` would add. Zero when the
 * operation only replaces or restores. */
static int select_extra(long long key, int *leaving_needs_park,
                        size_t *dest_parked) {
  int leaving = 0;
  size_t dest = parked_index(key);
  if (g_have_selected && key != g_selected_key && g_selected_has_store &&
      parked_index(g_selected_key) == k_parked_not_found)
    leaving = 1;
  if (leaving_needs_park)
    *leaving_needs_park = leaving;
  if (dest_parked)
    *dest_parked = dest;
  if (leaving && dest == k_parked_not_found)
    return 1;
  return 0;
}

static int ledger_would_exceed(long long key) {
  int extra;
  if (!g_limit_ready || g_agreed_limit == CPMDC_ORBITAL_SLOT_UNLIMITED)
    return 0;
  if (g_agreed_limit <= 0)
    return 1;
  extra = select_extra(key, NULL, NULL);
  if (extra <= 0)
    return 0;
  if (g_parked_n > (size_t)INT_MAX)
    return 1;
  if ((int)g_parked_n > g_agreed_limit - extra)
    return 1;
  return 0;
}

/* cpmdc_orbital_prepare reserves a slot before collective agreement, so this
 * commit contains no allocation and cannot fail on only one rank. */
static void ledger_commit(long long key) {
  int leaving = 0;
  size_t dest = k_parked_not_found;
  if (!g_have_selected) {
    g_have_selected = 1;
    g_selected_key = key;
    g_selected_has_store = g_live_store;
    return;
  }
  if (key == g_selected_key)
    return;
  (void)select_extra(key, &leaving, &dest);
  /* prepare reserved the slot, so this commit does not allocate. */
  if (leaving)
    g_parked[g_parked_n++] = g_selected_key;
  if (dest != k_parked_not_found) {
    parked_remove_at(dest);
    g_selected_has_store = 1;
  } else {
    g_selected_has_store = 0;
  }
  g_selected_key = key;
}
#endif

static void parse_orbital_slot_limit(int *limit_out, int *kind_out) {
  const char *raw = getenv("CPMDC_ORBITAL_SLOT_LIMIT");
  const char *begin;
  const char *end;
  int acc = 0;
  *limit_out = CPMDC_ORBITAL_SLOT_UNLIMITED;
  if (!raw) {
    *kind_out = CPMDC_ORBITAL_SLOT_KIND_UNSET;
    return;
  }
  begin = raw;
  while (*begin == ' ')
    begin++;
  end = raw + strlen(raw);
  while (end > begin && end[-1] == ' ')
    end--;
  *kind_out = CPMDC_ORBITAL_SLOT_KIND_INVALID;
  if (begin == end)
    return;
  if (*begin == '+')
    begin++;
  if (begin == end)
    return;
  for (const char *p = begin; p < end; ++p) {
    int digit;
    if (*p < '0' || *p > '9')
      return;
    digit = *p - '0';
    if (acc > (INT_MAX - digit) / CPMDC_ORBITAL_SLOT_DECIMAL_BASE)
      return;
    acc = acc * CPMDC_ORBITAL_SLOT_DECIMAL_BASE + digit;
  }
  if (acc <= 0)
    return;
  *limit_out = acc;
  *kind_out = CPMDC_ORBITAL_SLOT_KIND_SET;
}

#if defined(CPMDC_HAVE_MPI)
static int orbital_allreduce(int *values, int count) {
  int inited = 0;
  int finalized = 0;
  int size = 1;
  MPI_Comm comm;
  if (MPI_Initialized(&inited) != MPI_SUCCESS)
    return -1;
  if (!inited)
    return 0;
  if (MPI_Finalized(&finalized) != MPI_SUCCESS)
    return -1;
  if (finalized)
    return 0;
  if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0)
    comm = MPI_COMM_WORLD;
  if (MPI_Comm_size(comm, &size) != MPI_SUCCESS)
    return -1;
  if (size <= 1)
    return 0;
  if (MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_INT, MPI_MAX, comm) !=
      MPI_SUCCESS)
    return -1;
  return 0;
}
#else
static int orbital_allreduce(int *values, int count) {
  (void)values;
  (void)count;
  return 0;
}
#endif

int cpmdc_orbital_allreduce(int *values, int count) {
  return orbital_allreduce(values, count);
}

static double parse_cell_in_place(void) {
  const char *raw = getenv("CPMDC_CELL_IN_PLACE");
  char *end = NULL;
  double value;
  if (!raw || !*raw)
    return 0.0;
  value = strtod(raw, &end);
  if (end == raw)
    return 0.0;
  while (*end && isspace((unsigned char)*end))
    end++;
  if (*end || !isfinite(value) || value <= 0.0)
    return 0.0;
  return value;
}

/* The in-place threshold controls whether every parked key is discarded.
 * Read it on calculator rank zero and distribute one value before any rank
 * can make that decision. */
static int agree_cell_in_place(void) {
  double value = parse_cell_in_place();
#if defined(CPMDC_HAVE_MPI)
  int inited = 0;
  int finalized = 0;
  MPI_Comm comm;
  int rank = 0;
  if (MPI_Initialized(&inited) != MPI_SUCCESS)
    return -1;
  if (inited) {
    if (MPI_Finalized(&finalized) != MPI_SUCCESS)
      return -1;
    if (!finalized) {
      if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0)
        comm = MPI_COMM_WORLD;
      if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS)
        return -1;
      if (rank != 0)
        value = 0.0;
      if (MPI_Bcast(&value, 1, MPI_DOUBLE, 0, comm) != MPI_SUCCESS)
        return -1;
    }
  }
#endif
  cpmdc_embed_set_cell_in_place(value);
  return 0;
}

/* Agrees CPMDC_ORBITAL_SLOT_LIMIT across the calculator ranks, with
 * `local_failed` set when a check of this call already failed here. Every
 * rank calls this exactly once per evaluation, after its own checks, and
 * runs the collective even when its own check failed, so a rank whose
 * check failed alone does not leave the others waiting. A failure on any
 * rank fails the call on all of them. Returns 0 when no rank failed and
 * the limit agrees. Otherwise returns -1 and sets `*err` to the reason, or
 * to NULL on a rank whose own check failed, which reports its own reason. */
int cpmdc_orbital_agree(int local_failed, const char **err) {
  int limit = CPMDC_ORBITAL_SLOT_UNLIMITED;
  int kind = CPMDC_ORBITAL_SLOT_KIND_UNSET;
  int buf[6];
  int min_limit;
  g_limit_ready = 0;
  *err = NULL;
  parse_orbital_slot_limit(&limit, &kind);
  buf[0] = kind == CPMDC_ORBITAL_SLOT_KIND_INVALID;
  buf[1] = kind == CPMDC_ORBITAL_SLOT_KIND_SET ? limit : 0;
  buf[2] = kind == CPMDC_ORBITAL_SLOT_KIND_SET;
  buf[3] = kind == CPMDC_ORBITAL_SLOT_KIND_UNSET;
  buf[4] = kind == CPMDC_ORBITAL_SLOT_KIND_SET ? -limit : INT_MIN;
  buf[5] = local_failed != 0;
  if (orbital_allreduce(buf, 6) != 0 || agree_cell_in_place() != 0) {
    if (!local_failed)
      *err = k_orbital_limit_mpi;
    return -1;
  }
  if (buf[5]) {
    if (!local_failed)
      *err = k_orbital_limit_peer;
    return -1;
  }
  if (buf[0]) {
    *err = k_orbital_limit_invalid;
    return -1;
  }
  if (buf[2] && buf[3]) {
    *err = k_orbital_limit_differs;
    return -1;
  }
  if (buf[2]) {
    min_limit = -buf[4];
    if (min_limit != buf[1]) {
      *err = k_orbital_limit_differs;
      return -1;
    }
    g_agreed_limit = buf[1];
  } else {
    g_agreed_limit = CPMDC_ORBITAL_SLOT_UNLIMITED;
  }
  cpmdc_embed_set_orbital_slot_limit(g_agreed_limit);
  g_limit_ready = 1;
  return 0;
}

int cpmdc_orbital_key_agree(int *has_key, long long key, const char **err) {
  int flags[2] = {*has_key != 0, *has_key == 0};
#if defined(CPMDC_HAVE_MPI)
  int inited = 0;
  int finalized = 0;
  int size = 1;
  long long min_key = key;
  long long max_key = key;
  MPI_Comm comm;
  if (orbital_allreduce(flags, 2) != 0) {
    *err = k_orbital_limit_mpi;
    return -1;
  }
  if (flags[0] && flags[1]) {
    *err = k_orbital_key_differs;
    return -1;
  }
  if (MPI_Initialized(&inited) != MPI_SUCCESS ||
      (inited && MPI_Finalized(&finalized) != MPI_SUCCESS)) {
    *err = k_orbital_limit_mpi;
    return -1;
  }
  if (inited && !finalized && flags[0]) {
    if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0)
      comm = MPI_COMM_WORLD;
    if (MPI_Comm_size(comm, &size) != MPI_SUCCESS ||
        (size > 1 &&
         (MPI_Allreduce(MPI_IN_PLACE, &min_key, 1, MPI_LONG_LONG, MPI_MIN,
                        comm) != MPI_SUCCESS ||
          MPI_Allreduce(MPI_IN_PLACE, &max_key, 1, MPI_LONG_LONG, MPI_MAX,
                        comm) != MPI_SUCCESS))) {
      *err = k_orbital_limit_mpi;
      return -1;
    }
    if (min_key != max_key) {
      *err = k_orbital_key_differs;
      return -1;
    }
  }
#else
  (void)key;
#endif
  *has_key = flags[0];
  *err = NULL;
  return 0;
}

int cpmdc_orbital_prepare(long long key, const char **err) {
#if defined(CPMDC_HAS_CPMD)
  int code = CPMDC_ORBITAL_SLOT_OK;
  int fail_allocation = 0;
  const char *msg;
  if (!g_limit_ready)
    code = CPMDC_ORBITAL_SLOT_INVALID;
  else {
#if defined(CPMDC_ENABLE_TEST_FAULTS)
    {
      const char *fault = getenv("CPMDC_TEST_ORBITAL_FAULT");
      if (fault && strcmp(fault, "selector-allocation") == 0)
        fail_allocation = 1;
      if (fault && strcmp(fault, "prepare") == 0)
        code = CPMDC_ORBITAL_SLOT_NO_MEMORY;
    }
#endif
    if (code == CPMDC_ORBITAL_SLOT_OK)
      code = cpmdc_embed_prepare_orbitals(key, fail_allocation);
  }
  if (orbital_allreduce(&code, 1) != 0) {
    cpmdc_embed_abort_orbitals();
    *err = k_orbital_limit_mpi;
    return -1;
  }
  if (code != CPMDC_ORBITAL_SLOT_OK)
    cpmdc_embed_abort_orbitals();
  if (code == CPMDC_ORBITAL_SLOT_NO_MEMORY) {
    *err = k_orbital_limit_memory;
    return -1;
  }
  msg = cpmdc_orbital_limit_message(code);
  if (msg) {
    *err = msg;
    return -1;
  }
  cpmdc_embed_commit_orbitals();
  *err = NULL;
  return 0;
#else
  int fail;
  int leaving = 0;
  if (!g_limit_ready) {
    *err = k_orbital_limit_invalid;
    return -1;
  }
  /* The table grows before the collective, so a rank that cannot grow it
   * fails every rank and no rank is left evaluating alone. */
  fail = ledger_would_exceed(key) ? CPMDC_ORBITAL_SLOT_EXCEEDED : 0;
  (void)select_extra(key, &leaving, NULL);
  if (!fail && leaving && parked_reserve() != 0)
    fail = CPMDC_ORBITAL_SLOT_NO_MEMORY;
#if defined(CPMDC_ENABLE_TEST_FAULTS)
  {
    const char *fault = getenv("CPMDC_TEST_ORBITAL_FAULT");
    if (fault && strcmp(fault, "prepare") == 0)
      fail = CPMDC_ORBITAL_SLOT_NO_MEMORY;
  }
#endif
  if (orbital_allreduce(&fail, 1) != 0) {
    *err = k_orbital_limit_mpi;
    return -1;
  }
  if (fail == CPMDC_ORBITAL_SLOT_NO_MEMORY) {
    *err = k_orbital_limit_memory;
    return -1;
  }
  if (fail) {
    *err = k_orbital_limit_exceeded;
    return -1;
  }
  ledger_commit(key);
  *err = NULL;
  return 0;
#endif
}

void cpmdc_orbital_note_stored(void) {
#if !defined(CPMDC_HAS_CPMD)
  g_live_store = 1;
  if (g_have_selected)
    g_selected_has_store = 1;
#endif
}
