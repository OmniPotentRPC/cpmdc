/**
 * Bounded parked orbital storage under CPMDC_ORBITAL_SLOT_LIMIT:
 *   1) Invalid values (zero, negative, empty, non-integer, overflow) fail
 *      evaluation with a diagnostic naming CPMDC_ORBITAL_SLOT_LIMIT rather
 *      than silently falling back to unlimited storage.
 *   2) Visiting more bead geometries than CPMDC_ORBITAL_SLOT_LIMIT caps the
 *      parked slot count at the limit, evicts the oldest parked beads first,
 *      never evicts the active or newly selected bead, bounds orbital memory,
 *      and converges every geometry to the same energy and gradient as an
 *      unlimited session.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

#include "setup_child.h"

int cpmdc_embed_orbital_slots(void);
int cpmdc_embed_last_warm_source(void);

static const char *g_params = NULL;

enum {
  kAtoms = 3,
  kCoords = 3 * kAtoms,
  kBeadSteps = 7,
  kSlotLimit = 2
};

static const int kZ[kAtoms] = {8, 1, 1};

/* Five water bead geometries stretched along one O-H coordinate so each bead
 * sits more than 0.10 Angstrom RMSD per atom from every other bead, plus small
 * local displacements around bead 3 and bead 0. */
static const double kBead0[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.7572,
                                       -0.4692, 0.0, -0.7572, -0.4692};
static const double kBead1[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.9400,
                                       -0.5400, 0.0, -0.7572, -0.4692};
static const double kBead2[kCoords] = {0.0, 0.0,     0.1173,  0.0,    1.1300,
                                       -0.6150, 0.0, -0.7572, -0.4692};
static const double kBead3[kCoords] = {0.0, 0.0,     0.1173,  0.0,    1.3200,
                                       -0.6900, 0.0, -0.7572, -0.4692};
static const double kBead4[kCoords] = {0.0, 0.0,     0.1173,  0.0,    1.5100,
                                       -0.7650, 0.0, -0.7572, -0.4692};
static const double kBead3Near[kCoords] = {0.0, 0.0,     0.1233,  0.0,    1.3200,
                                           -0.6900, 0.0, -0.7572, -0.4692};
static const double kBead0Near[kCoords] = {0.0, 0.0,     0.1233,  0.0,    0.7572,
                                           -0.4692, 0.0, -0.7572, -0.4692};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
  if (buf && fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static long current_rss_bytes(void) {
  FILE *fp = fopen("/proc/self/statm", "r");
  unsigned long total_pages = 0;
  unsigned long rss_pages = 0;
  long page_size = sysconf(_SC_PAGESIZE);
  if (!fp)
    return 0;
  if (fscanf(fp, "%lu %lu", &total_pages, &rss_pages) != 2) {
    fclose(fp);
    return 0;
  }
  fclose(fp);
  if (page_size <= 0)
    return 0;
  return (long)(rss_pages * (unsigned long)page_size);
}

struct invalid_ctx {
  const unsigned char *params;
  size_t params_size;
};

struct invalid_out {
  int bad_count;
  int bad_rejected;
  int bad_named_var;
  int bad_last_error_matched;
  int unset_ok;
  char message[256];
};

static void invalid_limit_child(void *vctx, void *vout) {
  static const char *const kBadValues[] = {
      "0",   "+0",  "-1",  "-8",   "",
      "   ", "abc", "2.5", "2extra", "1 2",
      "999999999999999999999999999999"};
  const struct invalid_ctx *ctx = (const struct invalid_ctx *)vctx;
  struct invalid_out *o = (struct invalid_out *)vout;
  double grad[kCoords];
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  (void)cpmdc_session_select_orbitals(session, 0);
  o->bad_count = (int)(sizeof(kBadValues) / sizeof(kBadValues[0]));
  for (int i = 0; i < o->bad_count; ++i) {
    setenv("CPMDC_ORBITAL_SLOT_LIMIT", kBadValues[i], 1);
    CPMDCResult r =
        cpmdc_session_energy_gradient(session, kAtoms, kBead0, kZ, grad);
    const char *last = cpmdc_last_error();
    if (r.ok == 0)
      o->bad_rejected++;
    else if (o->message[0] == '\0')
      snprintf(o->message, sizeof(o->message), "accepted invalid '%s'",
               kBadValues[i]);
    if (strstr(r.message, "CPMDC_ORBITAL_SLOT_LIMIT") != NULL)
      o->bad_named_var++;
    else if (o->message[0] == '\0')
      snprintf(o->message, sizeof(o->message),
               "missing variable name for '%s': %s", kBadValues[i], r.message);
    if (last && strcmp(r.message, last) == 0)
      o->bad_last_error_matched++;
  }
  unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
  {
    CPMDCResult r =
        cpmdc_session_energy_gradient(session, kAtoms, kBead0, kZ, grad);
    o->unset_ok = r.ok;
    if (!r.ok && o->message[0] == '\0')
      snprintf(o->message, sizeof(o->message), "unset failed: %s", r.message);
  }
  cpmdc_session_destroy(session);
}

static void test_invalid_orbital_slot_limits_are_refused(void **state) {
  (void)state;
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);

  struct invalid_out out;
  struct invalid_ctx ctx = {params, params_size};
  RUN_SETUP_CHILD(invalid_limit_child, &ctx, &out);
  if (out.message[0] != '\0')
    print_message("invalid-limit check: %s\n", out.message);
  assert_true(out.bad_count > 0);
  assert_int_equal(out.bad_rejected, out.bad_count);
  assert_int_equal(out.bad_named_var, out.bad_count);
  assert_int_equal(out.bad_last_error_matched, out.bad_count);
  assert_int_equal(out.unset_ok, 1);
  free(params);
}

struct bead_ctx {
  const unsigned char *params;
  size_t params_size;
  int slot_limit;
};

struct bead_out {
  int ok[kBeadSteps];
  int slots[kBeadSteps];
  int source[kBeadSteps];
  long rss_bytes[kBeadSteps];
  double energy[kBeadSteps];
  double grad[kBeadSteps][kCoords];
  int select_rc;
  char message[256];
};

static void bead_sweep_child(void *vctx, void *vout) {
  const struct bead_ctx *ctx = (const struct bead_ctx *)vctx;
  struct bead_out *o = (struct bead_out *)vout;
  const double *geoms[kBeadSteps] = {kBead0, kBead1,     kBead2,    kBead3,
                                     kBead4, kBead3Near, kBead0Near};
  const long long keys[kBeadSteps] = {0, 1, 2, 3, 4, 3, 0};

  setenv("CPMDC_WARM_RMSD", "0.05", 1);
  if (ctx->slot_limit > 0) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%d", ctx->slot_limit);
    setenv("CPMDC_ORBITAL_SLOT_LIMIT", buf, 1);
  } else {
    unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
  }

  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  for (int i = 0; i < kBeadSteps; ++i) {
    o->select_rc |= cpmdc_session_select_orbitals(session, keys[i]);
    CPMDCResult r =
        cpmdc_session_energy_gradient(session, kAtoms, geoms[i], kZ, o->grad[i]);
    o->ok[i] = r.ok;
    o->energy[i] = r.energy_h;
    o->slots[i] = cpmdc_embed_orbital_slots();
    o->source[i] = cpmdc_embed_last_warm_source();
    o->rss_bytes[i] = current_rss_bytes();
    if (!r.ok) {
      snprintf(o->message, sizeof(o->message), "step %d: %s", i, r.message);
      break;
    }
  }
  cpmdc_session_destroy(session);
}

static void test_slot_limit_bounds_parked_orbitals_and_matches_unlimited(
    void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);

  struct bead_out unlimited;
  struct bead_ctx uctx = {params, params_size, 0};
  RUN_SETUP_CHILD(bead_sweep_child, &uctx, &unlimited);
  if (unlimited.message[0] != '\0')
    print_message("unlimited run: %s\n", unlimited.message);

  struct bead_out bounded;
  struct bead_ctx bctx = {params, params_size, kSlotLimit};
  RUN_SETUP_CHILD(bead_sweep_child, &bctx, &bounded);
  if (bounded.message[0] != '\0')
    print_message("bounded run: %s\n", bounded.message);

  assert_int_equal(unlimited.select_rc, 0);
  assert_int_equal(bounded.select_rc, 0);

  const int want_unlimited_slots[kBeadSteps] = {0, 1, 2, 3, 4, 4, 4};
  const int want_bounded_slots[kBeadSteps] = {0, 1, 2, 2, 2, 2, 2};
  for (int i = 0; i < kBeadSteps; ++i) {
    assert_int_equal(unlimited.ok[i], 1);
    assert_int_equal(bounded.ok[i], 1);
    assert_int_equal(unlimited.slots[i], want_unlimited_slots[i]);
    assert_true(bounded.slots[i] <= kSlotLimit);
    assert_int_equal(bounded.slots[i], want_bounded_slots[i]);
  }

  /* Unlimited run keeps both bead 3 and bead 0 warm when revisited.
   * Bounded run (limit 2) keeps recent bead 3 warm without evicting it on
   * selection, whereas oldest bead 0 was evicted when bead 3 parked. */
  assert_int_equal(unlimited.source[5], 0);
  assert_int_equal(unlimited.source[6], 0);
  assert_int_equal(bounded.source[5], 0);
  assert_int_equal(bounded.source[6], -1);

  print_message(
      "orbital-slot-limit summary: unlimited_max_slots=%d bounded_max_slots=%d "
      "unlimited_rss_end=%ld bounded_rss_end=%ld\n",
      unlimited.slots[kBeadSteps - 1], bounded.slots[kBeadSteps - 1],
      unlimited.rss_bytes[kBeadSteps - 1], bounded.rss_bytes[kBeadSteps - 1]);

  /* Eviction changes only which starting orbitals remain cached, never the
   * converged SCF energy or nuclear gradient. */
  for (int i = 0; i < kBeadSteps; ++i) {
    assert_float_equal(bounded.energy[i], unlimited.energy[i], 1e-7);
    for (int k = 0; k < kCoords; ++k)
      assert_float_equal(bounded.grad[i][k], unlimited.grad[i][k], 1e-4);
  }

  free(params);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_invalid_orbital_slot_limits_are_refused),
      cmocka_unit_test(
          test_slot_limit_bounds_parked_orbitals_and_matches_unlimited),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
