/**
 * Stored orbitals per key: one session evaluates two geometries in turn
 * under two keys, as a calculator does for two images of a band or two
 * beads of a ring polymer. Each key keeps the converged orbitals of its
 * own last evaluation, and every keyed result equals a cold evaluation of
 * the same geometry.
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

#include <cmocka.h>

#include "setup_child.h"

/* Keys holding parked orbitals in this process, the selected key excluded.
 * Exported by the embed layer (src/cpmd_embed_c_api.F90) for tests. */
int cpmdc_embed_orbital_slots(void);
int cpmdc_embed_copy_orbitals(long long from_key, long long to_key);
/* Where the last SCF started: 0 the selected key's orbitals, 1 another
 * key's, -1 the guess (every stored geometry beyond CPMDC_WARM_RMSD), -2
 * the guess (nothing stored). Exported by the same layer for tests. */
int cpmdc_embed_last_warm_source(void);

static const char *g_params = NULL;

enum { kAtoms = 3, kSteps = 4 };

static const int kZ[kAtoms] = {8, 1, 1};
/* Water, the same molecule with one O-H bond stretched, and both again
 * after a small step, as two images of a band move between iterations. */
static const double kA[3 * kAtoms] = {0.0, 0.0,    0.1173, 0.0, 0.7572,
                                      -0.4692, 0.0, -0.7572, -0.4692};
static const double kB[3 * kAtoms] = {0.0, 0.0,    0.1173, 0.0, 0.9500,
                                      -0.5500, 0.0, -0.7572, -0.4692};
static const double kA2[3 * kAtoms] = {0.0, 0.0,    0.1273, 0.0, 0.7572,
                                       -0.4692, 0.0, -0.7572, -0.4692};
static const double kB2[3 * kAtoms] = {0.0, 0.0,    0.1173, 0.0, 0.9600,
                                       -0.5550, 0.0, -0.7572, -0.4692};

struct keyed_ctx {
  const unsigned char *params;
  size_t params_size;
  /* -1: evaluate only the last geometry, with no key, as a cold call. */
  int keyed;
};

struct keyed_out {
  int ok[kSteps];
  double energy[kSteps];
  double grad[kSteps][3 * kAtoms];
  int slots[kSteps];
  int select_rc;
  char message[256];
};

static void keyed_child(void *vctx, void *vout) {
  const struct keyed_ctx *ctx = (const struct keyed_ctx *)vctx;
  struct keyed_out *o = (struct keyed_out *)vout;
  const double *steps[kSteps] = {kA, kB, kA2, kB2};
  const long long keys[kSteps] = {0, 1, 0, 1};
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  if (ctx->keyed < 0) {
    /* Cold references: one geometry per child, so nothing is stored. */
    const int which = -ctx->keyed - 1;
    CPMDCResult r = cpmdc_session_energy_gradient(
        session, kAtoms, steps[which], kZ, o->grad[which]);
    o->ok[which] = r.ok;
    o->energy[which] = r.energy_h;
    if (!r.ok)
      snprintf(o->message, sizeof(o->message), "%s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  for (int i = 0; i < kSteps; ++i) {
    o->select_rc |= cpmdc_session_select_orbitals(session, keys[i]);
    CPMDCResult r =
        cpmdc_session_energy_gradient(session, kAtoms, steps[i], kZ, o->grad[i]);
    o->ok[i] = r.ok;
    o->energy[i] = r.energy_h;
    o->slots[i] = cpmdc_embed_orbital_slots();
    if (!r.ok) {
      snprintf(o->message, sizeof(o->message), "step %d: %s", i, r.message);
      break;
    }
  }
  cpmdc_session_destroy(session);
}

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

static void test_keys_keep_their_orbitals(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);

  struct keyed_out keyed;
  struct keyed_ctx kctx = {params, params_size, 1};
  RUN_SETUP_CHILD(keyed_child, &kctx, &keyed);
  if (keyed.message[0] != '\0')
    print_message("keyed run: %s\n", keyed.message);
  assert_int_equal(keyed.select_rc, 0);
  for (int i = 0; i < kSteps; ++i)
    assert_int_equal(keyed.ok[i], 1);
  /* Key 0 is selected for the first call, so nothing is parked. Every
   * later call parks the key it leaves; a key returned to takes its
   * orbitals back out of the table. */
  assert_int_equal(keyed.slots[0], 0);
  assert_int_equal(keyed.slots[1], 1);
  assert_int_equal(keyed.slots[2], 1);
  assert_int_equal(keyed.slots[3], 1);

  /* A warm start changes where the SCF starts, not where it ends: each
   * keyed result is the cold result of its geometry to the orbital
   * convergence of the deck (1e-5). */
  for (int i = 0; i < kSteps; ++i) {
    struct keyed_out cold;
    struct keyed_ctx cctx = {params, params_size, -(i + 1)};
    RUN_SETUP_CHILD(keyed_child, &cctx, &cold);
    if (cold.message[0] != '\0')
      print_message("cold run %d: %s\n", i, cold.message);
    assert_int_equal(cold.ok[i], 1);
    assert_float_equal(keyed.energy[i], cold.energy[i], 1e-7);
    for (int k = 0; k < 3 * kAtoms; ++k)
      assert_float_equal(keyed.grad[i][k], cold.grad[i][k], 1e-4);
  }
  free(params);
}

enum { kFarSteps = 5 };

struct far_out {
  int ok[kFarSteps];
  int source[kFarSteps];
  double energy[kFarSteps];
  double grad[kFarSteps][3 * kAtoms];
  char message[256];
};

struct copy_out {
  int ok[3];
  int source;
  int slots_after_copy;
  int copy_rc;
  double energy[3];
  double grad[3][3 * kAtoms];
  char message[256];
};

/* Keys 0, 1, 0, 1, 2 at A, B, A2, B2, A with a 0.05 Angstrom limit. B
 * lies 0.12 Angstrom (root-mean-square per atom) from A and A2 0.006 from
 * A, so key 1 starts from the guess, keys 0 and 1 then start from their
 * own orbitals, and the new key 2 at A starts from key 0's. */
static void far_child(void *vctx, void *vout) {
  const struct keyed_ctx *ctx = (const struct keyed_ctx *)vctx;
  struct far_out *o = (struct far_out *)vout;
  const double *steps[kFarSteps] = {kA, kB, kA2, kB2, kA};
  const long long keys[kFarSteps] = {0, 1, 0, 1, 2};
  setenv("CPMDC_WARM_RMSD", "0.05", 1);
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  for (int i = 0; i < kFarSteps; ++i) {
    cpmdc_session_select_orbitals(session, keys[i]);
    CPMDCResult r =
        cpmdc_session_energy_gradient(session, kAtoms, steps[i], kZ, o->grad[i]);
    o->ok[i] = r.ok;
    o->energy[i] = r.energy_h;
    o->source[i] = cpmdc_embed_last_warm_source();
    if (!r.ok) {
      snprintf(o->message, sizeof(o->message), "step %d: %s", i, r.message);
      break;
    }
  }
  cpmdc_session_destroy(session);
}

static void test_far_orbitals_give_way_to_the_guess(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);
  struct far_out far;
  struct keyed_ctx fctx = {params, params_size, 1};
  RUN_SETUP_CHILD(far_child, &fctx, &far);
  if (far.message[0] != '\0')
    print_message("far run: %s\n", far.message);
  const int want[kFarSteps] = {-2, -1, 0, 0, 1};
  for (int i = 0; i < kFarSteps; ++i) {
    assert_int_equal(far.ok[i], 1);
    assert_int_equal(far.source[i], want[i]);
  }
  /* Every start ends at the cold result of its geometry. */
  const int cold_of[kFarSteps] = {0, 1, 2, 3, 0};
  for (int i = 0; i < kFarSteps; ++i) {
    struct keyed_out cold;
    struct keyed_ctx cctx = {params, params_size, -(cold_of[i] + 1)};
    RUN_SETUP_CHILD(keyed_child, &cctx, &cold);
    assert_int_equal(cold.ok[cold_of[i]], 1);
    assert_float_equal(far.energy[i], cold.energy[cold_of[i]], 1e-7);
    for (int k = 0; k < 3 * kAtoms; ++k)
      assert_float_equal(far.grad[i][k], cold.grad[cold_of[i]][k], 1e-4);
  }
  free(params);
}

static void copy_child(void *vctx, void *vout) {
  const struct keyed_ctx *ctx = (const struct keyed_ctx *)vctx;
  struct copy_out *o = (struct copy_out *)vout;
  memset(o, 0, sizeof(*o));
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  cpmdc_session_select_orbitals(session, 0);
  CPMDCResult seed =
      cpmdc_session_energy_gradient(session, kAtoms, kA, kZ, o->grad[0]);
  o->ok[0] = seed.ok;
  o->energy[0] = seed.energy_h;
  if (!seed.ok) {
    snprintf(o->message, sizeof(o->message), "seed: %s", seed.message);
    cpmdc_session_destroy(session);
    return;
  }
  o->copy_rc = cpmdc_session_copy_orbitals(session, 0, 2);
  o->slots_after_copy = cpmdc_embed_orbital_slots();
  if (o->copy_rc != 0) {
    snprintf(o->message, sizeof(o->message), "copy: %s", cpmdc_last_error());
    cpmdc_session_destroy(session);
    return;
  }
  cpmdc_session_select_orbitals(session, 2);
  CPMDCResult warm =
      cpmdc_session_energy_gradient(session, kAtoms, kA2, kZ, o->grad[1]);
  o->ok[1] = warm.ok;
  o->energy[1] = warm.energy_h;
  o->source = cpmdc_embed_last_warm_source();
  if (!warm.ok) {
    snprintf(o->message, sizeof(o->message), "warm: %s", warm.message);
    cpmdc_session_destroy(session);
    return;
  }
  cpmdc_session_select_orbitals(session, 0);
  CPMDCResult own = cpmdc_session_energy_gradient(session, kAtoms, kA2, kZ, o->grad[2]);
  o->ok[2] = own.ok;
  o->energy[2] = own.energy_h;
  if (!own.ok)
    snprintf(o->message, sizeof(o->message), "own: %s", own.message);
  cpmdc_session_destroy(session);
}

static void test_copy_orbitals_seeds_target_key(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);
  struct copy_out out;
  struct keyed_ctx ctx = {params, params_size, 1};
  RUN_SETUP_CHILD(copy_child, &ctx, &out);
  if (out.message[0] != '\0')
    print_message("copy run: %s\n", out.message);
  assert_int_equal(out.ok[0], 1);
  assert_int_equal(out.copy_rc, 0);
  assert_true(out.slots_after_copy >= 1);
  assert_int_equal(out.ok[1], 1);
  assert_int_equal(out.ok[2], 1);
  assert_int_equal(out.source, 0);
  assert_float_equal(out.energy[1], out.energy[2], 1e-7);
  for (int k = 0; k < 3 * kAtoms; ++k)
    assert_float_equal(out.grad[1][k], out.grad[2][k], 1e-4);
  free(params);
}

struct edge_out {
  int rc_missing;
  int rc_parked;
  int rc_current;
  int slots_seed;
  int slots_unused;
  int slots_missing;
  int slots_before_parked;
  int slots_after_parked;
  int slots_before_reselect;
  int slots_after_reselect;
  int source_store;
  int source_target;
  int source_donor;
  int source_current;
  int ok_store;
  int ok_target;
  int ok_donor;
  int ok_current;
  char error[256];
  char message[256];
};

/* A is the source geometry. B is 0.12 Angstrom (root-mean-square per atom)
 * from A, past a 0.05 limit. Later points stay inside that limit of A or of
 * B and never repeat a geometry the session has already cached. A point on
 * the A side is placed so the copy under test is nearer than every other
 * parked geometry: A+0.01 is already parked before the occupied-target
 * check, so the target sits 0.004 toward -z from A and the donor check sits
 * 0.004 toward +z. */
static void shift_z(double out[3 * kAtoms], const double *in, double dz) {
  memcpy(out, in, sizeof(double) * 3u * (size_t)kAtoms);
  out[2] += dz;
}

static void edge_child(void *vctx, void *vout) {
  const struct keyed_ctx *ctx = (const struct keyed_ctx *)vctx;
  struct edge_out *o = (struct edge_out *)vout;
  double a2[3 * kAtoms], a3[3 * kAtoms], a4[3 * kAtoms];
  double b2[3 * kAtoms], b3[3 * kAtoms], b4[3 * kAtoms];
  double grad[3 * kAtoms];
  CPMDCResult r;
  memset(o, 0, sizeof(*o));
  shift_z(a2, kA, 0.01);
  shift_z(a3, kA, -0.004);
  shift_z(a4, kA, 0.004);
  shift_z(b2, kB, 0.01);
  shift_z(b3, kB, 0.015);
  shift_z(b4, kB, 0.025);
  setenv("CPMDC_WARM_RMSD", "0.05", 1);
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "session: %s", cpmdc_last_error());
    return;
  }
  cpmdc_session_select_orbitals(session, 0);
  r = cpmdc_session_energy_gradient(session, kAtoms, kA, kZ, grad);
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "seed: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  o->slots_seed = cpmdc_embed_orbital_slots();
  cpmdc_session_select_orbitals(session, 99);
  cpmdc_session_select_orbitals(session, 0);
  o->slots_unused = cpmdc_embed_orbital_slots();
  o->rc_missing = cpmdc_session_copy_orbitals(session, 9, 8);
  snprintf(o->error, sizeof(o->error), "%s", cpmdc_last_error());
  o->slots_missing = cpmdc_embed_orbital_slots();

  cpmdc_session_select_orbitals(session, 1);
  r = cpmdc_session_energy_gradient(session, kAtoms, kB, kZ, grad);
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "occupy: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  /* Key 1 is selected and holds B. Copying A onto it replaces that copy. */
  o->rc_current = cpmdc_session_copy_orbitals(session, 0, 1);
  r = cpmdc_session_energy_gradient(session, kAtoms, a2, kZ, grad);
  o->ok_current = r.ok;
  o->source_current = cpmdc_embed_last_warm_source();
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "current: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }

  cpmdc_session_select_orbitals(session, 5);
  r = cpmdc_session_energy_gradient(session, kAtoms, b4, kZ, grad);
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "park seed: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  cpmdc_session_select_orbitals(session, 2);
  r = cpmdc_session_energy_gradient(session, kAtoms, b2, kZ, grad);
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "live: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  o->slots_before_parked = cpmdc_embed_orbital_slots();
  /* Key 5 is parked at B. Copying A onto it replaces that occupied slot
   * and leaves key 2's live store alone. */
  o->rc_parked = cpmdc_session_copy_orbitals(session, 0, 5);
  o->slots_after_parked = cpmdc_embed_orbital_slots();
  r = cpmdc_session_energy_gradient(session, kAtoms, b3, kZ, grad);
  o->ok_store = r.ok;
  o->source_store = cpmdc_embed_last_warm_source();
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "store: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  cpmdc_session_select_orbitals(session, 5);
  r = cpmdc_session_energy_gradient(session, kAtoms, a3, kZ, grad);
  o->ok_target = r.ok;
  o->source_target = cpmdc_embed_last_warm_source();
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "target: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  cpmdc_session_select_orbitals(session, 0);
  r = cpmdc_session_energy_gradient(session, kAtoms, a4, kZ, grad);
  o->ok_donor = r.ok;
  o->source_donor = cpmdc_embed_last_warm_source();
  if (!r.ok) {
    snprintf(o->message, sizeof(o->message), "donor: %s", r.message);
    cpmdc_session_destroy(session);
    return;
  }
  o->slots_before_reselect = cpmdc_embed_orbital_slots();
  cpmdc_session_select_orbitals(session, 77);
  cpmdc_session_select_orbitals(session, 0);
  o->slots_after_reselect = cpmdc_embed_orbital_slots();
  cpmdc_session_destroy(session);
}

static void test_copy_and_select_keep_their_contracts(void **state) {
  (void)state;
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  size_t params_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  assert_non_null(params);
  struct edge_out out;
  struct keyed_ctx ctx = {params, params_size, 1};
  RUN_SETUP_CHILD(edge_child, &ctx, &out);
  if (out.message[0] != '\0')
    print_message("edge run: %s\n", out.message);
  assert_int_equal(out.slots_seed, 0);
  assert_int_equal(out.slots_unused, 0);
  assert_int_equal(out.rc_missing, -1);
  assert_non_null(strstr(out.error, "converged"));
  assert_int_equal(out.slots_missing, 0);
  assert_int_equal(out.rc_current, 0);
  assert_int_equal(out.ok_current, 1);
  assert_int_equal(out.source_current, 0);
  assert_int_equal(out.rc_parked, 0);
  assert_int_equal(out.slots_after_parked, out.slots_before_parked);
  assert_true(out.slots_before_parked >= 1);
  assert_int_equal(out.ok_store, 1);
  assert_int_equal(out.source_store, 0);
  assert_int_equal(out.ok_target, 1);
  assert_int_equal(out.source_target, 0);
  assert_int_equal(out.ok_donor, 1);
  assert_int_equal(out.source_donor, 0);
  assert_int_equal(out.slots_after_reselect, out.slots_before_reselect);
  free(params);
}

static void test_null_session_is_refused(void **state) {
  (void)state;
  assert_int_equal(cpmdc_session_select_orbitals(NULL, 3), -1);
  assert_non_null(strstr(cpmdc_last_error(), "null session"));
  assert_int_equal(cpmdc_session_copy_orbitals(NULL, 3, 4), -1);
  assert_non_null(strstr(cpmdc_last_error(), "null session"));
}

static void test_embed_copy_without_orbitals_is_refused(void **state) {
  (void)state;
  assert_int_equal(cpmdc_embed_copy_orbitals(1, 2), -1);
  assert_int_equal(cpmdc_embed_copy_orbitals(1, 1), -1);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_null_session_is_refused),
      cmocka_unit_test(test_embed_copy_without_orbitals_is_refused),
      cmocka_unit_test(test_keys_keep_their_orbitals),
      cmocka_unit_test(test_far_orbitals_give_way_to_the_guess),
      cmocka_unit_test(test_copy_orbitals_seeds_target_key),
      cmocka_unit_test(test_copy_and_select_keep_their_contracts),
  };
  int failed = cmocka_run_group_tests(tests, NULL, NULL);
  /* `extrapolate-wfn` prints this line when the embed comparison ran.
   * The OpenCPMD job greps for it. A reference run skips the embed
   * cases and must not print it. */
  if (cpmdc_available())
    printf("PASS=%d\n", failed == 0 ? 1 : 0);
  return failed;
}
