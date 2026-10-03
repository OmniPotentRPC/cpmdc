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

static void test_null_session_is_refused(void **state) {
  (void)state;
  assert_int_equal(cpmdc_session_select_orbitals(NULL, 3), -1);
  assert_non_null(strstr(cpmdc_last_error(), "null session"));
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_null_session_is_refused),
      cmocka_unit_test(test_keys_keep_their_orbitals),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
