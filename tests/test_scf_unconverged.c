/**
 * An SCF held to MAXITER 2 cannot reach the orbital threshold. ODIIS runs
 * out, the embed continues with PCG MINIMIZE, and that runs out too. The
 * failed result, and cpmdc_last_error, carry the iteration each pass
 * stopped at, MAXITER, the final gemax, and the minimizers.
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

static const char *g_params = NULL;

enum { kAtoms = 3 };

static const int kZ[kAtoms] = {8, 1, 1};
static const double kWater[3 * kAtoms] = {0.0, 0.0,    0.1173, 0.0, 0.7572,
                                          -0.4692, 0.0, -0.7572, -0.4692};

struct scf_ctx {
  const unsigned char *params;
  size_t params_size;
};

struct scf_out {
  int created;
  int ok;
  char message[512];
  char last_error[512];
};

static void scf_child(void *vctx, void *vout) {
  const struct scf_ctx *ctx = (const struct scf_ctx *)vctx;
  struct scf_out *o = (struct scf_out *)vout;
  double grad[3 * kAtoms];
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(o->message, sizeof(o->message), "%s", cpmdc_last_error());
    return;
  }
  o->created = 1;
  CPMDCResult r = cpmdc_session_energy_gradient(session, kAtoms, kWater, kZ, grad);
  o->ok = r.ok;
  snprintf(o->message, sizeof(o->message), "%s", r.message);
  snprintf(o->last_error, sizeof(o->last_error), "%s", cpmdc_last_error());
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

static void test_unconverged_scf_message(void **state) {
  (void)state;
  if (!cpmdc_available())
    skip();
  struct scf_ctx ctx;
  struct scf_out out;
  char why[256];
  unsigned char *params = read_file(g_params, &ctx.params_size);
  assert_non_null(params);
  ctx.params = params;
  int rc = run_setup_child(scf_child, &ctx, &out, sizeof(out), why, sizeof(why));
  free(params);
  if (rc != 0)
    fail_msg("child: %s", why);
  if (!out.created)
    fail_msg("session: %s", out.message);
  printf("message: %s\n", out.message);
  assert_int_equal(out.ok, 0);
  assert_string_equal(out.message, out.last_error);
  assert_non_null(strstr(out.message, "CPMD SCF did not converge: "));
  assert_non_null(strstr(out.message, "PCG MINIMIZE stopped at iteration 2 of MAXITER 2"));
  assert_non_null(strstr(out.message, "after ODIIS stopped at iteration 2 of MAXITER 2"));
  assert_non_null(strstr(out.message, "CONVERGENCE ORBITALS 1.000000E-05"));
  const char *at = strstr(out.message, " with gemax ");
  assert_non_null(at);
  double gemax = 0.0;
  assert_int_equal(sscanf(at, " with gemax %lf", &gemax), 1);
  assert_true(isfinite(gemax));
  assert_true(gemax > 1.0e-5);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS_BIN\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_unconverged_scf_message),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
