/**
 * One process sets OpenCPMD up again: after a basis change (new parameters,
 * or a new composition in a new session) and after a stopgm during setup.
 * Every energy matches, to 1e-8 Ha, the energy of a fresh process that sets
 * up that system once. The fresh-process references run in child processes
 * before this process touches CPMD.
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

#define ENERGY_TOL_HA 1e-8

struct molecule {
  int n;
  double pos[3 * 8];
  int z[8];
};

/* Positions in Angstrom, the convention of tests/force_input_water. */
static const struct molecule water = {
    3,
    {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692},
    {8, 1, 1}};

/* H2O2 uses the water deck: same elements, another composition. */
static const struct molecule peroxide = {
    4,
    {0.0, 0.7375, -0.0528, 0.0, -0.7375, -0.0528, 0.8190, 0.8657, 0.4220,
     -0.8190, -0.8657, 0.4220},
    {8, 8, 1, 1}};

/* H3Si-NH2 in the middle of the 10 Angstrom box: N-H 1.01, Si-N 1.72,
 * Si-H 1.48 Angstrom, tetrahedral at Si. */
static const struct molecule silylamine = {
    7,
    {4.14, 5.0, 5.0, 5.86, 5.0, 5.0, 3.79, 5.95, 5.0, 3.79, 4.05, 5.0,
     6.353, 6.395, 5.0, 6.353, 4.302, 6.208, 6.353, 4.302, 3.792},
    {7, 14, 1, 1, 1, 1, 1}};

struct blob {
  unsigned char *data;
  size_t size;
};

static struct blob g_water_params, g_silylamine_params, g_stop_params;

struct ref_ctx {
  const struct blob *params;
  const struct molecule *mol;
  int session;
};

struct ref_out {
  CPMDCResult r;
};

static struct ref_out g_ref_water, g_ref_silylamine, g_ref_water_session,
    g_ref_peroxide_session;
static int g_have_refs = 0;

static struct blob read_file(const char *path) {
  struct blob b = {NULL, 0};
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return b;
  if (fseek(fp, 0, SEEK_END) == 0) {
    long n = ftell(fp);
    if (n > 0) {
      rewind(fp);
      b.data = (unsigned char *)malloc((size_t)n);
      if (b.data && fread(b.data, 1, (size_t)n, fp) == (size_t)n) {
        b.size = (size_t)n;
      } else {
        free(b.data);
        b.data = NULL;
      }
    }
  }
  fclose(fp);
  return b;
}

static CPMDCResult one_shot(const struct blob *params,
                            const struct molecule *mol) {
  double forces[3 * 8];
  return cpmdc_energy_forces(mol->n, mol->pos, mol->z, params->data,
                             params->size, forces);
}

static CPMDCResult in_session(CPMDCSession *session,
                              const struct molecule *mol) {
  double forces[3 * 8];
  return cpmdc_session_energy_forces(session, mol->n, mol->pos, mol->z,
                                     forces);
}

static void ref_child(void *vctx, void *vout) {
  struct ref_ctx *ctx = (struct ref_ctx *)vctx;
  struct ref_out *o = (struct ref_out *)vout;
  if (ctx->session) {
    CPMDCSession *s =
        cpmdc_session_create(ctx->params->data, ctx->params->size);
    if (!s)
      return;
    o->r = in_session(s, ctx->mol);
    cpmdc_session_destroy(s);
  } else {
    o->r = one_shot(ctx->params, ctx->mol);
  }
}

static int fresh_process_energy(const struct blob *params,
                                const struct molecule *mol, int session,
                                struct ref_out *out) {
  char why[256];
  struct ref_ctx ctx = {params, mol, session};
  if (run_setup_child(ref_child, &ctx, out, sizeof(*out), why, sizeof(why)) !=
      0) {
    fprintf(stderr, "reference child: %s\n", why);
    return -1;
  }
  if (!out->r.ok) {
    fprintf(stderr, "reference child: %s\n", out->r.message);
    return -1;
  }
  return 0;
}

static int group_setup(void **state) {
  (void)state;
  if (!cpmdc_available())
    return 0;
  if (fresh_process_energy(&g_water_params, &water, 0, &g_ref_water) != 0 ||
      fresh_process_energy(&g_silylamine_params, &silylamine, 0,
                           &g_ref_silylamine) != 0 ||
      fresh_process_energy(&g_water_params, &water, 1,
                           &g_ref_water_session) != 0 ||
      fresh_process_energy(&g_water_params, &peroxide, 1,
                           &g_ref_peroxide_session) != 0)
    return -1;
  g_have_refs = 1;
  return 0;
}

static void require_embed(void) {
  if (!cpmdc_available()) {
    print_message("[  SKIP   ] libcpmd embed not linked\n");
    skip();
  }
  assert_true(g_have_refs);
}

static void assert_energy(CPMDCResult r, const struct ref_out *ref,
                          const char *what) {
  if (!r.ok)
    fail_msg("%s: %s", what, r.message);
  print_message("%s: E = %.10f Ha, fresh process %.10f Ha, diff %.3e\n", what,
                r.energy_h, ref->r.energy_h, r.energy_h - ref->r.energy_h);
  if (fabs(r.energy_h - ref->r.energy_h) > ENERGY_TOL_HA)
    fail_msg("%s: energy %.12f differs from the fresh process %.12f", what,
             r.energy_h, ref->r.energy_h);
}

/* New parameters each call: cpmdc_embed_reset_state tears the setup down. */
static void test_new_params_set_up_again(void **state) {
  (void)state;
  require_embed();
  assert_energy(one_shot(&g_water_params, &water), &g_ref_water, "water");
  assert_energy(one_shot(&g_silylamine_params, &silylamine), &g_ref_silylamine,
                "silylamine after water");
  assert_energy(one_shot(&g_water_params, &water), &g_ref_water,
                "water after silylamine");
}

/* A new session with another composition: the cold path of its first call
 * tears the setup of the previous session down. A session itself refuses
 * a topology change. */
static CPMDCResult new_session(const struct molecule *mol) {
  CPMDCSession *s =
      cpmdc_session_create(g_water_params.data, g_water_params.size);
  if (!s)
    fail_msg("session: %s", cpmdc_last_error());
  CPMDCResult r = in_session(s, mol);
  cpmdc_session_destroy(s);
  return r;
}

static void test_new_composition_new_session(void **state) {
  (void)state;
  require_embed();
  assert_energy(new_session(&water), &g_ref_water_session, "session water");
  assert_energy(new_session(&peroxide), &g_ref_peroxide_session,
                "session H2O2 after water");
  assert_energy(new_session(&water), &g_ref_water_session,
                "session water after H2O2");
}

/* HF_INIT stops on exact exchange with Hockney, late in the setup. */
static void test_setup_after_stopgm(void **state) {
  (void)state;
  require_embed();
  CPMDCResult stop = one_shot(&g_stop_params, &water);
  assert_int_equal(stop.ok, 0);
  print_message("stopped setup: %s\n", stop.message);
  assert_non_null(strstr(stop.message, "stopgm"));
  assert_energy(one_shot(&g_water_params, &water), &g_ref_water,
                "water after a stopgm in setup");
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s water.bin silylamine.bin water_pbe0.bin\n",
            argv[0]);
    return 2;
  }
  g_water_params = read_file(argv[1]);
  g_silylamine_params = read_file(argv[2]);
  g_stop_params = read_file(argv[3]);
  if (!g_water_params.data || !g_silylamine_params.data ||
      !g_stop_params.data) {
    fprintf(stderr, "cannot read the parameter files\n");
    return 2;
  }
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_new_params_set_up_again),
      cmocka_unit_test(test_new_composition_new_session),
      cmocka_unit_test(test_setup_after_stopgm),
  };
  int rc = cmocka_run_group_tests(tests, group_setup, NULL);
  cpmdc_finalize();
  free(g_water_params.data);
  free(g_silylamine_params.data);
  free(g_stop_params.data);
  return rc;
}
