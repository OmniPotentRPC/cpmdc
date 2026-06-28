#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <setjmp.h>
#include <cmocka.h>

#include "schema/Potentials.capnp.h"
#include <capn.h>

static const char *g_params = NULL;
static const char *g_step = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  fseek(fp, 0, SEEK_END);
  long n = ftell(fp);
  rewind(fp);
  unsigned char *buf = malloc((size_t)n);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void test_unsupported_pbeoriginal_rejected(void **state) {
  (void)state;
  assert_int_equal(
      cpmdc_params_reject_unsupported_inputs("PBEoriginal", "&DFT\n FUNCTIONAL BLYP\n&END"),
      -1);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", " N_STREAMS\n  2\n"), -1);
  assert_int_equal(
      cpmdc_params_reject_unsupported_inputs("BLYP", " BLAS_N_STREAMS_PER_DEVICE\n  2\n"), 0);
}

static void test_pef_post_eval_full_surface(void **state) {
  (void)state;
  size_t ps = 0, ss = 0;
  unsigned char *params = read_file(g_params, &ps);
  unsigned char *step = read_file(g_step, &ss);
  assert_non_null(params);
  assert_non_null(step);

  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  size_t need = cpmdc_potential_result_size_for_force_input(step, ss);
  unsigned char *out = malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  CPMDCResult r =
      cpmdc_session_calculate_result(session, step, ss, out, need, &wrote);
  assert_int_equal(r.ok, 1);

  CPMDCEnergyComponents ec;
  assert_int_equal(cpmdc_last_energy_components(&ec), 0);
  assert_int_equal(ec.valid, 1);
  assert_true(fabs(ec.etot - r.energy_h) < 1e-12);

  CPMDCChargeIntegrals ch;
  assert_int_equal(cpmdc_last_charge_integrals(&ch), 0);
  assert_int_equal(ch.valid, 1);

  CPMDCMultiStateEnergies ms;
  assert_int_equal(cpmdc_last_multi_state_energies(&ms), 0);
  assert_int_equal(ms.valid, 1);
  assert_true(ms.count >= 1);
  assert_true(fabs(ms.values[0] - r.energy_h) < 1e-12);

  CPMDCMDTrajectoryRow md;
  assert_int_equal(cpmdc_last_md_trajectory_row(&md), 0);
  assert_int_equal(md.valid, 1);
  assert_true(md.count >= 12);
  assert_true(fabs(md.values[0] - r.energy_h) < 1e-12); /* etot */
  assert_true(md.values[11] == 0.0); /* EKINC off MD */
  printf("md_energy_row etot=%.12f ekin=%.12f ekinc=%.12f count=%zu\n", md.values[0],
         md.values[1], md.values[11], md.count);

  CPMDCPropertySnapshot prop;
  assert_int_equal(cpmdc_last_property_snapshot(&prop), 0);
  assert_int_equal(prop.valid, 1);
  assert_int_equal(prop.dipole_count, 3);
  assert_int_equal(prop.polarizability_count, 9);
  printf("prop_snapshot valid=1 dipole_count=%zu hess_count=%zu pol_count=%zu\n",
         prop.dipole_count, prop.hessian_count, prop.polarizability_count);

  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, out, (int)wrote, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);
  assert_int_equal(view.componentsValid, 1);
  assert_int_equal(capn_len(view.energyComponents), 24);
  double wire_etot = capn_to_f64(capn_get64(view.energyComponents, 0));
  assert_true(fabs(wire_etot - ec.etot) < 1e-12);
  /* Full POD mirror on wire (PEF: only etot non-zero among first components). */
  for (int i = 0; i < 24; ++i) {
    double w = capn_to_f64(capn_get64(view.energyComponents, i));
    double pod = (&ec.etot)[i];
    assert_true(fabs(w - pod) < 1e-12);
  }
  assert_true(cpmdc_feature_find("catalog.cpmd.QMMM") != NULL);
  printf("qmmm_catalog_discoverable=1\n");

  capn_free(&arena);
  free(out);
  cpmdc_session_destroy(session);
  free(params);
  free(step);
  cpmdc_finalize();
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin step.bin\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_unsupported_pbeoriginal_rejected),
      cmocka_unit_test(test_pef_post_eval_full_surface),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
