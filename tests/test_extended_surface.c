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

static void test_pef_post_eval_honesty_and_wire(void **state) {
  (void)state;
  size_t ps = 0, ss = 0;
  unsigned char *params = read_file(g_params, &ps);
  unsigned char *step = read_file(g_step, &ss);
  assert_non_null(params);
  assert_non_null(step);

  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  size_t need = cpmdc_potential_result_size_for_force_input(step, ss);
  assert_true(need >= 80);
  unsigned char *out = malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  CPMDCResult r =
      cpmdc_session_calculate_result(session, step, ss, out, need, &wrote);
  assert_int_equal(r.ok, 1);
  assert_true(isfinite(r.energy_h));

  CPMDCEnergyComponents ec;
  assert_int_equal(cpmdc_last_energy_components(&ec), 0);
  assert_int_equal(ec.valid, 1);
  assert_true(fabs(ec.etot - r.energy_h) < 1e-12);
  /* PEF: only etot non-zero among components */
  assert_true(ec.ekin == 0.0 && ec.exc == 0.0);

  CPMDCChargeIntegrals ch;
  assert_int_equal(cpmdc_last_charge_integrals(&ch), -1);
  assert_int_equal(ch.valid, 0);

  CPMDCMultiStateEnergies ms;
  assert_int_equal(cpmdc_last_multi_state_energies(&ms), -1);
  assert_int_equal(ms.valid, 0);

  CPMDCMDTrajectoryRow md;
  assert_int_equal(cpmdc_last_md_trajectory_row(&md), -1);
  assert_int_equal(md.valid, 0);
  printf("embed_md_props_skipped=1 (PEF SCF path; MD ENERGY/EKINC not harvested)\n");

  CPMDCPropertySnapshot prop;
  assert_int_equal(cpmdc_last_property_snapshot(&prop), -1);
  assert_int_equal(prop.valid, 0);
  printf("embed_prop_skipped=1 (PEF SCF path; PROP Hessian/dipole not harvested)\n");

  /* Cap'n Proto: energy matches; componentsValid only if embed snapshot live. */
  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, out, (int)wrote, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);
  assert_true(fabs(view.energy - r.energy_h) < 1e-9 ||
              fabs(view.energy - r.energy_h * 27.211386245988) < 1e-6);
  /* PEF sets snapshot_total_only => componentsValid with etot==POD etot, other zeros. */
  assert_int_equal(view.componentsValid, 1);
  assert_true(view.energyComponents.p.type != CAPN_NULL);
  assert_int_equal(capn_len(view.energyComponents), 24);
  double wire_etot = capn_to_f64(capn_get64(view.energyComponents, 0));
  double wire_ekin = capn_to_f64(capn_get64(view.energyComponents, 1));
  assert_true(fabs(wire_etot - ec.etot) < 1e-12);
  assert_true(wire_ekin == 0.0);
  assert_int_equal(view.embedMdPropsSkipped, 1);
  assert_true(cpmdc_feature_find("catalog.cpmd.QMMM") != NULL);
  printf("qmmm_catalog_discoverable=1 feature_id=catalog.cpmd.QMMM\n");

  /* QMMM keyword renders when set in deck path — reject not applied to QMMM */
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", " QMMM\n"), 0);
  printf("qmmm_host_input_accepted=1 (catalog keyword; PEF has no QMMM forces)\n");

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
      cmocka_unit_test(test_pef_post_eval_honesty_and_wire),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
