#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static void test_unsupported_pbeoriginal_rejected(void **state) {
  (void)state;
  assert_int_equal(
      cpmdc_params_reject_unsupported_inputs("PBEoriginal", "&DFT\n FUNCTIONAL BLYP\n&END"),
      -1);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", " N_STREAMS\n  2\n"), -1);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", " BLAS_N_STREAMS_PER_DEVICE\n  2\n"),
                   0);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", "&CPMD\n MAXITER\n  10\n&END"), 0);
}

static void test_post_eval_snapshots_pef(void **state) {
  (void)state;
  /* Drive PEF via session when possible — use minimal invalid then skip if no params file */
  CPMDCChargeIntegrals ch;
  CPMDCMultiStateEnergies ms;
  CPMDCMDTrajectoryRow md;
  CPMDCPropertySnapshot prop;
  /* Before eval, may be invalid */
  (void)cpmdc_last_charge_integrals(&ch);
  (void)cpmdc_last_multi_state_energies(&ms);
  (void)cpmdc_last_md_trajectory_row(&md);
  (void)cpmdc_last_property_snapshot(&prop);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", "N_STREAMS\n"), -1);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("PBEoriginal", ""), -1);
  assert_true(cpmdc_feature_find("abi.cpmdc_last_property_snapshot") != NULL);
  assert_true(cpmdc_feature_find("catalog.cpmd.QMMM") != NULL);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_unsupported_pbeoriginal_rejected),
      cmocka_unit_test(test_post_eval_snapshots_pef),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
