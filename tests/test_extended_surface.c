#include "cpmdc.h"
#include "cpmdc_params.h"

#include <math.h>
#include <stdarg.h>
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

static void test_unsupported_rejected(void **state) {
  (void)state;
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("PBEoriginal", ""), -1);
  assert_int_equal(cpmdc_params_reject_unsupported_inputs("BLYP", "N_STREAMS\n"), -1);
}

static void test_dipole_dynamics_with_kpoints_rejected(void **state) {
  (void)state;
  const char *both =
      "&CPMD\n OPTIMIZE WAVEFUNCTION\n DIPOLE DYNAMICS\n&END\n"
      "&SYSTEM\n KPOINTS MONKHORST-PACK\n 2 2 2\n CUTOFF\n 40.0\n&END\n";
  const char *lower =
      "&cpmd\n dipole dynamics sample wannier\n&end\n"
      "&system\n  kpoints\n 1\n 0.0 0.0 0.0 1.0\n&end\n";
  const char *gamma_only =
      "&CPMD\n DIPOLE DYNAMICS\n&END\n&SYSTEM\n CUTOFF\n 40.0\n&END\n";
  const char *kpoints_only =
      "&CPMD\n OPTIMIZE WAVEFUNCTION\n&END\n&SYSTEM\n KPOINTS\n 1\n"
      " 0.0 0.0 0.0 1.0\n&END\n";
  /* RESTART KPOINTS in &CPMD is a restart option, not a k-point mesh. */
  const char *restart =
      "&CPMD\n DIPOLE DYNAMICS\n RESTART WAVEFUNCTION KPOINTS\n&END\n"
      "&SYSTEM\n CUTOFF\n 40.0\n&END\n";
  const char *reason = cpmdc_params_deck_conflict(both);
  assert_non_null(reason);
  assert_non_null(strstr(reason, "DIPOLE DYNAMICS"));
  assert_non_null(strstr(reason, "KPOINTS"));
  assert_non_null(cpmdc_params_deck_conflict(lower));
  assert_null(cpmdc_params_deck_conflict(gamma_only));
  assert_null(cpmdc_params_deck_conflict(kpoints_only));
  assert_null(cpmdc_params_deck_conflict(restart));
  assert_null(cpmdc_params_deck_conflict(""));
}

static void test_commented_dipole_dynamics_is_not_a_conflict(void **state) {
  (void)state;
  const char *commented_dipole =
      "&CPMD\n OPTIMIZE WAVEFUNCTION\n ! DIPOLE DYNAMICS\n&END\n"
      "&SYSTEM\n KPOINTS MONKHORST-PACK\n 2 2 2\n&END\n";
  const char *bang_tight =
      "&CPMD\n!DIPOLE DYNAMICS\n&END\n"
      "&SYSTEM\n KPOINTS\n 1\n 0.0 0.0 0.0 1.0\n&END\n";
  const char *inline_on_another_keyword =
      "&CPMD\n OPTIMIZE WAVEFUNCTION ! DIPOLE DYNAMICS\n&END\n"
      "&SYSTEM\n KPOINTS\n 1\n&END\n";
  const char *commented_kpoints =
      "&CPMD\n DIPOLE DYNAMICS\n&END\n"
      "&SYSTEM\n ! KPOINTS MONKHORST-PACK\n 2 2 2\n CUTOFF\n 40.0\n&END\n";
  const char *trailing_comments_still_conflict =
      "&CPMD\n DIPOLE DYNAMICS ! berry\n&END\n"
      "&SYSTEM\n KPOINTS ! monkhorst\n 2 2 2\n&END\n";
  const char *comment_does_not_open_a_section =
      "!&CPMD\n DIPOLE DYNAMICS\n"
      "&SYSTEM\n KPOINTS\n 1\n&END\n";

  assert_null(cpmdc_params_deck_conflict(commented_dipole));
  assert_null(cpmdc_params_deck_conflict(bang_tight));
  assert_null(cpmdc_params_deck_conflict(inline_on_another_keyword));
  assert_null(cpmdc_params_deck_conflict(commented_kpoints));
  assert_null(cpmdc_params_deck_conflict(comment_does_not_open_a_section));
  assert_non_null(cpmdc_params_deck_conflict(trailing_comments_still_conflict));
}

static void test_all_results_through_c_and_wire(void **state) {
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
  printf("energy_h=%.12f\n", r.energy_h);

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

  CPMDCMDTrajectoryRow md;
  assert_int_equal(cpmdc_last_md_trajectory_row(&md), 0);
  assert_int_equal(md.valid, 1);
  assert_true(md.count >= 12);
  assert_true(fabs(md.values[0] - r.energy_h) < 1e-12);
  assert_true(md.values[11] == 0.0);
  printf("pod_md etot=%.12f ekinc=%.12f count=%zu\n", md.values[0], md.values[11],
         md.count);

  CPMDCPropertySnapshot prop;
  assert_int_equal(cpmdc_last_property_snapshot(&prop), 0);
  assert_int_equal(prop.valid, 1);
  assert_int_equal(prop.dipole_count, 0);
  assert_int_equal(prop.polarizability_count, 9);
  for (int i = 0; i < 9; ++i)
    assert_true(prop.polarizability[i] == 0.0); /* never computed */
  assert_true(prop.hessian_count >= 3); /* PEF 2 atoms => 6 */
  printf("pod_prop dipole=%zu hess=%zu pol=%zu\n", prop.dipole_count,
         prop.hessian_count, prop.polarizability_count);
  double gradient[4096];
  size_t gradient_count = 0;
  assert_int_equal(cpmdc_last_gradient(gradient, 4096, &gradient_count), 0);
  assert_int_equal(gradient_count, prop.hessian_count);
  for (size_t i = 0; i < gradient_count; ++i)
    assert_true(fabs(gradient[i] - prop.hessian[i]) < 1e-12);
  gradient_count = 0;
  assert_int_equal(
      cpmdc_session_last_gradient(session, gradient, 1, &gradient_count), -1);
  assert_int_equal(gradient_count, prop.hessian_count);
  CPMDCDipole dipole;
  memset(&dipole, 0xab, sizeof(dipole));
  assert_int_equal(cpmdc_last_dipole(&dipole), -1);
  assert_int_equal(dipole.valid, 0);
  assert_true(dipole.values[0] == 0.0);
  {
    double held[4];
    memset(held, 0x5a, sizeof(held));
    assert_int_equal(cpmdc_last_charges(NULL, 4), -1);
    assert_int_equal(cpmdc_session_last_charges(NULL, held, 4), -1);
    if (cpmdc_last_charges(held, 4) != 0) {
      assert_int_equal(cpmdc_session_last_charges(session, held, 4), -1);
      for (int i = 0; i < 4; ++i)
        assert_true(((unsigned char *)&held[i])[0] == 0x5a);
    }
  }

  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, out, (int)wrote, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);

  assert_int_equal(view.componentsValid, 1);
  assert_int_equal(capn_len(view.energyComponents), 24);
  for (int i = 0; i < 24; ++i) {
    double w = capn_to_f64(capn_get64(view.energyComponents, i));
    double pod = (&ec.etot)[i];
    assert_true(fabs(w - pod) < 1e-12);
  }

  assert_int_equal(view.chargeValid, 1);
  assert_int_equal(capn_len(view.chargeIntegrals), 4);
  assert_true(fabs(capn_to_f64(capn_get64(view.chargeIntegrals, 0)) - ch.csumg) <
              1e-12);

  assert_int_equal(view.multiStateValid, 1);
  assert_true(capn_len(view.multiStateEnergies) >= 1);
  assert_true(fabs(capn_to_f64(capn_get64(view.multiStateEnergies, 0)) -
                   ms.values[0]) < 1e-12);

  assert_int_equal(view.mdTrajectoryValid, 1);
  assert_true(capn_len(view.mdTrajectoryRow) >= 12);
  assert_true(fabs(capn_to_f64(capn_get64(view.mdTrajectoryRow, 0)) - md.values[0]) <
              1e-12);
  assert_true(capn_to_f64(capn_get64(view.mdTrajectoryRow, 11)) == 0.0);

  assert_int_equal(capn_len(view.dipole), 0);
  assert_true(view.polarizability.p.type != CAPN_NULL);
  assert_int_equal(capn_len(view.polarizability), 9);
  assert_true(view.gradient.p.type != CAPN_NULL);
  assert_true((size_t)capn_len(view.gradient) == prop.hessian_count);
  for (size_t i = 0; i < prop.hessian_count; ++i) {
    double w = capn_to_f64(capn_get64(view.gradient, (int)i));
    assert_true(fabs(w - prop.hessian[i]) < 1e-12);
  }
  printf("wire_full components=%d charge=%d md=%d ms=%d grad_len=%d "
         "embedMdPropsSkipped_raw=%d\n",
         view.componentsValid, view.chargeValid, view.mdTrajectoryValid,
         view.multiStateValid, capn_len(view.gradient),
         (int)view.embedMdPropsSkipped);
  assert_int_equal(view.embedMdPropsSkipped, 0);
  assert_int_equal(view.chargeValid, 1);
  assert_int_equal(view.mdTrajectoryValid, 1);
  assert_int_equal(view.multiStateValid, 1);

  assert_true(cpmdc_feature_find("catalog.cpmd.QMMM") != NULL);
  printf("qmmm_catalog_ok=1\n");

  /* Dual consistency: second eval same PEF energy and wire etot */
  size_t wrote2 = 0;
  CPMDCResult r2 =
      cpmdc_session_calculate_result(session, step, ss, out, need, &wrote2);
  assert_int_equal(r2.ok, 1);
  assert_true(fabs(r2.energy_h - r.energy_h) < 1e-12);
  CPMDCMDTrajectoryRow md2;
  assert_int_equal(cpmdc_last_md_trajectory_row(&md2), 0);
  assert_true(fabs(md2.values[0] - md.values[0]) < 1e-12);
  printf("dual_run energy_h=%.12f md_etot=%.12f\n", r2.energy_h, md2.values[0]);

  capn_free(&arena);
  free(out);
  cpmdc_session_destroy(session);
  free(params);
  free(step);
}

static void water_step(double *positions, int *numbers);

static void test_nonfinite_position_beats_gradient_snapshot_limit(void **state) {
  (void)state;
  enum { atoms = 1366, components = atoms * 3 };
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  double *positions = calloc((size_t)components, sizeof(*positions));
  double *gradient = calloc((size_t)components, sizeof(*gradient));
  int *numbers = malloc((size_t)atoms * sizeof(*numbers));
  assert_non_null(positions);
  assert_non_null(gradient);
  assert_non_null(numbers);
  for (int i = 0; i < atoms; ++i)
    numbers[i] = 1;
  positions[4] = nan("");

  CPMDCResult shot =
      cpmdc_energy_gradient(atoms, positions, numbers, params, ps, gradient);
  assert_int_equal(shot.ok, 0);
  assert_string_equal(shot.message, "positions contain a non-finite component");
  assert_string_equal(cpmdc_last_error(), shot.message);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
  assert_true(shot.energy_h == 0.0);

  CPMDCResult result = cpmdc_session_energy_gradient(session, atoms, positions,
                                                    numbers, gradient);
  assert_int_equal(result.ok, 0);
  assert_string_equal(result.message, "positions contain a non-finite component");
  assert_string_equal(cpmdc_last_error(), result.message);
  assert_int_equal(cpmdc_last_error_code(), CPMDC_ERROR_INVALID_ARGUMENT);
  assert_true(result.energy_h == 0.0);
  assert_null(strstr(result.message, "1365"));

  double water_positions[9];
  int water_numbers[3];
  water_step(water_positions, water_numbers);
  CPMDCResult water = cpmdc_session_energy_gradient(
      session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(water.ok, 1);
  assert_string_equal(water.message, "ok");

  free(numbers);
  free(gradient);
  free(positions);
  cpmdc_session_destroy(session);
  free(params);
}

static void test_gradient_limit_fails_before_evaluation(void **state) {
  (void)state;
  enum { atoms = 1366, components = atoms * 3 };
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  double *positions = calloc(components, sizeof(*positions));
  double *gradient = calloc(components, sizeof(*gradient));
  int *numbers = malloc(atoms * sizeof(*numbers));
  assert_non_null(positions);
  assert_non_null(gradient);
  assert_non_null(numbers);
  for (int i = 0; i < atoms; ++i)
    numbers[i] = 1;
  CPMDCResult result =
      cpmdc_session_energy_gradient(session, atoms, positions, numbers, gradient);
  assert_int_equal(result.ok, 0);
  assert_non_null(strstr(result.message, "at most 1365 atoms"));
  free(numbers);
  free(gradient);
  free(positions);
  cpmdc_session_destroy(session);
  free(params);
}

static void water_step(double *positions, int *numbers) {
  static const double water[9] = {0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692,
                                  0.0, -0.7572, -0.4692};
  memcpy(positions, water, sizeof(water));
  numbers[0] = 8;
  numbers[1] = 1;
  numbers[2] = 1;
}

static void test_topology_change_is_reported_before_the_gradient_limit(
    void **state) {
  (void)state;
  enum { atoms = 1366, components = atoms * 3 };
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, ps);
  assert_non_null(session);
  double *positions = calloc(components, sizeof(*positions));
  double *gradient = calloc(components, sizeof(*gradient));
  int *numbers = malloc(atoms * sizeof(*numbers));
  assert_non_null(positions);
  assert_non_null(gradient);
  assert_non_null(numbers);
  for (int i = 0; i < atoms; ++i)
    numbers[i] = 1;

  double water_positions[9];
  int water_numbers[3];
  water_step(water_positions, water_numbers);

  /* A first call that is too large is refused for its size and leaves the
   * session without a topology. Water is then accepted as the topology and
   * evaluates cleanly on both the reference and OpenCPMD builds. */
  CPMDCResult big = cpmdc_session_energy_gradient(session, atoms, positions,
                                                  numbers, gradient);
  assert_int_equal(big.ok, 0);
  assert_non_null(strstr(big.message, "at most 1365 atoms"));
  CPMDCResult first = cpmdc_session_energy_gradient(
      session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(first.ok, 1);
  assert_string_equal(first.message, "ok");
  assert_string_equal(cpmdc_last_error(), "");

  /* The session now holds water. Another system, even one that is too
   * large, is a topology change. */
  CPMDCResult changed = cpmdc_session_energy_gradient(session, atoms, positions,
                                                      numbers, gradient);
  assert_int_equal(changed.ok, 0);
  assert_non_null(strstr(changed.message, "topology change"));
  assert_null(strstr(changed.message, "at most 1365 atoms"));
  assert_non_null(strstr(cpmdc_last_error(), "topology change"));

  /* The refused call did not disturb the stored topology. */
  CPMDCResult again = cpmdc_session_energy_gradient(
      session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(again.ok, 1);
  assert_string_equal(cpmdc_last_error(), "");
  assert_null(strstr(again.message, "topology change"));

  free(numbers);
  free(gradient);
  free(positions);
  cpmdc_session_destroy(session);
  free(params);
}

static void test_finalize_leaves_open_topology_uncommitted_and_preserves_fixed_topology(
    void **state) {
  (void)state;
  enum { atoms = 1366, components = atoms * 3 };
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  CPMDCSession *fixed_session = cpmdc_session_create(params, ps);
  assert_non_null(fixed_session);
  double *positions = calloc(components, sizeof(*positions));
  double *gradient = calloc(components, sizeof(*gradient));
  int *numbers = malloc(atoms * sizeof(*numbers));
  assert_non_null(positions);
  assert_non_null(gradient);
  assert_non_null(numbers);
  for (int i = 0; i < atoms; ++i)
    numbers[i] = 1;

  double water_positions[9];
  int water_numbers[3];
  water_step(water_positions, water_numbers);

  CPMDCResult seeded = cpmdc_session_energy_gradient(
      fixed_session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(seeded.ok, 1);

  cpmdc_finalize();

  /* After shutdown, a fresh session fails embed configuration before
   * committing any topology, leaving it configurable and without a topology. */
  CPMDCSession *open_session = cpmdc_session_create(params, ps);
  assert_non_null(open_session);
  CPMDCResult unavailable = cpmdc_session_energy_gradient(
      open_session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(unavailable.ok, 0);
  assert_non_null(strstr(unavailable.message, "CPMD embed not available"));
  assert_string_equal(cpmdc_last_error(), unavailable.message);
  assert_int_equal(
      cpmdc_session_set_constant_cutoff(open_session, 30.0, 0.0, 0.0), 0);

  CPMDCResult still_open = cpmdc_session_energy_gradient(
      open_session, atoms, positions, numbers, gradient);
  assert_int_equal(still_open.ok, 0);
  assert_non_null(strstr(still_open.message, "at most 1365 atoms"));
  assert_null(strstr(still_open.message, "topology change"));

  /* A session that already fixed a topology before shutdown still reports a
   * topology change ahead of both the atom limit and embed availability. */
  CPMDCResult fixed_changed = cpmdc_session_energy_gradient(
      fixed_session, atoms, positions, numbers, gradient);
  assert_int_equal(fixed_changed.ok, 0);
  assert_non_null(strstr(fixed_changed.message, "topology change"));
  assert_null(strstr(fixed_changed.message, "at most 1365 atoms"));
  assert_null(strstr(fixed_changed.message, "CPMD embed not available"));

  water_positions[0] += 0.05;
  CPMDCResult fixed_same = cpmdc_session_energy_gradient(
      fixed_session, 3, water_positions, water_numbers, gradient);
  assert_int_equal(fixed_same.ok, 0);
  assert_non_null(strstr(fixed_same.message, "CPMD embed not available"));

  free(numbers);
  free(gradient);
  free(positions);
  cpmdc_session_destroy(open_session);
  cpmdc_session_destroy(fixed_session);
  free(params);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin step.bin\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_unsupported_rejected),
      cmocka_unit_test(test_dipole_dynamics_with_kpoints_rejected),
      cmocka_unit_test(test_commented_dipole_dynamics_is_not_a_conflict),
      cmocka_unit_test(test_all_results_through_c_and_wire),
      cmocka_unit_test(test_nonfinite_position_beats_gradient_snapshot_limit),
      cmocka_unit_test(test_gradient_limit_fails_before_evaluation),
      cmocka_unit_test(test_topology_change_is_reported_before_the_gradient_limit),
      cmocka_unit_test(
          test_finalize_leaves_open_topology_uncommitted_and_preserves_fixed_topology),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
