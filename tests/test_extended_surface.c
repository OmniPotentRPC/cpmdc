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
  cpmdc_finalize();
}

/*
 * The legacy property snapshot holds 4096 gradient components. Systems around
 * and above that size are evaluated by the reference evaluator, which fills
 * the same snapshot code as an OpenCPMD run but costs nothing per atom. An
 * OpenCPMD build cannot run an SCF of 1366 atoms, so these tests skip there.
 */
enum { k_capacity_atoms = 1365, k_snapshot_components = 4096 };

struct big_system {
  int n_atoms;
  double *positions_ang;
  int *numbers;
};

static void big_system_init(struct big_system *sys, int n_atoms) {
  sys->n_atoms = n_atoms;
  sys->positions_ang = calloc((size_t)n_atoms * 3, sizeof(double));
  sys->numbers = malloc((size_t)n_atoms * sizeof(int));
  assert_non_null(sys->positions_ang);
  assert_non_null(sys->numbers);
  for (int i = 0; i < n_atoms; ++i) {
    sys->numbers[i] = (i % 2) ? 8 : 1;
    sys->positions_ang[3 * i] = 0.011 * (double)(1 + i % 37);
    sys->positions_ang[3 * i + 1] = -0.017 * (double)(1 + i % 23);
    sys->positions_ang[3 * i + 2] = 0.005 * (double)(1 + i % 91);
  }
}

static void big_system_free(struct big_system *sys) {
  free(sys->positions_ang);
  free(sys->numbers);
}

/* The reference gradient is k * Z * x_bohr, so its ratio to Z * x_bohr is one
 * constant. A gradient cut off at the snapshot capacity would break it past
 * the 4096th component. */
static void check_reference_gradient(const struct big_system *sys,
                                     const double *gradient) {
  const double bohr = 0.52917721092;
  double k = 0.0;
  int found = 0;
  for (int i = 0; i < sys->n_atoms; ++i) {
    for (int j = 0; j < 3; ++j) {
      double x = sys->positions_ang[3 * i + j] / bohr * (double)sys->numbers[i];
      double g = gradient[3 * i + j];
      assert_true(isfinite(g));
      if (!found) {
        k = g / x;
        found = 1;
      }
      assert_true(fabs(g - k * x) < 1e-9 * fabs(k * x));
    }
  }
  assert_true(k != 0.0);
}

/* Both last_gradient forms after an evaluation of n_atoms. Up to the capacity
 * they copy the gradient; above it they return -1 with the required length
 * and the snapshot holds no legacy copy. */
static void check_gradient_accessors(const CPMDCSession *session, int n_atoms,
                                     const double *gradient) {
  const size_t need = (size_t)n_atoms * 3;
  double *copy = calloc(need, sizeof(double));
  assert_non_null(copy);
  CPMDCPropertySnapshot *prop = malloc(sizeof(*prop));
  assert_non_null(prop);
  assert_int_equal(cpmdc_session_last_property_snapshot(session, prop), 0);
  assert_int_equal(prop->valid, 1);

  size_t count = 123;
  assert_int_equal(cpmdc_session_last_gradient(session, NULL, 0, &count), -1);
  assert_int_equal(count, need);
  count = 123;
  assert_int_equal(cpmdc_session_last_gradient(session, copy, 1, &count), -1);
  assert_int_equal(count, need);

  if (need <= k_snapshot_components) {
    assert_int_equal(prop->hessian_count, need);
    count = 0;
    assert_int_equal(cpmdc_session_last_gradient(session, copy, need, &count), 0);
    assert_int_equal(count, need);
    for (size_t i = 0; i < need; ++i)
      assert_true(copy[i] == gradient[i]);
  } else {
    assert_int_equal(prop->hessian_count, 0);
    count = 0;
    assert_int_equal(cpmdc_session_last_gradient(session, copy, need, &count),
                     -1);
    assert_int_equal(count, need);
    for (size_t i = 0; i < need; ++i)
      assert_true(copy[i] == 0.0);
  }
  free(prop);
  free(copy);
}

static void check_plain_gradient_accessor(int n_atoms, const double *gradient) {
  const size_t need = (size_t)n_atoms * 3;
  double *copy = calloc(need, sizeof(double));
  assert_non_null(copy);
  size_t count = 0;
  if (need <= k_snapshot_components) {
    assert_int_equal(cpmdc_last_gradient(copy, need, &count), 0);
    assert_int_equal(count, need);
    for (size_t i = 0; i < need; ++i)
      assert_true(copy[i] == gradient[i]);
  } else {
    assert_int_equal(cpmdc_last_gradient(copy, need, &count), -1);
    assert_int_equal(count, need);
  }
  count = 0;
  assert_int_equal(cpmdc_last_gradient(copy, 0, &count), -1);
  assert_int_equal(count, need);
  free(copy);
}

static CPMDCSession *open_session_or_skip(void) {
  if (cpmdc_available())
    skip();
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  CPMDCSession *session = cpmdc_session_create(params, ps);
  free(params);
  assert_non_null(session);
  return session;
}

static size_t encode_big_force_input(const struct big_system *sys,
                                     unsigned char **out) {
  size_t cap = 4096 + (size_t)sys->n_atoms * 64;
  unsigned char *buf = malloc(cap);
  assert_non_null(buf);
  struct capn c;
  capn_init_malloc(&c);
  capn_ptr root = capn_root(&c);
  ForceInput_ptr fp = new_ForceInput(root.seg);
  struct ForceInput v;
  memset(&v, 0, sizeof(v));
  v.pos = capn_new_list64(root.seg, 3 * sys->n_atoms);
  v.atmnrs = capn_new_list32(root.seg, sys->n_atoms);
  v.box = capn_new_list64(root.seg, 0);
  for (int i = 0; i < 3 * sys->n_atoms; ++i)
    capn_set64(v.pos, i, capn_from_f64(sys->positions_ang[i]));
  for (int i = 0; i < sys->n_atoms; ++i)
    capn_set32(v.atmnrs, i, (uint32_t)sys->numbers[i]);
  v.lengthUnit = (capn_text){4, "bohr", NULL};
  v.energyUnit = (capn_text){7, "hartree", NULL};
  write_ForceInput(&v, fp);
  int len = capn_setp(root, 0, fp.p) == 0 ? capn_write_mem(&c, buf, cap, 0) : -1;
  capn_free(&c);
  assert_true(len > 0);
  *out = buf;
  return (size_t)len;
}

static const int k_system_sizes[] = {k_capacity_atoms, k_capacity_atoms + 1,
                                     2 * k_capacity_atoms + 7};
enum { k_system_count = 3 };

static void test_session_gradient_above_snapshot_capacity(void **state) {
  (void)state;
  for (int s = 0; s < k_system_count; ++s) {
    struct big_system sys;
    big_system_init(&sys, k_system_sizes[s]);
    CPMDCSession *session = open_session_or_skip();
    double *gradient = calloc((size_t)sys.n_atoms * 3, sizeof(double));
    assert_non_null(gradient);
    CPMDCResult r = cpmdc_session_energy_gradient(
        session, sys.n_atoms, sys.positions_ang, sys.numbers, gradient);
    assert_int_equal(r.ok, 1);
    check_reference_gradient(&sys, gradient);
    check_gradient_accessors(session, sys.n_atoms, gradient);
    check_plain_gradient_accessor(sys.n_atoms, gradient);
    free(gradient);
    cpmdc_session_destroy(session);
    big_system_free(&sys);
  }
}

static void test_global_gradient_above_snapshot_capacity(void **state) {
  (void)state;
  if (cpmdc_available())
    skip();
  size_t ps = 0;
  unsigned char *params = read_file(g_params, &ps);
  assert_non_null(params);
  for (int s = 0; s < k_system_count; ++s) {
    struct big_system sys;
    big_system_init(&sys, k_system_sizes[s]);
    double *gradient = calloc((size_t)sys.n_atoms * 3, sizeof(double));
    assert_non_null(gradient);
    CPMDCResult r = cpmdc_energy_gradient(sys.n_atoms, sys.positions_ang,
                                          sys.numbers, params, ps, gradient);
    assert_int_equal(r.ok, 1);
    check_reference_gradient(&sys, gradient);
    check_plain_gradient_accessor(sys.n_atoms, gradient);
    memset(gradient, 0, (size_t)sys.n_atoms * 3 * sizeof(double));
    r = cpmdc_energy_forces(sys.n_atoms, sys.positions_ang, sys.numbers,
                            params, ps, gradient);
    assert_int_equal(r.ok, 1);
    for (int i = 0; i < 3 * sys.n_atoms; ++i)
      assert_true(isfinite(gradient[i]));
    assert_true(gradient[3 * (sys.n_atoms - 1)] != 0.0);
    free(gradient);
    big_system_free(&sys);
  }
  free(params);
}

static void test_calculate_au_above_snapshot_capacity(void **state) {
  (void)state;
  const double bohr = 0.52917721092;
  for (int s = 0; s < k_system_count; ++s) {
    struct big_system sys;
    big_system_init(&sys, k_system_sizes[s]);
    CPMDCSession *session = open_session_or_skip();
    size_t n3 = (size_t)sys.n_atoms * 3;
    double *positions_bohr = malloc(n3 * sizeof(double));
    double *forces = calloc(n3, sizeof(double));
    double *gradient = calloc(n3, sizeof(double));
    assert_non_null(positions_bohr);
    assert_non_null(forces);
    assert_non_null(gradient);
    for (size_t i = 0; i < n3; ++i)
      positions_bohr[i] = sys.positions_ang[i] / bohr;
    double energy = 0.0;
    assert_int_equal(cpmdc_session_calculate_au(session, sys.n_atoms,
                                                positions_bohr, sys.numbers,
                                                NULL, &energy, forces, NULL),
                     0);
    assert_true(isfinite(energy));
    for (size_t i = 0; i < n3; ++i)
      gradient[i] = -forces[i];
    check_reference_gradient(&sys, gradient);
    check_gradient_accessors(session, sys.n_atoms, gradient);
    CPMDCDipole dipole;
    assert_int_equal(cpmdc_session_last_dipole(session, &dipole), -1);
    free(gradient);
    free(forces);
    free(positions_bohr);
    cpmdc_session_destroy(session);
    big_system_free(&sys);
  }
}

static void check_wire_gradient(const struct big_system *sys,
                                const unsigned char *out, size_t wrote) {
  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, out, (int)wrote, 0), 0);
  PotentialResult_ptr pr;
  pr.p = capn_getp(capn_root(&arena), 0, 1);
  struct PotentialResult view;
  read_PotentialResult(&view, pr);
  const int n3 = 3 * sys->n_atoms;
  assert_int_equal(capn_len(view.forces), n3);
  assert_int_equal(capn_len(view.gradient), n3);
  for (int i = 0; i < n3; ++i) {
    double f = capn_to_f64(capn_get64(view.forces, i));
    double g = capn_to_f64(capn_get64(view.gradient, i));
    assert_true(isfinite(g));
    assert_true(fabs(g + f) <= 1e-9 * fabs(f));
  }
  assert_true(capn_to_f64(capn_get64(view.gradient, n3 - 1)) != 0.0);
  capn_free(&arena);
}

static void test_result_calls_above_snapshot_capacity(void **state) {
  (void)state;
  for (int s = 0; s < k_system_count; ++s) {
    struct big_system sys;
    big_system_init(&sys, k_system_sizes[s]);
    unsigned char *step = NULL;
    size_t step_size = encode_big_force_input(&sys, &step);
    size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
    assert_true(need > 0);
    unsigned char *out = malloc(need);
    assert_non_null(out);

    CPMDCSession *session = open_session_or_skip();
    size_t wrote = 0;
    CPMDCResult r = cpmdc_session_calculate_result(session, step, step_size,
                                                   out, need, &wrote);
    assert_int_equal(r.ok, 1);
    check_wire_gradient(&sys, out, wrote);
    size_t count = 0;
    int have = cpmdc_session_last_gradient(session, NULL, 0, &count);
    assert_int_equal(have, -1);
    assert_int_equal(count, (size_t)sys.n_atoms * 3);
    cpmdc_session_destroy(session);

    size_t ps = 0;
    unsigned char *params = read_file(g_params, &ps);
    assert_non_null(params);
    memset(out, 0, need);
    wrote = 0;
    r = cpmdc_calculate_result(params, ps, step, step_size, out, need, &wrote);
    assert_int_equal(r.ok, 1);
    check_wire_gradient(&sys, out, wrote);
    free(params);
    free(out);
    free(step);
    big_system_free(&sys);
  }
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
      cmocka_unit_test(test_session_gradient_above_snapshot_capacity),
      cmocka_unit_test(test_global_gradient_above_snapshot_capacity),
      cmocka_unit_test(test_calculate_au_above_snapshot_capacity),
      cmocka_unit_test(test_result_calls_above_snapshot_capacity),
      cmocka_unit_test(test_all_results_through_c_and_wire),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
