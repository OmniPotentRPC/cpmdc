#include "cpmdc.h"

#include <stdio.h>

#ifndef CPMDC_VERSION_STRING
#define CPMDC_VERSION_STRING "unknown"
#endif

static _Thread_local char g_last_error[512];

const char *cpmdc_last_error(void) { return g_last_error; }

static void stub_store_error(void) {
  snprintf(g_last_error, sizeof(g_last_error),
           "OpenCPMD embed not available in cpmdc stub");
}
#include <string.h>

static CPMDCResult stub_fail(void) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  snprintf(r.message, sizeof(r.message),
           "CPMD embed not available in cpmdc stub");
  return r;
}

int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  stub_store_error();
  return -1;
}

int cpmdc_configure(const void *config_capnp,
                    size_t config_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  stub_store_error();
  return -1;
}

CPMDCSession *cpmdc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  stub_store_error();
  return NULL;
}

int cpmdc_session_configure(CPMDCSession *session, const void *config_capnp,
                            size_t config_capnp_size_bytes) {
  (void)session;
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  stub_store_error();
  return -1;
}

CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  double *grad_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)grad_h_bohr;
  return stub_fail();
}

CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                         const int *atomic_numbers, const void *params_capnp,
                         size_t params_capnp_size_bytes) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return stub_fail();
}

CPMDCResult cpmdc_energy_forces(int n_atoms, const double *positions_ang,
                                const int *atomic_numbers,
                                const void *params_capnp,
                                size_t params_capnp_size_bytes,
                                double *forces_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)forces_h_bohr;
  return stub_fail();
}

const char *cpmdc_version(void) { return "cpmdc-stub/" CPMDC_VERSION_STRING; }

int cpmdc_abi_version(void) { return CPMDC_ABI_VERSION; }

int cpmdc_available(void) { return 0; }

void cpmdc_finalize(void) {}

CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                   size_t params_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return NULL;
}

int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes) {
  (void)session;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return -1;
}

void cpmdc_session_destroy(CPMDCSession *session) { (void)session; }

CPMDCResult cpmdc_session_energy(CPMDCSession *session, int n_atoms,
                                 const double *positions_ang,
                                 const int *atomic_numbers) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  return stub_fail();
}

CPMDCResult cpmdc_session_energy_gradient(CPMDCSession *session, int n_atoms,
                                          const double *positions_ang,
                                          const int *atomic_numbers,
                                          double *grad_h_bohr) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)grad_h_bohr;
  return stub_fail();
}

CPMDCResult cpmdc_session_energy_forces(CPMDCSession *session, int n_atoms,
                                        const double *positions_ang,
                                        const int *atomic_numbers,
                                        double *forces_h_bohr) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)forces_h_bohr;
  return stub_fail();
}

CPMDCResult cpmdc_session_calculate_forces(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, double *forces_h_bohr,
    size_t forces_len) {
  (void)session;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)forces_h_bohr;
  (void)forces_len;
  return stub_fail();
}

CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  (void)session;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

CPMDCResult cpmdc_calculate_result(const void *params_capnp,
                                   size_t params_capnp_size_bytes,
                                   const void *force_input_capnp,
                                   size_t force_input_capnp_size_bytes,
                                   void *potential_result_capnp,
                                   size_t potential_result_capnp_capacity_bytes,
                                   size_t *potential_result_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

CPMDCResult cpmdc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  return 0;
}


int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out) {
  if (out)
    memset(out, 0, sizeof(*out));
  return -1;
}
int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out) {
  if (out)
    memset(out, 0, sizeof(*out));
  return -1;
}
int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out) {
  if (out)
    memset(out, 0, sizeof(*out));
  return -1;
}
int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out) {
  if (out)
    memset(out, 0, sizeof(*out));
  return -1;
}

int cpmdc_last_energy_components(CPMDCEnergyComponents *out) {
  if (out)
    memset(out, 0, sizeof(*out));
  return -1;
}
