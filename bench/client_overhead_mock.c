/* Copyright (c) cpmdc developers
 *
 * Stand-in for libcpmdc used by bench/client_overhead.py. It exports the
 * session entry points the Python package binds and does the least work
 * that still writes every output: one pass over the positions and one over
 * the forces. The time of a Python call into it is therefore the client
 * side of an evaluation, with no SCF and no atom count limit of the
 * reference evaluator.
 */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define MOCK_EXPORT __attribute__((visibility("default")))

#define MOCK_FORCE_SCALE 0.5
#define MOCK_STRESS_SCALE 0.25

typedef struct MockResult {
  int ok;
  double energy_h;
  char message[512];
} MockResult;

typedef struct MockStress {
  int valid;
  double values[9];
} MockStress;

typedef struct MockDipole {
  int valid;
  double values[3];
} MockDipole;

typedef struct MockSession {
  int calls;
} MockSession;

static MockSession g_session;
static MockStress g_stress;

static double fill_forces(size_t n3, const double *positions, double *forces) {
  double energy = 0.0;
  for (size_t i = 0; i < n3; ++i) {
    double p = positions ? positions[i] : (double)(i % 7u);
    forces[i] = -MOCK_FORCE_SCALE * p;
    energy += 0.5 * MOCK_FORCE_SCALE * p * p;
  }
  return energy;
}

static MockResult result_ok(double energy_h) {
  MockResult r;
  memset(&r, 0, sizeof(r));
  r.ok = 1;
  r.energy_h = energy_h;
  snprintf(r.message, sizeof(r.message), "ok");
  return r;
}

static void set_stress(int periodic) {
  g_stress.valid = periodic != 0;
  for (int i = 0; i < 9; ++i)
    g_stress.values[i] = periodic ? MOCK_STRESS_SCALE * (double)(i + 1) : 0.0;
}

MOCK_EXPORT void *cpmdc_session_create(const void *params, size_t size) {
  (void)params;
  (void)size;
  g_session.calls = 0;
  return &g_session;
}

MOCK_EXPORT void cpmdc_session_destroy(void *session) { (void)session; }

MOCK_EXPORT MockResult cpmdc_session_energy_forces(void *session, int n_atoms,
                                                   const double *positions,
                                                   const int *numbers,
                                                   double *forces) {
  (void)numbers;
  ((MockSession *)session)->calls++;
  set_stress(0);
  return result_ok(fill_forces((size_t)n_atoms * 3u, positions, forces));
}

MOCK_EXPORT MockResult cpmdc_session_calculate_forces(
    void *session, const void *force_input, size_t force_input_size,
    double *forces, size_t forces_len) {
  (void)force_input;
  (void)force_input_size;
  ((MockSession *)session)->calls++;
  set_stress(1);
  return result_ok(fill_forces(forces_len, NULL, forces));
}

MOCK_EXPORT int cpmdc_session_calculate_au(void *session, int n_atoms,
                                           const double *positions,
                                           const int *numbers,
                                           const double *cell, double *energy,
                                           double *forces, double *virial) {
  (void)numbers;
  ((MockSession *)session)->calls++;
  *energy = fill_forces((size_t)n_atoms * 3u, positions, forces);
  for (int i = 0; i < 9; ++i)
    virial[i] = cell ? MOCK_STRESS_SCALE * (double)(i + 1) : 0.0;
  set_stress(cell != NULL);
  return 0;
}

MOCK_EXPORT int cpmdc_session_select_orbitals(void *session, long long key) {
  (void)session;
  (void)key;
  return 0;
}

MOCK_EXPORT int cpmdc_adopt_calculator_comm(const void *comm, size_t bytes,
                                            int ranks) {
  (void)comm;
  (void)bytes;
  (void)ranks;
  return -1;
}

MOCK_EXPORT int cpmdc_session_last_stress(const void *session,
                                          MockStress *out) {
  (void)session;
  *out = g_stress;
  return g_stress.valid ? 0 : -1;
}

MOCK_EXPORT int cpmdc_session_last_dipole(const void *session,
                                          MockDipole *out) {
  (void)session;
  memset(out, 0, sizeof(*out));
  return -1;
}

MOCK_EXPORT const char *cpmdc_last_error(void) { return ""; }
MOCK_EXPORT const char *cpmdc_version(void) { return "cpmdc-bench-mock"; }
MOCK_EXPORT int cpmdc_abi_version(void) { return 0; }
MOCK_EXPORT int cpmdc_available(void) { return 0; }
