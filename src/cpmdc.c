/* C ABI frontend: Cap'n Proto sessions + unit carriers (nwchemc pattern).
 * Engine work is in Fortran bind(C) embed surface (cpmd_embed_c_api.F90). */
#include "cpmdc.h"
#include "cpmdc_params.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>

void cpmdc_stop_arm(void);
void cpmdc_stop_disarm(void);
int cpmdc_stop_code(void);
jmp_buf *cpmdc_stop_jmp(void);

int cpmdc_embed_init(void);
int cpmdc_embed_available(void);
int cpmdc_embed_set_config(const char *functional, int functional_len,
                           double cutoff_ry, int charge, int multiplicity,
                           const char *input_deck, int input_deck_len,
                           const char *cpmd_root, int cpmd_root_len);
int cpmdc_embed_set_deck(const char *deck, int deck_len);
int cpmdc_embed_energy_grad(int n_atoms, const double *positions_ang,
                            const int *atomic_numbers, const double *cell_ang,
                            int has_cell, double *energy_h,
                            double *grad_h_bohr);
void cpmdc_embed_finalize(void);

/* Last Cap'n Proto params bytes for geometry-aware deck render on eval. */
static unsigned char *g_params_bytes = NULL;
static size_t g_params_size = 0;

struct CPMDCSession {
  unsigned char *params_bytes;
  size_t params_size;
  char functional[64];
  double cutoff_ry;
  int charge;
  int multiplicity;
  char input_deck[CPMDC_BLOCKS];
  char cpmd_root[1024];
  int topology_fixed;
  size_t fixed_n_atoms;
  int *fixed_atomic_numbers;
  double *step_positions_ang;
  int *step_atomic_numbers;
  size_t step_atom_capacity;
  int embed_configured;
};

static int ensure_embed_init(void) {
  static int attempted = 0;
  static int ok = 0;
  if (!attempted) {
    attempted = 1;
    ok = cpmdc_embed_init() != 0;
  }
  return ok;
}

static CPMDCResult fail_msg(const char *msg) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  snprintf(r.message, sizeof(r.message), "%s", msg ? msg : "error");
  return r;
}

static int apply_params_buffer(const void *params_capnp, size_t params_size,
                               char *functional, size_t functional_size,
                               double *cutoff_ry, int *charge, int *multiplicity,
                               char *input_deck, size_t input_deck_size,
                               char *cpmd_root, size_t cpmd_root_size) {
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(params_capnp, params_size, &arena, &root) != 0)
    return -1;
  struct CPMDParams view;
  read_CPMDParams(&view, root);
  snprintf(functional, functional_size, "%s",
           cpmdc_params_text_or(view.functional, "BLYP"));
  *cutoff_ry = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  *charge = view.charge;
  *multiplicity = view.multiplicity > 0 ? view.multiplicity : 1;
  snprintf(cpmd_root, cpmd_root_size, "%s",
           cpmdc_params_text_or(view.cpmdRoot, ""));
  if (cpmdc_params_render_input_deck(root, input_deck, input_deck_size) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  cpmdc_params_release(&arena);
  return 0;
}

static int configure_embed_from_session(CPMDCSession *session) {
  if (!session)
    return -1;
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return -1;
  int ok = cpmdc_embed_set_config(
      session->functional, (int)strlen(session->functional), session->cutoff_ry,
      session->charge, session->multiplicity, session->input_deck,
      (int)strlen(session->input_deck), session->cpmd_root,
      (int)strlen(session->cpmd_root));
  session->embed_configured = ok != 0;
  return ok != 0 ? 0 : -1;
}

int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes) {
  char functional[64];
  double cutoff_ry = 70.0;
  int charge = 0;
  int multiplicity = 1;
  char input_deck[CPMDC_BLOCKS];
  char cpmd_root[1024];
  if (apply_params_buffer(params_capnp, params_capnp_size_bytes, functional,
                          sizeof(functional), &cutoff_ry, &charge, &multiplicity,
                          input_deck, sizeof(input_deck), cpmd_root,
                          sizeof(cpmd_root)) != 0)
    return -1;
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return -1;
  free(g_params_bytes);
  g_params_bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!g_params_bytes)
    return -1;
  memcpy(g_params_bytes, params_capnp, params_capnp_size_bytes);
  g_params_size = params_capnp_size_bytes;
  if (cpmdc_embed_set_config(functional, (int)strlen(functional), cutoff_ry,
                             charge, multiplicity, input_deck,
                             (int)strlen(input_deck), cpmd_root,
                             (int)strlen(cpmd_root)) == 0)
    return -1;
  /* Full rendered deck (method sections); geometry merged at energy_grad. */
  (void)cpmdc_embed_set_deck(input_deck, (int)strlen(input_deck));
  return 0;
}

/* Merge ForceInput geometry into Cap'n Proto–derived CPMD INPUT deck. */
static int push_geometry_deck(int n_atoms, const double *positions_ang,
                              const int *atomic_numbers, const double *cell_ang,
                              int has_cell) {
  if (!g_params_bytes || g_params_size == 0)
    return -1;
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(g_params_bytes, g_params_size, &arena, &root) != 0)
    return -1;
  char deck[CPMDC_BLOCKS];
  if (cpmdc_params_render_deck_with_geometry(root, n_atoms, positions_ang,
                                             atomic_numbers, cell_ang, has_cell,
                                             deck, sizeof(deck)) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  cpmdc_params_release(&arena);
  return cpmdc_embed_set_deck(deck, (int)strlen(deck)) != 0 ? 0 : -1;
}

static CPMDCResult energy_gradient_cell(int n_atoms, const double *positions_ang,
                                        const int *atomic_numbers,
                                        const double *cell_ang, int has_cell,
                                        double *grad_h_bohr) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  r.message[0] = '\0';
  if (n_atoms <= 0 || !positions_ang || !atomic_numbers || !grad_h_bohr) {
    snprintf(r.message, sizeof(r.message), "invalid arguments");
    return r;
  }
  if (!ensure_embed_init() || !cpmdc_embed_available()) {
    snprintf(r.message, sizeof(r.message), "CPMD embed not available");
    return r;
  }
  double cell[9] = {0};
  if (has_cell && cell_ang)
    memcpy(cell, cell_ang, sizeof(cell));
  double energy = 0.0;
  if (push_geometry_deck(n_atoms, positions_ang, atomic_numbers, cell,
                         has_cell ? 1 : 0) != 0) {
    snprintf(r.message, sizeof(r.message), "Cap'n Proto deck render failed");
    return r;
  }
  cpmdc_stop_arm();
  int ok;
  if (setjmp(*cpmdc_stop_jmp()) != 0) {
    cpmdc_stop_disarm();
    return fail_msg("CPMD stopgm during embed SCF");
  }
  ok = cpmdc_embed_energy_grad(n_atoms, positions_ang, atomic_numbers, cell,
                               has_cell ? 1 : 0, &energy, grad_h_bohr);
  cpmdc_stop_disarm();
  if (!ok) {
    snprintf(r.message, sizeof(r.message), "CPMD energy/gradient failed");
    return r;
  }
  r.ok = 1;
  r.energy_h = energy;
  snprintf(r.message, sizeof(r.message), "ok");
  return r;
}

CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  double *grad_h_bohr) {
  if (cpmdc_set_params(params_capnp, params_capnp_size_bytes) != 0)
    return fail_msg("embed config failed");
  return energy_gradient_cell(n_atoms, positions_ang, atomic_numbers, NULL, 0,
                              grad_h_bohr);
}

CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                         const int *atomic_numbers, const void *params_capnp,
                         size_t params_capnp_size_bytes) {
  double *grad = NULL;
  if (n_atoms > 0) {
    grad = (double *)calloc((size_t)n_atoms * 3u, sizeof(double));
    if (!grad)
      return fail_msg("out of memory");
  }
  CPMDCResult r = cpmdc_energy_gradient(n_atoms, positions_ang, atomic_numbers,
                                        params_capnp, params_capnp_size_bytes,
                                        grad);
  free(grad);
  return r;
}

CPMDCResult cpmdc_energy_forces(int n_atoms, const double *positions_ang,
                                const int *atomic_numbers,
                                const void *params_capnp,
                                size_t params_capnp_size_bytes,
                                double *forces_h_bohr) {
  CPMDCResult r = cpmdc_energy_gradient(n_atoms, positions_ang, atomic_numbers,
                                        params_capnp, params_capnp_size_bytes,
                                        forces_h_bohr);
  if (r.ok && forces_h_bohr) {
    for (int i = 0; i < n_atoms * 3; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
  }
  return r;
}

static int session_reserve_step_atoms(CPMDCSession *session, size_t n_atoms) {
  if (!session)
    return -1;
  if (session->step_atom_capacity >= n_atoms)
    return 0;
  double *pos = (double *)realloc(session->step_positions_ang,
                                  n_atoms * 3u * sizeof(double));
  int *z = (int *)realloc(session->step_atomic_numbers, n_atoms * sizeof(int));
  if (!pos || !z) {
    free(pos);
    free(z);
    return -1;
  }
  session->step_positions_ang = pos;
  session->step_atomic_numbers = z;
  session->step_atom_capacity = n_atoms;
  return 0;
}

static int session_accept_topology(CPMDCSession *session, size_t n_atoms,
                                   const int *atomic_numbers) {
  if (!session || !atomic_numbers || n_atoms == 0)
    return -1;
  if (!session->topology_fixed) {
    int *copy = (int *)malloc(n_atoms * sizeof(int));
    if (!copy)
      return -1;
    memcpy(copy, atomic_numbers, n_atoms * sizeof(int));
    session->fixed_atomic_numbers = copy;
    session->fixed_n_atoms = n_atoms;
    session->topology_fixed = 1;
    return 0;
  }
  if (session->fixed_n_atoms != n_atoms)
    return -1;
  for (size_t i = 0; i < n_atoms; ++i) {
    if (session->fixed_atomic_numbers[i] != atomic_numbers[i])
      return -1;
  }
  return 0;
}

CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                   size_t params_capnp_size_bytes) {
  if (!params_capnp || params_capnp_size_bytes == 0)
    return NULL;
  CPMDCSession *session = (CPMDCSession *)calloc(1, sizeof(*session));
  if (!session)
    return NULL;
  session->params_bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!session->params_bytes) {
    free(session);
    return NULL;
  }
  memcpy(session->params_bytes, params_capnp, params_capnp_size_bytes);
  session->params_size = params_capnp_size_bytes;
  if (apply_params_buffer(session->params_bytes, session->params_size,
                          session->functional, sizeof(session->functional),
                          &session->cutoff_ry, &session->charge,
                          &session->multiplicity, session->input_deck,
                          sizeof(session->input_deck), session->cpmd_root,
                          sizeof(session->cpmd_root)) != 0) {
    free(session->params_bytes);
    free(session);
    return NULL;
  }
  free(g_params_bytes);
  g_params_bytes = (unsigned char *)malloc(session->params_size);
  if (g_params_bytes) {
    memcpy(g_params_bytes, session->params_bytes, session->params_size);
    g_params_size = session->params_size;
  }
  (void)configure_embed_from_session(session);
  return session;
}

int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes) {
  if (!session || !params_capnp || params_capnp_size_bytes == 0)
    return -1;
  if (session->topology_fixed)
    return -1;
  unsigned char *bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!bytes)
    return -1;
  if (apply_params_buffer(params_capnp, params_capnp_size_bytes,
                          session->functional, sizeof(session->functional),
                          &session->cutoff_ry, &session->charge,
                          &session->multiplicity, session->input_deck,
                          sizeof(session->input_deck), session->cpmd_root,
                          sizeof(session->cpmd_root)) != 0) {
    free(bytes);
    return -1;
  }
  memcpy(bytes, params_capnp, params_capnp_size_bytes);
  free(session->params_bytes);
  session->params_bytes = bytes;
  session->params_size = params_capnp_size_bytes;
  return configure_embed_from_session(session);
}

void cpmdc_session_destroy(CPMDCSession *session) {
  if (!session)
    return;
  free(session->params_bytes);
  free(session->fixed_atomic_numbers);
  free(session->step_positions_ang);
  free(session->step_atomic_numbers);
  free(session);
}

static CPMDCResult session_energy_gradient_cell(
    CPMDCSession *session, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const double *cell_ang, int has_cell,
    double *grad_h_bohr) {
  if (!session)
    return fail_msg("null session");
  if (session_accept_topology(session, (size_t)n_atoms, atomic_numbers) != 0)
    return fail_msg("topology change requires a new session");
  if (!session->embed_configured && configure_embed_from_session(session) != 0)
    return fail_msg("CPMD embed not available");
  return energy_gradient_cell(n_atoms, positions_ang, atomic_numbers, cell_ang,
                              has_cell, grad_h_bohr);
}

CPMDCResult cpmdc_session_energy_gradient(CPMDCSession *session, int n_atoms,
                                          const double *positions_ang,
                                          const int *atomic_numbers,
                                          double *grad_h_bohr) {
  return session_energy_gradient_cell(session, n_atoms, positions_ang,
                                      atomic_numbers, NULL, 0, grad_h_bohr);
}

CPMDCResult cpmdc_session_energy(CPMDCSession *session, int n_atoms,
                                 const double *positions_ang,
                                 const int *atomic_numbers) {
  double *grad = NULL;
  if (n_atoms > 0) {
    grad = (double *)calloc((size_t)n_atoms * 3u, sizeof(double));
    if (!grad)
      return fail_msg("out of memory");
  }
  CPMDCResult r = cpmdc_session_energy_gradient(session, n_atoms, positions_ang,
                                                atomic_numbers, grad);
  free(grad);
  return r;
}

CPMDCResult cpmdc_session_energy_forces(CPMDCSession *session, int n_atoms,
                                        const double *positions_ang,
                                        const int *atomic_numbers,
                                        double *forces_h_bohr) {
  CPMDCResult r = cpmdc_session_energy_gradient(session, n_atoms, positions_ang,
                                                atomic_numbers, forces_h_bohr);
  if (r.ok && forces_h_bohr) {
    for (int i = 0; i < n_atoms * 3; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
  }
  return r;
}

CPMDCResult cpmdc_session_calculate_forces(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, double *forces_h_bohr,
    size_t forces_len) {
  CPMDCResult r = fail_msg("invalid arguments");
  if (!session || !force_input_capnp || force_input_capnp_size_bytes == 0 ||
      !forces_h_bohr)
    return r;
  struct capn arena;
  ForceInput_ptr force_input;
  if (cpmdc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                             &arena, &force_input) != 0)
    return fail_msg("invalid ForceInput message");
  size_t n_atoms = 0;
  int has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      n_atoms == 0 || forces_len < n_atoms * 3u) {
    cpmdc_params_release(&arena);
    return fail_msg("invalid ForceInput geometry");
  }
  if (session_reserve_step_atoms(session, n_atoms) != 0) {
    cpmdc_params_release(&arena);
    return fail_msg("out of memory");
  }
  double cell_ang[9];
  if (cpmdc_force_input_copy_geometry(
          force_input, session->step_positions_ang, session->step_atomic_numbers,
          session->step_atom_capacity, cell_ang, &has_cell) != 0) {
    cpmdc_params_release(&arena);
    return fail_msg("invalid ForceInput geometry");
  }
  cpmdc_params_release(&arena);
  r = session_energy_gradient_cell(session, (int)n_atoms,
                                   session->step_positions_ang,
                                   session->step_atomic_numbers, cell_ang,
                                   has_cell, forces_h_bohr);
  if (r.ok) {
    for (size_t i = 0; i < n_atoms * 3u; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
  }
  return r;
}

CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  r.message[0] = '\0';
  if (!session || !force_input_capnp || force_input_capnp_size_bytes == 0 ||
      !potential_result_capnp_size_bytes) {
    snprintf(r.message, sizeof(r.message), "invalid arguments");
    return r;
  }
  *potential_result_capnp_size_bytes = 0;

  struct capn arena;
  ForceInput_ptr force_input;
  if (cpmdc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                             &arena, &force_input) != 0) {
    snprintf(r.message, sizeof(r.message), "invalid ForceInput message");
    return r;
  }

  size_t n_atoms = 0;
  int has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      n_atoms > (size_t)INT_MAX || n_atoms > SIZE_MAX / 3u) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput geometry");
    return r;
  }
  size_t force_count = n_atoms * 3u;
  size_t required_size = cpmdc_potential_result_flat_size(force_count);
  *potential_result_capnp_size_bytes = required_size;
  if (required_size == 0 || force_count > (size_t)INT_MAX) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput geometry");
    return r;
  }

  double energy_factor = 1.0;
  double force_factor = 1.0;
  if (cpmdc_force_input_result_factors(force_input, &energy_factor,
                                       &force_factor) != 0) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput result units");
    return r;
  }
  if (!potential_result_capnp ||
      potential_result_capnp_capacity_bytes < required_size) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "PotentialResult buffer too small");
    return r;
  }
  if (session_reserve_step_atoms(session, n_atoms) != 0) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "out of memory");
    return r;
  }

  double cell_ang[9];
  if (cpmdc_force_input_copy_geometry(
          force_input, session->step_positions_ang, session->step_atomic_numbers,
          session->step_atom_capacity, cell_ang, &has_cell) != 0) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput geometry");
    return r;
  }
  cpmdc_params_release(&arena);

  double *forces = (double *)malloc(force_count * sizeof(*forces));
  if (!forces) {
    snprintf(r.message, sizeof(r.message), "out of memory");
    return r;
  }
  r = session_energy_gradient_cell(session, (int)n_atoms,
                                   session->step_positions_ang,
                                   session->step_atomic_numbers, cell_ang,
                                   has_cell, forces);
  if (r.ok) {
    for (size_t i = 0; i < force_count; ++i)
      forces[i] = -forces[i] * force_factor;
    if (cpmdc_potential_result_write(r.energy_h * energy_factor, forces,
                                     force_count, potential_result_capnp,
                                     potential_result_capnp_capacity_bytes,
                                     potential_result_capnp_size_bytes) != 0) {
      r.ok = 0;
      snprintf(r.message, sizeof(r.message), "PotentialResult write failed");
    }
  }
  free(forces);
  return r;
}

CPMDCResult cpmdc_calculate_result(const void *params_capnp,
                                   size_t params_capnp_size_bytes,
                                   const void *force_input_capnp,
                                   size_t force_input_capnp_size_bytes,
                                   void *potential_result_capnp,
                                   size_t potential_result_capnp_capacity_bytes,
                                   size_t *potential_result_capnp_size_bytes) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  r.message[0] = '\0';
  if (!params_capnp || params_capnp_size_bytes == 0 || !force_input_capnp ||
      force_input_capnp_size_bytes == 0 || !potential_result_capnp_size_bytes) {
    snprintf(r.message, sizeof(r.message), "invalid arguments");
    return r;
  }
  *potential_result_capnp_size_bytes = 0;

  size_t required_size = cpmdc_potential_result_size_for_force_input(
      force_input_capnp, force_input_capnp_size_bytes);
  *potential_result_capnp_size_bytes = required_size;
  if (required_size == 0) {
    snprintf(r.message, sizeof(r.message), "invalid ForceInput geometry");
    return r;
  }
  if (!potential_result_capnp ||
      potential_result_capnp_capacity_bytes < required_size) {
    snprintf(r.message, sizeof(r.message), "PotentialResult buffer too small");
    return r;
  }

  CPMDCSession *session =
      cpmdc_session_create(params_capnp, params_capnp_size_bytes);
  if (!session) {
    snprintf(r.message, sizeof(r.message), "embed config failed");
    return r;
  }
  r = cpmdc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  cpmdc_session_destroy(session);
  return r;
}

size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  struct capn arena;
  ForceInput_ptr force_input;
  if (cpmdc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                             &arena, &force_input) != 0)
    return 0;
  size_t n_atoms = 0;
  int has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      n_atoms > (size_t)INT_MAX || n_atoms > SIZE_MAX / 3u) {
    cpmdc_params_release(&arena);
    return 0;
  }
  (void)has_cell;
  size_t force_count = n_atoms * 3u;
  if (force_count > (size_t)INT_MAX) {
    cpmdc_params_release(&arena);
    return 0;
  }
  size_t result_size = cpmdc_potential_result_flat_size(force_count);
  cpmdc_params_release(&arena);
  return result_size;
}

const char *cpmdc_version(void) { return "cpmdc/0.1.0"; }

int cpmdc_available(void) {
  if (!ensure_embed_init())
    return 0;
  return cpmdc_embed_available() != 0;
}

void cpmdc_finalize(void) { cpmdc_embed_finalize(); }
