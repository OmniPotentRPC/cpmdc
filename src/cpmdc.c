/* C ABI frontend: Cap'n Proto sessions + unit carriers (nwchemc pattern).
 * Engine work is in Fortran bind(C) embed surface (cpmd_embed_c_api.F90). */
#include "cpmdc.h"
#include "cpmdc_adopt.h"
#include "cpmdc_cell_rank.h"
#include "cpmdc_embed_host.h"
#include "cpmdc_embed_image.h"
#include "cpmdc_params.h"
#include "cpmdc_scf_report.h"

#ifndef CPMDC_VERSION_STRING
#define CPMDC_VERSION_STRING "unknown"
#endif

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* LP64 layout callers already allocate. Charges are the next image field,
 * not bytes inside this struct. */
_Static_assert(sizeof(CPMDCPropertySnapshot) == 32896,
               "CPMDCPropertySnapshot stays the historical size");
_Static_assert(offsetof(CPMDCEmbedImage, charges) ==
                   offsetof(CPMDCEmbedImage, prop) + sizeof(CPMDCPropertySnapshot),
               "charges follow the public snapshot");
_Static_assert(offsetof(CPMDCEmbedImage, stress) ==
                   offsetof(CPMDCEmbedImage, charges) + sizeof(CPMDCEmbedCharges),
               "stress follows the charge record");
#if defined(CPMDC_HAVE_MPI)
#include <mpi.h>
#endif

static void cpmdc_store_error(const char *msg);
static void cpmdc_store_error_code(const char *msg, CPMDCErrorCode code);

void cpmdc_stop_arm(void);
void cpmdc_stop_disarm(void);
int cpmdc_stop_code(void);
void cpmdc_note_stop(int code);
int cpmdc_agree_cold_stop(void);

int cpmdc_embed_init(void);
int cpmdc_embed_available(void);
int cpmdc_embed_reset_state(CPMDCEmbedImage *image);
void cpmdc_orbital_clear_key(int drop_history);
void cpmdc_orbital_ledger_reset(void);
int cpmdc_embed_detach_image(CPMDCEmbedImage *image);
int cpmdc_embed_set_config(const char *functional, int functional_len,
                           double cutoff_ry, int charge, int multiplicity,
                           const char *input_deck, int input_deck_len,
                           const char *cpmd_root, int cpmd_root_len,
                           CPMDCEmbedImage *image);
int cpmdc_embed_set_deck(const char *deck, int deck_len, CPMDCEmbedImage *image);
/* Fortran capnp-fortran decode of CPMDParams into embed knobs. */
int cpmdc_embed_apply_params(const void *params_capnp,
                             size_t params_capnp_size_bytes,
                             const char *input_deck, int input_deck_len,
                             const char *functional_ov, int functional_ov_len,
                             double cutoff_ov, int has_charge_ov, int charge_ov,
                             int has_mult_ov, int mult_ov, char *functional_out,
                             int functional_cap, double *cutoff_out,
                             int *charge_out, int *mult_out, char *deck_out,
                             int deck_cap, char *root_out, int root_cap);
int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len, char *cpmd_root,
                           int cpmd_root_len);
int cpmdc_embed_compose_cold_deck_image(int n_atoms, const double *positions_ang,
                                        const int *atomic_numbers,
                                        const double *cell_ang, int has_cell,
                                        char *deck_out, int deck_cap,
                                        int *deck_len, CPMDCEmbedImage *image);
int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len);
int cpmdc_embed_energy_grad(int n_atoms, const double *positions_ang,
                            const int *atomic_numbers, const double *cell_ang,
                            int has_cell, double *energy_h,
                            double *grad_h_bohr, int n_external,
                            const double *external_positions_ang,
                            const double *external_charges_e,
                            const double *external_widths_ang,
                            double *external_forces_h_bohr,
                            CPMDCEmbedImage *image);
void cpmdc_embed_finalize(void);
int cpmdc_embed_copy_orbitals(long long from_key, long long to_key);
int cpmdc_embed_slot_limit_code(void);
int cpmdc_embed_drop_stale_orbitals(int n_atoms, const int *atomic_numbers,
                                    const double *cell_ang, int has_cell,
                                    CPMDCEmbedImage *image);
int cpmdc_orbital_agree(int local_failed, const char **err);
int cpmdc_orbital_allreduce(int *values, int count);
int cpmdc_orbital_key_agree(int *has_key, long long key, const char **err);
int cpmdc_orbital_prepare(long long key, const char **err);
void cpmdc_orbital_ledger_finalize(void);
void cpmdc_orbital_note_stored(void);
const char *cpmdc_orbital_limit_message(int code);
int cpmdc_embed_bind_calculator(int ranks_per_calc);
#if defined(CPMDC_HAS_CPMD)
int cpmdc_embed_store_comm(int fcomm);
int cpmdc_embed_adopted_fcomm(void);
#endif


/* Last Cap'n Proto params bytes for geometry-aware deck render on eval. */
static unsigned char *g_params_bytes = NULL;
static CPMDCScalarOverrides g_overrides;
static char g_overrides_functional[64];
static int g_has_overrides = 0;
static _Thread_local CPMDCErrorCode g_last_error_code = CPMDC_ERROR_NONE;
static _Thread_local char g_last_error[512];
static _Thread_local char g_embed_detail[512];
/* The code that belongs to g_embed_detail. A note from the Fortran side
 * replaces it when the failure has a more specific reason. */
static _Thread_local CPMDCErrorCode g_embed_code = CPMDC_ERROR_EVALUATION_FAILED;

static int cpmdc_set_params_ov(const void *params_capnp,
                               size_t params_capnp_size_bytes,
                               const CPMDCScalarOverrides *overrides);
static CPMDCSession *cpmdc_session_create_ov(
    const void *params_capnp, size_t params_capnp_size_bytes,
    const CPMDCScalarOverrides *overrides);
static int cpmdc_session_set_params_ov(CPMDCSession *session,
                                       const void *params_capnp,
                                       size_t params_capnp_size_bytes,
                                       const CPMDCScalarOverrides *overrides);
static size_t g_params_size = 0;

/** Persistent method state plus topology guards for session evaluations. */
struct CPMDCSession {
  /** Serialized unpacked flat `CPMDParams` bytes owned by the session. */
  unsigned char *params_bytes;
  /** Byte count of `params_bytes`. */
  size_t params_size;
  /** Effective DFT functional passed to the embed layer. */
  char functional[64];
  /** Effective plane-wave cutoff in Rydberg. */
  double cutoff_ry;
  /** Effective system charge. */
  int charge;
  /** Effective spin multiplicity. */
  int multiplicity;
  /** Rendered CPMD input deck for the session parameters. Owned. */
  char *input_deck;
  /** Common-overlay scalar overrides merged over the params (may be unset). */
  CPMDCScalarOverrides overrides;
  char overrides_functional[64];
  int has_overrides;
  /** cpmdc_session_set_constant_cutoff() E0, A, and sigma in Ry; E0 0 => off. */
  double constant_cutoff_ry;
  double constant_cutoff_height_ry;
  double constant_cutoff_width_ry;
  /** OpenCPMD source or build tree hint. */
  char cpmd_root[1024];
  /** Non-zero once the session has accepted one topology. */
  int topology_fixed;
  /** Atom count accepted by the first successful topology check. */
  size_t fixed_n_atoms;
  /** Atomic numbers from the latest accepted step. Identity is the element counts. */
  int *fixed_atomic_numbers;
  /** Scratch positions for the current step in Angstrom. */
  double *step_positions_ang;
  /** Scratch atomic numbers for the current step. */
  int *step_atomic_numbers;
  /** Allocated atom capacity for step scratch buffers. */
  size_t step_atom_capacity;
  /** Gaussian-smeared external MM charges for electrostatic embedding. */
  size_t external_count;
  double *external_positions_ang;
  double *external_charges_e;
  double *external_widths_ang;
  double *external_forces_h_bohr;
  int external_forces_valid;
  /** Non-zero when the last successful gradient is stored. */
  int result_cached;
  /** Positions, atomic numbers, and gradient of that call. Owned. */
  double *cached_positions_ang;
  int *cached_atomic_numbers;
  double *cached_grad_h_bohr;
  /** Hirshfeld charges of that call, caller order, when valid. Owned. */
  double *cached_charges_e;
  int cached_charges_valid;
  size_t cached_charge_count;
  size_t cached_n_atoms;
  double cached_energy_h;
  int cached_has_cell;
  double cached_cell_ang[9];
  /** Non-zero when the embed layer has accepted the effective config. */
  int embed_configured;
  /** Non-zero once the session has started an evaluation. */
  int evaluated;
  /** Non-zero once cpmdc_session_select_orbitals() named a key. */
  int has_orbital_key;
  /** The calculation whose stored orbitals the next SCF starts from. */
  long long orbital_key;
  /** Result of the last evaluation of this session. */
  CPMDCEmbedImage image;
};

static CPMDCSession *g_active_session = NULL;
/* The session whose evaluation last selected or cleared the orbital key.
 * It owns the selected key and its history, which configuring another
 * session does not change. */
static CPMDCSession *g_orbital_owner = NULL;
static int g_live_sessions = 0;
/* Live sessions that have started an evaluation. Only these tie CPMD state
 * to the stored communicator, so only these block an adopt. */
static int g_evaluated_sessions = 0;
static CPMDCEmbedImage g_image;

static CPMDCEmbedImage *active_image(void) {
  return g_active_session ? &g_active_session->image : &g_image;
}

#if defined(CPMDC_ENABLE_TEST_FAULTS)
static int test_fault(const char *name) {
  const char *fault = getenv("CPMDC_TEST_FAULT");
  return fault && strcmp(fault, name) == 0;
}
#else
#define test_fault(name) 0
#endif

static void copy_knob_text(char *dst, int cap, const char *src, size_t src_cap) {
  int i;
  int n;
  if (!dst || cap <= 0)
    return;
  n = cap - 1;
  if ((size_t)n > src_cap)
    n = (int)src_cap;
  for (i = 0; i < n && src[i] != '\0'; i++)
    dst[i] = src[i];
  dst[i] = '\0';
}

int cpmdc_embed_image_store_deck(CPMDCEmbedImage *image, const char *text,
                                 int len) {
  size_t n = 0;
  char *buf = NULL;
  if (!image || len < 0)
    return -1;
  if (len > 0 && !text)
    return -1;
  if (len > 0) {
    n = (size_t)len;
    for (size_t i = 0; i < n; ++i) {
      if (text[i] == '\0') {
        n = i;
        break;
      }
    }
  }
  if (n > 0) {
    buf = (char *)malloc(n + 1u);
    if (!buf)
      return -1;
    memcpy(buf, text, n);
    buf[n] = '\0';
  }
  free(image->input_deck);
  image->input_deck = buf;
  image->input_deck_len = (int)n;
  return 0;
}

int cpmdc_embed_get_config(char *functional, int functional_len, double *cutoff_ry,
                           int *charge, int *mult, char *input_deck,
                           int input_deck_len, char *cpmd_root, int cpmd_root_len) {
  const CPMDCEmbedImage *image = active_image();
  int deck_n;
  if (!image || !image->cfg_set || !cutoff_ry || !charge || !mult)
    return -1;
  deck_n = image->input_deck ? image->input_deck_len : 0;
  if (deck_n < 0)
    return -1;
  /* A short caller buffer is an error. The stored deck stays intact. */
  if (deck_n > 0 && (!input_deck || input_deck_len <= deck_n))
    return -1;
  if (input_deck && input_deck_len > 0) {
    if (deck_n > 0)
      memcpy(input_deck, image->input_deck, (size_t)deck_n);
    input_deck[deck_n] = '\0';
  }
  copy_knob_text(functional, functional_len, image->functional,
                 sizeof(image->functional));
  *cutoff_ry = image->cutoff_ry;
  *charge = image->cfg_charge;
  *mult = image->multiplicity;
  copy_knob_text(cpmd_root, cpmd_root_len, image->cpmd_root,
                 sizeof(image->cpmd_root));
  return 0;
}

static int ensure_embed_init(void) {
  static int attempted = 0;
  static int ok = 0;
  if (!attempted) {
    attempted = 1;
    ok = cpmdc_embed_init() != 0;
  }
  return ok;
}

static int cpmd_text_equals_ci(capn_text text, const char *lit);
static int cpmd_raw_has_section(capn_text block, const char *name);
static int cpmd_directives_have_kw(CPMDDirective_list directives,
                                   const char *keyword);
static int cpmd_set_matches_kw(capn_text key, const char *section,
                               const char *keyword);
static int cpmd_params_method_check_opaque_scalars(CPMDParams_ptr params);

static CPMDCResult fail_msg_code(const char *msg, CPMDCErrorCode code) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  snprintf(r.message, sizeof(r.message), "%s", msg ? msg : "error");
  cpmdc_store_error_code(r.message, code);
  return r;
}

static const char k_deck_external_potential_msg[] =
    "session external charges cannot be combined with a deck-defined "
    "EXTERNAL POTENTIAL";
static const char k_periodic_embedding_msg[] =
    "periodic electrostatic embedding is unsupported; use an isolated "
    "QM cell or clear the external charges";

/* Reasons that the evaluation preflight and its collective gates hand on as
 * plain text, so that every rank reports the same code for the same text. */
static CPMDCErrorCode code_for_gate_message(const char *msg) {
  static const struct {
    const char *prefix;
    CPMDCErrorCode code;
  } k_gate_codes[] = {
      {"null session", CPMDC_ERROR_INVALID_ARGUMENT},
      {"invalid arguments", CPMDC_ERROR_INVALID_ARGUMENT},
      {"atomic numbers must be", CPMDC_ERROR_INVALID_ARGUMENT},
      {"cell ", CPMDC_ERROR_INVALID_ARGUMENT},
      {"positions contain a non-finite", CPMDC_ERROR_INVALID_ARGUMENT},
      {"topology change requires", CPMDC_ERROR_INVALID_ARGUMENT},
      {"property gradient snapshot supports", CPMDC_ERROR_INVALID_ARGUMENT},
      {"CPMD embed not available", CPMDC_ERROR_UNAVAILABLE},
      {k_periodic_embedding_msg,
       CPMDC_ERROR_UNSUPPORTED_PERIODIC_EMBEDDING},
      {k_deck_external_potential_msg, CPMDC_ERROR_DECK_EXTERNAL_POTENTIAL},
  };
  if (!msg)
    return CPMDC_ERROR_EVALUATION_FAILED;
  for (size_t i = 0; i < sizeof(k_gate_codes) / sizeof(k_gate_codes[0]); ++i) {
    if (strncmp(msg, k_gate_codes[i].prefix, strlen(k_gate_codes[i].prefix)) ==
        0)
      return k_gate_codes[i].code;
  }
  return CPMDC_ERROR_EVALUATION_FAILED;
}

/* A message that an earlier call stored through cpmdc_store_error() keeps the
 * code it was stored with. Embed detail keeps its own code. Anything else is
 * an evaluation failure unless the gate text names a reason. */
static CPMDCResult fail_msg(const char *msg) {
  CPMDCErrorCode code = code_for_gate_message(msg);
  if (msg && g_last_error[0] != '\0' && strcmp(msg, g_last_error) == 0 &&
      g_last_error_code != CPMDC_ERROR_NONE)
    code = g_last_error_code;
  else if (msg && g_embed_detail[0] != '\0' &&
           strcmp(msg, g_embed_detail) == 0)
    code = g_embed_code;
  return fail_msg_code(msg, code);
}

/* CPMDC_DECK_OUT names a file that receives each method deck cpmdc renders
 * from a message, overwritten per render. The deck the OpenCPMD path composes
 * from it with the step's geometry is not written. */
static void write_deck_if_requested(const char *deck) {
  const char *path = getenv("CPMDC_DECK_OUT");
  if (!path || !*path || !deck)
    return;
  FILE *fp = fopen(path, "w");
  if (!fp)
    return;
  fputs(deck, fp);
  fclose(fp);
}

int cpmdc_params_method(const void *params_capnp, size_t params_size,
                        CPMDCMethodInfo *out) {
  struct capn arena;
  CPMDParams_ptr root;
  double cutoff_ry;
  if (!out) {
    cpmdc_store_error("CPMDParams method output is missing");
    return -1;
  }
  memset(out, 0, sizeof(*out));
  if (!params_capnp || params_size == 0) {
    cpmdc_store_error("CPMDParams buffer is empty");
    return -1;
  }
  if (cpmdc_params_root(params_capnp, params_size, &arena, &root) != 0) {
    cpmdc_store_error("CPMDParams parse failed");
    return -1;
  }
  if (cpmd_params_method_check_opaque_scalars(root) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  if (cpmdc_params_effective_config(root, out->functional,
                                    sizeof(out->functional), &cutoff_ry,
                                    &out->charge, &out->multiplicity) != 0) {
    cpmdc_params_release(&arena);
    cpmdc_store_error("CPMDParams configuration failed");
    return -1;
  }
  cpmdc_params_release(&arena);
  cpmdc_store_error("");
  return 0;
}

static int apply_params_buffer(const void *params_capnp, size_t params_size,
                               const CPMDCScalarOverrides *overrides,
                               char *functional, size_t functional_size,
                               double *cutoff_ry, int *charge, int *multiplicity,
                               char **input_deck, char *cpmd_root,
                               size_t cpmd_root_size) {
  struct capn arena;
  CPMDParams_ptr root;
  char *deck = NULL;
  if (input_deck)
    *input_deck = NULL;
  if (!input_deck) {
    cpmdc_store_error("CPMDParams deck buffer is missing");
    return -1;
  }
  if (!params_capnp || params_size == 0) {
    cpmdc_store_error("CPMDParams buffer is empty");
    return -1;
  }
  if (cpmdc_params_root(params_capnp, params_size, &arena, &root) != 0) {
    cpmdc_store_error("CPMDParams parse failed");
    return -1;
  }
  struct CPMDParams view;
  read_CPMDParams(&view, root);
  if (cpmdc_params_effective_config_ov(root, overrides, functional,
                                       functional_size, cutoff_ry, charge,
                                       multiplicity) != 0) {
    cpmdc_params_release(&arena);
    cpmdc_store_error("CPMDParams configuration failed");
    return -1;
  }
  snprintf(cpmd_root, cpmd_root_size, "%s",
           cpmdc_params_text_or(view.cpmdRoot, ""));
  if (cpmdc_params_render_input_deck_alloc(root, overrides, &deck) != 0) {
    cpmdc_params_release(&arena);
    const char *reason = cpmdc_params_render_reason();
    if (reason[0] != '\0') {
      char msg[512];
      snprintf(msg, sizeof(msg), "CPMDParams deck render failed: %s", reason);
      cpmdc_store_error(msg);
    } else {
      cpmdc_store_error("CPMDParams deck render failed");
    }
    return -1;
  }
  write_deck_if_requested(deck);
  if (cpmdc_params_reject_unsupported_inputs(functional, deck) != 0) {
    free(deck);
    cpmdc_params_release(&arena);
    cpmdc_store_error("CPMDParams input is not supported");
    return -1;
  }
  const char *conflict = cpmdc_params_deck_conflict(deck);
  if (conflict) {
    char msg[512];
    snprintf(msg, sizeof(msg), "CPMDParams input is not supported: %s",
             conflict);
    free(deck);
    cpmdc_params_release(&arena);
    cpmdc_store_error(msg);
    return -1;
  }
  cpmdc_params_release(&arena);
  *input_deck = deck;
  return 0;
}

static int fit_c_int(size_t n, int *out) {
  if (!out || n > (size_t)INT_MAX)
    return -1;
  *out = (int)n;
  return 0;
}

/* INT32-C: the gradient length is a signed atom count times three. */
static int grad_len_from_atoms(int n_atoms, size_t *out) {
  if (n_atoms <= 0 || n_atoms > INT_MAX / 3)
    return -1;
  if (out)
    *out = (size_t)n_atoms * 3u;
  return 0;
}

static int embed_apply_from_wire(const void *params_capnp, size_t params_size,
                                 const char *input_deck,
                                 const CPMDCScalarOverrides *overrides,
                                 CPMDCEmbedImage *image) {
  /* The embed records the configuration that cpmdc_params_effective_config_ov
   * resolves, the same one the deck was rendered from. The Fortran decode
   * receives it as the overlay, which it applies last, so it cannot record
   * a value the deck does not carry. */
  char fov[64];
  int fov_len = 0;
  double cutoff_ov = 0.0;
  int charge_ov = 0;
  int mult_ov = 1;
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(params_capnp, params_size, &arena, &root) != 0)
    return -1;
  int resolved = cpmdc_params_effective_config_ov(
      root, overrides, fov, sizeof(fov), &cutoff_ov, &charge_ov, &mult_ov);
  cpmdc_params_release(&arena);
  if (resolved != 0 || fit_c_int(strlen(fov), &fov_len) != 0)
    return -1;
  const int has_charge_ov = 1;
  const int has_mult_ov = 1;
  image->output_dir[0] = '\0';
  /* Same method on a new session keeps the process copy of c0. A
   * different basis drops it on the next force call. A stopgm still
   * clears that copy immediately. */
  if (cpmdc_embed_detach_image(image) == 0)
    return -1;
  int deck_len = 0;
  if (!input_deck || fit_c_int(strlen(input_deck), &deck_len) != 0)
    return -1;
  if (deck_len == INT_MAX)
    return -1;
  int deck_cap = deck_len + 1;
  char *deck_buf = (char *)malloc((size_t)deck_cap);
  if (!deck_buf)
    return -1;
  if (cpmdc_embed_apply_params(
          params_capnp, params_size, input_deck, deck_len, fov, fov_len,
          cutoff_ov, has_charge_ov, charge_ov, has_mult_ov, mult_ov,
          image->functional, 64, &image->cutoff_ry, &image->cfg_charge,
          &image->multiplicity, deck_buf, deck_cap, image->cpmd_root,
          1024) != 0) {
    free(deck_buf);
    return -1;
  }
  free(image->input_deck);
  image->input_deck = deck_buf;
  image->input_deck_len = (int)strlen(deck_buf);
  (void)cpmdc_params_copy_output_dir(params_capnp, params_size, image->output_dir,
                                     sizeof(image->output_dir));
  image->cfg_set = 1;
  return 0;
}

static int configure_embed_from_session(CPMDCSession *session) {
  if (!session)
    return -1;
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return -1;
  /* Serialized params: Fortran capnp-fortran decode; C still rendered deck. */
  if (embed_apply_from_wire(
          session->params_bytes, session->params_size, session->input_deck,
          session->has_overrides ? &session->overrides : NULL,
          &session->image) != 0)
    return -1;
  session->embed_configured = 1;
  g_active_session = session;
  return 0;
}

static void cpmdc_store_error_code(const char *msg, CPMDCErrorCode code) {
  snprintf(g_last_error, sizeof(g_last_error), "%s", msg ? msg : "");
  g_last_error_code = g_last_error[0] == '\0' ? CPMDC_ERROR_NONE : code;
}

static void cpmdc_store_error(const char *msg) {
  cpmdc_store_error_code(msg, CPMDC_ERROR_INVALID_ARGUMENT);
}

void cpmdc_note_embed_failure(const char *msg) {
  snprintf(g_embed_detail, sizeof(g_embed_detail), "%s", msg ? msg : "");
  g_embed_code = CPMDC_ERROR_EVALUATION_FAILED;
}

int cpmdc_embed_failure_pending(void) { return g_embed_detail[0] != '\0'; }

void cpmdc_note_deck_external_potential(void) {
  cpmdc_note_embed_failure(k_deck_external_potential_msg);
  g_embed_code = CPMDC_ERROR_DECK_EXTERNAL_POTENTIAL;
}

void cpmdc_note_periodic_embedding(void) {
  cpmdc_note_embed_failure(k_periodic_embedding_msg);
  g_embed_code = CPMDC_ERROR_UNSUPPORTED_PERIODIC_EMBEDDING;
}

void cpmdc_note_scf_unconverged(int minimizer, int iterations, int maxiter,
                                double gemax, double tolog,
                                int first_minimizer, int first_iterations) {
  (void)cpmdc_format_scf_unconverged(g_embed_detail, sizeof(g_embed_detail),
                                     minimizer, iterations, maxiter, gemax,
                                     tolog, first_minimizer, first_iterations);
  g_embed_code = CPMDC_ERROR_EVALUATION_FAILED;
}

const char *cpmdc_last_error(void) { return g_last_error; }
CPMDCErrorCode cpmdc_last_error_code(void) { return g_last_error_code; }

int cpmdc_embed_compose_cold_deck(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const double *cell_ang, int has_cell,
                                  char *deck_out, int deck_cap, int *deck_len) {
  g_embed_detail[0] = '\0';
  int ok = cpmdc_embed_compose_cold_deck_image(
      n_atoms, positions_ang, atomic_numbers, cell_ang, has_cell, deck_out,
      deck_cap, deck_len, active_image());
  if (!ok)
    cpmdc_store_error(g_embed_detail[0] ? g_embed_detail
                                        : "cold deck composition failed");
  return ok;
}

/* Map CommonMethodSpec.xcFunctionals to a CPMD functional token. Matched on
 * the sorted '+'-joined uppercase list. */
static const struct {
  const char *key;
  const char *token;
} k_cpmd_common_xc_map[] = {
    {"PBE", "PBE"},
    {"GGA_C_PBE+GGA_X_PBE", "PBE"},
    {"PBE0", "PBE0"},
    {"HYB_GGA_XC_PBEH", "PBE0"},
    {"BLYP", "BLYP"},
    {"GGA_C_LYP+GGA_X_B88", "BLYP"},
    {"B3LYP", "B3LYP"},
    {"HYB_GGA_XC_B3LYP", "B3LYP"},
    {"LDA", "LDA"},
    {"LDA_C_VWN+LDA_X", "LDA"},
};

static const capn_text k_empty_text = {0, "", 0};

static int cpmd_text_equals_ci(capn_text text, const char *lit) {
  int n = (int)strlen(lit);
  if (text.len != n || !text.str)
    return 0;
  for (int i = 0; i < n; ++i) {
    char a = text.str[i];
    char b = lit[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (b >= 'a' && b <= 'z')
      b = (char)(b - ('a' - 'A'));
    if (a != b)
      return 0;
  }
  return 1;
}

static int common_text_list_len(capn_ptr ptr) {
  capn_resolve(&ptr);
  if (ptr.type == CAPN_NULL)
    return 0;
  if (ptr.type != CAPN_PTR_LIST)
    return -1;
  return ptr.len;
}

/* Resolves the caller's list in place: capn_get32 reads 0 from a list
 * pointer that has not been resolved. */
static int common_list32_len(capn_list32 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST || list->p.datasz != 4)
    return -1;
  return list->p.len;
}

static int cpmd_common_xc_token(capn_ptr list, char *out, size_t out_size) {
  int n = common_text_list_len(list);
  if (n <= 0 || n > 4)
    return -1;
  char names[4][64];
  for (int i = 0; i < n; ++i) {
    capn_text entry = capn_get_text(list, i, k_empty_text);
    if (!entry.str || entry.len <= 0 || (size_t)entry.len >= sizeof(names[0]))
      return -1;
    for (int j = 0; j < entry.len; ++j) {
      char ch = entry.str[j];
      names[i][j] = (char)((ch >= 'a' && ch <= 'z') ? ch - ('a' - 'A') : ch);
    }
    names[i][entry.len] = '\0';
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (strcmp(names[i], names[j]) > 0) {
        char tmp[64];
        memcpy(tmp, names[i], sizeof(tmp));
        memcpy(names[i], names[j], sizeof(names[j]));
        memcpy(names[j], tmp, sizeof(tmp));
      }
    }
  }
  char key[280] = "";
  size_t used = 0;
  for (int i = 0; i < n; ++i) {
    int wrote = snprintf(key + used, sizeof(key) - used, "%s%s",
                         i > 0 ? "+" : "", names[i]);
    if (wrote < 0 || (size_t)wrote >= sizeof(key) - used)
      return -1;
    used += (size_t)wrote;
  }
  for (size_t i = 0;
       i < sizeof(k_cpmd_common_xc_map) / sizeof(k_cpmd_common_xc_map[0]);
       ++i) {
    if (strcmp(k_cpmd_common_xc_map[i].key, key) == 0) {
      snprintf(out, out_size, "%s", k_cpmd_common_xc_map[i].token);
      return 0;
    }
  }
  return -1;
}

#define CPMDC_RYDBERG_EV 13.605693122994
#define CPMDC_HARTREE_EV 27.211386245988
#define CPMDC_BOLTZMANN_EV_K 8.617333262e-5

static int cpmd_common_validate_basis_set(capn_text basis_set) {
  if (basis_set.len <= 0)
    return 0;
  if (!cpmd_text_equals_ci(basis_set, "planewave")) {
    cpmdc_store_error(
        "common overlay: CPMD is plane-wave only; basisSet must be "
        "\"planewave\" or unset");
    return -1;
  }
  return 0;
}

static int cpmd_common_validate_smearing(struct CommonMethodSpec *c,
                                         double *temperature_k) {
  *temperature_k = 0.0;
  capn_resolve(&c->smearing.p);
  if (c->smearing.p.type == CAPN_NULL)
    return 0;
  if (c->smearing.p.type != CAPN_STRUCT) {
    cpmdc_store_error("common overlay: invalid smearing");
    return -1;
  }
  struct CommonMethodSpec_Smearing smear;
  memset(&smear, 0, sizeof(smear));
  read_CommonMethodSpec_Smearing(&smear, c->smearing);
  if (smear.kind == CommonMethodSpec_Smearing_Kind_fermi) {
    /* CPMD free-energy functional with electron temperature; kT = width. */
    if (!isfinite(smear.widthEv) || smear.widthEv <= 0.0) {
      cpmdc_store_error("common overlay: fermi smearing needs widthEv > 0");
      return -1;
    }
    *temperature_k = smear.widthEv / CPMDC_BOLTZMANN_EV_K;
    return 0;
  }
  if (smear.kind == CommonMethodSpec_Smearing_Kind_none) {
    if (!isfinite(smear.widthEv) || smear.widthEv != 0.0) {
      cpmdc_store_error("common overlay: smearing widthEv requires kind = fermi");
      return -1;
    }
    return 0;
  }
  cpmdc_store_error(
      "common overlay: only fermi smearing has a CPMD lowering "
      "(free-energy functional)");
  return -1;
}

static int cpmd_common_validate_scalar_ranges(const struct CommonMethodSpec *c) {
  if (!isfinite(c->planewaveCutoffEv) || c->planewaveCutoffEv < 0.0) {
    cpmdc_store_error(
        "common overlay: planewaveCutoffEv must be finite and >= 0");
    return -1;
  }
  if (c->spinMultiplicity < 0) {
    cpmdc_store_error("common overlay: spinMultiplicity must be >= 0");
    return -1;
  }
  if (!isfinite(c->scfEnergyToleranceEv) || c->scfEnergyToleranceEv < 0.0) {
    cpmdc_store_error(
        "common overlay: scfEnergyToleranceEv must be finite and >= 0");
    return -1;
  }
  if (c->scfMaxIterations < 0) {
    cpmdc_store_error("common overlay: scfMaxIterations must be >= 0");
    return -1;
  }
  if (!isfinite(c->vanDerWaalsS6) || c->vanDerWaalsS6 < 0.0) {
    cpmdc_store_error("common overlay: vanDerWaalsS6 must be finite and >= 0");
    return -1;
  }
  return 0;
}

/* Lower the normalized overlay into a synthesized CPMDParams message.
 * Slice: functional, basis set, plane-wave cutoff, charge, multiplicity,
 * SCF tolerance/iterations, kMesh, fermi smearing, vdW method and s6; every
 * other set field is rejected loudly. Returns malloc'd flat bytes. */
static int cpmd_common_to_params(CommonMethodSpec_ptr common_root,
                                 unsigned char **out, size_t *out_size) {
  *out = NULL;
  *out_size = 0;
  struct CommonMethodSpec c;
  memset(&c, 0, sizeof(c));
  read_CommonMethodSpec(&c, common_root);

  if (cpmd_common_validate_scalar_ranges(&c) != 0)
    return -1;

  const char *grimme_token = NULL;
  if (c.vanDerWaalsMethod.len > 0) {
    if (cpmd_text_equals_ci(c.vanDerWaalsMethod, "DFT-D2"))
      grimme_token = "D2";
    else if (cpmd_text_equals_ci(c.vanDerWaalsMethod, "DFT-D3"))
      grimme_token = "D3";
    else if (cpmd_text_equals_ci(c.vanDerWaalsMethod, "DFT-D3(BJ)"))
      grimme_token = "D3BJ";
    if (!grimme_token) {
      cpmdc_store_error(
          "common overlay: unmapped vanDerWaalsMethod (use DFT-D2, DFT-D3, "
          "or DFT-D3(BJ))");
      return -1;
    }
  }
  if (c.vanDerWaalsS6 > 0.0 && !grimme_token) {
    cpmdc_store_error(
        "common overlay: vanDerWaalsS6 requires vanDerWaalsMethod");
    return -1;
  }
  if (c.relativityMethod.len > 0) {
    cpmdc_store_error(
        "common overlay: relativityMethod has no CPMD lowering yet");
    return -1;
  }
  double smear_temperature_k = 0.0;
  if (cpmd_common_validate_smearing(&c, &smear_temperature_k) != 0)
    return -1;
  if (cpmd_common_validate_basis_set(c.basisSet) != 0)
    return -1;

  char functional[64] = "";
  int xc_len = common_text_list_len(c.xcFunctionals);
  if (xc_len < 0) {
    cpmdc_store_error("common overlay: invalid xcFunctionals list");
    return -1;
  }
  if (xc_len > 0) {
    if (cpmd_common_xc_token(c.xcFunctionals, functional,
                             sizeof(functional)) != 0) {
      cpmdc_store_error("common overlay: unmapped xcFunctionals for CPMD");
      return -1;
    }
  }

  int kmesh_len = common_list32_len(&c.kMesh);
  if (kmesh_len != 0 && kmesh_len != 3) {
    cpmdc_store_error("common overlay: kMesh must have 3 Monkhorst divisions");
    return -1;
  }
  int want_cpmd_section = c.scfMaxIterations > 0 ||
                          c.scfEnergyToleranceEv > 0.0 ||
                          smear_temperature_k > 0.0;
  int nsections = (want_cpmd_section ? 1 : 0) + (kmesh_len == 3 ? 1 : 0) +
                  (grimme_token ? 1 : 0);

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  CPMDParams_ptr params = new_CPMDParams(root.seg);
  struct CPMDParams view;
  memset(&view, 0, sizeof(view));
  if (functional[0] != '\0') {
    view.functional.str = functional;
    view.functional.len = (int)strlen(functional);
  }
  view.cutOffRy = c.planewaveCutoffEv > 0.0
                      ? c.planewaveCutoffEv / CPMDC_RYDBERG_EV
                      : 70.0;
  view.charge = c.charge;
  view.hasCharge = (c.hasCharge || c.charge != 0) ? 1 : 0;
  view.multiplicity = c.spinMultiplicity > 0 ? c.spinMultiplicity : 1;
  view.task.str = "gradient";
  view.task.len = 8;
  if (nsections > 0) {
    CPMDInputSection_list sections =
        new_CPMDInputSection_list(root.seg, nsections);
    int idx = 0;
    if (want_cpmd_section) {
      struct CPMDCpmdSection cpmd_sec;
      memset(&cpmd_sec, 0, sizeof(cpmd_sec));
      cpmd_sec.optimizeWavefunction = 1;
      /* CPMD's primary SCF criterion is CONVERGENCE ORBITALS (a.u.). Its
       * optional second value is the total-energy change between SCF
       * steps; the free-energy (Lanczos) SCF stops on either one. The
       * energy tolerance goes to both, so a smeared SCF stops once the
       * energy settles even while the density still moves. */
      char conv_args[64];
      if (c.scfEnergyToleranceEv > 0.0) {
        double tol = c.scfEnergyToleranceEv / CPMDC_HARTREE_EV;
        snprintf(conv_args, sizeof(conv_args), "%.10g %.10g", tol, tol);
        capn_ptr args = capn_new_ptr_list(root.seg, 1);
        capn_text arg_text = {(int)strlen(conv_args), conv_args, NULL};
        capn_set_text(args, 0, arg_text);
        struct CPMDDirective conv;
        memset(&conv, 0, sizeof(conv));
        conv.keyword.str = "CONVERGENCE ORBITALS";
        conv.keyword.len = (int)strlen(conv.keyword.str);
        conv.args = args;
        CPMDDirective_list directives = new_CPMDDirective_list(root.seg, 1);
        set_CPMDDirective(&conv, directives, 0);
        cpmd_sec.directives = directives;
        cpmd_sec.convergenceOrbitals = 0.0;
      } else {
        cpmd_sec.convergenceOrbitals = 1.0e-6;
      }
      cpmd_sec.maxIter = c.scfMaxIterations;
      if (smear_temperature_k > 0.0) {
        cpmd_sec.freeEnergy = 1;
        cpmd_sec.temperatureElectron = smear_temperature_k;
      }
      CPMDCpmdSection_ptr cpmd_ptr = new_CPMDCpmdSection(root.seg);
      write_CPMDCpmdSection(&cpmd_sec, cpmd_ptr);
      struct CPMDInputSection sec;
      memset(&sec, 0, sizeof(sec));
      sec.which = CPMDInputSection_cpmd;
      sec.cpmd = cpmd_ptr;
      set_CPMDInputSection(&sec, sections, idx++);
    }
    if (grimme_token) {
      struct CPMDVdwSection vdw_sec;
      memset(&vdw_sec, 0, sizeof(vdw_sec));
      vdw_sec.empiricalCorrection = 1;
      vdw_sec.grimme.str = grimme_token;
      vdw_sec.grimme.len = (int)strlen(grimme_token);
      if (c.vanDerWaalsS6 > 0.0)
        vdw_sec.s6 = c.vanDerWaalsS6;
      CPMDVdwSection_ptr vdw_ptr = new_CPMDVdwSection(root.seg);
      write_CPMDVdwSection(&vdw_sec, vdw_ptr);
      struct CPMDInputSection sec;
      memset(&sec, 0, sizeof(sec));
      sec.which = CPMDInputSection_vdwParams;
      sec.vdwParams = vdw_ptr;
      set_CPMDInputSection(&sec, sections, idx++);
    }
    if (kmesh_len == 3) {
      struct CPMDSystemSection sys_sec;
      memset(&sys_sec, 0, sizeof(sys_sec));
      capn_list32 mesh = capn_new_list32(root.seg, 3);
      for (int i = 0; i < 3; ++i)
        capn_set32(mesh, i, (uint32_t)capn_get32(c.kMesh, i));
      sys_sec.kpointsMonkhorstPack = mesh;
      CPMDSystemSection_ptr sys_ptr = new_CPMDSystemSection(root.seg);
      write_CPMDSystemSection(&sys_sec, sys_ptr);
      struct CPMDInputSection sec;
      memset(&sec, 0, sizeof(sec));
      sec.which = CPMDInputSection_system;
      sec.system = sys_ptr;
      set_CPMDInputSection(&sec, sections, idx++);
    }
    view.inputSections = sections;
  }
  write_CPMDParams(&view, params);
  if (capn_setp(root, 0, params.p) != 0) {
    capn_free(&arena);
    cpmdc_store_error("common overlay: params synthesis failed");
    return -1;
  }
  size_t capacity = 4096u;
  unsigned char *buffer = NULL;
  int written = -1;
  for (int attempt = 0; attempt < 8 && written < 0; ++attempt) {
    unsigned char *next = (unsigned char *)realloc(buffer, capacity);
    if (!next) {
      free(buffer);
      capn_free(&arena);
      cpmdc_store_error("out of memory");
      return -1;
    }
    buffer = next;
    written = capn_write_mem(&arena, (uint8_t *)buffer, capacity, 0);
    if (written < 0)
      capacity *= 2u;
  }
  capn_free(&arena);
  if (written < 0) {
    free(buffer);
    cpmdc_store_error("common overlay: params serialization failed");
    return -1;
  }
  *out = buffer;
  *out_size = (size_t)written;
  return 0;
}

/* Flatten an in-arena struct pointer into a fresh malloc'd message. */
static int cpmd_write_ptr_flat(capn_ptr struct_ptr, unsigned char **out,
                               size_t *out_size) {
  *out = NULL;
  *out_size = 0;
  if (struct_ptr.type == CAPN_NULL) {
    cpmdc_store_error("CPMDParams serialization failed");
    return -1;
  }
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL || capn_setp(root, 0, struct_ptr) != 0) {
    capn_free(&arena);
    cpmdc_store_error("CPMDParams serialization failed");
    return -1;
  }
  size_t capacity = 4096u;
  unsigned char *buffer = NULL;
  int written = -1;
  for (int attempt = 0; attempt < 16 && written < 0; ++attempt) {
    unsigned char *next = (unsigned char *)realloc(buffer, capacity);
    if (!next) {
      free(buffer);
      capn_free(&arena);
      cpmdc_store_error("out of memory");
      return -1;
    }
    buffer = next;
    written = capn_write_mem(&arena, (uint8_t *)buffer, capacity, 0);
    if (written < 0)
      capacity *= 2u;
  }
  capn_free(&arena);
  if (written < 0) {
    free(buffer);
    cpmdc_store_error("CPMDParams serialization failed");
    return -1;
  }
  *out = buffer;
  *out_size = (size_t)written;
  return 0;
}

static int cpmd_raw_has_section(capn_text block, const char *name) {
  if (!block.str || block.len <= 0)
    return 0;
  size_t nlen = strlen(name);
  size_t n = (size_t)block.len;
  size_t i = 0;
  while (i < n) {
    size_t line = i;
    while (line < n && (block.str[line] == ' ' || block.str[line] == '\t'))
      line++;
    if (line + 1 + nlen <= n && block.str[line] == '&') {
      int match = 1;
      for (size_t j = 0; j < nlen; ++j) {
        char a = block.str[line + 1 + j];
        if (a >= 'a' && a <= 'z')
          a = (char)(a - ('a' - 'A'));
        if (a != name[j]) {
          match = 0;
          break;
        }
      }
      if (match && line + 1 + nlen < n) {
        char next = block.str[line + 1 + nlen];
        if ((next >= 'A' && next <= 'Z') || (next >= 'a' && next <= 'z'))
          match = 0;
      }
      if (match)
        return 1;
    }
    while (i < n && block.str[i] != '\n')
      i++;
    if (i < n)
      i++;
  }
  return 0;
}

static int cpmd_text_has_keyword_ci(const char *s, size_t len, const char *kw) {
  size_t klen = strlen(kw);
  while (len > 0 && (*s == ' ' || *s == '\t')) {
    s++;
    len--;
  }
  if (len < klen)
    return 0;
  for (size_t i = 0; i < klen; ++i) {
    char a = s[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (a != kw[i])
      return 0;
  }
  if (len > klen) {
    char next = s[klen];
    if (next != ' ' && next != '\t' && next != '\r' && next != '\n' &&
        next != '=')
      return 0;
  }
  return 1;
}

static int cpmd_text_is_kw_ci(const char *s, size_t len, const char *kw) {
  size_t klen = strlen(kw);
  while (len > 0 && (*s == ' ' || *s == '\t')) {
    s++;
    len--;
  }
  while (len > 0 && (s[len - 1] == ' ' || s[len - 1] == '\t' ||
                     s[len - 1] == '\r' || s[len - 1] == '\n')) {
    len--;
  }
  if (len != klen)
    return 0;
  for (size_t i = 0; i < klen; ++i) {
    char a = s[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (a != kw[i])
      return 0;
  }
  return 1;
}

static int cpmd_text_has_two_kw_ci(const char *s, size_t len, const char *kw1,
                                   const char *kw2) {
  size_t k1 = strlen(kw1);
  while (len > 0 && (*s == ' ' || *s == '\t')) {
    s++;
    len--;
  }
  if (len <= k1)
    return 0;
  for (size_t i = 0; i < k1; ++i) {
    char a = s[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (a != kw1[i])
      return 0;
  }
  if (s[k1] != ' ' && s[k1] != '\t' && s[k1] != '\r' && s[k1] != '\n')
    return 0;
  s += k1;
  len -= k1;
  while (len > 0 && (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n')) {
    s++;
    len--;
  }
  return cpmd_text_has_keyword_ci(s, len, kw2);
}

static int cpmd_directives_have_kw(CPMDDirective_list directives,
                                   const char *kw) {
  capn_resolve(&directives.p);
  if (directives.p.type == CAPN_NULL)
    return 0;
  if (directives.p.type != CAPN_LIST)
    return -1;
  for (int i = 0; i < directives.p.len; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (d.keyword.str && d.keyword.len > 0 &&
        cpmd_text_has_keyword_ci(d.keyword.str, (size_t)d.keyword.len, kw))
      return 1;
  }
  return 0;
}

static int cpmd_directives_have_temp_electron(CPMDDirective_list directives) {
  capn_resolve(&directives.p);
  if (directives.p.type == CAPN_NULL)
    return 0;
  if (directives.p.type != CAPN_LIST)
    return -1;
  for (int i = 0; i < directives.p.len; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (!d.keyword.str || d.keyword.len <= 0)
      continue;
    if (cpmd_text_has_two_kw_ci(d.keyword.str, (size_t)d.keyword.len,
                                "TEMPERATURE", "ELECTRON"))
      return 1;
    if (cpmd_text_is_kw_ci(d.keyword.str, (size_t)d.keyword.len,
                           "TEMPERATURE")) {
      int na = common_text_list_len(d.args);
      if (na < 0)
        return -1;
      for (int j = 0; j < na; ++j) {
        capn_text arg = capn_get_text(d.args, j, k_empty_text);
        if (!arg.str || arg.len <= 0)
          continue;
        const char *p = arg.str;
        size_t rem = (size_t)arg.len;
        while (rem > 0 &&
               (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) {
          p++;
          rem--;
        }
        if (rem == 0)
          continue;
        if (cpmd_text_has_keyword_ci(p, rem, "ELECTRON"))
          return 1;
        break;
      }
    }
  }
  return 0;
}

static int cpmd_set_matches_kw(capn_text key, const char *section,
                               const char *kw) {
  if (!key.str || key.len <= 2)
    return 0;
  const char *dot =
      (const char *)memchr(key.str, '.', (size_t)key.len);
  if (!dot || dot == key.str || dot == key.str + (size_t)key.len - 1)
    return 0;
  size_t sec_len = (size_t)(dot - key.str);
  if (sec_len != strlen(section))
    return 0;
  for (size_t i = 0; i < sec_len; ++i) {
    char a = key.str[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (a != section[i])
      return 0;
  }
  const char *k = dot + 1;
  size_t klen = (size_t)key.len - sec_len - 1u;
  return cpmd_text_has_keyword_ci(k, klen, kw);
}

static int cpmd_params_method_check_opaque_scalars(CPMDParams_ptr params) {
  struct CPMDParams view;
  read_CPMDParams(&view, params);

  int nblocks = common_text_list_len(view.inputBlocks);
  if (nblocks < 0) {
    cpmdc_store_error("CPMDParams method metadata: inputBlocks is invalid");
    return -1;
  }
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(view.inputBlocks, i, k_empty_text);
    if (cpmd_raw_has_section(block, "DFT")) {
      cpmdc_store_error(
          "CPMDParams method metadata is ambiguous: inputBlocks opens &DFT");
      return -1;
    }
    if (cpmd_raw_has_section(block, "SYSTEM")) {
      cpmdc_store_error(
          "CPMDParams method metadata is ambiguous: inputBlocks opens &SYSTEM");
      return -1;
    }
  }

  capn_resolve(&view.inputSections.p);
  if (view.inputSections.p.type == CAPN_NULL)
    return 0;
  if (view.inputSections.p.type != CAPN_LIST) {
    cpmdc_store_error("CPMDParams method metadata: inputSections is invalid");
    return -1;
  }
  for (int i = 0; i < view.inputSections.p.len; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    if (sec.which == CPMDInputSection_raw) {
      if (cpmd_raw_has_section(sec.raw, "DFT")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: a raw section opens "
            "&DFT");
        return -1;
      }
      if (cpmd_raw_has_section(sec.raw, "SYSTEM")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: a raw section opens "
            "&SYSTEM");
        return -1;
      }
      continue;
    }
    if (sec.which == CPMDInputSection_generic) {
      struct CPMDGenericSection gen;
      read_CPMDGenericSection(&gen, sec.generic);
      if (cpmd_text_equals_ci(gen.name, "DFT")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: a generic section opens "
            "&DFT");
        return -1;
      }
      if (cpmd_text_equals_ci(gen.name, "SYSTEM")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: a generic section opens "
            "&SYSTEM");
        return -1;
      }
      continue;
    }
    if (sec.which == CPMDInputSection_set) {
      struct CPMDSetDirective set_dir;
      read_CPMDSetDirective(&set_dir, sec.set);
      if (cpmd_set_matches_kw(set_dir.key, "DFT", "FUNCTIONAL")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: DFT.FUNCTIONAL set "
            "directive is present");
        return -1;
      }
      if (cpmd_set_matches_kw(set_dir.key, "SYSTEM", "CUTOFF") ||
          cpmd_set_matches_kw(set_dir.key, "SYSTEM", "CHARGE") ||
          cpmd_set_matches_kw(set_dir.key, "SYSTEM", "MULTIPLICITY")) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: SYSTEM scalar set "
            "directive is present");
        return -1;
      }
      continue;
    }
    if (sec.which == CPMDInputSection_dft) {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      int has_functional = cpmd_directives_have_kw(body.directives, "FUNCTIONAL");
      if (has_functional < 0) {
        cpmdc_store_error("CPMDParams method metadata: dft directives are invalid");
        return -1;
      }
      if (has_functional) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: typed dft section "
            "contains FUNCTIONAL");
        return -1;
      }
      continue;
    }
    if (sec.which == CPMDInputSection_system) {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      int has_cutoff = cpmd_directives_have_kw(body.directives, "CUTOFF");
      int has_charge = cpmd_directives_have_kw(body.directives, "CHARGE");
      int has_mult = cpmd_directives_have_kw(body.directives, "MULTIPLICITY");
      if (has_cutoff < 0 || has_charge < 0 || has_mult < 0) {
        cpmdc_store_error(
            "CPMDParams method metadata: system directives are invalid");
        return -1;
      }
      if (has_cutoff || has_charge || has_mult) {
        cpmdc_store_error(
            "CPMDParams method metadata is ambiguous: typed system section "
            "contains CUTOFF, CHARGE, or MULTIPLICITY");
        return -1;
      }
    }
  }
  return 0;
}

static int
cpmd_set_matches_temp_electron(const struct CPMDSetDirective *set_dir) {
  if (!set_dir->key.str || set_dir->key.len <= 2)
    return 0;
  const char *dot =
      (const char *)memchr(set_dir->key.str, '.', (size_t)set_dir->key.len);
  if (!dot || dot == set_dir->key.str ||
      dot == set_dir->key.str + (size_t)set_dir->key.len - 1)
    return 0;
  size_t sec_len = (size_t)(dot - set_dir->key.str);
  if (sec_len != 4)
    return 0;
  for (size_t i = 0; i < 4; ++i) {
    char a = set_dir->key.str[i];
    if (a >= 'a' && a <= 'z')
      a = (char)(a - ('a' - 'A'));
    if (a != "CPMD"[i])
      return 0;
  }
  const char *k = dot + 1;
  size_t klen = (size_t)set_dir->key.len - sec_len - 1u;
  if (cpmd_text_has_two_kw_ci(k, klen, "TEMPERATURE", "ELECTRON"))
    return 1;
  if (cpmd_text_is_kw_ci(k, klen, "TEMPERATURE") && set_dir->value.str &&
      set_dir->value.len > 0) {
    const char *v = set_dir->value.str;
    size_t vlen = (size_t)set_dir->value.len;
    while (vlen > 0 && (*v == ' ' || *v == '\t' || *v == '\r' || *v == '\n')) {
      v++;
      vlen--;
    }
    if (cpmd_text_has_keyword_ci(v, vlen, "ELECTRON"))
      return 1;
  }
  return 0;
}

static int cpmd_reject_invalid_list(const char *field, const char *list) {
  char msg[256];
  snprintf(msg, sizeof(msg),
           "common overlay: %s cannot merge: the cpmd arm has an invalid %s "
           "list",
           field, list);
  cpmdc_store_error(msg);
  return -1;
}

static int cpmd_arm_check_smearing(struct CPMDParams *arm_view) {
  int nblocks = common_text_list_len(arm_view->inputBlocks);
  if (nblocks < 0)
    return cpmd_reject_invalid_list("smearing", "inputBlocks");
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(arm_view->inputBlocks, i, k_empty_text);
    if (cpmd_raw_has_section(block, "CPMD")) {
      cpmdc_store_error(
          "common overlay: smearing cannot merge into a raw &CPMD block on "
          "the cpmd arm");
      return -1;
    }
  }
  capn_resolve(&arm_view->inputSections.p);
  if (arm_view->inputSections.p.type == CAPN_NULL)
    return 0;
  if (arm_view->inputSections.p.type != CAPN_LIST)
    return cpmd_reject_invalid_list("smearing", "inputSections");
  for (int i = 0; i < arm_view->inputSections.p.len; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, arm_view->inputSections, i);
    if (sec.which == CPMDInputSection_raw &&
        cpmd_raw_has_section(sec.raw, "CPMD")) {
      cpmdc_store_error(
          "common overlay: smearing cannot merge into a raw &CPMD block on "
          "the cpmd arm");
      return -1;
    }
    if (sec.which == CPMDInputSection_generic) {
      struct CPMDGenericSection gen;
      read_CPMDGenericSection(&gen, sec.generic);
      if (cpmd_text_equals_ci(gen.name, "CPMD")) {
        cpmdc_store_error(
            "common overlay: smearing cannot merge into a generic &CPMD "
            "section on the cpmd arm");
        return -1;
      }
    }
    if (sec.which == CPMDInputSection_set) {
      struct CPMDSetDirective set_dir;
      read_CPMDSetDirective(&set_dir, sec.set);
      if (cpmd_set_matches_kw(set_dir.key, "CPMD", "FREE-ENERGY") ||
          cpmd_set_matches_temp_electron(&set_dir)) {
        cpmdc_store_error(
            "common overlay: smearing conflicts with &CPMD free-energy or "
            "electron-temperature settings already on the cpmd arm");
        return -1;
      }
    }
    if (sec.which == CPMDInputSection_cpmd) {
      struct CPMDCpmdSection cpmd_sec;
      read_CPMDCpmdSection(&cpmd_sec, sec.cpmd);
      int has_fe = cpmd_directives_have_kw(cpmd_sec.directives, "FREE-ENERGY");
      int has_te = cpmd_directives_have_temp_electron(cpmd_sec.directives);
      if (has_fe < 0 || has_te < 0)
        return cpmd_reject_invalid_list("smearing", "&CPMD directives");
      if (cpmd_sec.freeEnergy || cpmd_sec.temperatureElectron != 0.0 ||
          has_fe > 0 || has_te > 0) {
        cpmdc_store_error(
            "common overlay: smearing conflicts with &CPMD free-energy or "
            "electron-temperature settings already on the cpmd arm");
        return -1;
      }
    }
  }
  return 0;
}

static int cpmd_arm_check_vdw_s6(struct CPMDParams *arm_view) {
  int nblocks = common_text_list_len(arm_view->inputBlocks);
  if (nblocks < 0)
    return cpmd_reject_invalid_list("vanDerWaalsS6", "inputBlocks");
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(arm_view->inputBlocks, i, k_empty_text);
    if (cpmd_raw_has_section(block, "VDW")) {
      cpmdc_store_error(
          "common overlay: vanDerWaalsS6 cannot merge into a raw &VDW block "
          "on the cpmd arm; use vdwParams without s6 or set s6 on the arm");
      return -1;
    }
  }
  int has_vdw_params = 0;
  capn_resolve(&arm_view->inputSections.p);
  if (arm_view->inputSections.p.type != CAPN_NULL) {
    if (arm_view->inputSections.p.type != CAPN_LIST)
      return cpmd_reject_invalid_list("vanDerWaalsS6", "inputSections");
    for (int i = 0; i < arm_view->inputSections.p.len; ++i) {
      struct CPMDInputSection sec;
      get_CPMDInputSection(&sec, arm_view->inputSections, i);
      if (sec.which == CPMDInputSection_raw &&
          cpmd_raw_has_section(sec.raw, "VDW")) {
        cpmdc_store_error(
            "common overlay: vanDerWaalsS6 cannot merge into a raw &VDW block "
            "on the cpmd arm; use vdwParams without s6 or set s6 on the arm");
        return -1;
      }
      if (sec.which == CPMDInputSection_vdw) {
        cpmdc_store_error(
            "common overlay: vanDerWaalsS6 cannot merge into a directive "
            "&VDW section on the cpmd arm; use vdwParams without s6 or set s6 "
            "on the arm");
        return -1;
      }
      if (sec.which == CPMDInputSection_generic) {
        struct CPMDGenericSection gen;
        read_CPMDGenericSection(&gen, sec.generic);
        if (cpmd_text_equals_ci(gen.name, "VDW")) {
          cpmdc_store_error(
              "common overlay: vanDerWaalsS6 cannot merge into a generic &VDW "
              "section on the cpmd arm; use vdwParams without s6 or set s6 on "
              "the arm");
          return -1;
        }
      }
      if (sec.which == CPMDInputSection_set) {
        struct CPMDSetDirective set_dir;
        read_CPMDSetDirective(&set_dir, sec.set);
        if (cpmd_set_matches_kw(set_dir.key, "VDW", "S6")) {
          cpmdc_store_error(
              "common overlay: vanDerWaalsS6 conflicts with VDW.S6 already "
              "set on the cpmd arm");
          return -1;
        }
      }
      if (sec.which == CPMDInputSection_vdwParams) {
        struct CPMDVdwSection vdw_sec;
        read_CPMDVdwSection(&vdw_sec, sec.vdwParams);
        int has_s6 = cpmd_directives_have_kw(vdw_sec.directives, "S6");
        if (has_s6 < 0)
          return cpmd_reject_invalid_list("vanDerWaalsS6", "&VDW directives");
        if (vdw_sec.s6 != 0.0 || has_s6) {
          cpmdc_store_error(
              "common overlay: vanDerWaalsS6 conflicts with s6 already set on "
              "the cpmd arm");
          return -1;
        }
        has_vdw_params = 1;
      }
    }
  }
  if (!has_vdw_params) {
    cpmdc_store_error(
        "common overlay: vanDerWaalsS6 requires a vdwParams section on the "
        "cpmd arm (or set s6 on the arm)");
    return -1;
  }
  return 0;
}

/* inputBlocks, a raw section, and a generic section are emitted as that
 * section and are not read back into the effective functional, cutoff,
 * charge, or multiplicity. A set directive or typed section directive for the
 * same keyword is appended after the generated line and is not read back
 * either. Every set overlay scalar is rejected in those cases, whether or not
 * a typed field also sets it, so the deck cannot carry two values for one
 * setting and the recorded settings cannot name a value the deck does not
 * render. */
static int cpmd_reject_opaque_scalar(struct CPMDParams *arm, const char *section,
                                     const char *set_kw, const char *field) {
  char msg[512];
  int nblocks = common_text_list_len(arm->inputBlocks);
  if (nblocks < 0)
    return cpmd_reject_invalid_list(field, "inputBlocks");
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(arm->inputBlocks, i, k_empty_text);
    if (cpmd_raw_has_section(block, section)) {
      snprintf(msg, sizeof(msg),
               "common overlay: %s conflicts with an inputBlocks &%s block on "
               "the cpmd arm",
               field, section);
      cpmdc_store_error(msg);
      return -1;
    }
  }
  capn_resolve(&arm->inputSections.p);
  if (arm->inputSections.p.type == CAPN_NULL)
    return 0;
  if (arm->inputSections.p.type != CAPN_LIST)
    return cpmd_reject_invalid_list(field, "inputSections");
  for (int i = 0; i < arm->inputSections.p.len; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, arm->inputSections, i);
    if (sec.which == CPMDInputSection_raw &&
        cpmd_raw_has_section(sec.raw, section)) {
      snprintf(msg, sizeof(msg),
               "common overlay: %s conflicts with a raw &%s block on the cpmd "
               "arm",
               field, section);
      cpmdc_store_error(msg);
      return -1;
    }
    if (sec.which == CPMDInputSection_generic) {
      struct CPMDGenericSection gen;
      read_CPMDGenericSection(&gen, sec.generic);
      if (cpmd_text_equals_ci(gen.name, section)) {
        snprintf(msg, sizeof(msg),
                 "common overlay: %s conflicts with a generic &%s section on "
                 "the cpmd arm",
                 field, section);
        cpmdc_store_error(msg);
        return -1;
      }
    }
    if (sec.which == CPMDInputSection_set) {
      struct CPMDSetDirective set_dir;
      read_CPMDSetDirective(&set_dir, sec.set);
      if (cpmd_set_matches_kw(set_dir.key, section, set_kw)) {
        snprintf(msg, sizeof(msg),
                 "common overlay: %s conflicts with %s.%s already set on the "
                 "cpmd arm",
                 field, section, set_kw);
        cpmdc_store_error(msg);
        return -1;
      }
    }
    CPMDDirective_list typed_directives;
    int has_typed_directives = 0;
    memset(&typed_directives, 0, sizeof(typed_directives));
    if (strcmp(section, "DFT") == 0 && sec.which == CPMDInputSection_dft) {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      typed_directives = body.directives;
      has_typed_directives = 1;
    } else if (strcmp(section, "SYSTEM") == 0 &&
               sec.which == CPMDInputSection_system) {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      typed_directives = body.directives;
      has_typed_directives = 1;
    }
    if (has_typed_directives) {
      int has_keyword = cpmd_directives_have_kw(typed_directives, set_kw);
      if (has_keyword < 0) {
        char list[32];
        snprintf(list, sizeof(list), "&%s directives", section);
        return cpmd_reject_invalid_list(field, list);
      }
      if (has_keyword) {
        snprintf(msg, sizeof(msg),
                 "common overlay: %s conflicts with a typed &%s directive "
                 "already on the cpmd arm",
                 field, section);
        cpmdc_store_error(msg);
        return -1;
      }
    }
  }
  return 0;
}

/* Merge a set common overlay over a native cpmd arm: the arm's bytes are the
 * params; overlay scalars fill fields no typed field sets, where a field
 * omitted or at its schema default (functional BLYP, cutoff 70 Ry, charge 0,
 * multiplicity 1) is unset, apply Fermi smearing when &CPMD has none, copy
 * vanDerWaalsS6 into vdwParams when s6 is unset, and require basisSet to be
 * unset or "planewave". An unmapped xcFunctionals token is rejected for any
 * arm functional. Every set overlay scalar is rejected when inputBlocks, a
 * raw section, a generic section, or a set or typed-section directive
 * already supplies that setting, even when a typed field wins over the
 * overlay. Conflicting or unmergeable
 * overlay fields are rejected through cpmdc_last_error(). */
static int cpmd_merge_common_into_arm(CPMDParams_ptr arm,
                                      CommonMethodSpec_ptr common_root,
                                      unsigned char **out, size_t *out_size,
                                      CPMDCScalarOverrides *ov,
                                      char *ov_functional, size_t ov_len,
                                      int *has_ov) {
  struct CommonMethodSpec c;
  memset(&c, 0, sizeof(c));
  read_CommonMethodSpec(&c, common_root);

  if (cpmd_common_validate_scalar_ranges(&c) != 0)
    return -1;
  if (common_list32_len(&c.kMesh) != 0 || c.scfMaxIterations != 0 ||
      c.scfEnergyToleranceEv != 0.0) {
    cpmdc_store_error(
        "common overlay: kMesh/SCF controls next to a native cpmd arm are "
        "not yet mergeable; set them on the arm's system/cpmd sections");
    return -1;
  }
  if (c.vanDerWaalsMethod.len > 0 || c.relativityMethod.len > 0) {
    cpmdc_store_error(
        "common overlay: vanDerWaals/relativity have no CPMD lowering yet");
    return -1;
  }
  if (cpmd_common_validate_basis_set(c.basisSet) != 0)
    return -1;

  double smear_temperature_k = 0.0;
  if (cpmd_common_validate_smearing(&c, &smear_temperature_k) != 0)
    return -1;

  int xc_len = common_text_list_len(c.xcFunctionals);
  if (xc_len < 0) {
    cpmdc_store_error("common overlay: invalid xcFunctionals list");
    return -1;
  }

  /* The arm is a struct here, so only an inputSections value that is not a
   * struct list makes this fail. The smearing and vanDerWaalsS6 checks
   * report that list themselves. */
  CPMDCTypedScalars typed;
  if (cpmdc_params_typed_scalars(arm, &typed) != 0) {
    const char *field = xc_len > 0                   ? "xcFunctionals"
                        : c.planewaveCutoffEv > 0.0 ? "planewaveCutoffEv"
                        : c.charge != 0             ? "charge"
                        : c.spinMultiplicity != 0   ? "spinMultiplicity"
                                                    : NULL;
    if (field)
      return cpmd_reject_invalid_list(field, "inputSections");
    memset(&typed, 0, sizeof(typed));
  }

  struct CPMDParams arm_view;
  memset(&arm_view, 0, sizeof(arm_view));
  read_CPMDParams(&arm_view, arm);

  memset(ov, 0, sizeof(*ov));
  *has_ov = 0;
  if (xc_len > 0) {
    if (cpmd_common_xc_token(c.xcFunctionals, ov_functional, ov_len) != 0) {
      cpmdc_store_error("common overlay: unmapped xcFunctionals for CPMD");
      return -1;
    }
    if (cpmd_reject_opaque_scalar(&arm_view, "DFT", "FUNCTIONAL",
                                  "xcFunctionals") != 0)
      return -1;
    if (!typed.functional) {
      ov->functional = ov_functional;
      *has_ov = 1;
    }
  }
  if (c.planewaveCutoffEv > 0.0) {
    if (cpmd_reject_opaque_scalar(&arm_view, "SYSTEM", "CUTOFF",
                                  "planewaveCutoffEv") != 0)
      return -1;
    if (!typed.cutoff) {
      ov->cutoff_ry = c.planewaveCutoffEv / CPMDC_RYDBERG_EV;
      *has_ov = 1;
    }
  }
  if (c.hasCharge || c.charge != 0) {
    if (cpmd_reject_opaque_scalar(&arm_view, "SYSTEM", "CHARGE", "charge") !=
        0)
      return -1;
    if (!typed.charge && c.charge != 0) {
      ov->has_charge = 1;
      ov->charge = c.charge;
      *has_ov = 1;
    }
  }
  if (c.spinMultiplicity != 0) {
    if (cpmd_reject_opaque_scalar(&arm_view, "SYSTEM", "MULTIPLICITY",
                                  "spinMultiplicity") != 0)
      return -1;
    if (!typed.multiplicity) {
      ov->has_multiplicity = 1;
      ov->multiplicity = c.spinMultiplicity;
      *has_ov = 1;
    }
  }
  if (smear_temperature_k > 0.0) {
    if (cpmd_arm_check_smearing(&arm_view) != 0)
      return -1;
    ov->free_energy = 1;
    ov->temperature_electron = smear_temperature_k;
    *has_ov = 1;
  }
  if (c.vanDerWaalsS6 > 0.0) {
    if (cpmd_arm_check_vdw_s6(&arm_view) != 0)
      return -1;
    ov->vdw_s6 = c.vanDerWaalsS6;
    *has_ov = 1;
  }
  return cpmd_write_ptr_flat(arm.p, out, out_size);
}

/* Resolve a PotentialConfig into effective CPMDParams bytes: the cpmd arm
 * wins where it differs from schema defaults; a set common overlay fills the
 * rest or synthesizes params when no arm is present. *out NULL with rc 0
 * means nothing to apply. */
static int cpmd_config_to_params_bytes(const void *config_capnp,
                                       size_t config_capnp_size_bytes,
                                       unsigned char **out, size_t *out_size,
                                       CPMDCScalarOverrides *ov,
                                       char *ov_functional, size_t ov_len,
                                       int *has_ov) {
  *out = NULL;
  *out_size = 0;
  memset(ov, 0, sizeof(*ov));
  *has_ov = 0;
  if (!config_capnp || config_capnp_size_bytes == 0) {
    cpmdc_store_error("PotentialConfig buffer is empty");
    return -1;
  }
  struct capn arena;
  if (capn_init_mem(&arena, (const uint8_t *)config_capnp,
                    config_capnp_size_bytes, 0) != 0) {
    cpmdc_store_error("PotentialConfig parse failed");
    return -1;
  }
  PotentialConfig_ptr config_root;
  config_root.p = capn_getp(capn_root(&arena), 0, 1);
  if (config_root.p.type != CAPN_STRUCT) {
    capn_free(&arena);
    cpmdc_store_error("PotentialConfig parse failed");
    return -1;
  }
  struct PotentialConfig config;
  memset(&config, 0, sizeof(config));
  read_PotentialConfig(&config, config_root);

  int rc = 0;
  if (config.which == PotentialConfig_cpmd) {
    capn_resolve(&config.cpmd.p);
    if (config.cpmd.p.type != CAPN_STRUCT) {
      cpmdc_store_error("PotentialConfig cpmd arm is empty");
      rc = -1;
    } else {
      capn_resolve(&config.common.p);
      if (config.common.p.type == CAPN_STRUCT) {
        rc = cpmd_merge_common_into_arm(config.cpmd, config.common, out,
                                        out_size, ov, ov_functional, ov_len,
                                        has_ov);
      } else {
        rc = cpmd_write_ptr_flat(config.cpmd.p, out, out_size);
        if (rc != 0)
          cpmdc_store_error("PotentialConfig cpmd arm serialization failed");
      }
    }
  } else if (config.which == PotentialConfig_none) {
    capn_resolve(&config.common.p);
    if (config.common.p.type == CAPN_STRUCT)
      rc = cpmd_common_to_params(config.common, out, out_size);
  } else {
    cpmdc_store_error(
        "PotentialConfig arm does not match the CPMD backend");
    rc = -1;
  }
  capn_free(&arena);
  return rc;
}

int cpmdc_configure(const void *config_capnp,
                    size_t config_capnp_size_bytes) {
  cpmdc_store_error("");
  unsigned char *params = NULL;
  size_t params_size = 0;
  CPMDCScalarOverrides ov;
  char ov_functional[64];
  int has_ov = 0;
  if (cpmd_config_to_params_bytes(config_capnp, config_capnp_size_bytes,
                                  &params, &params_size, &ov, ov_functional,
                                  sizeof(ov_functional), &has_ov) != 0)
    return -1;
  if (!params)
    return 0;
  int rc = cpmdc_set_params_ov(params, params_size, has_ov ? &ov : NULL);
  if (rc != 0 && g_last_error[0] == '\0')
    cpmdc_store_error("applying resolved CPMD params failed");
  free(params);
  return rc;
}

CPMDCSession *cpmdc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes) {
  cpmdc_store_error("");
  unsigned char *params = NULL;
  size_t params_size = 0;
  CPMDCScalarOverrides ov;
  char ov_functional[64];
  int has_ov = 0;
  if (cpmd_config_to_params_bytes(config_capnp, config_capnp_size_bytes,
                                  &params, &params_size, &ov, ov_functional,
                                  sizeof(ov_functional), &has_ov) != 0)
    return NULL;
  if (!params) {
    cpmdc_store_error("PotentialConfig carries no CPMD-applicable settings");
    return NULL;
  }
  CPMDCSession *session =
      cpmdc_session_create_ov(params, params_size, has_ov ? &ov : NULL);
  free(params);
  return session;
}

int cpmdc_session_configure(CPMDCSession *session, const void *config_capnp,
                            size_t config_capnp_size_bytes) {
  cpmdc_store_error("");
  if (!session) {
    cpmdc_store_error("invalid session");
    return -1;
  }
  unsigned char *params = NULL;
  size_t params_size = 0;
  CPMDCScalarOverrides ov;
  char ov_functional[64];
  int has_ov = 0;
  if (cpmd_config_to_params_bytes(config_capnp, config_capnp_size_bytes,
                                  &params, &params_size, &ov, ov_functional,
                                  sizeof(ov_functional), &has_ov) != 0)
    return -1;
  if (!params)
    return 0;
  int rc = cpmdc_session_set_params_ov(session, params, params_size,
                                       has_ov ? &ov : NULL);
  if (rc != 0 && g_last_error[0] == '\0')
    cpmdc_store_error("installing resolved CPMD params failed");
  free(params);
  return rc;
}

static int cpmdc_set_params_ov(const void *params_capnp,
                               size_t params_capnp_size_bytes,
                               const CPMDCScalarOverrides *overrides) {
  char functional[64];
  double cutoff_ry = 70.0;
  int charge = 0;
  int multiplicity = 1;
  char *input_deck = NULL;
  char cpmd_root[1024];
  if (apply_params_buffer(params_capnp, params_capnp_size_bytes, overrides,
                          functional, sizeof(functional), &cutoff_ry, &charge,
                          &multiplicity, &input_deck, cpmd_root,
                          sizeof(cpmd_root)) != 0) {
    free(input_deck);
    return -1;
  }
  if (overrides) {
    g_overrides = *overrides;
    if (overrides->functional) {
      snprintf(g_overrides_functional, sizeof(g_overrides_functional), "%s",
               overrides->functional);
      g_overrides.functional = g_overrides_functional;
    }
    g_has_overrides = 1;
  } else {
    g_has_overrides = 0;
  }
  if (!ensure_embed_init() || !cpmdc_embed_available()) {
    free(input_deck);
    cpmdc_store_error("CPMD embed not available");
    return -1;
  }
  free(g_params_bytes);
  g_params_bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!g_params_bytes) {
    free(input_deck);
    cpmdc_store_error("out of memory");
    return -1;
  }
  memcpy(g_params_bytes, params_capnp, params_capnp_size_bytes);
  g_params_size = params_capnp_size_bytes;
  /* Fortran decode of functional/cutoff/charge/mult/root + store C deck. */
  if (embed_apply_from_wire(params_capnp, params_capnp_size_bytes, input_deck,
                            overrides, &g_image) != 0) {
    free(input_deck);
    cpmdc_store_error("applying CPMD params failed");
    return -1;
  }
  free(input_deck);
  (void)functional;
  (void)cutoff_ry;
  (void)charge;
  (void)multiplicity;
  (void)cpmd_root;
  g_active_session = NULL;
  g_orbital_owner = NULL;
  cpmdc_orbital_clear_key(0);
  cpmdc_store_error("");
  return 0;
}

int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes) {
  return cpmdc_set_params_ov(params_capnp, params_capnp_size_bytes, NULL);
}

static CPMDCRankSpread g_rank_spread;

#if defined(CPMDC_HAVE_MPI)
/* The calculator communicator, or MPI_COMM_WORLD when none was bound or
 * adopted. -1 when consensus is off or MPI is not running. */
static int rank_consensus_comm(MPI_Comm *comm) {
  const char *v = getenv("CPMDC_RANK_CONSENSUS");
  int inited = 0;
  int finalized = 0;
  if (!v || strcmp(v, "1") != 0)
    return -1;
  if (MPI_Initialized(&inited) != MPI_SUCCESS || !inited ||
      MPI_Finalized(&finalized) != MPI_SUCCESS || finalized)
    return -1;
  if (cpmdc_adopted_comm(comm, sizeof(*comm)) != 0)
    *comm = MPI_COMM_WORLD;
  return 0;
}

static double rank_difference(double local, double root) {
  double d;
  if (local == root)
    return 0.0;
  d = fabs(local - root);
  return isnan(d) ? HUGE_VAL : d;
}
#endif

/* With consensus on, the flag holds only where it holds on every rank, so
 * every rank takes the same branch into the next collective. */
static int rank_consensus_all(int flag) {
  int failed = !flag;
  if (cpmdc_orbital_allreduce(&failed, 1) != 0)
    return 0;
  return !failed;
}

/* Public evaluation wrappers call this after local preflight and before any
 * rank can enter the evaluation collectives. */
static const char *collective_call_gate(const char *local_err) {
#if defined(CPMDC_HAVE_MPI)
  static char agreed[sizeof(((CPMDCResult *)0)->message)];
  int inited = 0;
  int finalized = 0;
  int rank = 0;
  int owner;
  MPI_Comm comm;
  if (MPI_Initialized(&inited) != MPI_SUCCESS)
    return "MPI collective failed";
  if (!inited)
    return local_err;
  if (MPI_Finalized(&finalized) != MPI_SUCCESS)
    return "MPI collective failed";
  if (finalized)
    return local_err;
  if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0)
    comm = MPI_COMM_WORLD;
  if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS)
    return "MPI collective failed";
  owner = local_err ? rank : INT_MAX;
  if (MPI_Allreduce(MPI_IN_PLACE, &owner, 1, MPI_INT, MPI_MIN, comm) !=
      MPI_SUCCESS)
    return "MPI collective failed";
  if (owner == INT_MAX)
    return NULL;
  if (rank == owner)
    snprintf(agreed, sizeof(agreed), "%s", local_err);
  if (MPI_Bcast(agreed, (int)sizeof(agreed), MPI_CHAR, owner, comm) !=
      MPI_SUCCESS)
    return "MPI collective failed";
  return agreed;
#else
  return local_err;
#endif
}

/* CPMD treats rank 0 of the calculator as holding the result. Every rank
 * fails when one failed and resets when one stopped. Otherwise every rank
 * takes rank 0's energy, gradient, and snapshots, and the spread records
 * how far the other ranks were from them. */
static int rank_consensus(int *ok, int *stop, size_t n3, double *energy,
                          double *grad, CPMDCEmbedImage *image) {
  memset(&g_rank_spread, 0, sizeof(g_rank_spread));
#if defined(CPMDC_HAVE_MPI)
  MPI_Comm comm;
  int size = 0;
  int flags[2];
  double dev[2] = {0.0, 0.0};
  double *root;
  if (rank_consensus_comm(&comm) != 0)
    return 0;
  if (MPI_Comm_size(comm, &size) != MPI_SUCCESS)
    return -1;
  root = (double *)malloc((n3 + 1u) * sizeof(*root));
  flags[0] = !*ok || !root;
  flags[1] = *stop != 0;
  if (MPI_Allreduce(MPI_IN_PLACE, flags, 2, MPI_INT, MPI_MAX, comm) !=
      MPI_SUCCESS) {
    free(root);
    return -1;
  }
  if (flags[1] && *stop == 0)
    *stop = 1;
  if (flags[0] && *ok) {
    *ok = 0;
    cpmdc_note_embed_failure(
        "CPMD energy/gradient failed on another rank of the calculator");
  }
  if (flags[0] || flags[1]) {
    free(root);
    return 0;
  }
  root[0] = *energy;
  memcpy(root + 1, grad, n3 * sizeof(*grad));
  if (MPI_Bcast(root, (int)(n3 + 1u), MPI_DOUBLE, 0, comm) != MPI_SUCCESS ||
      MPI_Bcast(image, (int)offsetof(CPMDCEmbedImage, warm_cell), MPI_BYTE, 0,
                comm) != MPI_SUCCESS) {
    free(root);
    return -1;
  }
  dev[0] = rank_difference(*energy, root[0]);
  for (size_t i = 0; i < n3; ++i) {
    double d = rank_difference(grad[i], root[i + 1u]);
    if (d > dev[1])
      dev[1] = d;
  }
  *energy = root[0];
  memcpy(grad, root + 1, n3 * sizeof(*grad));
  free(root);
  if (MPI_Allreduce(MPI_IN_PLACE, dev, 2, MPI_DOUBLE, MPI_MAX, comm) !=
      MPI_SUCCESS)
    return -1;
  g_rank_spread.valid = 1;
  g_rank_spread.ranks = size;
  g_rank_spread.energy_h = dev[0];
  g_rank_spread.gradient_h_bohr = dev[1];
#else
  (void)ok;
  (void)stop;
  (void)n3;
  (void)energy;
  (void)grad;
  (void)image;
#endif
  return 0;
}

static int line_defines_external_potential(const char *line, size_t len) {
  static const char key[] = "EXTERNAL POTENTIAL";
  const size_t key_len = sizeof(key) - 1u;
  size_t i;
  while (len > 0 && (*line == ' ' || *line == '\t')) {
    line++;
    len--;
  }
  if (len < key_len)
    return 0;
  for (i = 0; i < key_len; ++i) {
    char c = line[i];
    if (c >= 'a' && c <= 'z')
      c = (char)(c - 'a' + 'A');
    if (c != key[i])
      return 0;
  }
  if (len == key_len)
    return 1;
  return line[key_len] == ' ' || line[key_len] == '\t' || line[key_len] == '\r';
}

/* The rendered deck is what CPMD reads. EXTERNAL POTENTIAL and
 * EXTERNAL POTENTIAL ADD both start a line with that keyword and turn on
 * the deck-owned external field. */
int cpmdc_deck_defines_external_potential(const char *deck, size_t deck_len) {
  const char *p;
  size_t left;
  if (!deck || deck_len == 0)
    return 0;
  p = deck;
  left = deck_len;
  while (left > 0) {
    const char *nl = (const char *)memchr(p, '\n', left);
    size_t len = nl ? (size_t)(nl - p) : left;
    if (line_defines_external_potential(p, len))
      return 1;
    if (!nl)
      break;
    left -= len + 1u;
    p = nl + 1;
  }
  return 0;
}

static size_t trim_line_bounds(const char *line, size_t len,
                               const char **start_out) {
  while (len > 0 && (*line == ' ' || *line == '\t' || *line == '\r')) {
    line++;
    len--;
  }
  while (len > 0 &&
         (line[len - 1] == ' ' || line[len - 1] == '\t' ||
          line[len - 1] == '\r')) {
    len--;
  }
  *start_out = line;
  return len;
}

static int line_starts_keyword_ci(const char *line, size_t len,
                                  const char *key) {
  size_t key_len = strlen(key);
  if (len < key_len)
    return 0;
  for (size_t i = 0; i < key_len; ++i) {
    char c = line[i];
    if (c >= 'a' && c <= 'z')
      c = (char)(c - 'a' + 'A');
    if (c != key[i])
      return 0;
  }
  if (len == key_len)
    return 1;
  return line[key_len] == ' ' || line[key_len] == '\t' ||
         line[key_len] == '!';
}

static int token_is_isolated_symmetry_ci(const char *tok, size_t len) {
  static const char isolated[] = "ISOLATED";
  const size_t isolated_len = sizeof(isolated) - 1u;
  size_t tok_len = 0;
  while (tok_len < len && tok[tok_len] != ' ' && tok[tok_len] != '\t' &&
         tok[tok_len] != '!')
    tok_len++;
  if (tok_len == 0)
    return 0;
  if (tok_len == 1 && tok[0] == '0')
    return 1;
  if (tok_len >= isolated_len) {
    int match = 1;
    for (size_t i = 0; i < isolated_len; ++i) {
      char c = tok[i];
      if (c >= 'a' && c <= 'z')
        c = (char)(c - 'a' + 'A');
      if (c != isolated[i]) {
        match = 0;
        break;
      }
    }
    if (match)
      return 1;
  }
  return 0;
}

static int deck_or_cell_is_periodic(const char *deck, const double *cell_ang,
                                    int has_cell) {
  int in_sys = 0;
  int expect_sym_val = 0;
  int sym_isolated = 0;
  int sym_periodic = 0;
  int has_isolated_kw = 0;
  int has_periodic_kw = 0;
  int has_cell_vectors = 0;

  if (deck && deck[0] != '\0') {
    const char *p = deck;
    size_t left = strlen(deck);
    while (left > 0) {
      const char *nl = (const char *)memchr(p, '\n', left);
      size_t raw_len = nl ? (size_t)(nl - p) : left;
      const char *line = NULL;
      size_t len = trim_line_bounds(p, raw_len, &line);

      if (len > 0 && line[0] != '!') {
        if (line_starts_keyword_ci(line, len, "&SYSTEM")) {
          in_sys = 1;
          expect_sym_val = 0;
        } else if (in_sys && line_starts_keyword_ci(line, len, "&END")) {
          in_sys = 0;
          expect_sym_val = 0;
        } else if (in_sys && expect_sym_val) {
          if (token_is_isolated_symmetry_ci(line, len))
            sym_isolated = 1;
          else
            sym_periodic = 1;
          expect_sym_val = 0;
        } else if (in_sys) {
          if (line_starts_keyword_ci(line, len, "SYMMETRY")) {
            const char *arg = line + 8u;
            size_t arg_len = len - 8u;
            while (arg_len > 0 && (*arg == ' ' || *arg == '\t')) {
              arg++;
              arg_len--;
            }
            if (arg_len > 0 && *arg != '!') {
              if (token_is_isolated_symmetry_ci(arg, arg_len))
                sym_isolated = 1;
              else
                sym_periodic = 1;
            } else {
              expect_sym_val = 1;
            }
          } else if (line_starts_keyword_ci(line, len, "POISSON SOLVER") ||
                     line_starts_keyword_ci(line, len, "CLUSTER") ||
                     line_starts_keyword_ci(line, len, "ISOLATED MOLECULE") ||
                     line_starts_keyword_ci(line, len, "MOLECULE ISOLATED")) {
            has_isolated_kw = 1;
          } else if (line_starts_keyword_ci(line, len, "SURFACE") ||
                     line_starts_keyword_ci(line, len, "POLYMER")) {
            has_periodic_kw = 1;
          } else if (line_starts_keyword_ci(line, len, "CELL VECTORS")) {
            has_cell_vectors = 1;
          }
        }
      }
      if (!nl)
        break;
      left -= raw_len + 1u;
      p = nl + 1;
    }
  }

  if (sym_periodic || has_periodic_kw)
    return 1;
  if (sym_isolated || has_isolated_kw)
    return 0;
  if (has_cell_vectors)
    return 1;
  return has_cell && cell_ang &&
         cpmdc_cell_rank3_state(cell_ang) != CPMDC_CELL_ZERO;
}

static void session_drop_cache(CPMDCSession *session);

/* Drop energy, gradient, stress, charges, and the other per-evaluation
 * snapshots. A failed attempt must not leave the previous step readable.
 * The warm cell and the method knobs stay. */
static void image_invalidate_results(CPMDCEmbedImage *image) {
  if (!image)
    return;
  memset(image, 0, offsetof(CPMDCEmbedImage, warm_cell));
}

/* The cached step is the same result. Leaving it set makes the next
 * identical call succeed from the cache and restore only the charges. */
static void session_drop_results(CPMDCSession *session) {
  if (!session)
    return;
  session->external_forces_valid = 0;
  image_invalidate_results(&session->image);
  session_drop_cache(session);
}

static CPMDCResult fail_image(CPMDCEmbedImage *image, const char *msg) {
  image_invalidate_results(image);
  return fail_msg(msg);
}

/* Every rank of the calculator calls this after each cold-setup step.
 * stopgm can return on one rank only, and a close or directory failure
 * is noted on one rank only. The flag is agreed before that rank frees
 * the group or the others enter the next setup collective. A pending
 * embed failure is copied to every rank so each reports the same text. */
int cpmdc_agree_cold_stop(void) {
  int flags[2];
  flags[0] = cpmdc_stop_code() != 0 ? 1 : 0;
  flags[1] = g_embed_detail[0] != '\0' ? 1 : 0;
  if (cpmdc_orbital_allreduce(flags, 2) != 0) {
    cpmdc_note_stop(1);
    return 1;
  }
  if (flags[0] && cpmdc_stop_code() == 0)
    cpmdc_note_stop(1);
#if defined(CPMDC_HAVE_MPI)
  if (flags[1]) {
    int inited = 0;
    int finalized = 0;
    int rank = 0;
    int owner;
    int code;
    char detail[sizeof(g_embed_detail)];
    MPI_Comm comm;
    if (MPI_Initialized(&inited) != MPI_SUCCESS || !inited)
      return 1;
    if (MPI_Finalized(&finalized) != MPI_SUCCESS || finalized)
      return 1;
    if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0)
      comm = MPI_COMM_WORLD;
    if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
      cpmdc_note_stop(1);
      return 1;
    }
    owner = g_embed_detail[0] != '\0' ? rank : INT_MAX;
    if (MPI_Allreduce(MPI_IN_PLACE, &owner, 1, MPI_INT, MPI_MIN, comm) !=
        MPI_SUCCESS) {
      cpmdc_note_stop(1);
      return 1;
    }
    memcpy(detail, g_embed_detail, sizeof(detail));
    code = (int)g_embed_code;
    if (owner != INT_MAX) {
      if (MPI_Bcast(&code, 1, MPI_INT, owner, comm) != MPI_SUCCESS ||
          MPI_Bcast(detail, (int)sizeof(detail), MPI_CHAR, owner, comm) !=
              MPI_SUCCESS) {
        cpmdc_note_stop(1);
        return 1;
      }
      detail[sizeof(detail) - 1] = '\0';
      memcpy(g_embed_detail, detail, sizeof(g_embed_detail));
      g_embed_code = (CPMDCErrorCode)code;
    }
  }
#endif
  return flags[0] || flags[1];
}

#if defined(CPMDC_ENABLE_TEST_FAULTS)
/* The reference build has no setup list. One collective stands in for
 * the step that just returned, rank 1 notes a stop, and then every rank
 * either frees the stand-in group or enters the next collective. */
static int cpmdc_cold_setup_stop_drill(void) {
#if defined(CPMDC_HAVE_MPI)
  MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm setup = MPI_COMM_NULL;
  int rank = 0;
  int inited = 0;
  int agreed;
  if (MPI_Initialized(&inited) != MPI_SUCCESS || !inited)
    return -1;
  if (cpmdc_adopted_comm(&world, sizeof(world)) != 0)
    world = MPI_COMM_WORLD;
  if (MPI_Comm_rank(world, &rank) != MPI_SUCCESS)
    return -1;
  if (MPI_Comm_dup(world, &setup) != MPI_SUCCESS)
    return -1;
  if (MPI_Barrier(setup) != MPI_SUCCESS) {
    MPI_Comm_free(&setup);
    return -1;
  }
  if (rank == 1)
    cpmdc_note_stop(1);
  agreed = cpmdc_agree_cold_stop();
  if (agreed != 0) {
    if (MPI_Comm_free(&setup) != MPI_SUCCESS)
      return -1;
    return 1;
  }
  if (MPI_Barrier(setup) != MPI_SUCCESS) {
    MPI_Comm_free(&setup);
    return -1;
  }
  if (MPI_Comm_free(&setup) != MPI_SUCCESS)
    return -1;
  return 0;
#else
  return 0;
#endif
}
#endif

static CPMDCResult
energy_gradient_cell_with_params(const void *params_bytes, size_t params_size,
                                 const CPMDCScalarOverrides *overrides,
                                 int n_atoms, const double *positions_ang,
                                 const int *atomic_numbers,
                                 const double *cell_ang, int has_cell,
                                 double *grad_h_bohr, CPMDCEmbedImage *image,
                                 CPMDCSession *session) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  r.message[0] = '\0';
  g_embed_detail[0] = '\0';
  g_embed_code = CPMDC_ERROR_EVALUATION_FAILED;
  if (n_atoms > (int)(sizeof(image->prop.hessian) /
                      (3u * sizeof(image->prop.hessian[0]))))
    return fail_image(image,
                      "property gradient snapshot supports at most 1365 atoms");
  if (grad_len_from_atoms(n_atoms, NULL) != 0 || !positions_ang ||
      !atomic_numbers || !grad_h_bohr)
    return fail_image(image, "invalid arguments");
  {
    int unavailable = !ensure_embed_init() || !cpmdc_embed_available();
    if (cpmdc_orbital_allreduce(&unavailable, 1) != 0)
      return fail_image(image, "MPI collective failed");
    if (unavailable)
      return fail_image(image, "CPMD embed not available");
  }
  double cell[9] = {0};
  if (has_cell && cell_ang)
    memcpy(cell, cell_ang, sizeof(cell));
  double energy = 0.0;
  /* Geometry and forces travel only as C arrays into the embed (nwchemc
   * pattern). Method knobs live in embed module state from set_config /
   * session configure — do not materialize an INPUT deck per force. */
  (void)params_bytes;
  (void)params_size;
  (void)overrides;
  cpmdc_stop_arm();
#if defined(CPMDC_ENABLE_TEST_FAULTS)
  if (test_fault("cold-setup-stop")) {
    int drill = cpmdc_cold_setup_stop_drill();
    if (drill < 0) {
      cpmdc_stop_disarm();
      return fail_image(image, "CPMDC_RANK_CONSENSUS: an MPI collective failed");
    }
    if (drill > 0) {
      int stop = cpmdc_stop_code();
      cpmdc_stop_disarm();
      if (stop == 0)
        cpmdc_note_stop(1);
      (void)cpmdc_embed_reset_state(image);
      image_invalidate_results(image);
      return fail_msg("CPMD stopgm during embed SCF");
    }
  }
#endif
  int n_external =
      session && session->external_count <= (size_t)INT_MAX
          ? (int)session->external_count
          : 0;
  /* A non-optional Fortran assumed-size array must have an interoperable
   * C address even when its logical extent is zero. */
  double external_dummy = 0.0;
  const double *external_positions =
      n_external > 0 ? session->external_positions_ang : &external_dummy;
  const double *external_charges =
      n_external > 0 ? session->external_charges_e : &external_dummy;
  const double *external_widths =
      n_external > 0 ? session->external_widths_ang : &external_dummy;
  double *external_forces =
      n_external > 0 ? session->external_forces_h_bohr : &external_dummy;
  if (session)
    session->external_forces_valid = 0;
  int ok = cpmdc_embed_energy_grad(
      n_atoms, positions_ang, atomic_numbers, cell, has_cell ? 1 : 0, &energy,
      grad_h_bohr, n_external, external_positions, external_charges,
      external_widths, external_forces, image);
  int stop = cpmdc_stop_code();
  cpmdc_stop_disarm();
  if (rank_consensus(&ok, &stop, (size_t)n_atoms * 3u, &energy, grad_h_bohr,
                     image) != 0)
    return fail_image(image, "CPMDC_RANK_CONSENSUS: an MPI collective failed");
  if (stop != 0) {
    /* stopgm returned on this rank. The wavefunction, forces and module
     * state of this call are undefined. Drop them so the next call sets
     * CPMD up again. The caller handles the error. This rank does not
     * abort the process. */
    (void)cpmdc_embed_reset_state(image);
    image_invalidate_results(image);
    return fail_msg("CPMD stopgm during embed SCF");
  }
  if (!ok) {
    const char *detail =
        g_embed_detail[0]
            ? g_embed_detail
            : "CPMD energy/gradient failed: orbitals not converged within "
              "MAXITER, or no energy";
    return fail_image(image, detail);
  }
  /* CPMDC_TEST_FAULT=charge-snapshot records one charge per atom when the
   * evaluator stored none, and a dipole when none was recomputed.
   * OpenCPMD keeps the charges and dipole it already wrote. */
  if (test_fault("charge-snapshot") && n_atoms > 0 &&
      n_atoms <= CPMDC_EMBED_CHARGE_CAP) {
    if (!image->charges.valid) {
      int atom;
      for (atom = 0; atom < n_atoms; ++atom)
        image->charges.values[atom] = 0.01 * (double)atomic_numbers[atom];
      image->charges.count = (size_t)n_atoms;
      image->charges.valid = 1;
    }
    if (image->prop.valid && image->prop.dipole_count != 3) {
      image->prop.dipole[0] = 0.11;
      image->prop.dipole[1] = -0.22;
      image->prop.dipole[2] = 0.33;
      image->prop.dipole_count = 3;
    }
  }
  cpmdc_orbital_note_stored();
  r.ok = 1;
  r.energy_h = energy;
  if (session)
    session->external_forces_valid = 1;
  snprintf(r.message, sizeof(r.message), "ok");
  cpmdc_store_error("");
  return r;
}

#define CPMDC_ATOMIC_NUMBER_MAX 118

/* Zero when an entry is outside 1..CPMDC_ATOMIC_NUMBER_MAX. A NULL array is
 * left to the argument checks of the caller. */
static int atomic_numbers_in_range(int n_atoms, const int *atomic_numbers) {
  if (!atomic_numbers)
    return 1;
  for (int i = 0; i < n_atoms; ++i) {
    if (atomic_numbers[i] < 1 || atomic_numbers[i] > CPMDC_ATOMIC_NUMBER_MAX)
      return 0;
  }
  return 1;
}

/* Agrees the orbital slot limit across the calculator ranks once the checks
 * of this call have run, and always runs the collective so that a check
 * failing on one rank does not leave the others in it. `local_err` is the
 * reason a check of this rank failed, or NULL. Returns the reason the call
 * fails, or NULL when it goes on. A reason of this rank comes first. */
static const char *orbital_limit_gate(const char *local_err) {
  const char *agreed_err = collective_call_gate(local_err);
  const char *limit_err = NULL;
  if (cpmdc_orbital_agree(agreed_err != NULL, &limit_err) == 0)
    return NULL;
  if (agreed_err)
    return agreed_err;
  return limit_err ? limit_err : "CPMDC_ORBITAL_SLOT_LIMIT could not be agreed";
}

/* Geometry the array one-shots can reject before cpmdc_set_params().
 * Nothing here installs parameters, detaches the active session, clears
 * the orbital key, or drops a stored result. The reference evaluator
 * folds a non-finite coordinate into the energy and still reports
 * success, so that input is refused with the session preflight text.
 * These wrappers do not take a cell. */
static const char *array_geometry_error(int n_atoms, const double *positions_ang,
                                        const int *atomic_numbers,
                                        const double *grad_h_bohr) {
  if (grad_len_from_atoms(n_atoms, NULL) != 0 || !positions_ang ||
      !atomic_numbers || !grad_h_bohr)
    return "invalid arguments";
  if (!atomic_numbers_in_range(n_atoms, atomic_numbers))
    return "atomic numbers must be in 1..118";
  {
    int ncomp = n_atoms * 3;
    int i;
    for (i = 0; i < ncomp; ++i) {
      if (!isfinite(positions_ang[i]))
        return "positions contain a non-finite component";
    }
  }
  if (n_atoms > (int)(sizeof(g_image.prop.hessian) /
                      (3u * sizeof(g_image.prop.hessian[0]))))
    return "property gradient snapshot supports at most 1365 atoms";
  return NULL;
}

CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  double *grad_h_bohr) {
  const char *err = array_geometry_error(n_atoms, positions_ang, atomic_numbers,
                                         grad_h_bohr);
  /* Agree the geometry before installing parameters. A rank that rejects
   * the input must not leave the others inside cpmdc_set_params(), and a
   * refusal must not detach the session or drop a stored result. */
  err = orbital_limit_gate(err);
  if (err)
    return fail_msg(err);
  if (cpmdc_set_params(params_capnp, params_capnp_size_bytes) != 0) {
    err = cpmdc_last_error();
    if (!err || err[0] == '\0')
      err = "embed config failed";
  }
  err = orbital_limit_gate(err);
  if (err)
    return fail_msg(err);
  return energy_gradient_cell_with_params(
      g_params_bytes, g_params_size, g_has_overrides ? &g_overrides : NULL,
      n_atoms, positions_ang, atomic_numbers, NULL, 0, grad_h_bohr, &g_image,
      NULL);
}

CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                         const int *atomic_numbers, const void *params_capnp,
                         size_t params_capnp_size_bytes) {
  size_t n3 = 0;
  double *grad = NULL;
  const char *err = NULL;
  if (grad_len_from_atoms(n_atoms, &n3) != 0)
    err = "invalid arguments";
  else {
    grad = test_fault("energy-allocation")
               ? NULL
               : (double *)calloc(n3, sizeof(double));
    if (!grad)
      err = "out of memory";
  }
  err = collective_call_gate(err);
  if (err) {
    free(grad);
    /* The gradient was never evaluated, so the previous result stays. */
    return fail_msg(err);
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
    size_t n3 = 0;
    if (grad_len_from_atoms(n_atoms, &n3) != 0) {
      image_invalidate_results(&g_image);
      return fail_msg("invalid arguments");
    }
    for (size_t i = 0; i < n3; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
  }
  return r;
}

static int session_reserve_step_atoms(CPMDCSession *session, size_t n_atoms) {
  if (!session)
    return -1;
  if (session->step_atom_capacity >= n_atoms)
    return 0;
  if (n_atoms > SIZE_MAX / (3u * sizeof(double)))
    return -1;
  double *pos = (double *)malloc(n_atoms * 3u * sizeof(double));
  int *z = (int *)malloc(n_atoms * sizeof(int));
  if (!pos || !z) {
    free(pos);
    free(z);
    return -1;
  }
  if (session->step_atom_capacity > 0 && session->step_positions_ang &&
      session->step_atomic_numbers) {
    memcpy(pos, session->step_positions_ang,
           session->step_atom_capacity * 3u * sizeof(double));
    memcpy(z, session->step_atomic_numbers,
           session->step_atom_capacity * sizeof(int));
  }
  free(session->step_positions_ang);
  free(session->step_atomic_numbers);
  session->step_positions_ang = pos;
  session->step_atomic_numbers = z;
  session->step_atom_capacity = n_atoms;
  return 0;
}

static int element_counts_match(size_t n_atoms, const int *left,
                                const int *right) {
  int counts[119];
  if (!left || !right)
    return 0;
  memset(counts, 0, sizeof(counts));
  for (size_t i = 0; i < n_atoms; ++i) {
    if (left[i] < 1 || left[i] > 118)
      return 0;
    counts[left[i]] += 1;
  }
  for (size_t i = 0; i < n_atoms; ++i) {
    if (right[i] < 1 || right[i] > 118)
      return 0;
    counts[right[i]] -= 1;
    if (counts[right[i]] < 0)
      return 0;
  }
  return 1;
}

static int coords_differ(const double *left, const double *right) {
  for (int k = 0; k < 3; ++k) {
    double delta;
    /* NaN minus anything is not ordered, so a magnitude test would
     * treat a non-finite coordinate as equal. */
    if (!isfinite(left[k]) || !isfinite(right[k]))
      return 1;
    delta = left[k] - right[k];
    if (delta < 0.0)
      delta = -delta;
    if (delta > 1.0e-8)
      return 1;
  }
  return 0;
}

static int cell_differs(int left_has, const double *left, int right_has,
                        const double *right) {
  if (left_has != right_has)
    return 1;
  if (!left_has)
    return 0;
  if (!left || !right)
    return 1;
  for (int i = 0; i < 9; ++i) {
    double delta;
    /* NaN minus anything is not ordered, so a magnitude test would
     * treat a non-finite component as equal. */
    if (!isfinite(left[i]) || !isfinite(right[i]))
      return 1;
    delta = left[i] - right[i];
    if (delta < 0.0)
      delta = -delta;
    if (delta > 1.0e-8)
      return 1;
  }
  return 0;
}

static void session_drop_cache(CPMDCSession *session) {
  if (!session)
    return;
  free(session->cached_positions_ang);
  free(session->cached_atomic_numbers);
  free(session->cached_grad_h_bohr);
  free(session->cached_charges_e);
  session->cached_positions_ang = NULL;
  session->cached_atomic_numbers = NULL;
  session->cached_grad_h_bohr = NULL;
  session->cached_charges_e = NULL;
  session->cached_charges_valid = 0;
  session->cached_charge_count = 0;
  session->cached_n_atoms = 0;
  session->result_cached = 0;
}

static int session_store_cache(CPMDCSession *session, size_t n_atoms,
                               const double *positions_ang,
                               const int *atomic_numbers, const double *grad,
                               double energy_h, const double *cell_ang,
                               int has_cell) {
  if (!session || !positions_ang || !atomic_numbers || !grad || n_atoms == 0)
    return -1;
  if (!isfinite(energy_h))
    return -1;
  for (size_t i = 0; i < n_atoms * 3u; ++i) {
    if (!isfinite(positions_ang[i]) || !isfinite(grad[i]))
      return -1;
  }
  if (session->cached_n_atoms != n_atoms) {
    double *pos = (double *)malloc(n_atoms * 3u * sizeof(double));
    int *z = (int *)malloc(n_atoms * sizeof(int));
    double *g = (double *)malloc(n_atoms * 3u * sizeof(double));
    double *q = (double *)malloc(n_atoms * sizeof(double));
    if (!pos || !z || !g || !q) {
      free(pos);
      free(z);
      free(g);
      free(q);
      return -1;
    }
    free(session->cached_positions_ang);
    free(session->cached_atomic_numbers);
    free(session->cached_grad_h_bohr);
    free(session->cached_charges_e);
    session->cached_positions_ang = pos;
    session->cached_atomic_numbers = z;
    session->cached_grad_h_bohr = g;
    session->cached_charges_e = q;
    session->cached_n_atoms = n_atoms;
  }
  memcpy(session->cached_positions_ang, positions_ang,
         n_atoms * 3u * sizeof(double));
  memcpy(session->cached_atomic_numbers, atomic_numbers,
         n_atoms * sizeof(int));
  memcpy(session->cached_grad_h_bohr, grad, n_atoms * 3u * sizeof(double));
  if (session->image.charges.valid &&
      session->image.charges.count == n_atoms && session->cached_charges_e) {
    memcpy(session->cached_charges_e, session->image.charges.values,
           n_atoms * sizeof(double));
    session->cached_charges_valid = 1;
    session->cached_charge_count = n_atoms;
  } else {
    session->cached_charges_valid = 0;
    session->cached_charge_count = 0;
  }
  session->cached_energy_h = energy_h;
  session->cached_has_cell = has_cell;
  if (has_cell && cell_ang)
    memcpy(session->cached_cell_ang, cell_ang, sizeof(session->cached_cell_ang));
  else
    memset(session->cached_cell_ang, 0, sizeof(session->cached_cell_ang));
  session->result_cached = 1;
  return 0;
}

/* The stored gradient is copied into the caller's atom order, including
 * the property snapshot read by cpmdc_session_last_gradient. Two atoms
 * match when the element and the position agree within 1e-8 Angstrom. */
static int session_copy_cached_gradient(CPMDCSession *session, size_t n_atoms,
                                        const double *positions,
                                        const int *atomic_numbers,
                                        const double *cell_ang, int has_cell,
                                        double *grad) {
  unsigned char *used;
  size_t *map;
  if (!session || !session->result_cached || !grad || !positions ||
      !atomic_numbers || session->cached_n_atoms != n_atoms)
    return 0;
  if (cell_differs(session->cached_has_cell, session->cached_cell_ang, has_cell,
                   cell_ang))
    return 0;
  used = (unsigned char *)calloc(n_atoms, 1);
  map = (size_t *)malloc(n_atoms * sizeof(*map));
  if (!used || !map) {
    free(used);
    free(map);
    return 0;
  }
  for (size_t i = 0; i < n_atoms; ++i) {
    size_t found = n_atoms;
    for (size_t j = 0; j < n_atoms; ++j) {
      if (used[j])
        continue;
      if (session->cached_atomic_numbers[j] != atomic_numbers[i])
        continue;
      if (coords_differ(session->cached_positions_ang + 3u * j,
                        positions + 3u * i))
        continue;
      found = j;
      break;
    }
    if (found == n_atoms) {
      free(used);
      free(map);
      return 0;
    }
    used[found] = 1;
    map[i] = found;
    grad[3u * i] = session->cached_grad_h_bohr[3u * found];
    grad[3u * i + 1u] = session->cached_grad_h_bohr[3u * found + 1u];
    grad[3u * i + 2u] = session->cached_grad_h_bohr[3u * found + 2u];
  }
  if (session->cached_charges_valid &&
      session->cached_charge_count == n_atoms && session->cached_charges_e) {
    for (size_t i = 0; i < n_atoms; ++i)
      session->image.charges.values[i] = session->cached_charges_e[map[i]];
    session->image.charges.valid = 1;
    session->image.charges.count = n_atoms;
  } else {
    session->image.charges.valid = 0;
    session->image.charges.count = 0;
  }
  /* The snapshot is the gradient in the order of the call that stored
   * it. Rewrite it from that stored gradient so a later hit, which
   * still matches against the stored order, does not permute twice. */
  if (session->image.prop.valid &&
      session->image.prop.hessian_count == n_atoms * 3u) {
    for (size_t i = 0; i < n_atoms; ++i) {
      size_t src = map[i];
      session->image.prop.hessian[3u * i] =
          session->cached_grad_h_bohr[3u * src];
      session->image.prop.hessian[3u * i + 1u] =
          session->cached_grad_h_bohr[3u * src + 1u];
      session->image.prop.hessian[3u * i + 2u] =
          session->cached_grad_h_bohr[3u * src + 2u];
    }
  }
  free(used);
  free(map);
  return 1;
}

static int session_check_topology(const CPMDCSession *session, size_t n_atoms,
                                  const int *atomic_numbers) {
  if (!session || !atomic_numbers || n_atoms == 0)
    return -1;
  if (!session->topology_fixed)
    return 0;
  if (session->fixed_n_atoms != n_atoms)
    return -1;
  if (!element_counts_match(n_atoms, session->fixed_atomic_numbers,
                            atomic_numbers))
    return -1;
  return 0;
}

/* The overrides a deck is rendered with: an overlay and the session's
 * constant cutoff. NULL when there is neither. */
static const CPMDCScalarOverrides *
render_overrides_with(const CPMDCSession *session,
                      const CPMDCScalarOverrides *overlay,
                      CPMDCScalarOverrides *scratch) {
  if (session->constant_cutoff_ry <= 0.0)
    return overlay;
  if (overlay)
    *scratch = *overlay;
  else
    memset(scratch, 0, sizeof(*scratch));
  scratch->constant_cutoff_ry = session->constant_cutoff_ry;
  scratch->constant_cutoff_height_ry = session->constant_cutoff_height_ry;
  scratch->constant_cutoff_width_ry = session->constant_cutoff_width_ry;
  return scratch;
}

/* The overrides a session renders its deck with: the common overlay and the
 * constant cutoff. NULL when it has neither. */
static const CPMDCScalarOverrides *
session_render_overrides(const CPMDCSession *session,
                         CPMDCScalarOverrides *scratch) {
  return render_overrides_with(
      session, session->has_overrides ? &session->overrides : NULL, scratch);
}

static CPMDCSession *cpmdc_session_create_ov(
    const void *params_capnp, size_t params_capnp_size_bytes,
    const CPMDCScalarOverrides *overrides) {
  cpmdc_store_error("");
  if (!params_capnp || params_capnp_size_bytes == 0) {
    cpmdc_store_error("CPMDParams buffer is empty");
    return NULL;
  }
  CPMDCSession *session = (CPMDCSession *)calloc(1, sizeof(*session));
  if (!session) {
    cpmdc_store_error("out of memory");
    return NULL;
  }
  session->params_bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!session->params_bytes) {
    free(session);
    cpmdc_store_error("out of memory");
    return NULL;
  }
  memcpy(session->params_bytes, params_capnp, params_capnp_size_bytes);
  session->params_size = params_capnp_size_bytes;
  if (overrides) {
    session->overrides = *overrides;
    if (overrides->functional) {
      snprintf(session->overrides_functional,
               sizeof(session->overrides_functional), "%s",
               overrides->functional);
      session->overrides.functional = session->overrides_functional;
    }
    session->has_overrides = 1;
  }
  if (apply_params_buffer(session->params_bytes, session->params_size,
                          session->has_overrides ? &session->overrides : NULL,
                          session->functional, sizeof(session->functional),
                          &session->cutoff_ry, &session->charge,
                          &session->multiplicity, &session->input_deck,
                          session->cpmd_root, sizeof(session->cpmd_root)) != 0) {
    free(session->input_deck);
    free(session->params_bytes);
    free(session);
    return NULL;
  }
  (void)configure_embed_from_session(session);
  cpmdc_store_error("");
  g_live_sessions++;
  return session;
}

CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                   size_t params_capnp_size_bytes) {
  return cpmdc_session_create_ov(params_capnp, params_capnp_size_bytes, NULL);
}

/* Replaces the parameters and the overlay together. Everything is validated
 * into temporaries first, so a call that fails leaves the session as it was:
 * parameters, deck, effective scalars, overlay, and embed image. Only the
 * last error changes. */
static int cpmdc_session_set_params_ov(CPMDCSession *session,
                                       const void *params_capnp,
                                       size_t params_capnp_size_bytes,
                                       const CPMDCScalarOverrides *overrides) {
  if (!session) {
    cpmdc_store_error("invalid session");
    return -1;
  }
  if (!params_capnp || params_capnp_size_bytes == 0) {
    cpmdc_store_error("CPMDParams buffer is empty");
    return -1;
  }
  if (session->topology_fixed) {
    cpmdc_store_error("session topology is already fixed");
    return -1;
  }
  unsigned char *bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!bytes) {
    cpmdc_store_error("out of memory");
    return -1;
  }
  CPMDCScalarOverrides overlay;
  char overlay_functional[sizeof(session->overrides_functional)];
  memset(&overlay, 0, sizeof(overlay));
  if (overrides) {
    overlay = *overrides;
    if (overrides->functional) {
      snprintf(overlay_functional, sizeof(overlay_functional), "%s",
               overrides->functional);
      overlay.functional = overlay_functional;
    }
  }
  char functional[sizeof(session->functional)];
  char cpmd_root[sizeof(session->cpmd_root)];
  double cutoff_ry = session->cutoff_ry;
  int charge = session->charge;
  int multiplicity = session->multiplicity;
  memcpy(functional, session->functional, sizeof(functional));
  memcpy(cpmd_root, session->cpmd_root, sizeof(cpmd_root));
  char *deck = NULL;
  CPMDCScalarOverrides render;
  if (apply_params_buffer(
          params_capnp, params_capnp_size_bytes,
          render_overrides_with(session, overrides ? &overlay : NULL, &render),
          functional, sizeof(functional), &cutoff_ry, &charge, &multiplicity,
          &deck, cpmd_root, sizeof(cpmd_root)) != 0) {
    free(deck);
    free(bytes);
    return -1;
  }
  memcpy(bytes, params_capnp, params_capnp_size_bytes);

  const CPMDCSession before = *session;
  memcpy(session->functional, functional, sizeof(functional));
  memcpy(session->cpmd_root, cpmd_root, sizeof(cpmd_root));
  session->cutoff_ry = cutoff_ry;
  session->charge = charge;
  session->multiplicity = multiplicity;
  session->params_bytes = bytes;
  session->params_size = params_capnp_size_bytes;
  session->input_deck = deck;
  if (overrides) {
    session->overrides = overlay;
    if (overrides->functional) {
      snprintf(session->overrides_functional,
               sizeof(session->overrides_functional), "%s",
               overlay.functional);
      session->overrides.functional = session->overrides_functional;
    }
    session->has_overrides = 1;
  } else {
    session->has_overrides = 0;
  }
  if (configure_embed_from_session(session) != 0) {
    *session = before;
    free(deck);
    free(bytes);
    cpmdc_store_error("CPMD embed not available");
    return -1;
  }
  free(before.params_bytes);
  free(before.input_deck);
  cpmdc_store_error("");
  return 0;
}

int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes) {
  return cpmdc_session_set_params_ov(session, params_capnp,
                                     params_capnp_size_bytes, NULL);
}

int cpmdc_session_set_constant_cutoff(CPMDCSession *session, double cutoff_ry,
                                      double height_ry, double width_ry) {
  if (!session) {
    cpmdc_store_error("invalid session");
    return -1;
  }
  if (!isfinite(cutoff_ry) || !isfinite(height_ry) || !isfinite(width_ry)) {
    cpmdc_store_error("constant cutoff values must be finite");
    return -1;
  }
  if (session->topology_fixed) {
    cpmdc_store_error("session topology is already fixed");
    return -1;
  }
  double saved[3] = {session->constant_cutoff_ry,
                     session->constant_cutoff_height_ry,
                     session->constant_cutoff_width_ry};
  int on = cutoff_ry > 0.0;
  session->constant_cutoff_ry = on ? cutoff_ry : 0.0;
  session->constant_cutoff_height_ry = on && height_ry > 0.0 ? height_ry : 0.0;
  session->constant_cutoff_width_ry = on && width_ry > 0.0 ? width_ry : 0.0;
  char *deck = NULL;
  CPMDCScalarOverrides render;
  if (apply_params_buffer(session->params_bytes, session->params_size,
                          session_render_overrides(session, &render),
                          session->functional, sizeof(session->functional),
                          &session->cutoff_ry, &session->charge,
                          &session->multiplicity, &deck, session->cpmd_root,
                          sizeof(session->cpmd_root)) != 0) {
    free(deck);
    session->constant_cutoff_ry = saved[0];
    session->constant_cutoff_height_ry = saved[1];
    session->constant_cutoff_width_ry = saved[2];
    return -1;
  }
  free(session->input_deck);
  session->input_deck = deck;
  session->embed_configured = 0;
  cpmdc_store_error("");
  return 0;
}

/* A Gaussian narrower than this has 1/sigma^4 force terms that no grid or
 * double precision sum can resolve, and far below it they overflow to
 * infinity and then to NaN. The width is refused, not clamped. */
#define CPMDC_EXTERNAL_WIDTH_MIN_ANG 1.0e-12
#define CPMDC_EXTERNAL_WIDTH_DEFAULT_ANG 0.30

static double external_default_width_ang(void) {
  const char *text = getenv("CPMDC_EXTERNAL_CHARGE_WIDTH_ANG");
  char *end = NULL;
  double value = CPMDC_EXTERNAL_WIDTH_DEFAULT_ANG;
  if (!text || text[0] == '\0')
    return value;
  value = strtod(text, &end);
  return end != text && *end == '\0' && isfinite(value) &&
                 value >= CPMDC_EXTERNAL_WIDTH_MIN_ANG
             ? value
             : -1.0;
}

int cpmdc_session_set_external_charges(
    CPMDCSession *session, size_t n_charges, const double *positions_ang,
    const double *charges_e, const double *widths_ang) {
  cpmdc_store_error("");
  if (!session || n_charges > (size_t)(INT_MAX / 3) ||
      (n_charges > 0 && (!positions_ang || !charges_e))) {
    cpmdc_store_error_code("invalid external charge arguments",
                           CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  double default_width = external_default_width_ang();
  if (n_charges > 0 && !widths_ang && default_width < 0.0) {
    cpmdc_store_error_code(
        "CPMDC_EXTERNAL_CHARGE_WIDTH_ANG must be finite and at least "
        "1e-12 Angstrom",
        CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  double *positions = NULL, *charges = NULL, *widths = NULL, *forces = NULL;
  if (n_charges > 0) {
    positions = (double *)malloc(3u * n_charges * sizeof(*positions));
    charges = (double *)malloc(n_charges * sizeof(*charges));
    widths = (double *)malloc(n_charges * sizeof(*widths));
    forces = (double *)calloc(3u * n_charges, sizeof(*forces));
    if (!positions || !charges || !widths || !forces) {
      free(positions);
      free(charges);
      free(widths);
      free(forces);
      cpmdc_store_error_code("out of memory",
                             CPMDC_ERROR_EVALUATION_FAILED);
      return -1;
    }
    for (size_t i = 0; i < 3u * n_charges; ++i) {
      if (!isfinite(positions_ang[i])) {
        free(positions);
        free(charges);
        free(widths);
        free(forces);
        cpmdc_store_error_code("external charge positions must be finite",
                               CPMDC_ERROR_INVALID_ARGUMENT);
        return -1;
      }
      positions[i] = positions_ang[i];
    }
    for (size_t i = 0; i < n_charges; ++i) {
      double width = widths_ang ? widths_ang[i] : default_width;
      if (!isfinite(charges_e[i]) || !isfinite(width) ||
          width < CPMDC_EXTERNAL_WIDTH_MIN_ANG) {
        free(positions);
        free(charges);
        free(widths);
        free(forces);
        cpmdc_store_error_code(
            "external charges must be finite and smearing widths finite and "
            "at least 1e-12 Angstrom",
            CPMDC_ERROR_INVALID_ARGUMENT);
        return -1;
      }
      charges[i] = charges_e[i];
      widths[i] = width;
    }
  }
  free(session->external_positions_ang);
  free(session->external_charges_e);
  free(session->external_widths_ang);
  free(session->external_forces_h_bohr);
  session->external_positions_ang = positions;
  session->external_charges_e = charges;
  session->external_widths_ang = widths;
  session->external_forces_h_bohr = forces;
  session->external_count = n_charges;
  session->external_forces_valid = 0;
  session_drop_cache(session);
  cpmdc_store_error("");
  return 0;
}

int cpmdc_session_last_external_charge_forces(
    const CPMDCSession *session, double *forces_h_bohr, size_t forces_len) {
  if (!session) {
    cpmdc_store_error_code("null session", CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  if (!session->external_forces_valid) {
    cpmdc_store_error_code("external charge forces are not available",
                           CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  if (forces_len != 3u * session->external_count) {
    cpmdc_store_error_code("external charge force buffer has the wrong length",
                           CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  if (forces_len > 0 && !forces_h_bohr) {
    cpmdc_store_error_code("external charge force buffer is null",
                           CPMDC_ERROR_INVALID_ARGUMENT);
    return -1;
  }
  if (forces_len > 0)
    memcpy(forces_h_bohr, session->external_forces_h_bohr,
           forces_len * sizeof(*forces_h_bohr));
  cpmdc_store_error("");
  return 0;
}

void cpmdc_session_destroy(CPMDCSession *session) {
  if (!session)
    return;
  if (g_live_sessions > 0)
    g_live_sessions--;
  if (session->evaluated && g_evaluated_sessions > 0)
    g_evaluated_sessions--;
  if (g_orbital_owner == session) {
    cpmdc_orbital_clear_key(1);
    cpmdc_orbital_ledger_reset();
    g_orbital_owner = NULL;
  }
  if (g_active_session == session)
    g_active_session = NULL;
  free(session->params_bytes);
  if (session->image.input_deck != session->input_deck)
    free(session->image.input_deck);
  free(session->input_deck);
  free(session->fixed_atomic_numbers);
  free(session->step_positions_ang);
  free(session->step_atomic_numbers);
  free(session->external_positions_ang);
  free(session->external_charges_e);
  free(session->external_widths_ang);
  free(session->external_forces_h_bohr);
  session_drop_cache(session);
  free(session);
}

static CPMDCResult session_energy_gradient_cell(
    CPMDCSession *session, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const double *cell_ang, int has_cell,
    double *grad_h_bohr) {
  const char *err = NULL;
  cpmdc_store_error("");
  /* Every check that can fail runs before the collective that agrees the
   * orbital slot limit. A rank that fails one still joins that collective,
   * so the ranks that passed do not wait for it, and its own reason is the
   * one it reports. */
  if (!session)
    err = "null session";
  /* A missing buffer is a bad call, including on a session that already
   * stored a result. The cache lookup reads these pointers. */
  else if (grad_len_from_atoms(n_atoms, NULL) != 0 || !positions_ang ||
           !atomic_numbers || !grad_h_bohr || (has_cell && !cell_ang))
    err = "invalid arguments";
  else if (!atomic_numbers_in_range(n_atoms, atomic_numbers))
    err = "atomic numbers must be in 1..118";
  if (!err && has_cell) {
    int cell_state = cpmdc_cell_rank3_state(cell_ang);
    if (cell_state == CPMDC_CELL_NONFINITE)
      err = "cell contains a non-finite component";
    else if (cell_state == CPMDC_CELL_RANK_DEFICIENT)
      err = "cell does not have three independent vectors";
    else if (cell_state == CPMDC_CELL_VOLUME_UNDERFLOW)
      err = "cell volume is too small for binary64";
  }
  if (!err) {
    int ncomp = n_atoms * 3;
    for (int i = 0; i < ncomp; ++i) {
      if (!isfinite(positions_ang[i])) {
        err = "positions contain a non-finite component";
        break;
      }
    }
  }
  if (!err &&
      session_check_topology(session, (size_t)n_atoms, atomic_numbers) != 0)
    err = "topology change requires a new session";
  /* After the topology check, so a session that already holds another
   * system reports the topology change. The topology is latched only
   * after this preflight, so a first call that is too large leaves the
   * session without a topology and a valid call can follow. */
  else if (n_atoms > (int)(sizeof(session->image.prop.hessian) /
                           (3u * sizeof(session->image.prop.hessian[0]))))
    err = "property gradient snapshot supports at most 1365 atoms";
  /* Decided before the collective gates, from state every rank holds. */
  if (!err && session->external_count > 0) {
    if (deck_or_cell_is_periodic(session->input_deck, cell_ang, has_cell))
      err = k_periodic_embedding_msg;
#if !defined(CPMDC_HAS_CPMD)
    if (!err &&
        cpmdc_deck_defines_external_potential(
            session->input_deck,
            session->input_deck ? strlen(session->input_deck) : 0u))
      err = k_deck_external_potential_msg;
#endif
  }
  err = orbital_limit_gate(err);
  if (err) {
    session_drop_results(session);
    return fail_msg(err);
  }
  /* The copy is committed only after every preflight collective succeeds, so
   * a key, stale-reset, or prepare failure leaves the topology open. */
  int *pending_topology = NULL;
  if (!session->topology_fixed) {
    pending_topology = test_fault("topology-allocation")
                           ? NULL
                           : (int *)malloc((size_t)n_atoms * sizeof(*pending_topology));
    if (!pending_topology)
      err = "out of memory";
    else
      memcpy(pending_topology, atomic_numbers,
             (size_t)n_atoms * sizeof(*pending_topology));
  }
  err = collective_call_gate(err);
  if (err) {
    free(pending_topology);
    session_drop_results(session);
    return fail_msg(err);
  }
  if (rank_consensus_all(
          grad_h_bohr &&
          session_copy_cached_gradient(session, (size_t)n_atoms, positions_ang,
                                       atomic_numbers, cell_ang, has_cell,
                                       grad_h_bohr))) {
    if (pending_topology) {
      session->fixed_atomic_numbers = pending_topology;
      session->fixed_n_atoms = (size_t)n_atoms;
      session->topology_fixed = 1;
    } else {
      memcpy(session->fixed_atomic_numbers, atomic_numbers,
             (size_t)n_atoms * sizeof(*atomic_numbers));
    }
    /* The cache holds a successful step, and replacing the charges
     * drops it, so the force buffer still holds that step. A failed
     * call in between may have cleared the flag. */
    session->external_forces_valid = 1;
    CPMDCResult cached;
    cached.ok = 1;
    cached.energy_h = session->cached_energy_h;
    cached.message[0] = '\0';
    cpmdc_store_error("");
    return cached;
  }
  session->external_forces_valid = 0;
  if ((!session->embed_configured || g_active_session != session) &&
      configure_embed_from_session(session) != 0)
    err = "CPMD embed not available";
  err = collective_call_gate(err);
  if (err) {
    free(pending_topology);
    session_drop_results(session);
    return fail_msg(err);
  }
  {
    int has_orbital_key = session->has_orbital_key;
    const char *slot_err = NULL;
    if (cpmdc_orbital_key_agree(&has_orbital_key, session->orbital_key,
                                &slot_err) != 0) {
      free(pending_topology);
      session_drop_results(session);
      return fail_msg(slot_err);
    }
    if (has_orbital_key) {
      int stale_failed;
      double cell[9] = {0};
      if (has_cell && cell_ang)
        memcpy(cell, cell_ang, sizeof(cell));
      /* A new basis frees every key's orbitals before this key parks one. */
      stale_failed =
          !cpmdc_embed_drop_stale_orbitals(n_atoms, atomic_numbers, cell,
                                            has_cell ? 1 : 0, &session->image);
      if (test_fault("drop-stale"))
        stale_failed = 1;
      if (cpmdc_orbital_allreduce(&stale_failed, 1) != 0) {
        free(pending_topology);
        session_drop_results(session);
        return fail_msg("MPI collective failed");
      }
      if (stale_failed) {
        free(pending_topology);
        session_drop_results(session);
        return fail_msg("failed to discard stale orbitals");
      }
      if (cpmdc_orbital_prepare(session->orbital_key, &slot_err) != 0) {
        free(pending_topology);
        session_drop_results(session);
        return fail_msg(slot_err);
      }
    } else {
      cpmdc_orbital_clear_key(0);
    }
    g_orbital_owner = session;
  }
  if (pending_topology) {
    session->fixed_atomic_numbers = pending_topology;
    session->fixed_n_atoms = (size_t)n_atoms;
    session->topology_fixed = 1;
  } else {
    memcpy(session->fixed_atomic_numbers, atomic_numbers,
           (size_t)n_atoms * sizeof(*atomic_numbers));
  }
  if (!session->evaluated) {
    session->evaluated = 1;
    g_evaluated_sessions++;
  }
  {
    CPMDCResult evaluated = energy_gradient_cell_with_params(
        session->params_bytes, session->params_size,
        session->has_overrides ? &session->overrides : NULL, n_atoms,
        positions_ang, atomic_numbers, cell_ang, has_cell, grad_h_bohr,
        &session->image, session);
    if (evaluated.ok && grad_h_bohr) {
      if (session_store_cache(session, (size_t)n_atoms, positions_ang,
                              atomic_numbers, grad_h_bohr, evaluated.energy_h,
                              cell_ang, has_cell) != 0)
        session_drop_cache(session);
    } else {
      session_drop_cache(session);
    }
    return evaluated;
  }
}

int cpmdc_session_select_orbitals(CPMDCSession *session, long long key) {
  if (!session) {
    cpmdc_store_error("null session");
    return -1;
  }
  session->has_orbital_key = 1;
  session->orbital_key = key;
  /* Record only. The evaluation selects the key, so a key named and then
   * replaced before an evaluation does not join the warm-start donors. */
  cpmdc_store_error("");
  return 0;
}

int cpmdc_session_copy_orbitals(CPMDCSession *session, long long from_key,
                                long long to_key) {
  const char *err = NULL;
  int has_from_key = 1;
  int has_to_key = 1;
  int available;
  int availability[2];
  if (!session)
    err = "null session";
  else if (!session->evaluated)
    err = "no stored orbitals; evaluate at least once before copy";
  err = orbital_limit_gate(err);
  if (err) {
    cpmdc_store_error(err);
    return -1;
  }
  if (cpmdc_orbital_key_agree(&has_from_key, from_key, &err) != 0 ||
      cpmdc_orbital_key_agree(&has_to_key, to_key, &err) != 0) {
    cpmdc_store_error(err);
    return -1;
  }
  available = cpmdc_available() != 0 || test_fault("embed-available");
  availability[0] = available;
  availability[1] = !available;
  if (cpmdc_orbital_allreduce(availability, 2) != 0) {
    cpmdc_store_error("MPI collective failed");
    return -1;
  }
  if (availability[0] && availability[1]) {
    cpmdc_store_error("CPMD embed availability differs across calculator ranks");
    return -1;
  }
  /* The reference evaluator has no orbital store. An OpenCPMD build
   * duplicates the converged copy and refuses a key that never converged. */
  if (availability[0]) {
    if ((!session->embed_configured || g_active_session != session) &&
        configure_embed_from_session(session) != 0)
      err = "CPMD embed not available";
    err = collective_call_gate(err);
    if (err) {
      cpmdc_store_error(err);
      return -1;
    }
    if (cpmdc_embed_copy_orbitals(from_key, to_key) != 0) {
      const char *why = cpmdc_orbital_limit_message(cpmdc_embed_slot_limit_code());
      err = why ? why : "from_key has no converged orbitals";
    }
    err = collective_call_gate(err);
    if (err) {
      cpmdc_store_error(err);
      return -1;
    }
  }
  cpmdc_store_error("");
  return 0;
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
  size_t n3 = 0;
  double *grad = NULL;
  const char *err = NULL;
  if (grad_len_from_atoms(n_atoms, &n3) != 0)
    err = "invalid arguments";
  else {
    grad = test_fault("session-energy-allocation")
               ? NULL
               : (double *)calloc(n3, sizeof(double));
    if (!grad)
      err = "out of memory";
  }
  err = collective_call_gate(err);
  if (err) {
    free(grad);
    session_drop_results(session);
    return fail_msg(err);
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
    size_t n3 = 0;
    if (grad_len_from_atoms(n_atoms, &n3) != 0) {
      session_drop_results(session);
      return fail_msg("invalid arguments");
    }
    for (size_t i = 0; i < n3; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
  }
  return r;
}

CPMDCResult cpmdc_session_calculate_forces(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, double *forces_h_bohr,
    size_t forces_len) {
  CPMDCResult r;
  const char *err = NULL;
  int arena_open = 0;
  struct capn arena;
  ForceInput_ptr force_input;
  size_t n_atoms = 0;
  int has_cell = 0;
  double cell_ang[9] = {0};
  if (!session || !force_input_capnp || force_input_capnp_size_bytes == 0 ||
      !forces_h_bohr)
    err = "invalid arguments";
  else if (cpmdc_force_input_root(force_input_capnp,
                                  force_input_capnp_size_bytes, &arena,
                                  &force_input) != 0)
    err = "invalid ForceInput message";
  else {
    arena_open = 1;
  }
  if (!err &&
      (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
       n_atoms == 0 || n_atoms > (size_t)(INT_MAX / 3) ||
       forces_len < n_atoms * 3u))
    err = "invalid ForceInput geometry";
  if (!err &&
      (test_fault("calculate-forces-reserve") ||
       session_reserve_step_atoms(session, n_atoms) != 0))
    err = "out of memory";
  if (!err &&
      cpmdc_force_input_copy_geometry(
          force_input, session->step_positions_ang,
          session->step_atomic_numbers, session->step_atom_capacity, cell_ang,
          &has_cell) != 0)
    err = "invalid ForceInput geometry";
  if (arena_open)
    cpmdc_params_release(&arena);
  err = collective_call_gate(err);
  if (err) {
    session_drop_results(session);
    return fail_msg(err);
  }
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
  const char *err = NULL;
  int arena_open = 0;
  struct capn arena;
  ForceInput_ptr force_input;
  size_t n_atoms = 0;
  int has_cell = 0;
  size_t force_count = 0;
  size_t required_size = 0;
  double energy_factor = 1.0;
  double force_factor = 1.0;
  double stress_factor = 1.0;
  double cell_ang[9] = {0};
  double *forces = NULL;
  const char *write_err;
  if (!session || !force_input_capnp || force_input_capnp_size_bytes == 0 ||
      !potential_result_capnp_size_bytes)
    err = "invalid arguments";
  if (potential_result_capnp_size_bytes)
    *potential_result_capnp_size_bytes = 0;
  if (!err &&
      cpmdc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                             &arena, &force_input) != 0)
    err = "invalid ForceInput message";
  else if (!err)
    arena_open = 1;
  if (!err &&
      (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
       n_atoms > (size_t)INT_MAX || n_atoms > SIZE_MAX / 3u))
    err = "invalid ForceInput geometry";
  if (!err) {
    force_count = n_atoms * 3u;
    required_size = cpmdc_potential_result_flat_size(force_count);
    *potential_result_capnp_size_bytes = required_size;
    if (required_size == 0 || force_count > (size_t)INT_MAX)
      err = "invalid ForceInput geometry";
  }
  if (!err && cpmdc_force_input_result_factors(force_input, &energy_factor,
                                                &force_factor) != 0)
    err = "invalid ForceInput result units";
  if (!err &&
      cpmdc_force_input_stress_result_factors(force_input, &energy_factor,
                                              &stress_factor) != 0)
    err = "invalid ForceInput stress units";
  if (!err &&
      (!potential_result_capnp ||
       potential_result_capnp_capacity_bytes < required_size))
    err = "PotentialResult buffer too small";
  if (!err &&
      (test_fault("calculate-result-reserve") ||
       session_reserve_step_atoms(session, n_atoms) != 0))
    err = "out of memory";
  if (!err &&
      cpmdc_force_input_copy_geometry(
          force_input, session->step_positions_ang,
          session->step_atomic_numbers, session->step_atom_capacity, cell_ang,
          &has_cell) != 0)
    err = "invalid ForceInput geometry";
  if (arena_open)
    cpmdc_params_release(&arena);
  if (!err) {
    forces = test_fault("calculate-result-allocation")
                 ? NULL
                 : (double *)malloc(force_count * sizeof(*forces));
    if (!forces)
      err = "out of memory";
  }
  err = collective_call_gate(err);
  if (err) {
    free(forces);
    session_drop_results(session);
    return fail_msg(err);
  }
  r = session_energy_gradient_cell(session, (int)n_atoms,
                                   session->step_positions_ang,
                                   session->step_atomic_numbers, cell_ang,
                                   has_cell, forces);
  if (r.ok && test_fault("embed-energy-gradient"))
    r = fail_msg("embedded energy gradient failed");
  write_err = r.ok ? NULL : r.message;
  if (r.ok) {
    for (size_t i = 0; i < force_count; ++i)
      forces[i] = -forces[i] * force_factor;
    if (test_fault("potential-result-write") ||
        cpmdc_potential_result_write(
            r.energy_h * energy_factor, forces, force_count, stress_factor,
            &session->image, potential_result_capnp,
            potential_result_capnp_capacity_bytes,
            potential_result_capnp_size_bytes) != 0)
      write_err = "PotentialResult write failed";
  }
  write_err = collective_call_gate(write_err);
  if (r.ok && write_err) {
    r.ok = 0;
    snprintf(r.message, sizeof(r.message), "%s", write_err);
    cpmdc_store_error(r.message);
  }
  if (!r.ok)
    session_drop_results(session);
  free(forces);
  return r;
}

/* Positions and the cell go to Angstrom with the factor a ForceInput in
 * "bohr" gets, so this step and that message hand the embed the same bits. */
int cpmdc_session_calculate_au(CPMDCSession *session, int n_atoms,
                               const double *positions_bohr,
                               const int *atomic_numbers,
                               const double *cell_bohr, double *energy_h,
                               double *forces_h_bohr, double *virial_h) {
  size_t n3 = 0;
  const char *err = NULL;
  double cell_ang[9] = {0};
  int has_cell = cell_bohr != NULL;
  double volume = 0.0;
  if (!session)
    err = "null session";
  else if (grad_len_from_atoms(n_atoms, &n3) != 0 || !atomic_numbers ||
           !positions_bohr || !energy_h || !forces_h_bohr)
    err = "invalid arguments";
  else if (test_fault("calculate-au-reserve") ||
           session_reserve_step_atoms(session, (size_t)n_atoms) != 0)
    err = "out of memory";
  if (!err && has_cell) {
    int cell_finite = 1;
    for (int i = 0; i < 9; ++i) {
      cell_ang[i] = cell_bohr[i] * CPMDC_BOHR_TO_ANGSTROM;
      cell_finite = cell_finite && isfinite(cell_bohr[i]);
    }
    /* A non-finite cell is refused below with its own message. */
    if (virial_h && cell_finite && cpmdc_cell_volume3(cell_bohr, &volume) != 0)
      err = "cell volume is too large for binary64";
  }
  err = collective_call_gate(err);
  if (err) {
    session_drop_results(session);
    (void)fail_msg(err);
    return -1;
  }
  for (size_t i = 0; i < n3; ++i)
    session->step_positions_ang[i] = positions_bohr[i] * CPMDC_BOHR_TO_ANGSTROM;
  memcpy(session->step_atomic_numbers, atomic_numbers,
         (size_t)n_atoms * sizeof(int));
  CPMDCResult r = session_energy_gradient_cell(
      session, n_atoms, session->step_positions_ang,
      session->step_atomic_numbers, cell_ang, has_cell, forces_h_bohr);
  if (r.ok && test_fault("embed-energy-gradient"))
    r = fail_msg("embedded energy gradient failed");
  {
    /* A rank whose evaluation failed still joins the gate, so the ranks
     * that passed learn of it instead of waiting on the virial agreement. */
    const char *virial_err = r.ok ? NULL : r.message;
    double virial[9] = {0};
    if (r.ok && virial_h) {
      int have = has_cell && session->image.stress.valid;
      for (int i = 0; i < 9; ++i) {
        virial[i] = have ? volume * session->image.stress.values[i] : 0.0;
        if (!virial_err && !isfinite(virial[i]))
          virial_err = "virial is not finite";
      }
      if (!virial_err && test_fault("calculate-au-virial"))
        virial_err = "virial is not finite";
    }
    virial_err = collective_call_gate(virial_err);
    if (!r.ok) {
      session_drop_results(session);
      return -1;
    }
    if (virial_err) {
      session_drop_results(session);
      (void)fail_msg(virial_err);
      return -1;
    }
    for (size_t i = 0; i < n3; ++i)
      forces_h_bohr[i] = -forces_h_bohr[i];
    *energy_h = r.energy_h;
    if (virial_h)
      memcpy(virial_h, virial, sizeof(virial));
  }
  return 0;
}

CPMDCResult cpmdc_calculate_result(const void *params_capnp,
                                   size_t params_capnp_size_bytes,
                                   const void *force_input_capnp,
                                   size_t force_input_capnp_size_bytes,
                                   void *potential_result_capnp,
                                   size_t potential_result_capnp_capacity_bytes,
                                   size_t *potential_result_capnp_size_bytes) {
  const char *err = NULL;
  size_t required_size = 0;
  CPMDCSession *session = NULL;
  if (!params_capnp || params_capnp_size_bytes == 0 || !force_input_capnp ||
      force_input_capnp_size_bytes == 0 || !potential_result_capnp_size_bytes)
    err = "invalid arguments";
  if (potential_result_capnp_size_bytes)
    *potential_result_capnp_size_bytes = 0;
  if (!err) {
    required_size = cpmdc_potential_result_size_for_force_input(
        force_input_capnp, force_input_capnp_size_bytes);
    *potential_result_capnp_size_bytes = required_size;
  }
  if (!err && required_size == 0) {
    err = cpmdc_last_error();
    if (!err || err[0] == '\0')
      err = "invalid ForceInput geometry";
  }
  if (!err &&
      (!potential_result_capnp ||
       potential_result_capnp_capacity_bytes < required_size))
    err = "PotentialResult buffer too small";
  if (!err)
    session = cpmdc_session_create(params_capnp, params_capnp_size_bytes);
  if (!err && !session) {
    err = cpmdc_last_error();
    if (!err || err[0] == '\0')
      err = "embed config failed";
  }
  err = collective_call_gate(err);
  if (err) {
    cpmdc_session_destroy(session);
    return fail_msg(err);
  }
  CPMDCResult r = cpmdc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  cpmdc_session_destroy(session);
  return r;
}

CPMDCResult cpmdc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  const char *err = NULL;
  size_t required_size = 0;
  CPMDCSession *session = NULL;
  if (!config_capnp || config_capnp_size_bytes == 0 || !force_input_capnp ||
      force_input_capnp_size_bytes == 0 || !potential_result_capnp_size_bytes)
    err = "invalid arguments";
  if (potential_result_capnp_size_bytes)
    *potential_result_capnp_size_bytes = 0;
  if (!err) {
    required_size = cpmdc_potential_result_size_for_force_input(
        force_input_capnp, force_input_capnp_size_bytes);
    *potential_result_capnp_size_bytes = required_size;
  }
  if (!err && required_size == 0) {
    err = cpmdc_last_error();
    if (!err || err[0] == '\0')
      err = "invalid ForceInput geometry";
  }
  if (!err &&
      (!potential_result_capnp ||
       potential_result_capnp_capacity_bytes < required_size))
    err = "PotentialResult buffer too small";
  if (!err)
    session =
        cpmdc_session_create_from_config(config_capnp, config_capnp_size_bytes);
  if (!err && !session) {
    err = cpmdc_last_error();
    if (!err || err[0] == '\0')
      err = "embed config failed";
  }
  err = collective_call_gate(err);
  if (err) {
    cpmdc_session_destroy(session);
    return fail_msg(err);
  }
  CPMDCResult r = cpmdc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  cpmdc_session_destroy(session);
  return r;
}


int cpmdc_bind_calculator(int ranks_per_calc) {
  int calc = cpmdc_embed_bind_calculator(ranks_per_calc);
  if (calc < 0)
    cpmdc_store_error("CPMD calculator bind refused");
  else
    cpmdc_store_error("");
  return calc;
}

int cpmdc_adopt_calculator_comm(const void *comm, size_t comm_bytes,
                                int ranks_per_calc) {
  char err[512];
  int index = -1;
  err[0] = '\0';
  if (cpmdc_adopt_resolve(comm, comm_bytes, ranks_per_calc,
                          g_evaluated_sessions > 0, &index, err, sizeof(err)) !=
      0) {
    cpmdc_store_error(err[0] != '\0' ? err : "CPMD calculator bind refused");
    return -1;
  }
#if defined(CPMDC_HAS_CPMD)
  {
    MPI_Comm ccomm;
    memcpy(&ccomm, comm, sizeof(ccomm));
    if (cpmdc_embed_store_comm(MPI_Comm_c2f(ccomm)) != 0) {
      cpmdc_store_error("CPMD calculator bind refused");
      return -1;
    }
  }
#endif
  if (cpmdc_adopt_commit(comm, comm_bytes) != 0) {
    cpmdc_store_error("CPMD calculator bind refused");
    return -1;
  }
  cpmdc_store_error("");
  return index;
}

int cpmdc_adopted_comm(void *out, size_t nbytes) {
#if defined(CPMDC_HAS_CPMD)
  MPI_Fint fcomm;
  MPI_Comm ccomm;
  if (!out || nbytes != sizeof(MPI_Comm))
    return -1;
  fcomm = cpmdc_embed_adopted_fcomm();
  /* -1 means unset. 0 is MPI_COMM_WORLD's Fortran handle on Open MPI. */
  if (fcomm < 0)
    return -1;
  ccomm = MPI_Comm_f2c(fcomm);
  memcpy(out, &ccomm, sizeof(ccomm));
  return 0;
#else
  return cpmdc_adopt_copy(out, nbytes);
#endif
}

int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out) {
  if (!out)
    return -1;
  *out = active_image()->charge;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_charge_integrals(const CPMDCSession *session,
                                        CPMDCChargeIntegrals *out) {
  if (!session || !out)
    return -1;
  *out = session->image.charge;
  return out->valid ? 0 : -1;
}

int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out) {
  if (!out)
    return -1;
  *out = active_image()->multi;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_multi_state_energies(const CPMDCSession *session,
                                            CPMDCMultiStateEnergies *out) {
  if (!session || !out)
    return -1;
  *out = session->image.multi;
  return out->valid ? 0 : -1;
}

int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out) {
  if (!out)
    return -1;
  *out = active_image()->md;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_md_trajectory_row(const CPMDCSession *session,
                                         CPMDCMDTrajectoryRow *out) {
  if (!session || !out)
    return -1;
  *out = session->image.md;
  return out->valid ? 0 : -1;
}

int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out) {
  if (!out)
    return -1;
  *out = active_image()->prop;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_property_snapshot(const CPMDCSession *session,
                                         CPMDCPropertySnapshot *out) {
  if (!session || !out)
    return -1;
  *out = session->image.prop;
  return out->valid ? 0 : -1;
}

static int copy_gradient(const CPMDCEmbedImage *image, double *out,
                         size_t capacity, size_t *count) {
  size_t n;
  if (count)
    *count = 0;
  if (!image || !out || !count || !image->prop.valid)
    return -1;
  n = image->prop.hessian_count;
  if (n > sizeof(image->prop.hessian) / sizeof(image->prop.hessian[0]))
    return -1;
  *count = n;
  if (n == 0 || capacity < n)
    return -1;
  memcpy(out, image->prop.hessian, n * sizeof(*out));
  return 0;
}

int cpmdc_last_gradient(double *out, size_t capacity, size_t *count) {
  return copy_gradient(active_image(), out, capacity, count);
}

int cpmdc_session_last_gradient(const CPMDCSession *session, double *out,
                                size_t capacity, size_t *count) {
  if (!session) {
    if (count)
      *count = 0;
    return -1;
  }
  return copy_gradient(&session->image, out, capacity, count);
}

static int copy_dipole(const CPMDCEmbedImage *image, CPMDCDipole *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!image || !image->prop.valid || image->prop.dipole_count != 3)
    return -1;
  memcpy(out->values, image->prop.dipole, sizeof(out->values));
  out->valid = 1;
  return 0;
}

int cpmdc_last_dipole(CPMDCDipole *out) {
  return copy_dipole(active_image(), out);
}

int cpmdc_session_last_dipole(const CPMDCSession *session, CPMDCDipole *out) {
  return copy_dipole(session ? &session->image : NULL, out);
}

static int copy_charges(const CPMDCEmbedImage *image, double *charges,
                        size_t length) {
  size_t n;
  if (!image || !charges || !image->charges.valid)
    return -1;
  n = image->charges.count;
  if (n == 0 || n > CPMDC_EMBED_CHARGE_CAP || length < n)
    return -1;
  memcpy(charges, image->charges.values, n * sizeof(*charges));
  return 0;
}

int cpmdc_last_charges(double *charges, size_t length) {
  return copy_charges(active_image(), charges, length);
}

int cpmdc_session_last_charges(const CPMDCSession *session, double *charges,
                               size_t length) {
  return copy_charges(session ? &session->image : NULL, charges, length);
}

int cpmdc_last_stress(CPMDCStressTensor *out) {
  if (!out)
    return -1;
  *out = active_image()->stress;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_stress(const CPMDCSession *session, CPMDCStressTensor *out) {
  if (!session || !out)
    return -1;
  *out = session->image.stress;
  return out->valid ? 0 : -1;
}

int cpmdc_last_rank_spread(CPMDCRankSpread *out) {
  if (!out)
    return -1;
  *out = g_rank_spread;
  return out->valid ? 0 : -1;
}

int cpmdc_last_energy_components(CPMDCEnergyComponents *out) {
  if (!out)
    return -1;
  *out = active_image()->energy;
  return out->valid ? 0 : -1;
}

int cpmdc_session_last_energy_components(const CPMDCSession *session,
                                         CPMDCEnergyComponents *out) {
  if (!session || !out)
    return -1;
  *out = session->image.energy;
  return out->valid ? 0 : -1;
}

size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  struct capn arena;
  ForceInput_ptr force_input;
  if (cpmdc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                             &arena, &force_input) != 0) {
    cpmdc_store_error("invalid ForceInput message");
    return 0;
  }
  size_t n_atoms = 0;
  int has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      n_atoms > (size_t)INT_MAX || n_atoms > SIZE_MAX / 3u) {
    cpmdc_params_release(&arena);
    cpmdc_store_error("invalid ForceInput geometry");
    return 0;
  }
  (void)has_cell;
  size_t force_count = n_atoms * 3u;
  if (force_count > (size_t)INT_MAX) {
    cpmdc_params_release(&arena);
    cpmdc_store_error("invalid ForceInput geometry");
    return 0;
  }
  size_t result_size = cpmdc_potential_result_flat_size(force_count);
  cpmdc_params_release(&arena);
  if (result_size == 0) {
    cpmdc_store_error("invalid ForceInput geometry");
    return 0;
  }
  cpmdc_store_error("");
  return result_size;
}

const char *cpmdc_version(void) { return "cpmdc/" CPMDC_VERSION_STRING; }

int cpmdc_abi_version(void) { return CPMDC_ABI_VERSION; }

int cpmdc_available(void) {
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return 0;
  /* A ready reference evaluator is not an OpenCPMD archive. */
#ifndef CPMDC_HAS_CPMD
  return 0;
#else
  return 1;
#endif
}

void cpmdc_finalize(void) {
  g_active_session = NULL;
  g_orbital_owner = NULL;
  free(g_image.input_deck);
  g_image.input_deck = NULL;
  g_image.input_deck_len = 0;
  cpmdc_orbital_ledger_finalize();
  cpmdc_embed_finalize();
}
