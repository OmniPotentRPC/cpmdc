/* C ABI frontend: Cap'n Proto sessions + unit carriers (nwchemc pattern).
 * Engine work is in Fortran bind(C) embed surface (cpmd_embed_c_api.F90). */
#include "cpmdc.h"
#include "cpmdc_params.h"

#ifndef CPMDC_VERSION_STRING
#define CPMDC_VERSION_STRING "unknown"
#endif

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
/* Fortran capnp-fortran decode of CPMDParams into embed knobs. */
int cpmdc_embed_apply_params(const void *params_capnp,
                             size_t params_capnp_size_bytes,
                             const char *input_deck, int input_deck_len,
                             const char *functional_ov, int functional_ov_len,
                             double cutoff_ov, int has_charge_ov, int charge_ov,
                             int has_mult_ov, int mult_ov);
int cpmdc_embed_get_config(char *functional, int functional_len,
                           double *cutoff_ry, int *charge, int *mult,
                           char *input_deck, int input_deck_len, char *cpmd_root,
                           int cpmd_root_len);
int cpmdc_embed_energy_grad(int n_atoms, const double *positions_ang,
                            const int *atomic_numbers, const double *cell_ang,
                            int has_cell, double *energy_h,
                            double *grad_h_bohr);
int cpmdc_embed_last_energy_components(
    int *valid, double *etot, double *ekin, double *epseu, double *enl,
    double *eht, double *ehep, double *ehee, double *ehii, double *exc,
    double *vxc, double *egc, double *esr, double *eeig, double *eband,
    double *entropy, double *eself, double *ecnstr, double *amu, double *ebogo,
    double *eext, double *etddft, double *ehsic, double *erestr,
    double *eefield);
void cpmdc_embed_finalize(void);
int cpmdc_embed_last_charge_integrals(int *valid, double *csumg, double *csumr,
                                      double *csums, double *csumsabs);
int cpmdc_embed_last_multi_state(int *valid, int *count, double *values,
                                 int capacity);
int cpmdc_embed_last_md_row(int *valid, int *count, double *values, int capacity);
int cpmdc_embed_last_properties(int *valid, int *hess_count, double *hess,
                                int hess_cap, int *dip_count, double *dip,
                                int *pol_count, double *pol);
int cpmdc_embed_last_stress(int *valid, double *stress);


/* Last Cap'n Proto params bytes for geometry-aware deck render on eval. */
static unsigned char *g_params_bytes = NULL;
static CPMDCScalarOverrides g_overrides;
static char g_overrides_functional[64];
static int g_has_overrides = 0;

static int cpmdc_set_params_ov(const void *params_capnp,
                               size_t params_capnp_size_bytes,
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
  /** Rendered CPMD input deck for the session parameters. */
  char input_deck[CPMDC_BLOCKS];
  /** Common-overlay scalar overrides merged over the params (may be unset). */
  CPMDCScalarOverrides overrides;
  char overrides_functional[64];
  int has_overrides;
  /** OpenCPMD source or build tree hint. */
  char cpmd_root[1024];
  /** Non-zero once the session has accepted one topology. */
  int topology_fixed;
  /** Atom count accepted by the first successful topology check. */
  size_t fixed_n_atoms;
  /** Ordered atomic numbers accepted by the first successful topology check. */
  int *fixed_atomic_numbers;
  /** Scratch positions for the current step in Angstrom. */
  double *step_positions_ang;
  /** Scratch atomic numbers for the current step. */
  int *step_atomic_numbers;
  /** Allocated atom capacity for step scratch buffers. */
  size_t step_atom_capacity;
  /** Non-zero when the embed layer has accepted the effective config. */
  int embed_configured;
};

static CPMDCSession *g_active_session = NULL;

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
                               const CPMDCScalarOverrides *overrides,
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
  if (cpmdc_params_effective_config_ov(root, overrides, functional,
                                       functional_size, cutoff_ry, charge,
                                       multiplicity) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  snprintf(cpmd_root, cpmd_root_size, "%s",
           cpmdc_params_text_or(view.cpmdRoot, ""));
  if (cpmdc_params_render_input_deck_ov(root, overrides, input_deck,
                                        input_deck_size) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  if (cpmdc_params_reject_unsupported_inputs(functional, input_deck) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  cpmdc_params_release(&arena);
  return 0;
}

static int embed_apply_from_wire(const void *params_capnp, size_t params_size,
                                 const char *input_deck,
                                 const CPMDCScalarOverrides *overrides) {
  const char *fov = "";
  int fov_len = 0;
  double cutoff_ov = 0.0;
  int has_charge_ov = 0, charge_ov = 0;
  int has_mult_ov = 0, mult_ov = 1;
  if (overrides) {
    if (overrides->functional && overrides->functional[0] != '\0') {
      fov = overrides->functional;
      fov_len = (int)strlen(overrides->functional);
    }
    if (overrides->cutoff_ry > 0.0)
      cutoff_ov = overrides->cutoff_ry;
    if (overrides->has_charge) {
      has_charge_ov = 1;
      charge_ov = overrides->charge;
    }
    if (overrides->has_multiplicity) {
      has_mult_ov = 1;
      mult_ov = overrides->multiplicity > 0 ? overrides->multiplicity : 1;
    }
  }
  return cpmdc_embed_apply_params(
      params_capnp, params_size, input_deck, (int)strlen(input_deck), fov,
      fov_len, cutoff_ov, has_charge_ov, charge_ov, has_mult_ov, mult_ov);
}

static int configure_embed_from_session(CPMDCSession *session) {
  if (!session)
    return -1;
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return -1;
  /* Serialized params: Fortran capnp-fortran decode; C still rendered deck. */
  if (embed_apply_from_wire(
          session->params_bytes, session->params_size, session->input_deck,
          session->has_overrides ? &session->overrides : NULL) != 0)
    return -1;
  session->embed_configured = 1;
  g_active_session = session;
  return 0;
}

/* Diagnostic channel for the int-returning configuration entry points. */
static _Thread_local char g_last_error[512];

static void cpmdc_store_error(const char *msg) {
  snprintf(g_last_error, sizeof(g_last_error), "%s", msg ? msg : "");
}

const char *cpmdc_last_error(void) { return g_last_error; }

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

static int common_list32_len(capn_list32 list) {
  capn_resolve(&list.p);
  if (list.p.type == CAPN_NULL)
    return 0;
  if (list.p.type != CAPN_LIST || list.p.datasz != 4)
    return -1;
  return list.p.len;
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

/* Lower the normalized overlay into a synthesized CPMDParams message.
 * Slice: functional, plane-wave cutoff, charge, multiplicity; every other
 * set field is rejected loudly. Returns malloc'd flat bytes. */
static int cpmd_common_to_params(CommonMethodSpec_ptr common_root,
                                 unsigned char **out, size_t *out_size) {
  *out = NULL;
  *out_size = 0;
  struct CommonMethodSpec c;
  memset(&c, 0, sizeof(c));
  read_CommonMethodSpec(&c, common_root);

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
  if (c.relativityMethod.len > 0) {
    cpmdc_store_error(
        "common overlay: relativityMethod has no CPMD lowering yet");
    return -1;
  }
  double smear_temperature_k = 0.0;
  capn_resolve(&c.smearing.p);
  if (c.smearing.p.type == CAPN_STRUCT) {
    struct CommonMethodSpec_Smearing smear;
    read_CommonMethodSpec_Smearing(&smear, c.smearing);
    if (smear.kind == CommonMethodSpec_Smearing_Kind_fermi) {
      /* CPMD free-energy functional with electron temperature; kT = width. */
      if (smear.widthEv <= 0.0) {
        cpmdc_store_error("common overlay: fermi smearing needs widthEv > 0");
        return -1;
      }
      smear_temperature_k = smear.widthEv / 8.617333262e-5;
    } else if (smear.kind != CommonMethodSpec_Smearing_Kind_none) {
      cpmdc_store_error(
          "common overlay: only fermi smearing has a CPMD lowering "
          "(free-energy functional)");
      return -1;
    }
  }
  if (c.basisSet.len > 0) {
    const char *pw = "planewave";
    int is_pw = c.basisSet.len == (int)strlen(pw);
    for (int i = 0; is_pw && i < c.basisSet.len; ++i) {
      char ch = c.basisSet.str[i];
      if ((char)((ch >= 'A' && ch <= 'Z') ? ch + ('a' - 'A') : ch) != pw[i])
        is_pw = 0;
    }
    if (!is_pw) {
      cpmdc_store_error(
          "common overlay: CPMD is plane-wave only; basisSet must be "
          "\"planewave\" or unset");
      return -1;
    }
  }

  char functional[64] = "";
  if (common_text_list_len(c.xcFunctionals) > 0) {
    if (cpmd_common_xc_token(c.xcFunctionals, functional,
                             sizeof(functional)) != 0) {
      cpmdc_store_error("common overlay: unmapped xcFunctionals for CPMD");
      return -1;
    }
  }

  int kmesh_len = common_list32_len(c.kMesh);
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
      /* CPMD's primary SCF criterion is CONVERGENCE ORBITALS (a.u.). */
      cpmd_sec.convergenceOrbitals =
          c.scfEnergyToleranceEv > 0.0
              ? c.scfEnergyToleranceEv / CPMDC_HARTREE_EV
              : 1.0e-6;
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
  if (struct_ptr.type == CAPN_NULL)
    return -1;
  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL || capn_setp(root, 0, struct_ptr) != 0) {
    capn_free(&arena);
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
    return -1;
  }
  *out = buffer;
  *out_size = (size_t)written;
  return 0;
}

/* Merge a set common overlay over a native cpmd arm: the arm's bytes are the
 * params; overlay scalars fill fields the arm carries at schema defaults
 * (functional BLYP, cutoff 70 Ry, charge 0, multiplicity 1). Overlay fields
 * that would require synthesizing sections next to an arm are not yet
 * mergeable. */
static int cpmd_merge_common_into_arm(CPMDParams_ptr arm,
                                      CommonMethodSpec_ptr common_root,
                                      unsigned char **out, size_t *out_size,
                                      CPMDCScalarOverrides *ov,
                                      char *ov_functional, size_t ov_len,
                                      int *has_ov) {
  struct CommonMethodSpec c;
  memset(&c, 0, sizeof(c));
  read_CommonMethodSpec(&c, common_root);

  if (common_list32_len(c.kMesh) != 0 || c.scfMaxIterations > 0 ||
      c.scfEnergyToleranceEv > 0.0) {
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

  char base_functional[64];
  double base_cutoff = 70.0;
  int base_charge = 0;
  int base_mult = 1;
  if (cpmdc_params_effective_config(arm, base_functional,
                                    sizeof(base_functional), &base_cutoff,
                                    &base_charge, &base_mult) != 0) {
    cpmdc_store_error("PotentialConfig cpmd arm is not readable");
    return -1;
  }

  memset(ov, 0, sizeof(*ov));
  *has_ov = 0;
  if (common_text_list_len(c.xcFunctionals) > 0 &&
      strcmp(base_functional, "BLYP") == 0) {
    if (cpmd_common_xc_token(c.xcFunctionals, ov_functional, ov_len) != 0) {
      cpmdc_store_error("common overlay: unmapped xcFunctionals for CPMD");
      return -1;
    }
    ov->functional = ov_functional;
    *has_ov = 1;
  }
  if (c.planewaveCutoffEv > 0.0 && base_cutoff == 70.0) {
    ov->cutoff_ry = c.planewaveCutoffEv / CPMDC_RYDBERG_EV;
    *has_ov = 1;
  }
  if (c.charge != 0 && base_charge == 0) {
    ov->has_charge = 1;
    ov->charge = c.charge;
    *has_ov = 1;
  }
  if (c.spinMultiplicity > 1 && base_mult == 1) {
    ov->has_multiplicity = 1;
    ov->multiplicity = c.spinMultiplicity;
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
  CPMDCSession *session = cpmdc_session_create(params, params_size);
  if (session && has_ov) {
    /* Re-install with the overlay overrides so scalars and decks merge. */
    if (cpmdc_session_set_params_ov(session, params, params_size, &ov) != 0) {
      cpmdc_session_destroy(session);
      session = NULL;
    }
  }
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
  char input_deck[CPMDC_BLOCKS];
  char cpmd_root[1024];
  if (apply_params_buffer(params_capnp, params_capnp_size_bytes, overrides,
                          functional, sizeof(functional), &cutoff_ry, &charge,
                          &multiplicity, input_deck, sizeof(input_deck),
                          cpmd_root, sizeof(cpmd_root)) != 0)
    return -1;
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
  if (!ensure_embed_init() || !cpmdc_embed_available())
    return -1;
  free(g_params_bytes);
  g_params_bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!g_params_bytes)
    return -1;
  memcpy(g_params_bytes, params_capnp, params_capnp_size_bytes);
  g_params_size = params_capnp_size_bytes;
  /* Fortran decode of functional/cutoff/charge/mult/root + store C deck. */
  if (embed_apply_from_wire(params_capnp, params_capnp_size_bytes, input_deck,
                            overrides) != 0)
    return -1;
  (void)functional;
  (void)cutoff_ry;
  (void)charge;
  (void)multiplicity;
  (void)cpmd_root;
  g_active_session = NULL;
  return 0;
}

int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes) {
  return cpmdc_set_params_ov(params_capnp, params_capnp_size_bytes, NULL);
}

/* Merge ForceInput geometry into a Cap'n Proto-derived CPMD INPUT deck. */
static int push_geometry_deck_from_params(const void *params_bytes,
                                          size_t params_size,
                                          const CPMDCScalarOverrides *overrides,
                                          int n_atoms,
                                          const double *positions_ang,
                                          const int *atomic_numbers,
                                          const double *cell_ang,
                                          int has_cell) {
  if (!params_bytes || params_size == 0)
    return -1;
  struct capn arena;
  CPMDParams_ptr root;
  if (cpmdc_params_root(params_bytes, params_size, &arena, &root) != 0)
    return -1;
  char deck[CPMDC_BLOCKS];
  if (cpmdc_params_render_deck_with_geometry_ov(root, overrides, n_atoms,
                                                positions_ang, atomic_numbers,
                                                cell_ang, has_cell, deck,
                                                sizeof(deck)) != 0) {
    cpmdc_params_release(&arena);
    return -1;
  }
  cpmdc_params_release(&arena);
  return cpmdc_embed_set_deck(deck, (int)strlen(deck)) != 0 ? 0 : -1;
}

static CPMDCResult
energy_gradient_cell_with_params(const void *params_bytes, size_t params_size,
                                 const CPMDCScalarOverrides *overrides,
                                 int n_atoms, const double *positions_ang,
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
  /* Geometry and forces travel only as C arrays into the embed (nwchemc
   * pattern). Method knobs live in embed module state from set_config /
   * session configure — do not materialize an INPUT deck per force. */
  (void)params_bytes;
  (void)params_size;
  (void)overrides;
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

static CPMDCResult energy_gradient_cell(int n_atoms, const double *positions_ang,
                                        const int *atomic_numbers,
                                        const double *cell_ang, int has_cell,
                                        double *grad_h_bohr) {
  return energy_gradient_cell_with_params(g_params_bytes, g_params_size,
                                          g_has_overrides ? &g_overrides : NULL,
                                          n_atoms,
                                          positions_ang, atomic_numbers,
                                          cell_ang, has_cell, grad_h_bohr);
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
                          session->has_overrides ? &session->overrides : NULL,
                          session->functional, sizeof(session->functional),
                          &session->cutoff_ry, &session->charge,
                          &session->multiplicity, session->input_deck,
                          sizeof(session->input_deck), session->cpmd_root,
                          sizeof(session->cpmd_root)) != 0) {
    free(session->params_bytes);
    free(session);
    return NULL;
  }
  (void)configure_embed_from_session(session);
  return session;
}

static int cpmdc_session_set_params_ov(CPMDCSession *session,
                                       const void *params_capnp,
                                       size_t params_capnp_size_bytes,
                                       const CPMDCScalarOverrides *overrides) {
  if (!session || !params_capnp || params_capnp_size_bytes == 0)
    return -1;
  if (session->topology_fixed)
    return -1;
  unsigned char *bytes = (unsigned char *)malloc(params_capnp_size_bytes);
  if (!bytes)
    return -1;
  if (overrides) {
    session->overrides = *overrides;
    if (overrides->functional) {
      snprintf(session->overrides_functional,
               sizeof(session->overrides_functional), "%s",
               overrides->functional);
      session->overrides.functional = session->overrides_functional;
    }
    session->has_overrides = 1;
  } else {
    session->has_overrides = 0;
  }
  if (apply_params_buffer(params_capnp, params_capnp_size_bytes,
                          session->has_overrides ? &session->overrides : NULL,
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

int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes) {
  return cpmdc_session_set_params_ov(session, params_capnp,
                                     params_capnp_size_bytes, NULL);
}

void cpmdc_session_destroy(CPMDCSession *session) {
  if (!session)
    return;
  if (g_active_session == session)
    g_active_session = NULL;
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
  if ((!session->embed_configured || g_active_session != session) &&
      configure_embed_from_session(session) != 0)
    return fail_msg("CPMD embed not available");
  return energy_gradient_cell_with_params(session->params_bytes,
                                          session->params_size,
                                          session->has_overrides
                                              ? &session->overrides
                                              : NULL,
                                          n_atoms,
                                          positions_ang, atomic_numbers,
                                          cell_ang, has_cell, grad_h_bohr);
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
  double stress_factor = 1.0;
  if (cpmdc_force_input_result_factors(force_input, &energy_factor,
                                       &force_factor) != 0) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput result units");
    return r;
  }
  if (cpmdc_force_input_stress_result_factors(force_input, &energy_factor,
                                              &stress_factor) != 0) {
    cpmdc_params_release(&arena);
    snprintf(r.message, sizeof(r.message), "invalid ForceInput stress units");
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
                                     force_count, stress_factor,
                                     potential_result_capnp,
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

CPMDCResult cpmdc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  r.message[0] = '\0';
  if (!config_capnp || config_capnp_size_bytes == 0 || !force_input_capnp ||
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
      cpmdc_session_create_from_config(config_capnp, config_capnp_size_bytes);
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


int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0;
  int rc = cpmdc_embed_last_charge_integrals(&valid, &out->csumg, &out->csumr,
                                             &out->csums, &out->csumsabs);
  out->valid = valid;
  return rc;
}

int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0, count = 0;
  int rc = cpmdc_embed_last_multi_state(&valid, &count, out->values,
                                        (int)(sizeof(out->values) / sizeof(out->values[0])));
  out->valid = valid;
  out->count = count > 0 ? (size_t)count : 0;
  return rc;
}

int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0, count = 0;
  int rc = cpmdc_embed_last_md_row(&valid, &count, out->values,
                                   (int)(sizeof(out->values) / sizeof(out->values[0])));
  out->valid = valid;
  out->count = count > 0 ? (size_t)count : 0;
  return rc;
}

int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0, hc = 0, dc = 0, pc = 0;
  int rc = cpmdc_embed_last_properties(
      &valid, &hc, out->hessian,
      (int)(sizeof(out->hessian) / sizeof(out->hessian[0])), &dc, out->dipole,
      &pc, out->polarizability);
  out->valid = valid;
  out->hessian_count = hc > 0 ? (size_t)hc : 0;
  out->dipole_count = dc > 0 ? (size_t)dc : 0;
  out->polarizability_count = pc > 0 ? (size_t)pc : 0;
  return rc;
}

int cpmdc_last_stress(CPMDCStressTensor *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0;
  int rc = cpmdc_embed_last_stress(&valid, out->values);
  out->valid = valid;
  return rc;
}

int cpmdc_last_energy_components(CPMDCEnergyComponents *out) {
  if (!out)
    return -1;
  memset(out, 0, sizeof(*out));
  if (!ensure_embed_init())
    return -1;
  int valid = 0;
  int rc = cpmdc_embed_last_energy_components(
      &valid, &out->etot, &out->ekin, &out->epseu, &out->enl, &out->eht,
      &out->ehep, &out->ehee, &out->ehii, &out->exc, &out->vxc, &out->egc,
      &out->esr, &out->eeig, &out->eband, &out->entropy, &out->eself,
      &out->ecnstr, &out->amu, &out->ebogo, &out->eext, &out->etddft,
      &out->ehsic, &out->erestr, &out->eefield);
  out->valid = valid;
  return rc;
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

const char *cpmdc_version(void) { return "cpmdc/" CPMDC_VERSION_STRING; }

int cpmdc_abi_version(void) { return CPMDC_ABI_VERSION; }

int cpmdc_available(void) {
  if (!ensure_embed_init())
    return 0;
  return cpmdc_embed_available() != 0;
}

void cpmdc_finalize(void) {
  g_active_session = NULL;
  cpmdc_embed_finalize();
}
