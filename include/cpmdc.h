#pragma once

#include "cpmdc_features.h"

#include <stddef.h>

/**
 * @brief Numeric ABI generation of this header.
 *
 * Matches the shared-library soversion and cpmdc_abi_version(); bumps only on
 * an incompatible ABI change.
 */
#define CPMDC_ABI_VERSION 0

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc.h
 * @brief Stable C ABI for configuring and evaluating embedded OpenCPMD.
 *
 * The parameter buffer accepted by this API is an unpacked flat Cap'n Proto
 * message whose root type is `CPMDParams` from `schema/Potentials.capnp`.
 *
 * Direct-call "socket" style: create one `CPMDCSession` from `CPMDParams`,
 * then pass a serialized `ForceInput` for each geometry step and receive an
 * unpacked flat `PotentialResult` (same carrier as nwchemc / rgpot).
 */

/**
 * @brief Result returned by energy / gradient / forces entry points.
 */
typedef struct CPMDCResult {
  /** Non-zero when the calculation succeeds. */
  int ok;
  /** Total energy in Hartree (CPMD native a.u. energy). */
  double energy_h;
  /** Null-terminated status or error message. */
  char message[512];
} CPMDCResult;

/**
 * @brief In-process snapshot of OpenCPMD `ener_com` scalars (Hartree a.u.).
 *
 * Filled after a successful embed SCF (`wfopts`) or reference PEF evaluation.
 * Hosts read this via `cpmdc_last_energy_components()` without parsing CLI
 * ENERGY files or opening a network socket. Field names mirror `ener_com_t`
 * in OpenCPMD `ener.mod.F90`. Zero fields are valid (not set for that run).
 */
typedef struct CPMDCEnergyComponents {
  /** Non-zero when the snapshot was written by a successful evaluation. */
  int valid;
  double etot;
  double ekin;
  double epseu;
  double enl;
  double eht;
  double ehep;
  double ehee;
  double ehii;
  double exc;
  double vxc;
  double egc;
  double esr;
  double eeig;
  double eband;
  double entropy;
  double eself;
  double ecnstr;
  double amu;
  double ebogo;
  double eext;
  double etddft;
  double ehsic;
  double erestr;
  double eefield;
} CPMDCEnergyComponents;

/**
 * @brief OpenCPMD `chrg_t` density integrals (post-SCF module state).
 */
typedef struct CPMDCChargeIntegrals {
  int valid;
  double csumg;
  double csumr;
  double csums;
  double csumsabs;
} CPMDCChargeIntegrals;

/**
 * @brief Flattened CAS22-class multi-state catalog (`ener_c` + `ener_d`).
 *
 * Layout is backend-defined; `count` is the number of doubles copied into
 * `values` (caller provides capacity). Returns -1 when no snapshot.
 */
typedef struct CPMDCMultiStateEnergies {
  int valid;
  size_t count;
  double values[64];
} CPMDCMultiStateEnergies;

/**
 * @brief One ENERGY-file-equivalent trajectory row (Hartree).
 *
 * Layout (count >= 12 after a successful eval):
 *   [0] etot  [1] ekin  [2] epseu  [3] enl  [4] eht  [5] exc
 *   [6] ehep  [7] ehee  [8] ehii   [9] esr  [10] eself
 *   [11] EKINC (fictitious electronic KE; 0 for BO/SCF-only wfopt, filled in MD)
 */
typedef struct CPMDCMDTrajectoryRow {
  int valid;
  size_t count;
  double values[32];
} CPMDCMDTrajectoryRow;

/**
 * @brief PROP-style property snapshot after a successful evaluation.
 *
 * - dipole[3]: OpenCPMD `ddip%pdipole` (a.u.) when linked; PEF zeros
 * - polarizability[9]: filled when PROP/aoresponse available; else zeros with count 9
 * - hessian[]: nuclear gradient dE/dR packed as [natoms*3] (full Hessian needs
 *   dedicated PROP/Hessian run; gradient is always available after force eval)
 */
typedef struct CPMDCPropertySnapshot {
  int valid;
  size_t hessian_count;
  double hessian[4096];
  size_t dipole_count;
  double dipole[3];
  size_t polarizability_count;
  double polarizability[9];
} CPMDCPropertySnapshot;

/** Opaque handle for repeated evaluations with one Cap'n Proto parameter set. */
typedef struct CPMDCSession CPMDCSession;

/**
 * @brief Apply CPMD method parameters from a Cap'n Proto message.
 *
 * Callers do not need C setter functions for individual CPMD keywords. Build
 * one `CPMDParams` message with top-level fields, structured `inputSections`,
 * and literal `inputBlocks`, then pass its bytes to this function or to
 * `cpmdc_session_create()`.
 *
 * Feature discovery mirrors the schema carriers: typed fields such as
 * `params.inputSections.cpmd.maxIter`, `params.inputSections.system.cell`,
 * `params.inputSections.dft.hfxScreening`, and
 * `params.inputSections.atoms.pseudopotentials`; catalog sections such as
 * `catalog.section.VDW`; and escape hatches such as
 * `params.inputSections.raw`. The same serialized params buffer is accepted by
 * `cpmdc_calculate_result()` for one-shot calls.
 *
 * @param params_capnp Pointer to an unpacked flat `CPMDParams` message.
 * @param params_capnp_size_bytes Size of `params_capnp` in bytes.
 * @return 0 on success, -1 on parse or configuration failure.
 */
int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes);

/**
 * @brief Apply configuration from a `PotentialConfig` message.
 *
 * The `cpmd` union arm carries `CPMDParams` and wins wholesale when present.
 * With the arm unset, a set `common` overlay (`CommonMethodSpec`) lowers to
 * synthesized `CPMDParams`: functional, plane-wave cutoff, charge,
 * multiplicity, SCF controls (CONVERGENCE ORBITALS / MAXITER), and the
 * Monkhorst-Pack kMesh. Overlay fields without a CPMD lowering are rejected
 * and reported through `cpmdc_last_error()`.
 *
 * @return 0 on success, -1 on parse, lowering, or apply failure.
 */
int cpmdc_configure(const void *config_capnp, size_t config_capnp_size_bytes);

/**
 * @brief Compute energy and nuclear gradient for an atomic configuration.
 *
 * Positions are Angstrom; gradient is Hartree/Bohr (CPMD ionic forces are
 * negated into a nuclear gradient for API symmetry with nwchemc).
 */
CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  double *grad_h_bohr);

/**
 * @brief Compute total energy only (no gradient allocation).
 */
CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                         const int *atomic_numbers, const void *params_capnp,
                         size_t params_capnp_size_bytes);

/**
 * @brief Compute energy and nuclear forces (negative gradient, Hartree/Bohr).
 */
CPMDCResult cpmdc_energy_forces(int n_atoms, const double *positions_ang,
                                const int *atomic_numbers,
                                const void *params_capnp,
                                size_t params_capnp_size_bytes,
                                double *forces_h_bohr);

/**
 * @brief Create a persistent evaluation session from a Cap'n Proto message.
 *
 * The session owns a copy of the serialized message so callers may release the
 * input buffer after this call returns.
 */
CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                   size_t params_capnp_size_bytes);

/**
 * @brief Replace Cap'n Proto parameters before the session accepts topology.
 */
int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes);

/**
 * @brief Create a persistent session from a `PotentialConfig` message.
 *
 * Resolves the config exactly like `cpmdc_configure()` and installs the
 * effective `CPMDParams` on a new session.
 */
CPMDCSession *cpmdc_session_create_from_config(const void *config_capnp,
                                               size_t config_capnp_size_bytes);

/**
 * @brief Configure an existing session from a `PotentialConfig` message.
 *
 * Accepted only before the session evaluates; resolution matches
 * `cpmdc_configure()`.
 */
int cpmdc_session_configure(CPMDCSession *session, const void *config_capnp,
                            size_t config_capnp_size_bytes);

/** @brief Release a persistent evaluation session. */
void cpmdc_session_destroy(CPMDCSession *session);

/**
 * @brief Compute energy and nuclear gradient with session-owned parameters.
 *
 * Positions are Angstrom. The gradient buffer must have `n_atoms * 3` doubles
 * and is filled in Hartree/Bohr.
 */
CPMDCResult cpmdc_session_energy_gradient(CPMDCSession *session, int n_atoms,
                                          const double *positions_ang,
                                          const int *atomic_numbers,
                                          double *grad_h_bohr);

/**
 * @brief Compute total energy with session-owned parameters.
 *
 * Positions are Angstrom and the returned energy is Hartree.
 */
CPMDCResult cpmdc_session_energy(CPMDCSession *session, int n_atoms,
                                 const double *positions_ang,
                                 const int *atomic_numbers);

/**
 * @brief Compute energy and nuclear forces with session-owned parameters.
 *
 * Positions are Angstrom. The forces buffer must have `n_atoms * 3` doubles and
 * is filled in Hartree/Bohr.
 */
CPMDCResult cpmdc_session_energy_forces(CPMDCSession *session, int n_atoms,
                                        const double *positions_ang,
                                        const int *atomic_numbers,
                                        double *forces_h_bohr);

/**
 * @brief Compute energy and forces for one Cap'n Proto `ForceInput` step.
 *
 * Session keeps persistent `CPMDParams`; each call supplies geometry. Returned
 * energy/forces use CPMD native units: Hartree and Hartree/Bohr.
 */
CPMDCResult cpmdc_session_calculate_forces(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, double *forces_h_bohr,
    size_t forces_len);

/**
 * @brief Compute forces for one `ForceInput` step and write `PotentialResult`.
 *
 * Direct-call socket entry point: method state in the session, geometry in
 * `ForceInput`, output energy/forces converted to `ForceInput.energyUnit` and
 * `energyUnit / lengthUnit`.
 *
 * When `potential_result_capnp_capacity_bytes` is too small, returns `ok == 0`,
 * writes the required byte count to `potential_result_capnp_size_bytes`, and
 * does not evaluate CPMD.
 */
CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);

/**
 * @brief One-shot Cap'n Proto entry point (params + ForceInput -> PotentialResult).
 *
 * Multi-step callers should create one session and call
 * `cpmdc_session_calculate_result()` per step.
 */
CPMDCResult cpmdc_calculate_result(const void *params_capnp,
                                   size_t params_capnp_size_bytes,
                                   const void *force_input_capnp,
                                   size_t force_input_capnp_size_bytes,
                                   void *potential_result_capnp,
                                   size_t potential_result_capnp_capacity_bytes,
                                   size_t *potential_result_capnp_size_bytes);

/**
 * @brief Byte count needed for a `PotentialResult` for the given `ForceInput`.
 *
 * Parses geometry only; does not initialize or evaluate CPMD. Returns 0 when
 * the message is invalid or too large for the C ABI.
 */
size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes);

/** @brief Compiled library version string. */
const char *cpmdc_version(void);

/**
 * @brief Diagnostic message for the most recent int-returning configuration
 *        call on this thread.
 *
 * Covers `cpmdc_set_params()`, `cpmdc_configure()`, and the session setup
 * entry points. Returns an empty string when the last such call succeeded.
 */
const char *cpmdc_last_error(void);

/**
 * @brief Numeric ABI generation of the compiled library.
 *
 * Compare against the CPMDC_ABI_VERSION the consumer compiled with.
 */
int cpmdc_abi_version(void);

/** @brief 1 when the embedded OpenCPMD runtime is available. */
int cpmdc_available(void);

/** @brief Finalize an owned embedded CPMD runtime. */
void cpmdc_finalize(void);

/**
 * @brief Copy the last in-process `ener_com` energy decomposition.
 *
 * Call after a successful `cpmdc_energy*`, `cpmdc_session_energy*`, or
 * `cpmdc_*_calculate_*` evaluation in the same process. Returns 0 when
 * `out->valid` is set; -1 when no successful evaluation has run yet (stub
 * builds always return -1 and leave `out` zeroed). Values are Hartree.
 */
int cpmdc_last_energy_components(CPMDCEnergyComponents *out);

int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out);
int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out);
int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out);
int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out);

#ifdef __cplusplus
}
#endif
