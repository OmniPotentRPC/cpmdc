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
  /** `etot`, Hartree. */
  double etot;
  /** `ekin`, Hartree. */
  double ekin;
  /** `epseu`, Hartree. */
  double epseu;
  /** `enl`, Hartree. */
  double enl;
  /** `eht`, Hartree. */
  double eht;
  /** `ehep`, Hartree. */
  double ehep;
  /** `ehee`, Hartree. */
  double ehee;
  /** `ehii`, Hartree. */
  double ehii;
  /** `exc`, Hartree. */
  double exc;
  /** `vxc`, Hartree. */
  double vxc;
  /** `egc`, Hartree. */
  double egc;
  /** `esr`, Hartree. */
  double esr;
  /** `eeig`, Hartree. */
  double eeig;
  /** `eband`, Hartree. */
  double eband;
  /** `entropy`, Hartree. */
  double entropy;
  /** `eself`, Hartree. */
  double eself;
  /** `ecnstr`, Hartree. */
  double ecnstr;
  /** `amu`, Hartree. */
  double amu;
  /** `ebogo`, Hartree. */
  double ebogo;
  /** `eext`, Hartree. */
  double eext;
  /** `etddft`, Hartree. */
  double etddft;
  /** `ehsic`, Hartree. */
  double ehsic;
  /** `erestr`, Hartree. */
  double erestr;
  /** `eefield`, Hartree. */
  double eefield;
} CPMDCEnergyComponents;

/**
 * @brief OpenCPMD `chrg_t` density integrals (post-SCF module state).
 */
typedef struct CPMDCChargeIntegrals {
  /** Non-zero when the snapshot was written. */
  int valid;
  /** `csumg` from OpenCPMD `chrg`, electrons. */
  double csumg;
  /** `csumr` from OpenCPMD `chrg`, electrons. */
  double csumr;
  /** `csums` from OpenCPMD `chrg`, electrons. */
  double csums;
  /** `csumsabs` from OpenCPMD `chrg`, electrons. */
  double csumsabs;
} CPMDCChargeIntegrals;

/**
 * @brief Flattened CAS22-class multi-state catalog (`ener_c` + `ener_d`).
 *
 * Layout is backend-defined; `count` is the number of doubles copied into
 * `values` (caller provides capacity). Returns -1 when no snapshot.
 */
typedef struct CPMDCMultiStateEnergies {
  /** Non-zero when the snapshot was written. */
  int valid;
  /** Number of doubles copied into `values`. */
  size_t count;
  /** `ener_c` and `ener_d` values, Hartree. */
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
  /** Non-zero when the snapshot was written. */
  int valid;
  /** Number of doubles copied into `values`. */
  size_t count;
  /** Trajectory terms in the order given above, Hartree. */
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
  /** Non-zero when the snapshot was written. */
  int valid;
  /** Number of doubles copied into `hessian`. */
  size_t hessian_count;
  /** Nuclear gradient dE/dR, packed as `natoms * 3`. */
  double hessian[4096];
  /** Number of doubles copied into `dipole`. */
  size_t dipole_count;
  /** OpenCPMD `ddip%pdipole` in atomic units, or zeros for a reference PEF. */
  double dipole[3];
  /** Number of doubles copied into `polarizability`. */
  size_t polarizability_count;
  /** Polarizability when PROP is available; otherwise zeros. */
  double polarizability[9];
} CPMDCPropertySnapshot;

/**
 * @brief Cartesian stress tensor after a successful PEF evaluation.
 *
 * Layout is row-major [xx, xy, xz, yx, yy, yz, zx, zy, zz] in Hartree/Bohr^3
 * (OpenCPMD `paiu/omega` after `totstr` when `cntl%tpres`). `valid` is set
 * only when that tensor was computed for a periodic cell. A positive cell
 * volume is not enough: an isolated (cluster/Hockney) box leaves `valid`
 * unset. Returns 0 when `out->valid` is set; -1 when stress was not computed.
 */
typedef struct CPMDCStressTensor {
  /** Non-zero only when stress was computed for a periodic cell. */
  int valid;
  /** Row-major `[xx, xy, xz, yx, yy, yz, zx, zy, zz]`, Hartree/Bohr^3. */
  double values[9];
} CPMDCStressTensor;

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
 *         Failure text is `cpmdc_last_error()`.
 */
int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes);

/**
 * @brief Apply configuration from a `PotentialConfig` message.
 *
 * The `cpmd` union arm carries `CPMDParams` and wins wholesale when present.
 * With the arm unset, a set `common` overlay (`CommonMethodSpec`) lowers to
 * synthesized `CPMDParams`: functional, plane-wave cutoff, charge,
 * multiplicity, MAXITER, and the Monkhorst-Pack kMesh. Setting both the
 * cpmd arm and the overlay is rejected (capnp cannot distinguish unset arm
 * fields from defaults). Overlay fields without a CPMD lowering are rejected
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
 * @brief Bind this rank to one CPMD calculator of @p ranks_per_calc ranks.
 *
 * Collective on MPI_COMM_WORLD, once, before the first energy call.
 * World size must divide into groups of ranks_per_calc. Each group is one
 * CPMD calculator: its own communicator, its own wavefunction, the same deck.
 * A value of zero or less takes the whole world as one group.
 * Returns the group index, or -1 when the split is refused or the library
 * carries no CPMD backend. Failure text is `cpmdc_last_error()`.
 * A second call returns the same index and does not split again.
 */
int cpmdc_bind_calculator(int ranks_per_calc);

/**
 * @brief Install a communicator the caller already split.
 *
 * `comm` points at an `MPI_Comm` of `comm_bytes` bytes. The call does not
 * call `MPI_Init` or `MPI_Comm_split`. An OpenCPMD build stores that
 * communicator as `mp_comm_world` and returns this rank's calculator
 * index, `world_rank / ranks_per_calc`. A value of zero or less means one
 * calculator. Returns -1 when `comm` is null, MPI is not initialized, the
 * world does not divide, or this library has no CPMD backend. Failure
 * text is `cpmdc_last_error()`. A second call returns the same index and
 * does not replace the communicator.
 */
int cpmdc_adopt_calculator_comm(const void *comm, size_t comm_bytes,
                                int ranks_per_calc);

/**
 * @brief Copy the adopted `mp_comm_world` into `out`.
 *
 * `out` receives an `MPI_Comm` of `nbytes` bytes. Returns 0 when a
 * calculator communicator is stored, and -1 when `out` is null, `nbytes`
 * is not that size, or no communicator has been adopted. Does not write
 * `cpmdc_last_error()`.
 */
int cpmdc_adopted_comm(void *out, size_t nbytes);

/**
 * @brief Create a persistent evaluation session from a Cap'n Proto message.
 *
 * The session owns a copy of the serialized message so callers may release the
 * input buffer after this call returns. Returns NULL on failure; the reason
 * is `cpmdc_last_error()`.
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
 * @brief One-shot Cap'n Proto entry point resolving a full `PotentialConfig`
 * (native arm plus common overlay, exactly like `cpmdc_configure()`).
 *
 * Multi-step callers should create one session via
 * `cpmdc_session_create_from_config()` and call
 * `cpmdc_session_calculate_result()` per step.
 */
CPMDCResult cpmdc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);

/**
 * @brief Byte count needed for a `PotentialResult` for the given `ForceInput`.
 *
 * Parses geometry only; does not initialize or evaluate CPMD. Returns 0 when
 * the message is invalid or too large for the C ABI. That failure is
 * reported through `cpmdc_last_error()`.
 */
size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes);

/**
 * @brief Write a Cap'n Proto `Capabilities` message describing this backend.
 *
 * Loaders negotiate against the message before dispatch: backend name and
 * version, ABI generation, availability, the calculate operations the ABI
 * serves, the `CommonMethodSpec` fields the overlay lowers, and the
 * `PotentialConfig` arms accepted. A stub build reports the same operation
 * surface with `available = false`.
 *
 * Returns 0 on success. On a too-small buffer (including the pure size query
 * `capabilities_capnp == NULL`, `capabilities_capnp_capacity_bytes == 0`)
 * returns -1 with `*capabilities_capnp_size_bytes` set to the required size.
 * That -1 does not write `cpmdc_last_error()`.
 */
int cpmdc_capabilities_result(void *capabilities_capnp,
                              size_t capabilities_capnp_capacity_bytes,
                              size_t *capabilities_capnp_size_bytes);

/** @brief Compiled library version string. */
const char *cpmdc_version(void);

/**
 * @brief Diagnostic for the most recent public call on this thread that
 *        reports failure through this string.
 *
 * Written by `cpmdc_set_params()`, `cpmdc_configure()`,
 * `cpmdc_bind_calculator()`, `cpmdc_adopt_calculator_comm()`,
 * `cpmdc_session_create()`,
 * `cpmdc_session_set_params()`, `cpmdc_session_create_from_config()`,
 * `cpmdc_session_configure()`, `cpmdc_potential_result_size_for_force_input()`,
 * and every evaluation entry point (`cpmdc_energy()`,
 * `cpmdc_energy_gradient()`, `cpmdc_energy_forces()`, the session variants,
 * and the calculate-result calls). When the call returns `CPMDCResult`, the
 * text matches `message`, including a missing pseudopotential directory.
 * Empty after a successful call.
 *
 * Snapshot readers (`cpmdc_last_stress()` and the other `cpmdc_last_*` /
 * `cpmdc_session_last_*` getters) do not write it: -1 means the snapshot is
 * absent or the output pointer is null. `cpmdc_adopted_comm()` and
 * `cpmdc_capabilities_result()` do not write it either. A capabilities
 * return of -1 is the size query.
 */
const char *cpmdc_last_error(void);

/**
 * @brief Numeric ABI generation of the compiled library.
 *
 * Compare against the CPMDC_ABI_VERSION the consumer compiled with.
 */
int cpmdc_abi_version(void);

/** @brief 1 when OpenCPMD is linked, ready, and not finalized; 0 on the reference build; 0 after finalize and from the link stub. */
int cpmdc_available(void);

/** @brief Finalize an owned embedded CPMD runtime. */
void cpmdc_finalize(void);

/**
 * @brief Copy the last in-process `ener_com` energy decomposition.
 *
 * The no-session entry points read the active calculator. A session entry
 * point reads that session, including after another session has run.
 * Returns 0 when `out->valid` is set. Values are Hartree.
 */
int cpmdc_last_energy_components(CPMDCEnergyComponents *out);

/** @brief Copy the last energy decomposition stored on @p session. */
int cpmdc_session_last_energy_components(const CPMDCSession *session,
                                         CPMDCEnergyComponents *out);

/**
 * @brief Copy the last in-process density integrals.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out);

/** @brief Copy the last density integrals stored on @p session. */
int cpmdc_session_last_charge_integrals(const CPMDCSession *session,
                                        CPMDCChargeIntegrals *out);

/**
 * @brief Copy the last multi-state energy catalog.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out);

/** @brief Copy the last multi-state energy catalog stored on @p session. */
int cpmdc_session_last_multi_state_energies(const CPMDCSession *session,
                                            CPMDCMultiStateEnergies *out);

/**
 * @brief Copy the last trajectory row.
 *
 * Values are Hartree. Returns 0 when `out->valid` is set. Does not write
 * `cpmdc_last_error()`.
 */
int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out);

/** @brief Copy the last trajectory row stored on @p session. */
int cpmdc_session_last_md_trajectory_row(const CPMDCSession *session,
                                         CPMDCMDTrajectoryRow *out);

/**
 * @brief Copy the last property snapshot.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out);

/** @brief Copy the last property snapshot stored on @p session. */
int cpmdc_session_last_property_snapshot(const CPMDCSession *session,
                                         CPMDCPropertySnapshot *out);

/**
 * @brief Copy the last stress tensor.
 *
 * Returns 0 when `out->valid` is set, which is only after stress was computed
 * for a periodic cell. Does not write `cpmdc_last_error()`.
 */
int cpmdc_last_stress(CPMDCStressTensor *out);

/** @brief Copy the last stress tensor stored on @p session. */
int cpmdc_session_last_stress(const CPMDCSession *session,
                              CPMDCStressTensor *out);

#ifdef __cplusplus
}
#endif
