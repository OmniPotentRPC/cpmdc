#pragma once

#include "cpmdc_export.h"
#include "cpmdc_features.h"

#include <stddef.h>

/**
 * @brief Numeric ABI generation of this header.
 *
 * Matches the shared-library soversion and cpmdc_abi_version(); bumps only on
 * an incompatible ABI change.
 */
/* [[[cog
major = open("version.txt").read().strip().split(".")[0]
cog.outl(f"#define CPMDC_ABI_VERSION {major}")
]]] */
#define CPMDC_ABI_VERSION 0
/* [[[end]]] */

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
 * - dipole[3]: OpenCPMD molecular dipole (electron Bohr) for the current
 *   geometry when `dipole_count == 3`
 * - polarizability[9]: always zeros, with count 9; no build computes it
 * - hessian[]: deprecated name of the nuclear gradient dE/dR packed as
 *   [n_atoms*3]. Use cpmdc_last_gradient() or
 *   cpmdc_session_last_gradient() instead. A full Hessian needs a dedicated
 *   PROP/Hessian run.
 *
 * The fixed legacy array limits evaluations to 1365 atoms. libcpmdc rejects a
 * larger system instead of truncating its gradient.
 */
typedef struct CPMDCPropertySnapshot {
  /** Non-zero when the snapshot was written. */
  int valid;
  /** Deprecated: number of doubles copied into `hessian`. */
  size_t hessian_count;
  /** Deprecated: nuclear gradient dE/dR. Use the gradient accessors. */
  double hessian[4096];
  /** Number of doubles copied into `dipole`. */
  size_t dipole_count;
  /** Molecular dipole in electron Bohr; meaningful only when count is 3. */
  double dipole[3];
  /** Number of doubles copied into `polarizability`. */
  size_t polarizability_count;
  /** Always zeros: no build computes the polarizability. */
  double polarizability[9];
} CPMDCPropertySnapshot;

/**
 * @brief Molecular dipole after a successful evaluation.
 *
 * `valid` is set only when OpenCPMD recomputed the dipole from the converged
 * wavefunction for the current geometry. Values are electron Bohr.
 */
typedef struct CPMDCDipole {
  /** Non-zero only when this evaluation recomputed the dipole. */
  int valid;
  /** Cartesian molecular dipole in electron Bohr. */
  double values[3];
} CPMDCDipole;

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

/**
 * @brief How far the ranks of one calculator were from its rank 0.
 *
 * Written by an evaluation that ran with `CPMDC_RANK_CONSENSUS=1` under a
 * running MPI. The differences are taken before every rank receives rank 0's
 * energy and gradient. A NaN difference is reported as infinity. Both are 0
 * when every rank held the same values.
 */
typedef struct CPMDCRankSpread {
  /** Non-zero when the last evaluation ran the consensus step. */
  int valid;
  /** Ranks in the calculator communicator. */
  int ranks;
  /** Largest `|E_r - E_0|` over the ranks `r`, Hartree. */
  double energy_h;
  /** Largest nuclear-gradient component difference from rank 0, Hartree/Bohr. */
  double gradient_h_bohr;
} CPMDCRankSpread;

/**
 * @brief Effective scalar method settings decoded from `CPMDParams`.
 *
 * Structured section values take precedence over top-level values, as they
 * do when rendering the CPMD input deck.
 */
typedef struct CPMDCMethodInfo {
  /** Effective DFT functional, null-terminated. */
  char functional[64];
  /** Effective total system charge. */
  int charge;
  /** Effective spin multiplicity, `2S + 1`. */
  int multiplicity;
} CPMDCMethodInfo;

/** Opaque handle for repeated evaluations with one Cap'n Proto parameter set. */
typedef struct CPMDCSession CPMDCSession;

/**
 * @brief Decode the effective functional, charge, and multiplicity.
 *
 * This validates and reads the message without applying it or initializing
 * OpenCPMD.
 *
 * @return 0 on success, -1 on a null output or an invalid message.
 *         Failure text is `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_params_method(const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  CPMDCMethodInfo *out);

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
 * @return 0 on success, -1 on a NULL or empty buffer, a parse, render, or
 *         configuration failure, an embed that is not available, or an
 *         allocation failure. Failure text is `cpmdc_last_error()`.
 *
 * The call replaces the process-wide parameters that `cpmdc_energy()`,
 * `cpmdc_energy_gradient()`, and `cpmdc_energy_forces()` evaluate with, and
 * makes the process-wide calculator the active one, which the no-session
 * `cpmdc_last_*` getters read. No `CPMDCSession` is modified.
 */
CPMDC_API int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes);

/**
 * @brief Apply configuration from a `PotentialConfig` message.
 *
 * The `cpmd` union arm carries `CPMDParams`. With `common` unset, the arm is
 * applied as given. A set `common` overlay (`CommonMethodSpec`) is merged
 * over that arm. The arm bytes stay, and the deck renderers and the recorded
 * settings apply the same precedence to the functional, cutoff, charge, and
 * multiplicity lines cpmdc writes. A value a typed field sets beats the overlay,
 * the overlay beats a schema default, and a typed `dft` or `system` section
 * beats the top-level field. A field that is omitted or holds its schema
 * default (functional `BLYP`, cutoff 70 Ry, charge 0 with `hasCharge` false,
 * multiplicity 1) is unset, so it replaces neither the top-level arm value nor
 * a merged overlay value; a typed `system` charge of 0 with `hasCharge` false
 * is unset, and a non-zero charge or zero `system.charge` with
 * `system.hasCharge` true is explicit. `planewaveCutoffEv` is in eV and is
 * converted to Ry. `spinMultiplicity` 0 is unset, and any other value, an
 * explicit singlet included, is set. An `xcFunctionals` list is mapped before
 * that precedence, and an unmapped token is rejected for every arm functional.
 * Those opaque carriers are rendered as given and are not read back into
 * the recorded settings, so the carrier check runs for every set overlay
 * field, whether or not a typed field also sets it, including zero `charge`
 * with `hasCharge` true and an explicit singlet `spinMultiplicity` of 1: the
 * overlay field is rejected when `inputBlocks`, a raw section, or a generic
 * section opens the owning section (`&DFT` for `xcFunctionals`, `&SYSTEM` for
 * the cutoff, charge, and multiplicity), when a set directive names
 * `DFT.FUNCTIONAL`, `SYSTEM.CUTOFF`, `SYSTEM.CHARGE`, or
 * `SYSTEM.MULTIPLICITY`, or when the `directives` list of the typed `dft`
 * section holds `FUNCTIONAL` or that of the typed `system` section holds
 * `CUTOFF`, `CHARGE`, or `MULTIPLICITY`. Without an overlay, a `FUNCTIONAL`
 * directive in the typed `dft` section or a `DFT.FUNCTIONAL` set directive
 * is rejected when the deck is rendered, because the rendered `&DFT`
 * section already carries the functional. The same render rejects a
 * `CUTOFF`, `CHARGE`, or `MULTIPLICITY` directive in the typed `system`
 * section, and a `SYSTEM.CUTOFF`, `SYSTEM.CHARGE`, or `SYSTEM.MULTIPLICITY`
 * set directive. `CUTOFF` is written on every deck, including the default
 * 70 Ry. `CHARGE` is written only when it is not 0 and `MULTIPLICITY` only
 * when it is greater than 1; an omitted line is the schema default the
 * recorded configuration carries. Those directives are rejected for a
 * default value and for any other value alike. A default leaves the
 * directive as a value the recorded configuration does not carry, and any
 * other value leaves two lines. The render does not keep one and drop the
 * other. Fermi `smearing` is copied as
 * `FREE-ENERGY` and `TEMPERATURE ELECTRON` when `&CPMD` sets no free-energy
 * or electron-temperature control (including split `TEMPERATURE` and
 * `ELECTRON` arguments in `cpmd.directives` or `CPMD.TEMPERATURE` set
 * directives) and the arm has no raw or generic `&CPMD` block.
 * `vanDerWaalsS6` is copied into a `vdwParams` section whose `s6` is
 * unset. `basisSet` must be unset or `"planewave"`. Next to the arm,
 * `kMesh`, `scfMaxIterations`, `scfEnergyToleranceEv`, `vanDerWaalsMethod`,
 * and `relativityMethod` are rejected, as are already-set smearing or `s6`,
 * `vanDerWaalsS6` without `vdwParams`, a non-plane-wave `basisSet`, and a
 * negative or non-finite overlay scalar. With the arm unset, the overlay
 * lowers to synthesized `CPMDParams`: functional, basis set, plane-wave
 * cutoff, charge and `hasCharge`, multiplicity, SCF energy tolerance,
 * MAXITER, the Monkhorst-Pack kMesh, smearing, and the van der Waals method
 * and S6, the fields `cpmdc_capabilities_result()` lists. That overlay-only
 * lowering
 * also rejects `vanDerWaalsS6` without `vanDerWaalsMethod` and a non-zero
 * `smearing.widthEv` with `kind = none`. Overlay fields without a CPMD
 * lowering are rejected. Every rejection leaves a `cpmdc_last_error()` that
 * names its cause. An arm `inputBlocks`, `inputSections`, or `directives`
 * list the merge cannot read is reported with the overlay field that
 * needed it and the list, as in
 * `charge cannot merge: the cpmd arm has an invalid inputBlocks list`.
 *
 * A message with neither the arm nor an overlay applies nothing and returns
 * 0. An arm other than `none` and `cpmd` is rejected.
 *
 * @return 0 on success, -1 on a NULL or empty buffer, a parse, lowering, or
 *         apply failure, or an embed that is not available.
 */
CPMDC_API int cpmdc_configure(const void *config_capnp, size_t config_capnp_size_bytes);

/**
 * @brief Compute energy and nuclear gradient for an atomic configuration.
 *
 * Positions are Angstrom, `3 * n_atoms` doubles, atom-major. The gradient
 * buffer holds `3 * n_atoms` doubles and receives Hartree/Bohr (CPMD ionic
 * forces are negated into a nuclear gradient for API symmetry with nwchemc).
 * Atomic numbers must be in 1..118.
 *
 * Each call applies `params_capnp` like `cpmdc_set_params()` and then
 * evaluates. `cpmdc_set_params()` detaches the process image without clearing
 * the stored orbitals: a call whose functional, deck, cutoff, charge,
 * multiplicity, atom count, and elemental composition match the latched basis,
 * and whose step without a cell matches a latched state that also has no cell
 * (`CPMDC_CELL_IN_PLACE` does not apply without a cell, because
 * `embed_cell_moves_in_place` returns false when `has_cell` is 0), skips full
 * CPMD setup. It then starts the SCF from the stored orbitals when the
 * geometry is within `CPMDC_WARM_RMSD` (default 0.15 Angstrom) of a stored
 * geometry, and from the deck's guess without full setup beyond that limit. Any
 * other call runs full CPMD setup and starts from a cold SCF. On failure `ok`
 * is 0, `message` holds the reason, and `cpmdc_last_error()` repeats it.
 */
CPMDC_API CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                            const int *atomic_numbers,
                                            const void *params_capnp,
                                            size_t params_capnp_size_bytes,
                                            double *grad_h_bohr);

/**
 * @brief Compute the total energy. The gradient is evaluated into a temporary
 * array that the call frees before it returns.
 */
CPMDC_API CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                                   const int *atomic_numbers, const void *params_capnp,
                                   size_t params_capnp_size_bytes);

/**
 * @brief Compute energy and nuclear forces (negative gradient, Hartree/Bohr).
 */
CPMDC_API CPMDCResult cpmdc_energy_forces(int n_atoms, const double *positions_ang,
                                          const int *atomic_numbers,
                                          const void *params_capnp,
                                          size_t params_capnp_size_bytes,
                                          double *forces_h_bohr);

/**
 * @brief Bind this rank to one CPMD calculator of @p ranks_per_calc ranks.
 *
 * Collective on MPI_COMM_WORLD, before the first energy call. Each group is
 * one CPMD calculator: its own communicator, its own wavefunction, the same
 * deck.
 *
 * On an OpenCPMD build the call starts MPI with `MPI_Init` when MPI is not
 * running, before it looks at @p ranks_per_calc. No other build starts MPI:
 * the reference build and the link stub return -1 without calling MPI.
 *
 * How ranks are divided. A value of zero or less means the world size, which
 * makes one group of every rank. Otherwise the world size must be a multiple
 * of @p ranks_per_calc. There is no remainder group and no smaller last
 * group: a world size that is not a multiple, or a value above the world
 * size, is refused on every rank and nothing is split. The groups are
 * consecutive blocks of world ranks. Group `g` holds the world ranks
 * `g * ranks_per_calc` to `(g + 1) * ranks_per_calc - 1`, and there are
 * `world_size / ranks_per_calc` groups. A rank is in group
 * `world_rank / ranks_per_calc` (integer division) and has rank
 * `world_rank % ranks_per_calc` inside it; the group's first rank is its
 * parent rank.
 *
 * Returns the group index, or -1 when the world size is not a multiple of
 * @p ranks_per_calc, `MPI_Comm_split` returns an error, or the library
 * carries no CPMD backend (the reference build and the stub). Failure text is
 * `cpmdc_last_error()`.
 *
 * The first call that is not refused splits `MPI_COMM_WORLD` with
 * `MPI_Comm_split`. A later call does not split again, whether the
 * communicator in use came from an earlier bind or from
 * `cpmdc_adopt_calculator_comm()`. It checks its own @p ranks_per_calc as
 * above and returns `world_rank / ranks_per_calc` for that value, which is
 * the index of the group in use only when the value is the same. A later
 * call whose value does not divide the world size, or is above it, returns -1
 * on every rank like the first call, and the split in use stays as it was.
 */
CPMDC_API int cpmdc_bind_calculator(int ranks_per_calc);

/**
 * @brief Install a communicator the caller already split.
 *
 * `comm` points at an `MPI_Comm` of `comm_bytes` bytes. Collective on that
 * communicator. Does not call `MPI_Init` or `MPI_Comm_split`. Rank and size
 * come from `comm`. When `ranks_per_calc` is greater than zero, that size
 * must equal it. A mismatch returns -1 and `cpmdc_last_error()` names both
 * numbers. Zero or less skips the size check.
 * The group must be one aligned block of world ranks, starting at a multiple
 * of the adopted size. Any other layout returns -1 on every rank of the
 * group. The calculator index is the lowest world rank in that block divided
 * by the size. Indices are unique when every group has the same size, which
 * the same `ranks_per_calc > 0` on every rank enforces.
 * The handle is borrowed: the call does not use `MPI_Comm_dup`. It stores
 * the caller's communicator, and an OpenCPMD build also stores it as
 * `mp_comm_world`. The caller keeps ownership and must not free the
 * communicator until the last session is destroyed or a later adopt replaces
 * it. An adopt after `cpmdc_bind_calculator` replaces the communicator that
 * call split, and the caller of bind cannot free that one.
 * A later call replaces the communicator when no live `CPMDCSession` has
 * started an evaluation. A session that was created but has not evaluated
 * does not block the call. If any rank of the group has a live session that
 * has evaluated, every rank returns -1 and the stored communicator stays as
 * it was.
 * Returns the calculator index, or -1 when `comm` is null, `comm_bytes` is
 * not the size of an `MPI_Comm`, `comm` is `MPI_COMM_NULL`, MPI is not
 * initialized, or this build has no MPI. Failure text is
 * `cpmdc_last_error()`.
 * Call adopt, `cpmdc_session_create`, and `cpmdc_session_destroy` from one
 * thread: `g_evaluated_sessions` and `g_active_session` are process globals, and
 * OpenCPMD state is process-global too.
 */
CPMDC_API int cpmdc_adopt_calculator_comm(const void *comm, size_t comm_bytes,
                                          int ranks_per_calc);

/**
 * @brief Copy the stored calculator communicator into `out`.
 *
 * `out` receives an `MPI_Comm` of `nbytes` bytes. On an OpenCPMD build that
 * handle is `mp_comm_world`. Returns 0 when a calculator communicator is
 * stored, and -1 when `out` is null, `nbytes` is not that size, or none is
 * stored. `cpmdc_adopt_calculator_comm()` stores one on every build with
 * MPI; `cpmdc_bind_calculator()` stores its split on an OpenCPMD build only.
 * Does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_adopted_comm(void *out, size_t nbytes);

/**
 * @brief Create a persistent evaluation session from a Cap'n Proto message.
 *
 * The session owns a copy of the serialized message so callers may release the
 * input buffer after this call returns. Returns NULL on failure; the reason
 * is `cpmdc_last_error()`.
 */
CPMDC_API CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                             size_t params_capnp_size_bytes);

/**
 * @brief Replace Cap'n Proto parameters before the session accepts topology.
 *
 * The first evaluation call that passes its argument checks fixes the
 * topology, whether or not the evaluation then succeeds. The session keeps a
 * copy of the new bytes, drops any `common` overlay that
 * `cpmdc_session_create_from_config()` or `cpmdc_session_configure()`
 * installed, and keeps a constant cutoff set by
 * `cpmdc_session_set_constant_cutoff()`. The call is all or nothing: when it
 * returns -1 the session keeps its previous parameters, rendered deck,
 * effective scalars, overlay, and embed state, and only `cpmdc_last_error()`
 * changes.
 *
 * @return 0, or -1 when `session` is `NULL`, the buffer is NULL or empty,
 * the topology is fixed, the message does not parse or render, the embed
 * is not available, or an allocation fails. Every failure returns -1.
 * Failure text is `cpmdc_last_error()`; it is `out of memory` when the copy
 * of the buffer cannot be allocated.
 */
CPMDC_API int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                                       size_t params_capnp_size_bytes);

/**
 * @brief Evaluate the session with a constant-cutoff kinetic functional.
 *
 * The session deck gains `&SYSTEM CONSTANT CUTOFF`, which replaces G^2 in the
 * kinetic energy by G^2 + A [1 + erf((G^2 - E0) / sigma)] (Bernasconi et al.,
 * J. Phys. Chem. Solids 56, 501 (1995)). Energy and stress then belong to the
 * constant cutoff E0 while the plane-wave basis keeps the deck's cutoff, so
 * the stress carries no Pulay stress and follows the energy across a change
 * of cell. All values are Rydberg. `cutoff_ry` is E0; zero or less turns the
 * functional off. `height_ry` is A; zero or less selects 50 times the deck's
 * cutoff. `width_ry` is sigma; zero or less selects E0 / 10. E0 + 2 sigma must
 * not exceed the deck's cutoff.
 *
 * @return 0, or -1 when `session` is `NULL`, a value is not finite, the
 * session has accepted a topology, E0 + 2 sigma exceeds the deck's cutoff,
 * the message also sets `system.constantCutoff`, or the message gives
 * `&SYSTEM` as raw text. The session keeps its previous setting on failure,
 * and the reason is `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_session_set_constant_cutoff(CPMDCSession *session, double cutoff_ry,
                                                double height_ry, double width_ry);

/**
 * @brief Create a persistent session from a `PotentialConfig` message.
 *
 * Resolves the config exactly like `cpmdc_configure()` and installs the
 * effective `CPMDParams` on a new session.
 *
 * @return The session, or `NULL` on a failure `cpmdc_configure()` reports,
 * and when the message carries nothing CPMD can apply (neither the `cpmd`
 * arm nor an overlay), unlike `cpmdc_configure()`, which returns 0 then.
 * The reason is `cpmdc_last_error()`.
 */
CPMDC_API CPMDCSession *cpmdc_session_create_from_config(const void *config_capnp,
                                                         size_t config_capnp_size_bytes);

/**
 * @brief Configure an existing session from a `PotentialConfig` message.
 *
 * Accepted only before the session accepts a topology, as for
 * `cpmdc_session_set_params()`; resolution matches `cpmdc_configure()`. A
 * message with nothing CPMD can apply changes nothing and returns 0. The call
 * is all or nothing: when it returns -1 the session keeps its previous
 * parameters, overlay, and effective settings, and only `cpmdc_last_error()`
 * changes.
 *
 * @return 0, or -1 when `session` is `NULL`, the config fails as in
 * `cpmdc_configure()`, or it has something to apply and the topology is
 * fixed. The reason is `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_session_configure(CPMDCSession *session, const void *config_capnp,
                                      size_t config_capnp_size_bytes);

/** @brief Release a persistent evaluation session. */
CPMDC_API void cpmdc_session_destroy(CPMDCSession *session);

/**
 * @brief Name the calculation the next evaluations of this session belong to.
 *
 * The process stores the converged orbitals of one SCF. With a key named,
 * those orbitals are kept per key: an evaluation under a key starts from
 * the orbitals the last converged evaluation under the same key left, and
 * a key that has none yet starts from the orbitals already stored. Of the
 * stored orbitals of every key, the ones that converged nearest to the
 * new geometry start the SCF, and only within CPMDC_WARM_RMSD (root-mean-
 * square distance per atom, default 0.15 Angstrom); farther, the deck's
 * guess starts it. A host
 * that runs several calculations on one calculator, the images of a band
 * or the beads of a ring polymer, names each by a stable key, so each
 * starts from its own previous geometry instead of from its neighbour's.
 * Every rank of a calculator names the same key before the same
 * evaluation. The call records the key. The evaluation applies it, so a
 * key that is named and then replaced before an evaluation does not park
 * orbitals. A key also keeps its own `EXTRAPOLATE WFN` ring, so switching
 * keys preserves each trajectory's predictor history.
 * CPMDC_ORBITAL_SLOT_LIMIT, a positive integer, is the most parked orbital
 * copies plus populated extrapolation-ring entries, including the selected
 * key's entries. An evaluation that would add another copy past that limit
 * fails, cpmdc_last_error() names the variable, and every stored copy and the
 * selected key are left unchanged. Unset storage is unbounded. A basis change
 * clears every key's orbitals and history. A session that never names a key
 * keeps the single process history exactly as the unkeyed route does.
 *
 * @return 0, or -1 when `session` is `NULL`.
 */
CPMDC_API int cpmdc_session_select_orbitals(CPMDCSession *session, long long key);

/**
 * @brief Copy stored orbitals from one key to another in this session.
 *
 * The call duplicates the converged orbitals `from_key` held at copy time,
 * and the geometry they converged at, into `to_key`, replacing whatever
 * `to_key` already held. The selected key is unchanged. If `to_key` is the
 * selected key, its live store is replaced immediately; otherwise the live
 * store is unchanged. Selecting `to_key` afterwards starts from that copy.
 * A donor whose last pass did not converge is refused, including a parked
 * copy and a copy onto the same key. A copy that would park a new key past
 * CPMDC_ORBITAL_SLOT_LIMIT returns -1, cpmdc_last_error() names that
 * variable, and the copies already stored stay. The reference build has no
 * orbital store and returns 0, unless the limit is set to a value that is
 * not a positive integer, which fails the call. The link stub returns -1,
 * including a same key copy. Every rank of a calculator names the same
 * keys, the same rule `cpmdc_session_select_orbitals()` has.
 *
 * @return 0, or -1 when `session` is `NULL`, no evaluation has run yet,
 * `from_key` has no converged orbitals, or the copy cannot be completed.
 * The reason is `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_session_copy_orbitals(CPMDCSession *session, long long from_key,
                                          long long to_key);

/**
 * @brief Compute energy and nuclear gradient with session-owned parameters.
 *
 * Positions are Angstrom. The gradient buffer must have `n_atoms * 3` doubles
 * and is filled in Hartree/Bohr.
 */
CPMDC_API CPMDCResult cpmdc_session_energy_gradient(CPMDCSession *session, int n_atoms,
                                                    const double *positions_ang,
                                                    const int *atomic_numbers,
                                                    double *grad_h_bohr);

/**
 * @brief Compute total energy with session-owned parameters.
 *
 * Positions are Angstrom and the returned energy is Hartree.
 */
CPMDC_API CPMDCResult cpmdc_session_energy(CPMDCSession *session, int n_atoms,
                                           const double *positions_ang,
                                           const int *atomic_numbers);

/**
 * @brief Compute energy and nuclear forces with session-owned parameters.
 *
 * Positions are Angstrom. The forces buffer must have `n_atoms * 3` doubles and
 * is filled in Hartree/Bohr.
 */
CPMDC_API CPMDCResult cpmdc_session_energy_forces(CPMDCSession *session, int n_atoms,
                                                  const double *positions_ang,
                                                  const int *atomic_numbers,
                                                  double *forces_h_bohr);

/**
 * @brief Compute energy and forces for one Cap'n Proto `ForceInput` step.
 *
 * Session keeps persistent `CPMDParams`; each call supplies geometry. The
 * positions and cell are read in `ForceInput.lengthUnit`. The returned
 * energy and the forces use CPMD native units whatever the message names:
 * Hartree and Hartree/Bohr, forces being the negative nuclear gradient.
 * `forces_h_bohr` receives `3 * n_atoms` doubles and `forces_len` is its
 * capacity in doubles; a `forces_len` below `3 * n_atoms` fails with
 * `ok == 0`. The same topology and stored-result rules as
 * `cpmdc_session_energy_forces()` apply.
 */
CPMDC_API CPMDCResult cpmdc_session_calculate_forces(
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
 * When `potential_result_capnp_capacity_bytes` is too small, or
 * `potential_result_capnp` is `NULL`, returns `ok == 0`, writes the required
 * byte count to `potential_result_capnp_size_bytes`, and does not evaluate
 * CPMD. `potential_result_capnp_size_bytes` is required: a `NULL` there
 * fails with `ok == 0`. It is set to 0 when the `ForceInput` does not parse
 * or its geometry is invalid.
 * Size the buffer with `cpmdc_potential_result_size_for_force_input()`.
 */
CPMDC_API CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);

/**
 * @brief Compute energy, forces, and virial for one step in atomic units.
 *
 * The plain-array form of `cpmdc_session_calculate_result()` for callers
 * without a Cap'n Proto encoder (C, Fortran `BIND(C)`, LAMMPS plugins, i-PI
 * clients). The arguments come in the order of
 * `cpmdc_session_energy_forces()`: atom count, positions, atomic numbers.
 * Every quantity is atomic units.
 *
 * `session` comes from `cpmdc_session_create()` or
 * `cpmdc_session_create_from_config()` and is not yet destroyed. Every
 * pointer is borrowed for the call: the session copies what it keeps, and
 * the caller allocates and owns every buffer. The call is not thread-safe,
 * and with OpenCPMD linked it is collective on the calculator communicator,
 * so every rank of a calculator passes the same arguments in the same order.
 * The stub build returns -1.
 *
 * - `positions_bohr`: `3 * n_atoms` doubles, Bohr, atom-major row-major
 *   (`x0 y0 z0 x1 ...`).
 * - `atomic_numbers`: `n_atoms` values in 1..118. A value outside that range
 *   fails the call.
 * - `cell_bohr`: `NULL` for a step without a cell, or 9 doubles, Bohr,
 *   row-major with one lattice vector per row (`a_x a_y a_z b_x ...`), the
 *   layout of `ForceInput.box`. A cell with a non-zero entry replaces the
 *   deck's `CELL`. `NULL` and an all-zero cell give the deck the same
 *   treatment: its `CELL` stays, or a 12 Angstrom cube when it has none.
 *   The call still keeps an all-zero cell as a step with a cell, which is
 *   not a step without one. A `NULL` cell and an all-zero cell are different
 *   steps for the repeat match below. With OpenCPMD linked, stored orbitals
 *   follow the same rule: a step with an all-zero cell can start from them
 *   only when the latched step had a cell equal to it, never when it had
 *   none, and a `NULL` step only when the latched step had none.
 * - `energy_h`: receives the total energy, Hartree.
 * - `forces_h_bohr`: receives `3 * n_atoms` forces, Hartree/Bohr, laid out as
 *   the positions. Forces are the negative nuclear gradient, `-dE/dR`. On
 *   an isolated OpenCPMD deck the atoms are first moved rigidly so that their
 *   centre of mass is the middle of the deck's cell along every non-periodic
 *   lattice vector, and the net force of the grid is removed
 *   (`CPMDC_RECENTER`, default `mass`), so the energy and forces do not
 *   depend on where the caller puts the atoms. The caller's positions are
 *   not changed.
 * - `virial_h`: `NULL`, or receives 9 doubles, Hartree, row-major, the
 *   virial `-dE/d(strain)`. It is `|det(cell)|` times the stress of
 *   `cpmdc_session_last_stress()`, whose sign is `-(1/V) dE/d(strain)`.
 *   ASE's stress is `-virial / V`; i-PI takes this virial as it is. The
 *   virial is zero when `cell_bohr` is `NULL` or all zero, and when no
 *   stress was computed (an isolated deck or `CPMDC_STRESS=0`);
 *   `cpmdc_session_last_stress()` tells those apart.
 *
 * The step is the `cpmdc_session_calculate_result()` step of a `ForceInput`
 * with `lengthUnit = "bohr"`, `energyUnit = "hartree"`, and the same
 * positions, atomic numbers, and box (empty for `cell_bohr == NULL`), and
 * gives the same bits. `session_accept_topology` fixes the topology as soon
 * as an evaluation call passes its argument checks, even when the evaluation
 * then fails, and it is never rolled back. A repeat of the stored positions
 * and cell returns the stored result.
 *
 * A repeat matches the stored step atom by atom, in any atom order, when
 * the element agrees and every position component and every cell component
 * agrees within 1e-8 Angstrom, and when both steps have a cell or both have
 * none. A `NULL` cell and an all-zero cell are different steps for that
 * match. A repeat does not evaluate: it copies the stored forces into the
 * caller's atom order, returns the stored energy, and forms the virial from
 * the stress that evaluation stored.
 *
 * @return 0 on success; -1 on a `NULL` `session`, `positions_bohr`,
 *         `atomic_numbers`, `energy_h`, or `forces_h_bohr`; `n_atoms`
 *         outside 1..`INT_MAX / 3`; an atomic number outside 1..118; a
 *         change of the atom count or of the element counts after the
 *         topology is fixed; an allocation failure; or a failed evaluation.
 *         Failure text is `cpmdc_last_error()`, empty after a success.
 *         On -1 `energy_h` and `virial_h` are not written, and the contents
 *         of `forces_h_bohr` are unspecified.
 */
CPMDC_API int cpmdc_session_calculate_au(CPMDCSession *session, int n_atoms,
                                         const double *positions_bohr,
                                         const int *atomic_numbers,
                                         const double *cell_bohr, double *energy_h,
                                         double *forces_h_bohr, double *virial_h);

/**
 * @brief One-shot Cap'n Proto entry point (params + ForceInput -> PotentialResult).
 *
 * Multi-step callers should create one session and call
 * `cpmdc_session_calculate_result()` per step.
 */
CPMDC_API CPMDCResult cpmdc_calculate_result(const void *params_capnp,
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
CPMDC_API CPMDCResult cpmdc_calculate_result_from_config(
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
CPMDC_API size_t cpmdc_potential_result_size_for_force_input(
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
 * A `NULL` `capabilities_capnp_size_bytes` returns -1 and writes nothing.
 * That -1 does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_capabilities_result(void *capabilities_capnp,
                                        size_t capabilities_capnp_capacity_bytes,
                                        size_t *capabilities_capnp_size_bytes);

/**
 * @brief Compiled library version string: `cpmdc/` followed by the version,
 *        or `cpmdc-stub/` followed by the version in the link stub.
 */
CPMDC_API const char *cpmdc_version(void);

/**
 * @brief Diagnostic for the most recent public call on this thread that
 *        reports failure through this string.
 *
 * Written by `cpmdc_params_method()`, `cpmdc_set_params()`,
 * `cpmdc_configure()`,
 * `cpmdc_bind_calculator()`, `cpmdc_adopt_calculator_comm()`,
 * `cpmdc_session_create()`,
 * `cpmdc_session_set_params()`, `cpmdc_session_set_constant_cutoff()`,
 * `cpmdc_session_select_orbitals()`, `cpmdc_session_copy_orbitals()`,
 * `cpmdc_session_create_from_config()`,
 * `cpmdc_session_configure()`, `cpmdc_potential_result_size_for_force_input()`,
 * and every evaluation entry point (`cpmdc_energy()`,
 * `cpmdc_energy_gradient()`, `cpmdc_energy_forces()`, the session variants,
 * `cpmdc_session_calculate_au()`, and the calculate-result calls). When the call returns `CPMDCResult`, the
 * text matches `message`, including a missing pseudopotential directory.
 * Empty after a successful call.
 *
 * Snapshot readers (`cpmdc_last_stress()` and the other `cpmdc_last_*` /
 * `cpmdc_session_last_*` getters) do not write it: -1 means the snapshot is
 * absent or the output pointer is null. `cpmdc_adopted_comm()` and
 * `cpmdc_capabilities_result()` do not write it either. A capabilities
 * return of -1 is the size query.
 */
CPMDC_API const char *cpmdc_last_error(void);

/**
 * @brief Numeric ABI generation of the compiled library.
 *
 * Compare against the CPMDC_ABI_VERSION the consumer compiled with.
 */
CPMDC_API int cpmdc_abi_version(void);

/** @brief 1 when OpenCPMD is linked, ready, and not finalized; 0 on the reference build; 0 after finalize and from the link stub. */
CPMDC_API int cpmdc_available(void);

/**
 * @brief Mark the embedded runtime finalized.
 *
 * The reference build and an OpenCPMD build forget the active session, free
 * the stored input deck, and mark the runtime finalized for the rest of the
 * process. `cpmdc_available()` returns 0 from then on. Every call that needs
 * the runtime then fails with "CPMD embed not available": the no-session
 * evaluations, `cpmdc_set_params()`, `cpmdc_configure()`, and a session
 * evaluation of any step other than the one that session stored last.
 * `cpmdc_session_create()` still returns a session, but that session can only
 * answer a repeat. A session evaluation that repeats its stored step, by the
 * match rule of `cpmdc_session_calculate_au()`, still returns 0: it copies
 * the stored result from the session's cache without the runtime and leaves
 * `cpmdc_last_error()` empty. A step that fails does not clear that cache.
 *
 * The link stub has no runtime, so there it does nothing.
 * `cpmdc_available()` is 0 before and after, and every evaluation fails with
 * "OpenCPMD embed not available in cpmdc stub (cpmdc_energy)", with the
 * name of the function that was called in the parentheses.
 *
 * Neither build calls `MPI_Finalize`, destroys sessions, or frees OpenCPMD
 * module state.
 */
CPMDC_API void cpmdc_finalize(void);

/**
 * @brief Copy the last in-process `ener_com` energy decomposition.
 *
 * The no-session entry points read the active calculator. A session entry
 * point reads that session, including after another session has run.
 * Returns 0 when `out->valid` is set. Values are Hartree.
 */
CPMDC_API int cpmdc_last_energy_components(CPMDCEnergyComponents *out);

/** @brief Copy the last energy decomposition stored on @p session. */
CPMDC_API int cpmdc_session_last_energy_components(const CPMDCSession *session,
                                                   CPMDCEnergyComponents *out);

/**
 * @brief Copy the last in-process density integrals.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_charge_integrals(CPMDCChargeIntegrals *out);

/** @brief Copy the last density integrals stored on @p session. */
CPMDC_API int cpmdc_session_last_charge_integrals(const CPMDCSession *session,
                                                  CPMDCChargeIntegrals *out);

/**
 * @brief Copy the last multi-state energy catalog.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_multi_state_energies(CPMDCMultiStateEnergies *out);

/** @brief Copy the last multi-state energy catalog stored on @p session. */
CPMDC_API int cpmdc_session_last_multi_state_energies(const CPMDCSession *session,
                                                      CPMDCMultiStateEnergies *out);

/**
 * @brief Copy the last trajectory row.
 *
 * Values are Hartree. Returns 0 when `out->valid` is set. Does not write
 * `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_md_trajectory_row(CPMDCMDTrajectoryRow *out);

/** @brief Copy the last trajectory row stored on @p session. */
CPMDC_API int cpmdc_session_last_md_trajectory_row(const CPMDCSession *session,
                                                   CPMDCMDTrajectoryRow *out);

/**
 * @brief Copy the last property snapshot.
 *
 * Returns 0 when `out->valid` is set. Does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_property_snapshot(CPMDCPropertySnapshot *out);

/** @brief Copy the last property snapshot stored on @p session. */
CPMDC_API int cpmdc_session_last_property_snapshot(const CPMDCSession *session,
                                                   CPMDCPropertySnapshot *out);

/**
 * @brief Copy the nuclear gradient of the last evaluation.
 *
 * `capacity` is the number of doubles available in `out`. `count` receives
 * the required number. Returns 0 on success and -1 when no gradient is
 * available, an argument is null, or the buffer is too short. The legacy
 * property snapshot limits successful evaluations to 4096 components.
 */
CPMDC_API int cpmdc_last_gradient(double *out, size_t capacity, size_t *count);

/** @brief Copy the last nuclear gradient stored on @p session. */
CPMDC_API int cpmdc_session_last_gradient(const CPMDCSession *session,
                                          double *out, size_t capacity,
                                          size_t *count);

/**
 * @brief Copy the molecular dipole of the last evaluation.
 *
 * Returns 0 only when the dipole was recomputed for that geometry. Returns -1
 * without writing `cpmdc_last_error()` when it was unavailable.
 */
CPMDC_API int cpmdc_last_dipole(CPMDCDipole *out);

/** @brief Copy the last molecular dipole stored on @p session. */
CPMDC_API int cpmdc_session_last_dipole(const CPMDCSession *session,
                                        CPMDCDipole *out);

/**
 * @brief Copy the last stress tensor.
 *
 * Returns 0 when `out->valid` is set, which is only after stress was computed
 * for a periodic cell. Does not write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_stress(CPMDCStressTensor *out);

/** @brief Copy the last stress tensor stored on @p session. */
CPMDC_API int cpmdc_session_last_stress(const CPMDCSession *session,
                                        CPMDCStressTensor *out);

/**
 * @brief Copy the rank spread of the last evaluation in this process.
 *
 * Returns 0 when `out->valid` is set: the last evaluation ran with
 * `CPMDC_RANK_CONSENSUS=1` under a running MPI. A repeat of the stored
 * positions keeps the spread of the evaluation that stored them. Does not
 * write `cpmdc_last_error()`.
 */
CPMDC_API int cpmdc_last_rank_spread(CPMDCRankSpread *out);

#ifdef __cplusplus
}
#endif
