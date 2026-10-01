# Changelog

<!-- towncrier release notes start -->

## [0.3.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.3.0) - 2026-10-01

### Added

- Wire public **capnp-fortran** into the embed apply path: Fortran decodes
  serialized `CPMDParams` (functional / cutOffRy / charge / multiplicity /
  cpmdRoot, including system and dft section overrides) via
  `cpmdc_embed_apply_params`. Exploded `cpmdc_embed_set_config` is unused on the
  serialized-params path; C still renders the input deck and geometry merges
  still use `set_deck`. ([#capnp-fortran-apply](https://github.com/OmniPotentRPC/cpmdc/issues/capnp-fortran-apply))
- PEF embed enables OpenCPMD `cntl%tpres`, snapshots Cartesian stress
  (`paiu/omega` as Ha/Bohr^3, row-major 9-vector), exports via
  `cpmdc_last_stress` and `PotentialResult.stress` with ForceInput unit conversion. ([#embed-pef-stress](https://github.com/OmniPotentRPC/cpmdc/issues/embed-pef-stress))
- Add `tests/test_cpmd_inscan_token_fidelity.py` to check typed-section render
  tokens against OpenCPMD input sources (CONSTRAINTS/ISOTOPE/VELOCITIES/…). ([#inscan-token-fidelity](https://github.com/OmniPotentRPC/cpmdc/issues/inscan-token-fidelity))
- `cpmdc-restart` and `libcpmdc_restart` read and write a CPMD `RESTART` file,
  including an in-place coordinate patch that keeps the wavefunction records. ([#restart-file](https://github.com/OmniPotentRPC/cpmdc/issues/restart-file))
- `examples/strasbourg_blyp_ext_pot` runs one BLYP wavefunction optimization for an eOn image. ([#strasbourg-blyp-ext-pot](https://github.com/OmniPotentRPC/cpmdc/issues/strasbourg-blyp-ext-pot))
- `cpmdc_adopt_calculator_comm` stores a communicator the host already split and does not call `MPI_Comm_split`.

### Changed

- Cold OpenCPMD SCF prefers Cap'n-rendered `applied_input_deck` when it includes
  real `&ATOMS` PP lines. Method-only decks (empty C-render `&ATOMS` placeholder)
  keep typed `&CPMD`/`&SYSTEM`/`&DFT` text and merge geometry atoms from the C
  arrays instead of rebuilding a minimal BLYP deck. ([#embed-applied-deck-cold](https://github.com/OmniPotentRPC/cpmdc/issues/embed-applied-deck-cold))
- Linked and cold-deck OpenCPMD embed paths honor the applied DFT functional,
  system charge, and multiplicity from Cap'n Proto apply state (no longer
  hard-wire FUNCTIONAL BLYP only). ([#embed-honor-applied-functional](https://github.com/OmniPotentRPC/cpmdc/issues/embed-honor-applied-functional))
- A warm in-process call that holds a converged orbital copy returns from the wavefunction setup after the phase factors and does not build a starting guess. It writes no RESTART.1 or LATEST. The cold call still builds the guess. ([#embed-warm-guess-restart](https://github.com/OmniPotentRPC/cpmdc/issues/embed-warm-guess-restart))
- The MPI model states that cpmdc has no `with_mpi` option and takes MPI from the CPMD archive it is linked into. ([#mpi-from-archive](https://github.com/OmniPotentRPC/cpmdc/issues/mpi-from-archive))
- A library force call writes no RESTART.1, LATEST, GEOMETRY, or GEOMETRY.xyz. The orbitals for the next call stay in memory.

### Fixed

- The OpenCPMD path takes each listed element's pseudopotential file, LMAX, LOC, and KLEINMAN-BYLANDER from the CPMDParams message. KLEINMAN-BYLANDER is written on that element's *file line. Elements the message does not list still use the built-in table (H, C, N, O, Si, Ge). Any other element fails with its atomic number. ([#atoms-message-pseudopotentials](https://github.com/OmniPotentRPC/cpmdc/issues/atoms-message-pseudopotentials))
- A warm start keeps `c0` only after the SCF converges. An unconverged ODIIS
  pass leaves the previous orbitals, and a PCG MINIMIZE continuation restarts
  from that copy. A cutoff, cell, or state-count change drops the saved
  orbitals before the next SCF. ([#converged-warm-orbitals](https://github.com/OmniPotentRPC/cpmdc/issues/converged-warm-orbitals))
- cpmd.x writes RESTART.1, LATEST, GEOMETRY, and GEOMETRY.xyz in permanentDir, or in scratchDir when permanentDir is empty, or in the host working directory. The pseudopotential directory is left unchanged. ([#cpmd-output-directory](https://github.com/OmniPotentRPC/cpmdc/issues/cpmd-output-directory))
- Live OpenCPMD PEF multi-force now requests ionic forces via
  `embed_set_need_forces` (vendored OpenCPMD `tfor` patch), documents
  `CPMDC_PSEUDO_DIR` cold memfd resolution, and guards the force-export path in
  tests. ([#embed-bomd-forces](https://github.com/OmniPotentRPC/cpmdc/issues/embed-bomd-forces))
- Live OpenCPMD evaluations retain the converged orbital state that produced the
  reported energy and forces instead of applying and saving an extra optimizer
  update after the convergence test is satisfied. ([#embed-converged-state](https://github.com/OmniPotentRPC/cpmdc/issues/embed-converged-state))
- The stress snapshot is valid only after the tensor was computed for a periodic cell. An isolated cell leaves it unset even when the box volume is positive, and `CPMDC_STRESS=0` leaves it unset on a periodic cell too. ([#isolated-stress-invalid](https://github.com/OmniPotentRPC/cpmdc/issues/isolated-stress-invalid))
- `cpmdc_last_error()` is set by every configuration, session-setup, bind, size, and evaluation call that fails, and cleared when that call succeeds. Snapshot readers and the capabilities size query do not write it. ([#last-error-on-failure](https://github.com/OmniPotentRPC/cpmdc/issues/last-error-on-failure))
- The embed path keeps a rendered method deck at the length of the text. A preview or config buffer shorter than that text returns an error instead of a shortened deck. ([#long-method-deck](https://github.com/OmniPotentRPC/cpmdc/issues/long-method-deck))
- Rebuild patched OpenCPMD module objects and `libcpmd.a` through an explicit
  dependency-target helper, avoiding a `make lib` no-op against the existing
  `lib/` directory. ([#opencpmd-rebuild](https://github.com/OmniPotentRPC/cpmdc/issues/opencpmd-rebuild))
- A call with neither CPMDC_PSEUDO_DIR nor CPMD_PP_LIBRARY_PATH set now fails with its own message on the result and on cpmdc_last_error. ([#pseudopotential-directory](https://github.com/OmniPotentRPC/cpmdc/issues/pseudopotential-directory))
- Atoms of one element may be split in ForceInput. Positions are placed in CPMD species blocks and forces are returned in ForceInput order. ([#species-order-forces](https://github.com/OmniPotentRPC/cpmdc/issues/species-order-forces))
- The Sphinx version and release are the project version in meson.build. ([#sphinx-version](https://github.com/OmniPotentRPC/cpmdc/issues/sphinx-version))
- The embed link includes `obj/timetag.o`. OpenCPMD keeps that compilation
  stamp out of `libcpmd.a`, and `header` calls it. ([#timetag-embed-link](https://github.com/OmniPotentRPC/cpmdc/issues/timetag-embed-link))
- The OpenCPMD archive patches follow
  [pull request 9](https://github.com/OpenCPMD/CPMD/pull/9). `tistopgm` prints
  the call stack only when `trace_depth` lies inside `trace_names`, so a stop
  before `tistart` writes `LocalError` instead of faulting. ([#tistopgm-before-tistart](https://github.com/OmniPotentRPC/cpmdc/issues/tistopgm-before-tistart))
- A warm OpenCPMD step keeps its orbitals only while the cell matches the cold step. A different cell rebuilds the plane-wave basis. Silicon with an omitted projector channel is rendered as `LMAX=D LOC=D`. ([#warm-cell-stress](https://github.com/OmniPotentRPC/cpmdc/issues/warm-cell-stress))
- A `common` overlay's `scfEnergyToleranceEv` reaches `CONVERGENCE ORBITALS` as the total-energy change as well as the orbital threshold, so a smeared SCF stops once its energy settles within the tolerance.
- A call that needs a second OpenCPMD setup in one process, for a new basis or after a CPMD stop, fails with an error instead of ending the host. OpenCPMD's `tistart` executes `STOP` on its second call, which exited the process with status 0 in the middle of the call, and its setup routines allocate module arrays that nothing frees.
- A new session keeps the stored orbitals when the cutoff, cell, charge, multiplicity, functional, deck, and elemental composition match. Reordering the atomic numbers keeps them.
- A raw input block that already opens `&CPMD`, `&SYSTEM`, `&DFT`, or `&ATOMS` is that section. The renderer does not append a second copy, so a periodic `SYMMETRY 1` block is not followed by the isolated Hockney system.
- A session keeps a reordering of the same elements. A repeat of those positions returns the stored energy and forces.
- An OpenCPMD stop in the deck setup of an energy and gradient call, such as `HFX not implemented for HOCKNEY`, ends the call with an error before the SCF. The SCF ran on the half-initialized state and could crash the host.
- An energy and gradient call whose orbitals do not converge within MAXITER now fails with a message. OpenCPMD computes no ionic forces then, and the call used to report success with a zero gradient.
- An isolated (`SYMMETRY 0`) deck with a hybrid functional and no Poisson solver gets `POISSON SOLVER TUCKERMAN` in place of the `HOCKNEY` default, which OpenCPMD's exact exchange does not support.
- The Cap'n Proto renderer writes `MULTIPLICITY` into `&SYSTEM`, from the system section or else from the top-level `multiplicity`, and a multiplicity above 1 named only in the system section turns on `LSD`. The deck had no `MULTIPLICITY` line, so a triplet request ran as the OpenCPMD default spin state.
- The archive how-to applies `opencpmd_embed_geometry.patch` before `opencpmd_embed_rwfopt.patch`, then `opencpmd_c_mem_addrs.patch`, and copies `tools/LINUX-X86_64-GFORTRAN-MPI-PIC` into the clone before `configure.sh`. That configuration adds `-fPIC`, links `-lopenblas`, and adds `-fallow-argument-mismatch` when `mpif90` accepts it. The address patch makes `cGetMemAddrs` return the pointer as `size_t`.
- The installed-package consumer test configures its nested build with the C and Fortran compilers of the outer build, so a compiler launcher in the host's `CC` no longer breaks the nested `meson setup`.
- The stopgm hook patch only adds lines. The hook still runs before the call stack is printed, and a handled stop closes the log and returns.
- `cpmdc_available()` returns 1 only when OpenCPMD is linked and the runtime is ready. The reference build, a finalized runtime, and the link stub return 0.
- `tools/rebuild_opencpmd_embed.sh` rebuilds and re-archives one object per vendored OpenCPMD patch: `embed_ctrl.mod.o`, `initrun_driver.mod.o`, and `c_mem_utils.o` join the list beside `rinitwf_driver.mod.o`.

### Developer

- The effective-config XC fixtures ask for a doublet. Their HO molecule with charge +2 has 5 valence electrons, so the triplet they named cannot exist; it ran only because the renderer never wrote `MULTIPLICITY`.
- The effective-config fixtures run PBE0 through the XC driver (`USE_XC_DRIVER`, `FUNCTIONAL HYB_GGA_XC_PBE0`). OpenCPMD refuses old-code PBE0 with LSD, since its VWN correlation has no LSD form; cpmd.x converges the driver deck at -14.72424 Ha.
- The session-socket, e2e-single-point, e2e-set-directives-session, e2e-optimizer-session, stress-isolated and effective-config tests run each OpenCPMD setup in a child process and assert on the results it hands back. They exited with status 0 at the second in-process setup, through OpenCPMD's `tistart` STOP, so their later assertions never ran against OpenCPMD.
- The water parity test compares against -17.13449267 Ha, the energy that cpmd.x and the embed both give for the deck the water fixture renders (agreement to 1e-8). The old reference, -17.04926427 Ha, matched no rendering of that deck.


## [0.2.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.2.0) - 2026-07-03

### Added

- `cpmdc_abi_version()` and `CPMDC_ABI_VERSION` for numeric ABI gating.
- pkg-config (`cpmdc.pc`) and CMake package config (`find_package(cpmdc)`)
  for the installed library, with namespaced `include/cpmdc/` headers and a
  versioned SONAME.
- Release scaffolding: towncrier changelog, cog bump wiring, and
  `scripts/release_assert.py` version-sync gate.

### Changed

- `schema/Potentials.capnp` is now the canonical OmniPotentRPC
  `potentials-schema` release (v1.0.1), pinned via meson wrap with a
  byte-identical sync test; `ForceInput` gains the shared per-step
  charge/multiplicity override fields.
- `cpmdc_version()` reports the project version from the build instead of a
  hardcoded string.

### Fixed

- The wire `PotentialResult.hessian` field no longer aliases the nuclear
  gradient; gradient data stays on `gradient`.
- `cpmdc_potential_result_flat_size()` now budgets for both natoms*3 lists
  (forces and gradient), fixing result writes for systems past a few dozen
  atoms.

## [0.1.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.1.0) - 2026-06-28

- OpenCPMD driver in the ISO_C embed shell for working E2E single-point and
  multi-step session tests without OpenCPMD archives.
- Docs tree with Antics analytics, pixi/prek/CI wiring.
- Initial OmniPotentRPC `cpmdc` package modeled on `nwchemc`.
- Cap'n Proto `ForceInput` / `PotentialResult` / `CPMDParams` / `PotentialConfig`.
- Stable C ABI with session direct-call socket entry points.
- Fortran `iso_c_binding` embed shell (`cpmd_embed_c_api.f90`).
- Stub ABI and cmocka tests for parser + result sizing without OpenCPMD.
