# Changelog

<!-- towncrier release notes start -->

## [0.4.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.4.0) - 2026-11-13

### Added

- A job TOML compiles to the CPMDParams text. The in-process route encodes that text, and the file route renders the same message.
- An optional plane-wave scaffold times strong and weak scaling of a structure-factor force, and keeps a CUDA backend behind a build flag that stays off unless requested.
  The scaffold stays out of the default build, so an ordinary build does not require OpenMP for it.
  It does not change the embed force, the reference potential, or an SCF.
- Python `Session` is a handle over libcpmdc. `Potential` matches the rgpot evaluate shape, in eV and eV/Angstrom. `CPMDCCalculator` is an ASE calculator for energy, forces, and a fully periodic stress. The narrative tutorials and both Python examples run in the test suite.
- The README and the docs index now name the license and link `LICENSE`.
- The README now points questions and bug reports to the GitHub issue tracker.
  The docs index and the troubleshooting guide link it too.
  The troubleshooting guide also lists what a new report needs, with a command that prints the version and build options.
- The contributing guide now opens with the contributions that are welcome and how to send them.
  The README links that guide from Start Here and from a new Contributing section.
- The exported `cpmdc_*` symbols are compared with `scripts/abi_baseline.txt` on every test run.
- `CITATION.cff` gives the citation record for cpmdc.
  It names the author, the version, the license, and the repository.
  `pixi run version-sync` regenerates its version from `version.txt`.
  The release check fails when the two differ.
  The README and the docs index gain a Citation section.
  That section also asks papers with OpenCPMD energies to cite OpenCPMD.
- `CPMDC_RANK_CONSENSUS=1` gives every rank of a calculator the energy, gradient, and result snapshots of its rank 0 after each evaluation, so host optimizers that run on every rank step from the same numbers. Every rank fails when one rank failed, and a stored-result repeat counts only when every rank repeats. `cpmdc_last_rank_spread()` reports the largest energy and gradient difference from rank 0 before the broadcast, with a NaN difference reported as infinity. Unset, every call behaves as before.
- `cpmdc-fd-check` compares the forces of a session with central differences of its energy, after the CP2K `DEBUG_FORCES` check with `DX = 1e-3` Bohr.
  `--stress` adds the symmetric-strain check of `PotentialResult.stress` after `DEBUG_STRESS_TENSOR`.
  The tool is built with the tests and is not installed.
  On the reference evaluator the forces and the periodic stress agree to rounding.
- `cpmdc-restart export-con` and `import-con` are available when Meson is configured with `-Dwith_readcon=true`. The bridge multiplies by `CPMDC_RESTART_ANGSTROM_PER_BOHR`, writes a RoundTrip CON, and calls `cpmdc_restart_set_cell` on import. The default build does not link readcon. The CON file is the geometry that is resumed, and the binary RESTART stays the wavefunction.
- `cpmdc_session_select_orbitals()` keeps the converged orbitals per key, so a calculator that evaluates several images of a band or beads of a ring polymer in turn starts each SCF from that image's or bead's own previous orbitals. `tools/opencpmd_embed_rwfopt.patch` adds the per-key table (`embed_select_orbital_slot`) to OpenCPMD; a basis change frees every key's orbitals.
- `tools/opencpmd_embed_extrapolate.patch` honours `EXTRAPOLATE WFN` between in-process calls.
  OpenCPMD turns the keyword off outside Born-Oppenheimer molecular dynamics.
  A host-driven trajectory therefore started each self-consistent field (SCF) from the unchanged orbitals of the last call.
  With the patch, `rwfopt` keeps a ring of the converged orbitals of the previous calls and starts the SCF from the `extrapwf` prediction of Car-Parrinello molecular dynamics (CPMD), with the always stable predictor-corrector (ASPC) or the polynomial coefficients.
  The in-process route reached an SCF iteration count of 297 on a 40-step water trajectory with ASPC of order 4, against 292 for `cpmd.x` with the same deck and 441 for `cpmd.x` without extrapolation.
  `CPMDC_SCF_STEPS=1` prints `cpmdc_wfn_extrapolated K` for a call that started from a prediction over `K` earlier calls.
  Decks without the keyword run as before.
- `tools/opencpmd_embed_teardown.patch` adds `embed_teardown_state` to OpenCPMD. It frees the module arrays, the setup communicators, and the saved orbitals of a setup. It also resets the first-call state of 37 routines on the setup and wavefunction optimisation path. The rebuild helper and the archive how-to apply it after the other patches.

### Changed

- Stored orbitals start an SCF only when they converged within `CPMDC_WARM_RMSD` (root-mean-square distance per atom, default 0.15 Angstrom) of the new geometry; the nearest stored copy of any key is taken, and beyond the limit the deck's guess starts the SCF. On the Si3N4 OCI-NEB a far key had taken 143 to 530 SCF steps from another image's orbitals where the guess takes 37 to 63. `CPMDC_SCF_STEPS=1` also prints where each SCF started. A line that carries a distance also names the donor key and whether the distance used the periodic boundary.
- The `ForceInput` box is the cell: a non-zero box replaces the `CELL` and `SYMMETRY` of the deck (`SYMMETRY 1` for a cubic deck and box, `SYMMETRY 8` for another orthorhombic box, `CELL VECTORS` otherwise), as a file-based caller that writes its own cell does. An isolated deck takes only an orthorhombic box. An all-zero box keeps the deck's cell.
- The opt-in scaffold's `cpmdc_pwscale_eval_mpi` runs 2 collectives per
  call, down from 12 for the plane-wave split and 11 for the atom split.
  `cpmdc_pwscale_mpi_bind_create` checks a problem once, and
  `cpmdc_pwscale_eval_mpi_bound` then runs one `MPI_Allgather` per call.
  Every rank sums the partials in rank order, so energy and forces keep
  their bytes, match on every rank, and repeat from run to run.
  `bench/run_pwscale.sh` takes `BIND=core`.
  Its default, `none` with `OMP_PROC_BIND=close`, puts every rank on
  CPU 0; on an 8-core host 2, 4, and 8 ranks took 0.0960, 0.180, and
  0.416 s against 0.0747 s for one thread, and pinned ranks took 0.0380,
  0.0211, and 0.0143 s.
  At 16 atoms and 64 plane waves an 8-rank strong call went from 20.1 to
  7.6 us; at the default size the change stayed within the spread of two
  equal builds.
  `bench/run_pwscale_abba.sh` compares two builds in ABBA order, and
  `tests/test_pwscale_collectives.c` counts collectives through a PMPI
  shim.
- The reference evaluator's periodic stress changed sign. This is a behaviour change. The tensor was the energy divided by the volume. It is now `-dE/d(eps) / volume`, the OpenCPMD sign, so `cpmdc-fd-check --stress` agrees with it. The Python provider flips that sign for ASE and rgpot, which want `(1/V) dE/d(strain)`.
- The x86_64 OpenCPMD PIC file passes `-march=x86-64`, and the aarch64 file passes `-march=armv8-a`, on debug and release `FFLAGS` and `CFLAGS`.
  GCC 11 on RHEL 9 defaults to x86-64-v2, and GCC 14 on RHEL 10 defaults to x86-64-v3 with AVX2 and FMA.
  Both files also pass `-fPIC` and `-ffp-contract=off`.
  cpmdc C and Fortran compile with `-ffp-contract=off` so a multiply and an add stay separate on every build target.

### Fixed

- A positive-volume cell with no isolated keyword is written with cubic or orthorhombic symmetry, or with cell vectors when the box is tilted, while a zero box and a deck that already asks for an isolated system stay on symmetry 0 with the Hockney solver. ([#periodic-cell-symmetry](https://github.com/OmniPotentRPC/cpmdc/issues/periodic-cell-symmetry))
- A cell with off-diagonal components no longer stops the host with a Fortran "End of file" while the cold deck is written: the `CELL VECTORS` rows are joined by newlines inside one record.
- A session that already stored a result rejects a null position or atomic-number buffer with `invalid arguments` instead of reading the pointer.
- An opt-in MPI setup failure now returns on every rank together, and an FFT setup failure releases the buffers it already took.
  The optional CUDA compile line passes the scaffold include directory, and this change does not run CUDA.
- One process can set OpenCPMD up again, for a new basis or after a CPMD stop. The embed calls `cpmdc_embed_teardown` before a second cold setup and from `cpmdc_embed_reset_state`. It frees what the first setup allocated, and the second setup gives the energies of a fresh process. A second setup used to stop in `ratom` with `allocation problem`, and the embed refused it.
- Temporary CON files, parameter messages, and the quickstart encode step use the working directory. The quickstart writes encoded messages under `tutorial-work/`. A tracked-file check rejects host-local directory names.
- The Version section of the README described 0.1.0.
  `version.txt` held 0.3.0.
  The section now names the version from `version.txt`.
  `pixi run version-sync` regenerates that file.
  It also points to `CHANGELOG.md` for the changes in each release.
- The `cpmdc_warm_start` and `cpmdc_scf_steps` lines of `CPMDC_SCF_STEPS=1` now reach a standard error that is redirected to a file at the end of each call.
  The gfortran runtime buffers that unit when it is a regular file, so the lines used to appear only at exit, or not at all after an abort.
- The cold deck writes the caller's cell with twelve decimals. The CELL ratios b/a and c/a carried six, which moved b or c by up to 7e-6 A in a 14 A cell; `CELL VECTORS` carried eight.
- The first energy tutorial now compiles its host program with `mpicc` and links it with `mpif90`.
  The wrappers are needed because the default `libcpmdc.so` links the Message Passing Interface (MPI) library whenever Meson finds one.
  The tutorial transcripts now come from a run of version 0.3.0 against both the default build and a patched OpenCPMD archive.
- `cpmdc_adopt_calculator_comm` reads rank and size from the communicator the host passed.
  When `ranks_per_calc` is positive that size must match, and the index is the lowest world rank of one aligned block divided by that size.
  A later call replaces the communicator while no session exists.
- `tools/opencpmd_omp_info.patch` makes every rank of a calculator use the I/O node's OpenMP report size. A calculator that spans nodes with unequal rank counts otherwise calls `MPI_Allreduce` with different counts and can stop before the first electronic iteration.
- `tools/opencpmd_pcg_step.patch` divides both PCG line-search derivatives by their largest magnitude before finding the root, so a small orbital residual keeps its interpolated step. Equal derivatives keep the trial step. The convergence target is unchanged.

### Developer

- CI calls the rgpot ISA audit on unstripped `libcpmdc` and the test binaries.
  The allowlist names no cpmdc symbol.
- Sphinx treats a warning as a failed book. `docs/render_reference_well.py` draws the reference-well figure from the evaluator formula.
- The geometry-aware native orbital selector is the versioned symbol `embed_choose_warm_start_metric_v2`. The unversioned name stays undefined.
- The instruction set architecture (ISA) allowlist names the minimage kernels built for Advanced Vector Extensions (AVX), AVX2 or AVX-512.
  minimage calls each kernel only after the Rust standard library reports its features at run time.
  The allowlist also names three functions that minimage links but does not define.
  One rule matches `find_impl` and `find_prefilter_impl` of the AVX2 packed-pair searcher in the `memchr` crate, which the standard library depends on through the `object` crate.
  The standard library build links both functions and never calls them.
  The other rule names `_xgetbv` from the Rust core library.
  It runs only after the processor identification instruction (CPUID) reports that the processor and the operating system support extended state saving.
  The portable ISA test rejects an allow rule that is not anchored, uses a wildcard, matches a cpmdc symbol, or states no run-time check.
- The meson suite runs every stub-library documentation example. Each fenced block is executed or named with the tool it needs. The Python how-to, the finite-difference check, the Cap'n Proto messages, and the first-energy and ASE lessons run there. CI installs ASE before that suite.
- The version is written by hand only in `version.txt`.
  Meson and Sphinx read it, and `pixi run version-sync` regenerates the other copies with cog.
  CI fails when a copy differs.
- `pixi run lychee` reads the new `.lychee.toml`.
  It no longer stops on a missing file.
  The configuration accepts 429 responses from GitHub.
  It retries three times.
  It also checks the anchors of links such as `README.md#citation`.
  The contributing guide names the task.

### Miscellaneous

- `docs/AUDIT.md` grades the One Good Tutorial checklist again.
  The grade sits on top of `c631ed9`.
  Eight items pass, and Acknowledgments does not apply.
  The open gaps are the quickstart placement and the missing eOn output.
  The inventory tests still require the README tables.


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
