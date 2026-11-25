# Changelog

<!-- towncrier release notes start -->

## [0.7.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.7.0) - 2026-11-24

### Deprecated

- `CPMDCPropertySnapshot.hessian` and `hessian_count` are deprecated. The array
  contains the nuclear gradient dE/dR, not a Hessian. Existing binaries retain
  the field; new callers should use `cpmdc_last_gradient()` or
  `cpmdc_session_last_gradient()`.

### Added

- GROMACS-generated PBE and BLYP QM/MM inputs can load a matching
  `CPMDParams` message from `CPMDC_GROMACS_PARAMS`.

  The C ABI exposes `cpmdc_params_method()` and Python exposes
  `params_method()` for reading a message's effective functional, charge, and
  multiplicity without applying it.
- The GROMACS interface caps a covalent bond between a QM atom and an MM atom
  with a link hydrogen, so a QM group no longer has to hold whole molecules.
  The library reads each `&LINK` section of the generated input, or of a
  `-qmi` file, places a hydrogen at the IMOMM position of CP2K with
  `ALPHA_IMOMM` 1.38 unless given, adds it to the QM calculation, and returns
  its force to the two atoms of the bond in the shares of the chain rule. The
  virial, MPI runs, and mechanical embedding follow, and the MM atom of a cut
  bond leaves the QM-MM Coulomb sum.

  `cpmdc.gromacs.qmi_input()` and `write_qmi()` take `links` and `link_alpha`.
- `cpmdc.ase` registers `CPMDCCalculator` under `"cpmdc"` with `ase.calculators.calculator.register_calculator_class` (idempotently, so `get_calculator_class("cpmdc")` returns `CPMDCCalculator` and repeated registration does not raise), and `pyproject.toml` exposes `cpmdc = "cpmdc.ase"` in the `ase.plugins` entry-point group. The module builds `__ase_plugins__`, a `CalculatorPlugin` set, on first read, so `ase._4.plugins.plugins.calculators` lists `cpmdc` whichever of `cpmdc.ase` and `ase._4.plugins` is imported first, and the set is empty on ASE releases without the plugin API.
- `cpmdc_last_dipole()` and its session form return a molecular dipole only when
  linked OpenCPMD recomputed it for the current geometry. Python `Session` and
  `Potential` expose that value, and `CPMDCCalculator` provides ASE's `dipole`
  property in e*Angstrom. Reference, isolated, k-point, and unrequested
  evaluations report the property as unavailable instead of returning a zero or
  stale value.

  The dipole is a Berry-phase quantity of the periodic cell. It is defined only
  for a neutral cell, and only modulo the lattice quantum: it can jump by an
  integer number of electron charges times a lattice vector between geometries.
  `cpmdc_session_calculate_au()` has no dipole argument, because its argument
  list is fixed for the LAMMPS, GROMACS, and i-PI callers that bind to it; read
  `cpmdc_session_last_dipole()` after it returns 0.
  `Potential.caps()` reports `dipole` as true only when the evaluator can produce
  one: never for the reference evaluator, and for OpenCPMD until an evaluation
  comes back without a dipole. `Session.dipole_possible()` gives the same answer.
- `libcpmdc_gromacs` supports `&CPMDC EPSILON_R` and
  `cpmdc.gromacs.qmi_input(..., epsilon_r=...)` to divide the
  mechanical-embedding QM-MM Coulomb energy, forces, and Ewald virial by a
  positive relative dielectric constant while documenting GROMACS'
  unforwarded `epsilon-r` and six-decimal PDB charge formatting.

### Changed

- A deck that sets `DIPOLE DYNAMICS` and `KPOINTS` is now rejected when the parameters are installed, with an error that names both keywords, instead of stopping OpenCPMD in `localize2` during the first evaluation.
  The check reads the rendered deck, so it also covers parameters given as raw input blocks.
  Use a Gamma-point deck to get a dipole.
- For rgpot and other consumers of `PotentialResult` and `CPMDCPropertySnapshot`, the dipole fields now say whether a dipole exists.
  `CPMDCPropertySnapshot.dipole_count` is 3 only when OpenCPMD recomputed the dipole for the current geometry and 0 otherwise; it was always 3, with zeros from the reference evaluator and the stale `ddip%pdipole` from OpenCPMD.
  `PotentialResult.dipole` is an empty list in the same cases and was three zeros.
  When present, it is the molecular dipole `dipomod%moment%dmom` in electron Bohr, with the ionic term, in place of `ddip%pdipole`.
  A consumer must read an empty `dipole` list or `dipole_count == 0` as "no dipole", not as a dipole of zero.
  `CPMDCDipole` and `cpmdc_last_dipole()` carry the same `valid` flag, and Python `Potential.evaluate()` leaves out the `dipole` key.
  The dipole is a Berry-phase quantity: defined only for a neutral periodic cell, and only modulo the lattice quantum.
  The same release makes `PotentialResult.gradient` the negated `forces` in the units of the `ForceInput`.
  An evaluation of more than 1365 atoms is rejected before it runs, so no wire result carries a gradient for such a system.
- The README links to the C ABI reference, the CPMD option mapping, and the feature inventory reference instead of copying their session sample, carrier table, mapping rules, and inventory test list. The option mapping gains the manual lines of the old sample. The build-and-test how-to moved from `tutorials/quickstart.org` to `howto/quickstart.org`, its list of test suites moved to `reference/build-options.org`, and the old address keeps a page that links to the new one. The Re-check section of `docs/AUDIT.md` records gaps six and eight as closed.
- The eOn tutorial `tutorials/eon-rgpot.org` now shows the output of both builds of `libcpmdc`. The default build runs through `eonclient` in about a second and prints 39 iterations, 79 force calls, and a final energy of 10.63288405405 eV; the OpenCPMD build shows eOn's `[Matter]` table with the force norm of each iteration, the SCF iterations of each force call, and its `results.dat`. The page states which tolerances to expect, and `tests/test_tutorials_run.py` runs the page's `pos.con` and `config.ini` through `eonclient` and compares the numbers, or skips that check when no `eonclient` is on `PATH`.
- The embedding how-to, `docs/orgmode/howto/embedding.org`, is now organised by task: pick the entry point, check what the build supports, send the message, step a session, handle status and errors, keep orbitals between steps, choose a `CPMDParams` carrier, keep the wavefunction from a RESTART file, adopt an MPI communicator, set the OpenCPMD runtime inputs, drive cpmdc through rgpot or a socket, and check forces. It links to `docs/orgmode/reference/c-abi.org` for the call contracts, argument rules, return codes, and size rules that it used to repeat. The C ABI reference gained the sources a `CPMDParams` buffer can come from, the `NULL` and zero-size cases and the `nwchemc` parity of `cpmdc_session_calculate_result()`, the `forces_len` failure of `cpmdc_session_calculate_forces()`, the `cpmdc_last_energy_components()` note, and the borrowed communicator of `cpmdc_adopt_calculator_comm()`. The call routing of `cpmdc-ipi` moved into `docs/orgmode/howto/socket-drivers.org`, and the host helpers of the embed moved into the layer table of `docs/orgmode/explanation/architecture.org`. The README now sends call contracts to the C ABI reference, `tests/test_docs_examples.py` follows the new order of the code blocks, and the Re-check section of `docs/AUDIT.md` records the embedding item as closed.
- The i-PI documentation now uses rgpot's generic `rgpot-ipi-driver` CPMD route
  as the primary path and keeps `cpmdc-ipi` and `cpmdc-ipi-c` as compatibility
  fallbacks.
- `CPMDCCalculator` reads `Atoms.initial_charges` and `Atoms.initial_magmoms` instead of ignoring them. The total charge is the sum of the charges and the multiplicity is `|M| + 1` for the sum `M` of the collinear moments; each sum must be an integer within the float64 rounding bound of its terms, and non-collinear moments raise `ValueError`. A total that `CPMDParams` leaves unset (charge 0 with `hasCharge` false, multiplicity 1) is applied to the session, one that equals the message value passes silently, and one that differs from a set value, including charge 0 with `hasCharge` true, raises `ValueError` that names both values and both sources. The set value is kept. An `Atoms` with zero totals behaves as before. A change of `initial_charges` or `initial_magmoms` now evaluates the structure again, and `ignored_changes` is empty. See `howto/ase.org`.
- `Session.energy_forces`, `Session.calculate_au`, `Potential.evaluate`, `CPMDCCalculator`, and the LAMMPS adapter of `cpmdc.lammps` spend less time in Python per step. Positions, atomic numbers, and the `ForceInput` message are packed with `array` and `memoryview` instead of per element checks and argument unpacking, forces come back as rows from one `memoryview.tolist()`, and `CPMDCCalculator` converts each result once and compares compositions as sorted bytes. On the mock library of `bench/client_overhead.py` one `Session.energy_forces` step with a NumPy `(N, 3)` array takes 81 microseconds for 512 atoms where it took 1100, and `CPMDCCalculator` takes 565 where it took 1649; the tables are in `bench/results/client_overhead.org`. Results, units, and errors are unchanged: a coordinate or atomic number that is not a plain `float`, plain `int`, or C-contiguous `float64` or `int32` buffer still goes through the checked loops and keeps its message. `bench/client_overhead.py` measures the client side of an evaluation per call and compares two sets of results.
- `cpmdc.lammps.FixExternal` now keeps its existing API while delegating LAMMPS callback, unit, cell, virial, MPI and error handling to the generic `rgpot.lammps.FixExternal` frontend.
- `docs/orgmode/tutorials/quickstart.org` drops its stray stub-suite prose, sets `:nosearch:`, and is excluded from `sitemap.xml` and `llms-full.txt` so the moved URL redirect stub no longer duplicates `docs/orgmode/howto/quickstart.org` in search, sitemaps, or the LLM export. `docs/orgmode/howto/quickstart.org`, `docs/orgmode/reference/build-options.org`, `docs/orgmode/tutorials/ase-energy.org`, `docs/orgmode/tutorials/first-energy.org`, and `docs/orgmode/index.org` now match the `meson.build` toolchain and test suites (the `cpmdc` suite is the project-name suite that Meson gives every registered test, so it selects the whole default build and not only the geometry and warm-selector checks), the `examples/host_step.c` output lines, the `water_params_bin` target prerequisite and origin-shifted geometry of `examples/python/ase_energy.py`, the `CPMDC_WARM_RMSD` orbital reuse rule, and one-line page descriptions in `Start Here`.

### Fixed

- ASE and GROMACS adapter ergonomics were tightened for drop-in use:

  - `CPMDCCalculator.todict()` now keeps in-memory `params` bytes by emitting a
    JSON-safe `{"cpmdparams_base64": ...}` mapping, and the calculator accepts
    that mapping in `params` for a lossless rebuild path.
  - `cpmdc.gromacs.qmi_input()` and `write_qmi()` now validate `params` and output
    paths up front, with typed errors for binary or empty path inputs.
- An older `libcpmdc` without `cpmdc_params_method` still loads in Python;
  calling `params_method()` against it raises `CPMDCError` reporting the
  missing symbol.
- An older `libcpmdc` without `cpmdc_session_last_dipole` still loads at the
  same ABI. Energy and forces are unchanged, and the dipole stays unavailable.
  Asking ASE for the dipole names the missing symbol.
- Both LAMMPS routes now accept positive, unique atom IDs with gaps, reopen their session after atom deletion, and report multi-rank input errors without deadlocking.
- Documentation audit against the code. The `cpmdc_session_calculate_au()` header comment now states the session it needs, buffer ownership, collective use with OpenCPMD, the stub result, the centre of mass recentring, the topology fixed on the first evaluation call that passes its argument checks, the repeat match, and every `-1` cause with what is left unwritten. For `cell_bohr`, a NULL cell and an all-zero cell are read the same way by the deck, but the call keeps an all-zero cell as a step that has a cell, so the two are different steps for the repeat match and for stored OpenCPMD orbitals. The other public header comments gain the missing return codes, side effects, and ranges: `cpmdc_set_params()`, `cpmdc_configure()` (with a `cpmd` arm and a `common` overlay, `xcFunctionals`, `planewaveCutoffEv`, `charge`, and `spinMultiplicity` merge, each only while the arm still has its schema default, which is `BLYP`, 70 Ry, charge 0, and multiplicity 1; `kMesh`, `scfMaxIterations`, `scfEnergyToleranceEv`, `vanDerWaalsMethod`, and `relativityMethod` are rejected with -1, and so is an unmapped `xcFunctionals` token while the functional is still `BLYP`; `basisSet`, `smearing`, and `vanDerWaalsS6` are not applied), `cpmdc_energy_gradient()` (a step without a cell skips full setup only when the latched basis matches and the latched state also has no cell, since `embed_cell_moves_in_place` returns false when `has_cell` is 0; it then starts from the stored orbitals within `CPMDC_WARM_RMSD`, default 0.15 Angstrom, or from the deck's guess without full setup beyond that limit, and any other call runs full setup for a cold SCF), `cpmdc_adopt_calculator_comm()`, `cpmdc_adopted_comm()`, the session `set_params` and `configure` calls and `create_from_config`, `cpmdc_session_calculate_forces()`, `cpmdc_session_calculate_result()`, `cpmdc_capabilities_result()`, `cpmdc_version()` (including the `cpmdc-stub/` prefix of the link stub), `cpmdc_finalize()` (in the reference and OpenCPMD builds a session evaluation that repeats the step the session stored last still returns 0 from the session's cache after the call, while the no-session evaluations, `cpmdc_set_params()`, `cpmdc_configure()`, and any other step fail with `CPMD embed not available`; the link stub does nothing and its evaluations fail with `OpenCPMD embed not available in cpmdc stub` and the function name), `cpmdc_last_error()`, `cpmdc_restart_write_path()`, and `cpmdc_restart_expand_species_symbols()`.

  `cpmdc_bind_calculator()` and the MPI how-to now describe how ranks are divided. `ranks_per_calc` is the size of every group, and zero or less means the world size. The world size must be a multiple of it, and there is no remainder group. A value that does not divide the world size, or that is above it, returns `-1` on every rank and splits nothing. Group `g` holds world ranks `g * ranks_per_calc` to `(g + 1) * ranks_per_calc - 1`, and a later call does not split again, also after `cpmdc_adopt_calculator_comm()`: it returns `world_rank / ranks_per_calc` for its own value, or `-1` on every rank when that value does not divide the world size or is above it, and the split already in use stays as it was. `cpmdc` calls `MPI_Init` only in OpenCPMD builds, and it does so before it looks at the value. The reference build and the link stub return `-1` from `cpmdc_bind_calculator()` without calling MPI, and `cpmdc_adopt_calculator_comm()` returns `-1` there until the host has started MPI. The `cpmdc_session_set_params()` comment lists allocation failure among its errors; it returns `-1` with the text `out of memory`.

  The `cpmdc_pwscale.h` comments give the limits behind `CPMDC_PWSCALE_ERR_ARG`, state that `CPMDC_PWSCALE_ERR_NOFFTW` is returned by `bench/pwscale_fft.c` when built without FFTW or when thread initialisation or plan creation fails, and say that the rank timer ends after the rank-ordered sum. The pwscale rank how-to now says that the library timer of the bound call starts at the entry of that call, so the binding made by `cpmdc_pwscale_eval_mpi()` is not timed. The C ABI reference, environment reference, wavefunction-state explanation, install how-to, Python how-to, ASE how-to (`constant_cutoff`), socket how-to (`cpmdc-ipi` status 2), and quickstart (every Meson test suite from `meson.build`, including the stub suite with `adopt-comm-stub` and `bind-calculator-stub`) match the code.

  New tests cover the claims that had none: `test_session_calculate_au` checks the NULL and all-zero cell repeat match and, on an OpenCPMD build, the warm start source for each cell kind, and `bind-calculator-stub` and `bind-calculator-shell` run `cpmdc_bind_calculator()` on four ranks.
- Evaluations above the legacy snapshot capacity now fail before evaluation.
  They no longer return a silently truncated gradient. The ABI-compatible
  `hessian[4096]` field currently limits an evaluation to 1365 atoms.
  On a session that already holds another system, the call reports the topology change instead of the limit.
- Finalization now releases the reference orbital ledger and resets all of its selection and limit state before reinitialization.
- Installed libraries no longer include environment-controlled test fault hooks when tests are enabled.
- LAMMPS tests skip the `rgpot.lammps` route only when that optional module is absent. A present but broken module now fails collection, and pair-style plugin tests remain runnable without the Python frontend. That module is not in a released rgpot wheel; CI builds it from the pinned commit that adds the frontend. The `lammps` extra installs the LAMMPS Python module and is separate from the `test` extra.
- Meson tests now use a default setup that pins `OPENBLAS_NUM_THREADS`,
  `OMP_NUM_THREADS`, and `MKL_NUM_THREADS` to 1 so parallel `meson test`
  workers do not oversubscribe OpenBLAS. OpenCPMD live targets
  `extended-surface`, `session-cached-geometry`, and `configure-overlay`
  get 300 s timeouts when `with_cpmd` is enabled. Meson replaces those
  three variables on test children; parent-shell exports do not apply (see
  `docs/orgmode/reference/build-options.org`).
- Non-MPI builds no longer define an unused orbital-key mismatch message, and environment documentation no longer promises an unreachable peer-rejection diagnostic.
- Orbital copies now agree OpenCPMD availability across calculator ranks, MPI state-check failures are reported, and cell-in-place documentation identifies rank 0 as authoritative.
- Pin `potentials-schema` `1.16.0`, sync the vendored `schema/Potentials.capnp` and feature inventory, and honor `CPMDSystemSection.hasCharge`, `CPMDParams.hasCharge`, and `CommonMethodSpec.hasCharge` across effective config, Fortran embed apply, `&SYSTEM` deck rendering, `Capabilities.loweredCommonFields`, and `PotentialConfig` overlay lowering and conflict checks so an explicit zero charge is distinguished from an unset charge.
- Prepared orbital-slot moves are now discarded on every rank when any calculator rank cannot prepare the selection.
- Rank-local embedded evaluation failures now propagate through result serialization without leaving successful MPI ranks blocked in a collective.
- Rendering a deck with no overlay now rejects a `CUTOFF`, `CHARGE`, or `MULTIPLICITY` directive in the typed `system` section, and a `SYSTEM.CUTOFF`, `SYSTEM.CHARGE`, or `SYSTEM.MULTIPLICITY` set directive. Those lines are the resolved scalars the deck and the recorded configuration share. `CUTOFF` is written on every deck, including the default 70 Ry. `CHARGE` is written only when it is not 0 and `MULTIPLICITY` only when it is greater than 1; an omitted line is the schema default the recorded configuration carries. The directives are rejected for a default value and for any other value alike. A default leaves the directive as a value the recorded configuration does not carry, and any other value leaves two lines. The render names the directive and does not keep one line and drop the other.
- Session evaluations now validate topology read-only and commit `topology_fixed` only after all orbital preflight collectives succeed, leaving the session unchanged on any preflight failure. Keyed state, stale-orbital cleanup, and orbital preparation outcomes are agreed across all calculator ranks on both the reference and OpenCPMD paths so one rank's divergence or failure fails every rank with the same error instead of hanging.
- The README, the start-here index, the route overview, and the C surface map now describe `rgpot-ipi-driver` selecting `CPMDPot` in an rgpot build with the RPC stack (`--pot "CPMD params FILE"` or `-m CPMD -o params=FILE`, with `--preserve-open-cell`). `cpmdc-ipi` and `cpmdc-ipi-c` remain compatibility fallbacks. Wrapped Org code spans in the docs tree stay on one line so the plain-text export does not leave a literal tilde.
- The deck conflict check ignores a CPMD comment, so a commented
  `DIPOLE DYNAMICS` line is not treated as a request.
- The eOn tutorial check in `tests/test_tutorials_run.py` now compares the final energy within an absolute 1e-6 eV for the default build and 1e-4 eV for the OpenCPMD build, as `tutorials/eon-rgpot.org` states. It used to scale the tolerance by the size of the expected energy, which let errors of about 1e-5 eV and 0.05 eV pass.
- The installed `cpmdc.pc`, `cpmdc-restart.pc`, and `libcp2k.pc` files now derive `includedir` and `libdir` from `${pcfiledir}` before setting `prefix`, so consumers compile and link after an installation prefix is moved under both plain `pkg-config` and `pkgconf --define-prefix`. Installed shared libraries and plugins record an `$ORIGIN` install run path (`@loader_path` on macOS), or an empty one when they do not need another library from the prefix, so `libcpmdc_gromacs` still opens `libcpmdc` after that move. `cpmdc-installed-package-consumer` moves the installed prefix, rejects a `.so` or `.dylib` whose `RUNPATH` or `RPATH` points at the build tree or the original prefix, and runs the `cpmdc`, `cpmdc-restart`, `libcp2k`, `cpmdc-ipi-c`, and CMake consumer checks against the relocated prefix under both `pkg-config` modes.
- The pixi environment now sets `PKG_CONFIG` to its own `pkgconf`, so a host `pkg-config` and a host libcmocka can no longer pair with the conda `cmocka.h` and trip the configure-time cmocka probe in the readcon restart job.
- The pytest suite passes on ASE 3.23, the declared floor, through the current ASE 3 release. ASE 3.23 names the finite difference helpers `numeric_forces` and `numeric_stress` and calls the step `d`. ASE 3.24 moved them to `ase.calculators.fd` and renamed the step `eps`. ASE 3.26 added `force_consistent`. ASE 3.24 and 3.25 always differentiate `free_energy` and reject that keyword, which failed the six GROMACS stress comparisons. The shear formula is unchanged: off diagonal strain `d` divided by `4 d V` is the same tensor as strain `0.5 eps` divided by `2 eps V`. The GROMACS shim reports `free_energy` equal to `energy`, so that older derivative is the energy derivative. ASE 3.23 still stores `mu0 = 4 pi 1e-7` in its 2018 table, and ASE 3.24 stores the CODATA 2018 permeability. The calculator keeps cpmdc's own constants. The session treats geometries within 1e-8 Angstrom as identical, so the ASE gradient check stops at 1e-7 Angstrom. ASE 3.23 updates the displaced coordinate in place and lands on that identity, while ASE 3.24 rebuilds the coordinate and can slip past it. CI runs pytest on ASE 3.23.0 and on the newest ASE 3 release.
- The reference evaluator now computes periodic stress using the scaled LU cell volume so accepted near-singular cells receive a valid stress tensor instead of cancelling to zero, and fails explicitly if a non-zero periodic cell cannot represent a positive finite stress volume. Generated diagonal `CELL` lines stay within 80 columns, `cpmdc-ipi` enforces the ASCII integer grammar on `INIT` `batch_size`, and the C ABI and contract references document the cell validation and failure modes.
- Two `DIPOLE DYNAMICS` evaluations in one process no longer fault in `opeigr_para` when the second uses a larger basis, such as 40 Ry followed by 50 Ry.
  `localize2` kept the orbital maps and dipole matrices of the first basis; the OpenCPMD teardown patch now resets them like the other routines that keep first-call state.
  Rebuild `libcpmd.a` with the updated patch; `tools/rebuild_opencpmd_embed.sh` recompiles `localize_utils`.
- `CPMDCCalculator.read(label)` raises ASE's `ReadError` instead of setting the label and returning as if results had been restored. A step stores no file that holds results, so there is nothing to read; `parameters=` still loads a file written by `Parameters.write`.
- `CPMDCCalculator` keeps an explicit `CPMDParams` charge of 0 (`hasCharge` true) when `Atoms.initial_charges` sum to another integer. The sum is written only when `hasCharge` is false and the effective charge is still 0. A disagreeing total raises `ValueError` that names both values and both sources, the same refusal as any other set charge, and the step does not run. A typed `system` charge with `hasCharge` true is kept the same way, because the edited message is read back before the session opens. Initial charges that are all zero still leave the message charge unchanged. The embed records that message charge. It does not read `initial_charges`.
- `CPMDCCalculator` refuses a dipole that is not finite, the same way it
  refuses a non-finite energy, force, or stress.
- `CPMDCCalculator`, `cpmdc-ipi`, `cpmdc-ipi-c`, and native deck rendering now use the same float64 rank cutoff for periodic cells and preserve every representable cell component in generated CPMD input; the C library and the Python package decide the rank with the same iteration, so they agree on cells within rounding of the cutoff. Cells whose volume would underflow, duplicate or invalid i-PI `batch_size` values, non-UTF-8 `INIT` data, and negative or non-finite connection waits now report their cause instead of being treated as an ordinary request. Both i-PI drivers judge the cell volume on the Bohr values they receive, `cpmdc-ipi-c` stops on a non-finite position or cell value as `cpmdc-ipi` does, and the virial of `cpmdc_session_calculate_au()` and of both drivers no longer cancels to zero for an accepted near-singular cell; a cell volume above the float64 range or a non-finite virial is reported instead of being sent, and session cell, volume, and virial checks agree rank-local failures across MPI calculator ranks before returning.
- `cpmdc-ipi-c` records a normalised `$ORIGIN` run path from the resolved `bindir` to the resolved `libdir` (`@loader_path` on macOS). `-Dbindir=.`, absolute `bindir` and `libdir` under the prefix, a nested `libdir`, and paths that contain `..` or a trailing slash install a run path that opens `libcpmdc`. A layout whose relative path cannot be stored as one run path entry fails at configure time.
- `cpmdc.pc`, `cpmdc-restart.pc`, and `libcp2k.pc` now derive `prefix` and `includedir` from `${pcfiledir}` with the relative path from the pkg-config directory, so `-Dlibdir=.`, a nested `libdir`, and an absolute `includedir` under the prefix stay inside that prefix.
- `libcpmdc_gromacs` again ignores multi-value `&SECTION` headers outside
  `&CPMDC` on explicit `&CPMDC` inputs while keeping `&XC_FUNCTIONAL` arity
  checks for GROMACS-generated `&DFT` inputs.
- `libcpmdc_gromacs` refuses a `-qmi` input whose first `CHARGE` or
  `MULTIPLICITY` anywhere in the file differs from the charge or
  multiplicity of the `CPMDParams` message in `&CPMDC PARAMS`, matching
  GROMACS' own token-stream reader just as it does for the input that
  GROMACS generates. The explicit route ran the message's system and
  ignored the values that GROMACS used for the QM charge. Later
  `CHARGE` or `MULTIPLICITY` lines do not override the first occurrence,
  and a line without a single value or a non-integer first value is
  refused too. `qmi_input()` and `write_qmi()` take the deck's `charge`
  and `multiplicity` for a deck that is not neutral and singlet. A link
  hydrogen is not counted into either value.
- `test_extended_surface` now calls `cpmdc_finalize()` only in its final shutdown regression test rather than mid-suite, so `test_topology_change_is_reported_before_the_gradient_limit` evaluates on a live embed runtime and post-finalize checks confirm that failed embed preflight leaves open topologies uncommitted while preserving topology checks on already-fixed sessions.

### Developer

- AddressSanitizer leak detection stays on for serial C and Fortran tests
  that do not preload the sanitizer runtime or run under
  `tests/run_mpi_test.py`. `params-cp-dft-render` keeps its fixture paths
  in `main` instead of leaving a heap copy for the process exit, and
  `dlopen-leak` verifies that `dlopen` constructor leaks are reported.
  `meson.build` sets `detect_leaks=0` for the four ASan preload tests
  (`python-api`, `python-pytest`, `docs-examples`, and `tutorials-run`)
  because CPython and compilers started from those tests do not free
  their pools, and `tests/run_mpi_test.py` sets `detect_leaks=0` for the
  MPI tests (`adopt-comm-stub`, `adopt-comm-shell`,
  `orbital-limit-agree-mpi`, and `rank-consensus`, which is registered
  only for the reference build) and the `pwscale` tests
  (`pwscale-agreement`, `pwscale-logger`, `pwscale-fft-fail`,
  `pwscale-agreement-mpi-*`, `pwscale-mpi-fault`,
  `pwscale-collectives-*`). The pixi `test-stub-sanitize` task and the
  GCC 13 CI job use that same split.

  The `shared-dlopen-abi` test now pins `libibverbs.so.1` with
  `RTLD_NOLOAD|RTLD_NODELETE` before `dlclose`. `RTLD_NODELETE` keeps
  `libibverbs` mapped, so LeakSanitizer sees the allocation as reachable
  without a suppression file.
- The `opencpmd-embed` CI job now runs `tests/python/test_session.py` against `cpmd.x`, and fails on any skip other than the three that name the reference evaluator's stress formula. The test `test_session_preserves_an_accepted_small_triclinic_cell` had only run on the reference evaluator. Its cell has a third vector 1e-13 Angstrom long, so OpenCPMD holds 8 electrons in about 1e-13 cubic bohr, a plane wave basis of 5 waves for 4 orbitals and an electronic gradient of 1e13 from the first step. No SCF converges there, and the lattice vectors, the row convention and the mesh that OpenCPMD printed for it were correct: rotating the cell and the atoms together leaves the energy unchanged to 6e-14 hartree, and transposing the cell moves it by 6e-4 hartree. The test now puts the tilt components of 2.58e-12 and 1e-13 Angstrom beside 6 Angstrom edges and requires the energy, forces and stress of the untilted cell on both backends. The cell with the 1e-13 Angstrom vector stays as its own test, which expects the reference evaluator's stress and, on OpenCPMD, an error from CPMD that is not the rank guard's.
- The `rcut` and `alpha` assertions of `test_mechanical_embedding_matches_numpy_ewald`
  now check `1e-14` relative with no absolute slack, where `rel=1e-15` had been an
  absolute `1e-12` in effect, and the other tolerances of that test carry their
  error bounds.
- The `sdist` workflow now installs the `test` extra (pytest and ASE) into the installed-wheel check venv, matching cibuildwheel's `test-extras`.
- The installed-wheel test changes the working directory before loading the bundled shared library.
- The pytest session now checks NumPy's matrix product against a plain loop
  before collecting any test, and stops with a message naming the remedies when
  the two disagree. On aarch64 under `qemu-user` 8.2 with `-cpu max`, OpenBLAS
  selects its `armv9sme` kernel and `@` returns wrong values, which made the
  NumPy Ewald reference of `test_mechanical_embedding_matches_numpy_ewald` sum
  to zero and fail on its energy assertion while the library was correct.
- The sdist job now sets `CC`, `CXX` and `FC` to the GCC toolchain it installs, so Meson and its cmake module no longer pick a clang `c++` that cannot link.
- The two-rank LAMMPS cases in `tests/python/test_lammps.py` start their ranks with `OMPI_MCA_pml=ob1` and `OMPI_MCA_btl=self,sm` unless the caller sets them. On GitHub runners with the Azure MANA NIC, UCX finds `mana_0` but cannot open it, and the Open MPI ranks of the `lammps` Pixi environment hung until `test_two_ranks_match_one` and `test_two_rank_duplicate_id_error_does_not_deadlock` timed out. The cases run on one host and need only shared memory.
- `cpmdc.lammps.FixExternal` evaluates each LAMMPS step through the session adapter of the `rgpot.lammps` frontend, and no longer carries a second copy of that step. The convention tests build the fix with its real constructor against a stand-in LAMMPS object, so they need `rgpot.lammps` and the `opencpmd-embed` job installs it. `mypy --strict` passes whether or not `rgpot.lammps` is installed.
- `scripts/generate_llms_txt.py` now sorts Org inputs by repository-relative path, rejects a missing feature inventory explicitly, omits absolute paths from diagnostics, and reports stale `llms.txt` or `llms-full.txt` files with the exact regeneration command and a unified diff head; `tests/test_llms_txt.py` verifies repeatability across runs and checks the failure output.
- `tests/python/test_ase_conformance.py` also checks `CPMDCCalculator` against `SumCalculator` and `VelocityVerlet`, the parameters round trip through `todict` and a `Parameters` file with a `constant_cutoff` mapping, and that a change of tags, momenta, info, or a constraint, or an equal copy of the `Atoms`, reuses the stored step.


## [0.6.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.6.0) - 2026-11-22

### Changed

- The `lychee` Pixi task now runs in the `docs` environment with conda-forge lychee 0.24.2, matching the documentation workflow. Run it with `pixi run -e docs lychee`.
- The i-PI documentation now treats protocol drivers as optional rgpot frontends, keeps direct `Session.calculate_au()` calls and rgpot's Cap'n Proto RPC as the default routes, and records the rgpot command, open-cell, init-less handshake, and batched-`INIT` (`batch_size > 1`) gaps that currently keep the transitional cpmdc drivers in place. Element-list parsing used by the Python i-PI and LAMMPS adapters now lives in a backend-neutral internal module.
- `pair_style cpmdc` now refuses a run in which a compute asks for a per-atom energy or virial (`pe/atom`, `stress/atom`, `centroid/stress/atom`) in `init_style`, before the first evaluation, and the error names the compute. cpmdc returns a total energy and a cell derivative virial, neither of which has a unique per-atom split, so the style gives none; the LAMMPS interface reference explains why and lists the refused computes, and tests pin the refusal.

### Fixed

- A session evaluation and an orbital copy check their arguments and session before they agree `CPMDC_ORBITAL_SLOT_LIMIT` across the calculator ranks, so a null session or a bad argument is reported as that and not as the limit. A rank whose own check fails still joins the agreement and every other rank fails the call too, instead of waiting for it.
- An OpenCPMD evaluation under a session key frees the parked orbitals of an old basis before it selects its key, so a `CPMDC_ORBITAL_SLOT_LIMIT` that the old basis filled no longer fails the first evaluation on a new basis.
- Common-overlay merging now treats a field in a typed `dft` or `system` section, or a top-level `CPMDParams` field, that is omitted or holds its schema default (functional `BLYP`, cutoff 70 Ry, charge 0, multiplicity 1) as unset, preserves explicit arm values, and rejects every set overlay scalar, whether or not a typed field also sets it, next to the typed-section directive that renders the same setting: `FUNCTIONAL` in a typed `dft` section, or `CUTOFF`, `CHARGE`, or `MULTIPLICITY` in a typed `system` section. `cpmdc_params_render_input_deck_ov()`, `cpmdc_params_render_deck_with_geometry_ov()`, `cpmdc_params_effective_config_ov()`, and the configuration the embed records apply the same precedence to the functional, cutoff, charge, and multiplicity lines cpmdc writes, so a typed `dft` section naming `BLYP` next to an overlay `xcFunctionals` no longer writes `FUNCTIONAL BLYP` while recording the overlay functional. A value a typed field sets beats the overlay, the overlay beats a schema default, a typed section beats the top-level field, and the overlay cutoff is converted from eV to Ry before it is written. Without an overlay, a `FUNCTIONAL` directive in a typed `dft` section or a `DFT.FUNCTIONAL` set directive is rejected when the deck is rendered, instead of following the rendered `FUNCTIONAL` line in the same `&DFT` section.
- Fixed two heap buffer overflows in OpenCPMD's k-point atomic
  wavefunction initialization (`atomwf` and `gs_ortho_c`) when
  `nstate > atwp%nattot` via `tools/opencpmd_atomwf_kpts_ortho.patch`:
  `gs_ortho_c` now passes `np*n0` rather than `2*np*n0` to `mp_sum` on
  the `COMPLEX(real_8)` overlap matrix, and `atomwf` allocates a local
  `smat_c(nleft, MAX(atwp%nattot, nleft))` scratch array instead of
  reusing `zmatat(atwp%nattot, atwp%nattot)` while fixing the spin-up
  `c0(:, atwp%nattot+1:, ikind)` slice.
  `meson.build` raises the timeouts of `extended-surface`,
  `session-calculate-au`, `configure-overlay`, and `docs-examples` for
  sanitizer builds only.
- MPI evaluation preflights now reject rank-local errors collectively before any rank enters the orbital or energy collectives.
- Merging a `common` overlay (`CommonMethodSpec`) into a native `cpmd` arm in `cpmdc_configure()`, `cpmdc_session_create_from_config()`, `cpmdc_session_configure()`, and `cpmdc_calculate_result_from_config()` now validates `basisSet` (`"planewave"` or unset), lowers Fermi `smearing` into `&CPMD` when free-energy controls are absent, copies `vanDerWaalsS6` into a `vdwParams` section when `s6` is unset, and rejects conflicting, unmappable, or out-of-range overlay values through `cpmdc_last_error()`. Overlay-only lowering also rejects `vanDerWaalsS6` without `vanDerWaalsMethod` and non-zero `smearing.widthEv` with `kind = none`. An `xcFunctionals` list is mapped before arm-wins precedence, so an unmapped token is rejected for every arm functional. A typed `system` charge of 0 is the schema default and does not replace `CPMDParams.charge`; a non-zero system charge does. The carrier check runs for every set overlay field, whether or not a typed field also sets it: the overlay field is rejected when `inputBlocks`, a raw section, or a generic section opens the owning section (`&DFT` for `xcFunctionals`, `&SYSTEM` for the cutoff, charge, and multiplicity) or a set directive names `DFT.FUNCTIONAL`, `SYSTEM.CUTOFF`, `SYSTEM.CHARGE`, or `SYSTEM.MULTIPLICITY`. Those carriers stay in the rendered deck and are not read back into the recorded functional, cutoff, charge, or multiplicity. A `spinMultiplicity` of 0 is unset, and any other value, an explicit singlet included, is merged and checked like the other overlay scalars, so a singlet next to an `&SYSTEM` carrier is rejected. Every rejection names its cause in `cpmdc_last_error()`, and an arm `inputBlocks`, `inputSections`, or `directives` list that cannot be read is reported as `<field> cannot merge: the cpmd arm has an invalid <list> list`. Smearing conflict detection also handles split `TEMPERATURE` and `ELECTRON` arguments in `cpmd.directives` and `CPMD.TEMPERATURE` set directives.
- Meson now performs a configure-time cmocka compile+link probe that catches header/library mismatches early and aborts with a clear diagnostic instead of deferring failure to the test link step.
- On POSIX, RESTART publication now refuses unrecorded directory symlinks in the destination path instead of following them outside the checked symlink chain.
- On POSIX, `cpmdc_restart_write_path` probes the final component of the destination with `fstatat` and `AT_SYMLINK_NOFOLLOW` on the already opened no-follow directory, instead of `lstat` on the whole path, so no stat walks through a directory symlink that the recorded `RESTART` chain does not hold. A path through such a link now fails with "cannot open the RESTART directory".
- On POSIX, publishing a `RESTART` file to an absent final name uses `renameat2` with `RENAME_NOREPLACE` only where the headers and kernel provide it. If that operation is unavailable, plain `renameat` leaves a gap in which a newly appeared name can be overwritten. Replacing an existing file has the same gap because POSIX has no rename that succeeds only when the directory entry is the same inode.
- Parked orbital copies now stop at CPMDC_ORBITAL_SLOT_LIMIT. A further park fails the call and leaves the copies already stored in place.
- The Python API reference page `docs/orgmode/reference/python.org` lists the full signatures, parameter names, defaults, units, library and schema search orders, and raised exceptions of `cpmdc`, and `tests/python/test_reference_python.py` checks the page against `cpmdc.__all__` and the inspected signatures of `Session`, `Potential`, `Evaluation`, and `ConstantCutoff`.
- The table of parked orbital keys in the reference evaluator grows only when it is full, and not on every new key, so two keys that alternate no longer exhaust memory. A growth that cannot be sized or allocated fails the call on every rank and leaves the table unchanged.
- With OpenCPMD, a step whose box is triclinic, so that the deck gets `CELL VECTORS`, used to keep the minimum-image half lengths of an earlier orthorhombic box set up in the same process with `SYMMETRY 8`. OpenCPMD's `setsc` writes those lengths only for a `CELL` deck, so `pbc3` folded the triclinic cell as a monoclinic one of the old lengths instead of through its inverse cell matrix, and `CPMDC_CELL_IN_PLACE` refused every strain of that cell and set CPMD up again. `cpmdc_embed_teardown` now clears the lengths, as at the first setup of a process.
- `CPMDCCalculator` stores a relative `params` or `library` path as the absolute path it was read from. `todict()` used to return the path as given, so a calculator rebuilt from a trajectory, an `ase.db` row, or `CPMDCCalculator(**calc.todict())` in another directory either failed to find the file or evaluated another file of the same name without an error, and a copy that reopened its session there lost a relative `library`. `tests/python/test_ase_conformance.py` now also checks the calculator against ASE's own calculator tests: the cache and `set()`, `get_property(..., allow_calculation=False)`, the copy of the caller's `Atoms`, extended XYZ through `SinglePointCalculator`, a central difference of the energy over strain for every stress component, and a `FrechetCellFilter` relaxation.
- `CPMDCCalculator` stores a stress only for a step that sent a periodic cell. With OpenCPMD and a periodic CPMDParams deck, a structure with every `pbc` flag false used to get the stress of the deck's own cell, so `get_stress()` after `get_potential_energy()` returned a tensor where it raises `PropertyNotImplementedError` as the first request.
- `cpmdc_session_set_params` and `cpmdc_session_configure` now leave the session exactly as it was when they fail. A call that failed validation used to clear the overlay of an earlier `cpmdc_session_configure` and overwrite the effective functional, cutoff, charge, multiplicity, and `cpmd_root` of the session. The new parameters and overlay are validated into temporaries and installed only on success, and a failure while configuring the embed layer restores the previous state. Only `cpmdc_last_error()` changes.

### Developer

- A new pytest module holds the species grammar of `cpmdc._elements.species_numbers` fixed for the i-PI driver and the LAMMPS front end: the symbol table is checked against ASE's element table, and the documented forms, the `Co` against `CO` rule, zero counts, and every refusal path are covered.
- On OpenCPMD the LAMMPS harness checks the forces against a central difference of the energy at 1e-3 Angstrom steps (2.6e-5 eV/Angstrom), the plain-deck virial against strains of 1e-5 (2.5e-7 of the largest component), the virial on a basis held in place against a Richardson extrapolation of strains of 2.5e-3 and 5e-3 (4.8e-6), the constant-cutoff virial against new setups, NVE drift and the dt^2 scaling of its fluctuation, and an iso NPT run at a constant cutoff, on the periodic water deck at an orbital threshold of 1e-7. The unwrapped path runs there with forces held to 1e-9 eV/Angstrom. A plain deck with `REFERENCE CELL` gives a strain difference of +3.28 eV where the virial is -32.42 eV; that test is a strict xfail.
- Tests cover the growth of the parked key table under a low address space limit, the order of the argument checks and the limit agreement, an invalid and an empty `CPMDC_ORBITAL_SLOT_LIMIT` in `cpmdc_session_copy_orbitals()`, the reset on a new basis, and a two rank agreement. `tools/opencpmd_embed_extrapolate.patch` is regenerated against the current rwfopt patch and applies with fuzz 0. The tests use geometries that the real SCF converges from a cold start, so they also pass on the OpenCPMD build.
- The CI job `opencpmd-embed` also runs `tests/python/test_ase_conformance.py` and `tests/python/test_lammps_conventions.py` against `cpmd.x` through the `python-pytest` Meson test, and fails on a skip. Every other job runs them on the reference well, which returns no stress without a cell and an infinite energy for an absurd geometry, so five cases that fail with OpenCPMD went unnoticed.
- The LAMMPS harness, `tests/python/test_lammps.py`, collects a `cpmd_backend` test only for the library it describes. The `opencpmd` fixture and collection both read `cpmdc_available()`, so a run of either library reports no skips from the other. The reference well keeps four tests that use its quadratic energy: the unwrapped path at rtol 1e-10, NVE within 2e-3 of the fall into the well, NPT toward -20000 bar, and strains of 1e-5 against the virial. `docs/orgmode/howto/lammps.org` lists what each library collects with the measured numbers, and says how to build LAMMPS `stable_22Jul2025` and the patched OpenCPMD archive and run the harness on each.
- The LAMMPS tests, `tests/python/lammps_case.py`, and `bench/lammps_routes.py` hand every LAMMPS instance a fresh copy of the shared command-line list, because the `lammps` Python module of the 2Aug2023 stable series inserts its program name into the list it is given, so the second instance of a process started with `lammps lammps -log none` and stopped before the first command. The standalone `lammps_case.py` driver also closes its instance and then calls `lammps_kokkos_finalize` and `lammps_mpi_finalize` once, so that each rank of an `mpiexec` run leaves MPI cleanly with that module. It does not call `lammps.finalize()`, which also runs `lammps_python_finalize` and so finalizes the interpreter that hosts the module, a crash under the 22Jul2025 release. The two rank tests also check a single launched rank, and fail instead of skipping when `mpiexec -n 2` does not give LAMMPS two ranks.
- The POSIX build of `src/cpmdc_restart.c` no longer defines the unused `publish_open` and `publish_unlink` macros, and the `restart-unused-macros` meson test compiles the file with `-Werror=unused-macros` so an unexpanded platform macro fails the suite.
- The `opencpmd-embed` job also runs `tests/python/test_cold_start.py`. OpenCPMD keeps the orbitals of the last step after a `Session` closes, so the `cold_start` fixture parks them on another composition before a case that needs a step set up from its deck. The cases that scale the water positions by `1e200` now pin the `CPMD stopgm during embed SCF` reason on OpenCPMD instead of accepting any error.
- The installed-package test suite clears `CPMDC_LIBRARY`, `RGPOT_CPMDC_ENGINE`, `RGPOT_CPMD_ENGINE`, and `CPMDC_SCHEMA` before each test and checks that the Python package still finds its bundled shared library and schema after the working directory changes.


## [0.5.0](https://github.com/OmniPotentRPC/cpmdc/tree/v0.5.0) - 2026-11-19

### Removed

- `libcpmdc` exports only the 45 functions of `include/cpmdc.h` and `include/cpmdc_features.h`.
  The 56 internal `cpmdc_` glue symbols (`cpmdc_embed_*`, `cpmdc_params_*`, `cpmdc_stop_*`, and the rest) and the symbols of the static archives linked into the library (c-capnproto, capnp-fortran, minimage, OpenCPMD) are no longer in its dynamic symbol table.
  No public header declared any of them, so no public function changes, and `CPMDC_ABI_VERSION` and the soversion stay as they are.
  A program that looked one of them up with `dlsym` no longer finds it.
  `scripts/abi_baseline.txt` now lists the 45 public names.

### Deprecated

- The `schema` argument of `Session`, `Potential`, and `CPMDCCalculator`, and the `restart` argument of `CPMDCCalculator`, are ignored and warn with `DeprecationWarning`.
  cpmdc no longer needs the schema file to encode a cell and reads no restart file.
- `cpmdc.CpmdcError` is deprecated in favour of `cpmdc.CPMDCError`, which spells the acronym as `CPMDCCalculator` and the C types (`CPMDCSession`, `CPMDCResult`) do.
  Both names are the same class, so `except cpmdc.CpmdcError` still catches every error the package raises; reading the old name warns with `DeprecationWarning`, and it is no longer in `cpmdc.__all__`.
- `include/cpmdc_restart.h` names its opaque image `CPMDCRestart`, as `cpmdc.h` names `CPMDCSession` and `CPMDCResult`, and its prototypes use that name.
  `cpmdc_restart`, the 0.4.0 name, is deprecated; it stays a typedef of the same `struct cpmdc_restart`, so existing code compiles and links unchanged.

### Added

- A new reference page, `reference/contracts.org`, gives every public C function, type and constant and every public Python name its units, shape and layout, ownership, thread and MPI behaviour, return and error codes with the exception each maps to, and defaults, including the ASE property and parameter tables and the deprecated names.
- GROMACS can integrate on cpmdc forces in process, for a QM group or the whole system, through its CP2K QM/MM interface.
  `libcpmdc_gromacs`, built by Meson from `tools/gromacs/` (`-Dwith_gromacs_shim`), implements the libcp2k C API on `cpmdc_session_calculate_au()`, and `build/tools/gromacs/pkgconfig/libcp2k.pc` lets a double-precision GROMACS 2024 to 2026 configured with `-DGMX_CP2K=ON -DCP2K_DIR=...` link it in place of CP2K; the energy appears as `Quantum En.`.
  The cpmdc settings sit in a `&CPMDC` section of the grompp `-qmi` input, which `cpmdc.gromacs.qmi_input()` writes: the `CPMDParams` path, whether the GROMACS box is sent as the cell, and mechanical embedding, an Ewald (or, without a box, direct) QM-MM Coulomb sum of fixed QM charges with the MM charges.
  Electrostatic embedding and link atoms are not provided, and inputs that need them are refused with a `cpmdc-gromacs:` message.
  Run GROMACS on cpmdc forces (`docs/orgmode/howto/gromacs.org`) explains the choice of route over MiMiC, the NNP interface, PLUMED, and i-PI; `tests/python/test_gromacs.py` checks one-step energies and forces against `CPMDCCalculator` to 1e-10 and NVE conservation with GROMACS 2024.6, 2025.5, and 2026.4, the new `gromacs-qmmm` CI job builds GROMACS for it, and `bench/gromacs_routes.py` measures the per-step overhead beside `pair_style cpmdc`.
  The C ABI is unchanged.
- GROMACS can run NPT with cpmdc as the QM engine.
  `libcpmdc_gromacs` exports `cp2k_get_stress_tensor()`, the strain derivative of the energy it returned: cpmdc's stress with `CELL PERIODIC`, the symmetric part of the QM atoms' `sum x (x) f` with `CELL NONE`, plus the strain derivative of the mechanical-embedding Ewald sum, in cubic and triclinic boxes.
  GROMACS' QM/MM module never adds a virial, so `tools/gromacs/patches/qmmm-virial-2024.patch` (GROMACS 2024.6 and 2025.5) and `qmmm-virial-2026.patch` (2026.4) make it add `-0.5 sigma V` in kJ/mol from the main rank; `scripts/ci/build_gromacs_cpmdc.sh` applies them.
  The library now reads GROMACS 2026's box, sent one lattice vector per row, as well as the transposed box of 2024 and 2025.
  `tests/python/test_gromacs.py` checks the stress against ASE's finite-difference stress, GROMACS' virial and pressure against the stress times the volume and against a strain difference of GROMACS' potential energy, and the virial across domain decompositions; `tests/python/test_gromacs_opencpmd.py` checks the OpenCPMD stress the same way, and `bench/gromacs_opencpmd_npt.py` runs C-rescale or Parrinello-Rahman NPT with an OpenCPMD QM water beside a classical control.
  Run GROMACS on cpmdc forces (`docs/orgmode/howto/gromacs.org`) describes the settings.
  The C ABI is unchanged.
- LAMMPS can integrate on cpmdc forces in process, through two routes on `cpmdc_session_calculate_au()`.
  `cpmdc.lammps.FixExternal` attaches one session to a `fix external pf/callback` of a `lammps` Python object, and `cpmdcplugin.so`, built by Meson from `tools/lammps/` when the LAMMPS headers are found (`-Dwith_lammps_plugin`, `-Dlammps_include`), adds `pair_style cpmdc` to any LAMMPS with the `PLUGIN` package.
  Both convert `metal`, `real`, and `electron` units and refuse other styles, send the restricted triclinic box as a row-major cell with positions unwrapped by image flags, and return the energy, forces, and the virial `-dE/d(strain)` so `fix npt` sees the cpmdc pressure; under `mpirun`, world rank 0 owns the session and gathers the atoms by ID.
  `Session.calculate_au()` exposes the C call to Python, and `HARTREE_TO_KCAL_PER_MOL` joins the unit constants.
  Run LAMMPS on cpmdc forces (`docs/orgmode/howto/lammps.org`) covers both routes, units, the cell and virial conventions, and the MPI rule; `tests/python/test_lammps.py` checks them against `CPMDCCalculator` with the lammps pip wheel, and `bench/lammps_routes.py` measures their per-step overhead against the i-PI socket routes.
  The C ABI is unchanged.
- The environment reference lists every variable the library, the tests, the benchmarks, and the docs export read, including `CPMDC_PARAMS` of the ASE example and the `CPMDC_PWSCALE_MPI_FAULT` test hook.
- The plane-wave scaffold has a reference page that gives the units, array layout, ownership, return codes, and MPI and CUDA availability of every `cpmdc_pwscale_*` call, type, and constant.
- Variable-cell relaxation and NPT get a stress that follows the energy: `CPMDCCalculator(constant_cutoff=...)`, `Session(constant_cutoff=...)`, `Potential(constant_cutoff=...)`, the frozen `cpmdc.ConstantCutoff` dataclass, the C call `cpmdc_session_set_constant_cutoff()`, and a `[constant_cutoff]` table in the job TOML of `tools/job_toml.py` add OpenCPMD's `CONSTANT CUTOFF` to the `&SYSTEM` section of the deck, the modified kinetic functional of Bernasconi et al., J. Phys. Chem. Solids 56, 501 (1995).
  E0 is required; A defaults to 50 times the deck cutoff and sigma to E0 / 10, and E0 + 2 sigma must not exceed the deck cutoff.
  LAMMPS, GROMACS, and the socket drivers read the same setting from `system.constantCutoff = [A, sigma, E0]` in the message.
  The default is off, so existing messages render the same deck; `libcpmdc` now exports 46 functions, and the addition is backward compatible, so `CPMDC_ABI_VERSION` and the soversion stay the same.
  The optional `tools/opencpmd_constant_cutoff_precon.patch` adds the penalty to the ODIIS preconditioner, which otherwise stops at its step limit and falls back to PCG; `scripts/ci/build_opencpmd.sh` applies it.
  `bench/constant_cutoff_scan.py` measures the stress, the relaxed volume, and the time per stress call of the water deck against the cutoff.
- `CPMDCCalculator` takes every keyword through ASE's `Calculator.__init__`, as the tblite, xtb-python, and MACE calculators do, so `parameters=` reads a file written by `Parameters.write` and `atoms=` attaches the calculator.
  A step that the library rejects raises `cpmdc.ase.CalculationFailed`, which is both ASE's `CalculationFailed` and `CPMDCError`, with the working directory of the step and the library's reason in the message; the calculator stays usable.
  Use cpmdc as an ASE calculator (`docs/orgmode/howto/ase.org`) covers a single point, stress in a periodic cell, relaxing the positions and the cell, finite-difference checks, caching, the method settings, the working directory, and the errors.
  Its examples run in the `docs-examples` test, and the Python how-to, the Python reference, the ASE tutorial, the docs index, and the README point at it.
- `CPMDC_CELL_IN_PLACE` moves the cell of a periodic `CONSTANT CUTOFF` deck in place: a cell within that principal strain of the cell CPMD was set up at keeps the setup's plane waves, the stored orbitals, and the `EXTRAPOLATE WFN` history, so variable-cell dynamics and cell relaxations stay warm between setups and their energies lie on one surface whose derivative is the stress.
  The move also needs the strained E0 + 2 sigma sphere to fit inside the cutoff, the Gamma point, norm-conserving pseudopotentials, and the lattice class of the setup; any other cell, and every cell when the variable is unset, sets CPMD up again as before.
- `bench/constant_cutoff_npt.py` runs hydrostatic NPT dynamics of liquid water through `CPMDCCalculator` and ASE's `LangevinBAOAB`, with and without the constant-cutoff functional, evaluates the atomic and molecular virial pressure of one configuration at several bases, and block-averages the density; a run can continue from the state of another with a different piston, keeping the strain rate.
  The new explanation page `explanation/constant-cutoff-npt.org` gives the measured runs of 32 PBE waters at 350 K and 1 bar from `bench/results/constant_cutoff_npt.org` and `constant_cutoff_pressure.org`: at E0 = 85 Ry on a 100 Ry deck, the mean over 0.65 ps is 826 ± 6 kg/m³ at 352 K, the error a lower bound, 5% below the PBE density of Miceli et al. (2015), and the molecular pressure of the last state lies within 0.35 kbar of that at E0 = 136 Ry; from that state, fixed cutoffs of 100 and 120 Ry carry Pulay stresses of 17.4 and 2.9 kbar and stand 35% and 6% denser after 40 fs than the constant cutoff on the same deck.
  Its figure comes from `docs/render_constant_cutoff_npt.py`, which the `constant-cutoff-npt-plot` Meson test runs on `bench/results/constant_cutoff_npt.json`, and `tests/python/test_constant_cutoff_npt_table.py` requires the page to carry the generated tables.
  `tests/python/test_ase_calculator.py` checks the constant-cutoff stress against -dE/dV over cells moved in place, and `tests/python/test_lammps.py` and `tests/python/test_gromacs_opencpmd.py` check that `system.constantCutoff` gives LAMMPS and GROMACS the stress of `CPMDCCalculator` for the same configuration.
- `bench/npt_block_average.py` discards the detected equilibration of `bench/gromacs_opencpmd_npt.py` runs and gives Flyvbjerg-Petersen block-averaged volume, density and pressure with their blocking curves; `--classical-only` runs the classical control alone.
  Run GROMACS on cpmdc forces (`docs/orgmode/howto/gromacs.org`) now reports 15 ps of the 64-water C-rescale runs and a 1 ns classical control from `bench/results/gromacs_npt_crescale_64w.org`: the classical box settles at 918.5 +- 3.4 kg/m^3, and the box with the recentred BLYP water at 931.4 +- 7.3 kg/m^3.
- `cpmdc-ipi` (also `python -m cpmdc.ipi`) is a client for the i-PI socket protocol over UNIX and TCP sockets, so i-PI, ASE's `SocketIOCalculator` and `SocketServer`, and through i-PI a LAMMPS `fix ipi` run, can drive one cpmdc session.
  It returns Hartree energies, Hartree/Bohr forces, and the virial `-dE/d(strain)` in Hartree, keeps one `Session` across steps so wavefunction extrapolation still works, and stops with a clear error when the atom count changes; the species come from `--species` or `--structure` because the protocol carries none.
  Drive cpmdc over the i-PI socket protocol (`docs/orgmode/howto/socket-drivers.org`) covers i-PI, ASE, LAMMPS (untested), units, and limits.
  `cpmdc-ipi-c`, built by Meson from `tools/cpmdc_ipi.c` and installed beside `libcpmdc`, is the same client in C on `cpmdc_session_calculate_au()`, for hosts without Python; the pytest suite checks it against ASE, i-PI, and `cpmdc-ipi` to 1e-10.
- `cpmdc.ATOMIC_NUMBER_MIN` and `cpmdc.ATOMIC_NUMBER_MAX` (1 and 118) name the accepted range of atomic numbers, and `cpmdc.ipi.HEADER_BYTES` and `cpmdc.ipi.DEFAULT_SOCKET_PREFIX` are listed in `__all__`, so `from cpmdc.ipi import *` and the reference audit see them.
  The names were already importable.
- `cpmdc_session_calculate_au()` evaluates one session step from plain arrays in atomic units: Bohr positions, atomic numbers, and an optional row-major 3 by 3 cell in Bohr, returning the energy in Hartree, forces in Hartree/Bohr, and optionally the virial `-dE/d(strain)` in Hartree.
  Its arguments come in the order of `cpmdc_session_energy_forces()` and `Session.calculate_au`: session, `n_atoms`, `positions_bohr`, `atomic_numbers`, `cell_bohr`, then the outputs.
  A periodic cell no longer needs a Cap'n Proto `ForceInput`; the results are bit-for-bit those of `cpmdc_session_calculate_result()` with Bohr and Hartree units, which a C test checks for clusters and periodic cells.
  `libcpmdc` now exports 45 functions; the addition is backward compatible, so `CPMDC_ABI_VERSION` and the soversion stay the same.
- `include/cpmdc_export.h` defines `CPMDC_API`, which marks every public function: default visibility with GCC and Clang, and `dllexport` or `dllimport` on Windows (`CPMDC_STATIC` makes it empty for a static link).
  `libcpmdc` is built with hidden visibility, and on ELF and macOS it also links with an export list that `scripts/abi_inventory.py` writes from the same declarations, because gfortran does not apply `-fvisibility` to the Fortran `BIND(C)` glue.
- `libcpmdc_gromacs` runs under a GROMACS built with `-DGMX_MPI=ON`: it duplicates the simulation communicator GROMACS passes to `cp2k_create_force_env_comm`, broadcasts rank 0's input, positions, and results, and hands the communicator to OpenCPMD, with every rank as one calculator (`QM_RANKS ALL`, the default) or rank 0 alone (`QM_RANKS 1`, also `cpmdc.gromacs.qmi_input(qm_ranks=1)`).
  `tests/python/test_gromacs.py` checks `mpirun -np 2` and `-np 4` with domain decomposition and a QM water across domains against one rank to 1e-10, the `gromacs-qmmm` CI job builds GROMACS with MPI and against an installed cpmdc, and `bench/gromacs_opencpmd.py` runs one QM water in MM water on OpenCPMD and reports the NVE drift and the agreement between rank counts.
- `pip install cpmdc` installs the Python package with its own `libcpmdc` and schema; wheels cover Linux x86_64 and aarch64 and macOS arm64, without OpenCPMD or MPI.
  A periodic cell needs no other tool, because the package encodes each step as `ForceInput` itself; the `capnp` command is needed only to encode a `CPMDParams` text file.
- `tests/python` is a pytest suite for the ASE calculator on the reference evaluator: forces against finite differences, stress against numerical stress, units and sign, translation, caching, boundaries, `set()`, copies, and BFGS with a trajectory.
  Meson runs it as `python-pytest`, the wheel builds run it against the installed package, and the `test` extra installs pytest.
  The Python how-to gains a property table, a section on units, boundaries, and caching, and examples that run as written.
- cpmdc adds cpmdc_session_copy_orbitals to copy live session orbital state from one key to another key.
  The Python Session API now exposes select_orbitals and copy_orbitals for keyed orbital workflows.

### Changed

- An embed SCF that stops short of `CONVERGENCE ORBITALS`, after the ODIIS pass and the `PCG MINIMIZE` continuation, now fails with a message that names the minimizer of each pass, the iteration each stopped at, `MAXITER`, and the final `gemax`, for example `CPMD SCF did not converge: PCG MINIMIZE stopped at iteration 2 of MAXITER 2 with gemax ... above CONVERGENCE ORBITALS 1.000000E-05, after ODIIS stopped at iteration 2 of MAXITER 2`. `CPMDCResult.message` and `cpmdc_last_error()` carry it, so rgpot's `CPMD engine failed:` line shows it as it is. Before, the message was the fixed "orbitals not converged within MAXITER, or no energy".
- Each `CPMDCCalculator` step runs with `calc.directory` as the working directory, so a relative `CPMDC_PSEUDO_DIR` or `CPMD_PP_LIBRARY_PATH` resolves there, and so does the `RESTART.1` that a deck with `RESTART WAVEFUNCTION` reads when the message sets neither `permanentDir` nor `scratchDir`.
  The default directory, `.`, changes nothing.
  `charge`, `multiplicity`, `xc`, or any other keyword that is not `params`, `library`, or one of ASE's raises a `TypeError` that says the `CPMDParams` message holds the method, and a calculator without `params` raises `TypeError` with "params is required".
- The GROMACS, LAMMPS, and OpenCPMD archive how-tos keep their steps; their reference material moves to the new pages `reference/gromacs.org` (the `&CPMDC` keys, units, cell, embedding, QM virial, and rank layout) and `reference/lammps.org` (units, cell, and virial), and their explanations and measurements to `explanation/md-engines.org` (why the CP2K interface, the OpenCPMD runs in GROMACS, and the overhead per step of every route) and `explanation/opencpmd-patches.org` (what each patch changes inside CPMD, and why the MPI sentinels stay exported). The GROMACS and LAMMPS test descriptions move to the contributing guide. The GROMACS how-to names `write_qmi()` as the call its example makes.
- The Python modules new since 0.4.0 take their arguments as `Session` does: `cpmdc.ipi.Driver(params, atomic_numbers, ...)` names the atomic numbers `atomic_numbers` (also its attribute), `Driver.compute(positions, cell)` takes the positions before the i-PI cell, and `cpmdc.lammps.type_numbers(elements, n_types)` names the type count `n_types`.
- The README no longer repeats the ABI symbol list, the option-mapping sample, the feature ID namespaces, or the inventory test table; it links to the C ABI reference, the CPMD option mapping, and the feature inventory reference. `feature-inventory` now requires every ABI symbol in `reference/c-abi.org`, matched as a whole word, and the namespace terms in `reference/cpmd-options.org`. `inventory-docs` requires the inventory and docs tests in `reference/feature-inventory.org`, which now lists all of them. `cpmd-public-abi-inventory` plants an omission in a copy of the C ABI reference and expects the check to fail.
- The documentation landing page no longer repeats the README's input-line and runtime-contract tables; it links to the C ABI reference and the CPMD option mapping, which hold the complete versions. The landing page's SIMD row now says the force path has no hot kernel, matching the SIMD explanation, instead of saying the library has none.
- The embedding how-to replaces its copies of the result-buffer, session, and unit contracts with a short list of the rules a host loop follows, linked to the C ABI reference, and its account of force export on a live OpenCPMD link moves to `explanation/opencpmd-patches.org`. The environment reference keeps the rules of `CPMDC_RECENTER`; why the library moves the atoms, and the drift measurements, move to the new `explanation/isolated-cells.org`.
- The quickstart is now the how-to "Build and test cpmdc", filed under the how-to guides at the same address, with one test-selection table, its OpenCPMD steps left to the install how-to, and its failure table merged into troubleshooting. A new reference page, `reference/build-options.org`, lists every Meson option, installed artifact, and Pixi task, including `with_readcon`, `with_pwscale`, `with_lammps_plugin`, `lammps_include`, `with_gromacs_shim`, `test-readcon`, `cpmdc-ipi-c`, the LAMMPS plugin, and `libcpmdc_gromacs`. The README Documentation section links every page, and `readme-navigation` checks that it does.
- `cpmdc.Evaluation` is a dataclass, so it has a `repr` that lists the energy, forces and stress, and the keyword form `Evaluation(energy_h=..., forces_h_bohr=..., stress_h_bohr3=...)`.
  Instances still compare and hash by identity, as before.

### Fixed

- A second OpenCPMD setup in one process, after a cell change, no longer crashes in `stress`: `tools/opencpmd_embed_teardown.patch` now resets its first-call flag, so it allocates again the stress arrays the teardown frees. Before, the stress of any cell after the first (an ASE strain difference or cell relaxation, or GROMACS pressure coupling with `CELL PERIODIC`) faulted.
- An OpenCPMD build of `libcpmdc` no longer hides the MPI Fortran sentinel common blocks (`MPI_IN_PLACE` and the rest, of Open MPI and MPICH); with private copies, OpenCPMD's in-place collectives overwrote their buffers and the first calculation crashed or hung.
  `scripts/ci/build_gromacs_cpmdc.sh` builds GROMACS against the BLAS `libcpmdc` links instead of GROMACS' bundled BLAS, which took the place of OpenCPMD's and gave wrong forces on atoms with nonlocal pseudopotentials; `libcpmdc_gromacs` stops `mdrun` when it finds that BLAS.
- Documentation corrections: `CPMDC_DECK_OUT` holds only the method deck, not the deck the OpenCPMD path composes; `CPMDC_PSEUDO_DIR` and `CPMD_PP_LIBRARY_PATH` are read in `src/cpmdc_embed_host.c`; `CPMDC_RECENTER` also takes `0`, `MASS`, `GEOMETRY`, and `OFF`; a cold call is the first call of the process, not of a session; the `cpmdc_configure()` header comment lists every `common` overlay field it lowers; `cpmdc_session_calculate_au()` rejects `n_atoms` above `INT_MAX / 3` as well; an empty unit string reads as Angstrom or Hartree; the archive how-to lists `opencpmd_omp_info.patch` and the optional `opencpmd_pcg_step.patch`; the release checklist adds news fragments instead of editing `CHANGELOG.md`, and its version-sync list includes `CITATION.cff` and `README.md`; the surface map describes its `@since` column as it is; the README stress unit renders as Bohr³.
  The C ABI reference says that releases up to and including 0.4.0 export the 56 internal glue functions, instead of describing the current release in the past tense.
- On POSIX, a dangling `RESTART` symlink, such as `RESTART` pointing at a missing `RESTART.1`, is published to that target, following a chain of more than one dangling link to its actual end. A relative target is resolved from the directory that holds the link, and the new file takes mode 0666 under the umask, as `fopen` did. A symlink loop, and any other failure to read the link, is left in place and reported. Windows does not resolve a `RESTART` symlink at all.
- On POSIX, publishing through a `RESTART` symlink rechecks the chain immediately before `renameat` and fails if a followed symlink or the final component is no longer what was resolved, so a link swapped for a regular file in that window is not overwritten. An absent final name uses `renameat2` with `RENAME_NOREPLACE` where that call exists. Replacing an existing file still has a gap between the check and the rename, because POSIX has no rename that succeeds only when the directory entry is the same inode.
- On an isolated OpenCPMD deck each evaluation moves the atoms rigidly to the middle of the deck's cell, by their centre of mass, and removes the net force the fixed cell leaves, so energy, forces and `sum x (x) f` no longer depend on where a host puts the atoms.
  A BLYP water 0.5 Angstrom from a face of its 10 Angstrom Hockney cell had been 5.0e-4 Hartree off, and one 0.5 Angstrom past a face 8.4e-3 Hartree off with a net force of 0.17 Hartree/Bohr, which entered the `CELL NONE` virial of `libcpmdc_gromacs` as the QM region diffused in NPT runs.
  `CPMDC_RECENTER` chooses `mass`, `geometry` or `off`, and `CPMDC_RECENTER_MARGIN` makes a call fail when an atom ends up closer than that many Angstrom to a face.
  The 64-water C-rescale QM/MM run whose box had grown to 2.5 nm^3 in 15 ps, a mean of 844.6 kg/m^3, now holds 931.4 +- 7.3 kg/m^3, against 918.5 +- 3.4 kg/m^3 for the 1 ns classical control.
- Reference contract checking now recognizes struct field tables from the table shape alone, so a nonstandard line before the table still associates fields with the right struct and rejects bogus fields. Meson logs when the pwscale FFT failure test is not registered because FFTW is missing, and the Pixi environment pins FFTW so that test is available there.
- The C ABI reference gives every function of `cpmdc.h`, `cpmdc_features.h`, and `cpmdc_restart.h` its full parameter list and contract (the array and per-session snapshot calls had `...` or a bare name), lists `CPMDCRankSpread` with the other snapshots, and states the return, naming, argument-order, and unit-suffix conventions in one table.
  The OpenCPMD archive how-to and the wavefunction-state explanation no longer present the hidden `cpmdc_embed_*` and stop-hook functions as if they were callable.
- The Python reference lists `ATOMIC_NUMBER_MIN`, `ATOMIC_NUMBER_MAX`, the `Evaluation` fields, and the `cpmdc.ipi` header and socket-prefix constants, and its error table gives the `TypeError` and `ValueError` cases of the boundary checks.
- The README's session example allocates the `PotentialResult` buffer it passes, checks the result, and frees the buffer; the LAMMPS how-to marks its two `meson setup` lines as alternatives.
- The RST export keeps the space or line break after a superscript or subscript, which ox-rst ended with an escaping backslash, so `explanation/md-engines.org` read "g/cm³in a 1.24 nm box".
- The Sphinx build is warning-free again: the first-energy tutorial no longer puts `CPMDC_RECENTER` in a link description, where the RST export read the underscore as a subscript, and the `cpmdc_pwscale.h` comments write the tolerance with `fabs()` instead of `|x|`, which Doxyrest exported as RST substitutions.
- The eOn tutorial now builds `eonclient` against rgpot `5fec533f` in cpmdc's Pixi environment, says to start in the directory that holds `cpmdc-tutorial` (the copy of `water.params.bin` failed from where the first tutorial ends), and shows the `results.dat` and final energy of a real run: 7 force calls and -464.1301944953 eV on one rank, and the same minimisation on four ranks. It notes that `rgpot_version` in `results.dat` is a constant compiled into eOn, not the rgpot it loaded.
- The evaluation calls of the C library return a failed `CPMDCResult` with the message `atomic numbers must be in 1..118` for an atomic number outside 1..118.
  A session used to accept such a number on its first step, because the range was only checked when a later step was compared with the first, and then refused the identical step with `topology change requires a new session`.
  The check runs before a session fixes its topology, so a refused step leaves the session usable.
- The header and the C ABI reference now say what `CPMDCPropertySnapshot` holds: `polarizability` is always zeros with a count of 9 (no build computes it), and `hessian` is the nuclear gradient in Hartree/Bohr, cut at 4096 values, so a system of more than 1365 atoms is truncated and `hessian_count` gives the number kept.
  They claimed that polarizability was filled when PROP is available and did not mention the limit.
  The `cpmdc_energy` comment no longer says that no gradient is allocated; the call evaluates the gradient into a temporary array.
- The macOS `libcpmdc` link succeeds on `ld64`. The build documentation now explains the Darwin linker constraint behind the serial wheel fix: `libcpmdc` references `cpmd_stopgm_hook` only when `with_cpmd=true` and `cpmd_root` is non-empty (an OpenCPMD build), where the archive defines it, and the macOS export list adds MPI Fortran sentinel common block wildcard patterns to `-exported_symbols_list` only in that build, only when those symbols can be provided, so neither serial builds nor Open MPI builds fail on absent literal names.
- `CPMDCCalculator.calculate()` called directly replaces every stored result, so the stress of an earlier periodic step is not left beside the energy of a cluster, and a step that fails stores nothing instead of leaving the previous step under the new structure.
  A step whose energy, forces, or stress is NaN or infinite raises `cpmdc.ase.CalculationFailed`, and so does a failure to open the session.
  The Voigt stress averages the two off-diagonal elements of each shear component, as ASE's `full_3x3_to_voigt_6_stress` does, instead of taking the upper triangle.
- `CPMDCCalculator` follows ASE's `Calculator` contract.
  `params` and `library` are `parameters` with `default_parameters`, and `todict()` returns them for a trajectory, leaving out a `params` buffer given as bytes.
  `set()` rejects unknown keys, opens the new session before it drops the old one, clears the results, and returns the changed keys.
  `copy.copy`, `copy.deepcopy`, and `pickle` work, and a copy opens its own session; so does a calculator after `close()`.
  A fully periodic structure without three independent cell vectors raises `ValueError` instead of running as a cluster.
  Stress on free boundaries or on an isolated deck raises `PropertyNotImplementedError` with the reason, and asking again does not evaluate the step again.
  A new composition opens a new session instead of failing; a reordering of the same atoms keeps the session.
  A change of `initial_charges` or `initial_magmoms` no longer evaluates again, because `CPMDParams` sets the charge and multiplicity.
  `Session` encodes a periodic `ForceInput` in process instead of running `capnp encode` on every step, and a closed `Session` raises `CpmdcError`.
- `EXTRAPOLATE WFN` on more than one QM rank no longer hangs at the second call: `tools/opencpmd_embed_extrapolate.patch` broadcasts the request, which only the I/O rank read, so every rank takes the prediction.
- `Session.calculate_au` no longer returns NaN for non-finite positions or cell, and `Session.energy_forces` no longer truncates `14.9` to the atomic number 14.
  Atomic numbers must be integers (a whole float such as `8.0` is still accepted) in 1..118, and orbital keys must be integers that fit a signed 64-bit value: `Session.select_orbitals(2**70)` used to wrap to key 0 and `select_orbitals(1.5)` to key 1.
  Strings, `None` and `bool` raise `TypeError` and the other bad values raise `ValueError`, with a message that names the argument, instead of a bare ctypes error or a misleading "topology change" error on the next step.
- `cpmdc.gromacs.qmi_input` and `write_qmi` reject a non-finite or non-positive `ewald_rcut` or `ewald_precision`, a non-finite `qm_charges` entry, and a `charge` or `multiplicity` that is not an integer, with a message that names the option.
  A `multiplicity` below 1 raises `ValueError`.
  They used to write `nan` or `inf` into the file and truncate `charge=1.7` to 1.
- `cpmdc.lammps.FixExternal` refuses a step whose energy, forces, or virial is NaN or infinite, as `CPMDCCalculator` does, and opens a new session after a change of composition even when the previous attempt to open one failed.
- `cpmdc_adopt_calculator_comm()` no longer refuses a communicator while a session exists that has not evaluated yet. Only a live session that has started an evaluation blocks the adopt, and the error now says "a session has evaluated". rgpot's `CPMDPot` creates its session before it binds calculator groups, so `mpirun -np 4 eonclient` with the cpmdc backend stopped on every rank with "ranks_per_image=0 does not divide the MPI world into calculator groups".
- `cpmdc_restart.h` and the C ABI reference state the failure return of every RESTART function (-1, or `NULL` with the reason in `err`), that `cpmdc_restart_write_mem` returns a `malloc` buffer for `free`, and that velocities are Bohr per atomic unit of time, not Bohr; the reference gives each function its own contract instead of a list of names.
- `cpmdc_restart_write_path()`, and with it `cpmdc-restart patch-positions`, replaces a `RESTART` file in one step. The image goes to a short-named temporary file in the effective destination's directory, which is flushed and synced with `fsync`, renamed over the destination, and followed by an `fsync` of the directory; on Windows `_commit` and `MoveFileEx` with `MOVEFILE_REPLACE_EXISTING` and `MOVEFILE_WRITE_THROUGH` do the same. A reader sees the previous file or the finished one, a failure before the rename removes the temporary file and leaves the old `RESTART` in place, an existing file keeps its mode, and a new file takes mode 0666 under the umask as a plain `fopen` would. On POSIX, a symbolic-link destination updates its target without replacing the link; on Windows the destination is written as given, with no symlink resolution. Before, the writer truncated the destination first, so a reader or a crash during the write saw a partial file.
- `cpmdc_session_copy_orbitals` refuses a donor whose last pass did not converge, including a parked copy and a copy onto the same key. An unconverged pass that already held a converged copy also empties the `EXTRAPOLATE WFN` history: `tools/opencpmd_embed_rwfopt.patch` marks the donor stale and `tools/opencpmd_embed_extrapolate.patch` breaks that history.
- `test_bfgs_moves_the_cell_through_a_cell_filter` again requires the cell relaxation of the periodic water to end below its starting energy: with a 30 Ry constant cutoff on the 40 Ry deck, OpenCPMD converges within 10% of the starting volume, where the fixed-basis Pulay stress had collapsed the cell. `test_a_fixed_basis_cell_filter_collapses_the_cell` and `test_stress_follows_the_energy_over_a_finite_strain` keep the fixed-basis artifact in view, and the ASE and GROMACS how-tos point to the constant cutoff instead of warning about the cutoff.

### Developer

- AddressSanitizer leak detection stays on for the C and Fortran tests.
  `params-cp-dft-render` keeps its fixture paths in `main` instead of
  leaving a heap copy for the process exit. Tests that preload `libasan`
  into Python set `detect_leaks=0`, because CPython and a compiler started
  from those tests do not free their pools. The pixi sanitizer task and
  the GCC 13 CI job use that same split. `shared-dlopen-abi` suppresses
  the 24 byte `ld-linux` allocation in `_dl_catch_exception` that GCC 14
  LeakSanitizer reports during `dlopen`.
- The Meson suite passes against an OpenCPMD build: `session-cached-geometry` uses the cluster water deck, the OpenCPMD tests that run SCFs on one core get longer timeouts, `docs-examples` and `tutorials-run` link and run against the build directory under test and compare an OpenCPMD transcript by its energy and forces, and `python-api` uses finite-difference tolerances and a strain that suit an SCF and OpenCPMD's fixed-basis stress.
  `scripts/ci/fetch_pseudopotentials.sh` assembles the pseudopotentials the tests read, `Si_MT_BLYP.psp` included, from the pinned Regtests checkout.
  `tests/python/test_recentre_opencpmd.py` checks that a cluster translated across its isolated cell keeps its energy, forces and virial, and `bench/cluster_drift_scan.py` prints the same scan.
- The `reference-audit` test, also run by `pixi run test-stub` with the rest of the Meson suite, fails when a public name has no contract entry in `reference/contracts.org`, or when an entry names a symbol that no longer exists.
  It reads the C names from the headers that `meson install` places (the install plan of the build directory) and the Python names from the imported package.
  A C function needs entries for its units, array layout, ownership, threads and MPI ranks, return and error codes, Python mapping, and defaults.
  A Python name needs its signature, units, shape and dtype, ownership and threads, the exceptions it raises, and its defaults; the signature must equal what `inspect` reports, so a changed default or a renamed parameter fails until the page follows.
- The `reference-contract` test, also run as `pixi run check-api-reference`, fails when a public function, type, enum constant, macro, or struct field of a header in `include/` is missing from its reference page, `reference/c-abi.org` or `reference/pwscale.org` (a function needs its parameter list in a table row, and a struct field, enum constant, or macro needs a table row that starts with it), when a header has no page assigned or its page is not in the index toctree, when a module of `python/cpmdc` has no `__all__` or leaves a public function or class out of it, when a name in `__all__` or a public member of an exported class, such as `CPMDCCalculator.ignored_changes`, is missing from the Python reference, or when the README or an Org page names a `cpmdc_` identifier that no public header declares.
- The sdist workflow installs the `test` extra before the installed-wheel check, so pytest is in that environment. The `ase` extra does not install pytest, and the check runs `python -m pytest`.
- `docs-examples` covers the README as well as the Org pages, accepts Org blocks indented under list items, compiles the C fragments against `include/`, lists the blocks that need OpenCPMD separately from the other skips and requires their text to say so, and checks every shell block's meson options, pixi tasks, meson test names, repository paths, build outputs, and `--species` strings. It now runs the deck checks of the CPMDParams and debug-deck how-tos, the troubleshooting introspection, the eOn input files, the README `cpmdc-restart` commands, and the example host run of the first-energy tutorial.
- `python/cpmdc` passes `mypy --strict` with the `[tool.mypy]` table of `pyproject.toml`, and the `python-mypy-strict` Meson test runs it (it is skipped when mypy is not installed; the sanitizer and wheel jobs of the workflow install it).
  `cpmdc.gromacs.write_qmi` now lists the keyword options of `qmi_input` in its signature instead of accepting `**options`; the accepted keywords are unchanged.
- `schema/inventory/surface_map_abi.txt` is generated by `scripts/abi_inventory.py` from the symbols `libcpmdc` exports and the declarations in `include/cpmdc.h` and `include/cpmdc_features.h`; it lists the 45 public functions, which are all the library exports.
  The `c-abi-inventory` Meson test and the CI stub job fail when the file is stale, when a public declaration lacks `CPMDC_API` or is not exported, or when the library exports any other symbol; `meson compile -C build abi-inventory` regenerates it.
- `test_cpmd_restart` makes its scratch directory under `TMPDIR`, falling back to the C library's `P_tmpdir`, and removes it at the end.
- `tests/python/test_ase_conformance.py` checks `CPMDCCalculator` with ASE's own tools on the default build: honest `implemented_properties`, the cache and `system_changes`, one session for a run of steps, `PropertyNotImplementedError` and `CalculationFailed`, `BFGS` with `FixAtoms`, copies and pickles, and the eV and eV/Angstrom units against `ase.units`.
- `tests/test_extrapolate_wfn.c` compiles again in an OpenCPMD build tree with GCC 14, which stopped on the implicit declaration of `setenv` under `-std=c11`; it defines `_POSIX_C_SOURCE 200809L` like the other tests that call POSIX.
  The `extrapolate-wfn` Meson test is now defined in every build and skips with exit 77 when `cpmdc_available()` is 0, as in the reference build.
  `scripts/ci/build_opencpmd.sh` builds the patched OpenCPMD archive, and the new CI job `opencpmd-embed` compiles the `with_cpmd` build against it and fails unless `extrapolate-wfn` passes against the `cpmd.x` run in `tests/golden`.

### Miscellaneous

- The Re-check section of `docs/AUDIT.md` records gaps seven and eight as closed, with the eOn run, the steps it fixed, and the checks that moved to the reference pages.


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
