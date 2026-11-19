# cpmdc

`cpmdc` is the C ABI layer for embedding
[OpenCPMD](https://github.com/OpenCPMD/CPMD) in OmniPotentRPC tools. It gives
non-Fortran callers one stable boundary for CPMD setup, geometry steps, feature
discovery, and result messages.

Documentation: [cpmdc.rgoswami.me](https://cpmdc.rgoswami.me).

Use this repository when you need to:

- pass CPMD method setup as a serialized Cap'n Proto `CPMDParams` message
- pass each geometry step as a serialized `ForceInput` message
- receive native C results or a serialized `PotentialResult` message
- load the backend from C, C++, Rust, Python, or an RPC host without exposing
  Fortran or C++ internals

`CPMDParams` owns method and backend setup. `ForceInput` owns per-step
coordinates, species, cell, and requested units. `PotentialResult` is the
portable result carrier used by RPC-facing callers.

## Start Here

| Need | Do this first |
| --- | --- |
| Choose the API path first | Pick a row from [Main API Choices](#main-api-choices), then read the matching call contract in [Embedding cpmdc](docs/orgmode/howto/embedding.org) |
| Run LAMMPS on cpmdc forces | Load `pair_style cpmdc` from the plugin, or attach `cpmdc.lammps.FixExternal`, as [Run LAMMPS on cpmdc forces](docs/orgmode/howto/lammps.org) shows |
| Call it from ASE or Python | Run `pip install "cpmdc[ase]"` as in [Python Quickstart](#python-quickstart), then read [Use cpmdc as an ASE calculator](docs/orgmode/howto/ase.org), or [Call cpmdc from Python](docs/orgmode/howto/python.org) without ASE |
| Build and test the stub library | Run the default Meson commands in [Build and Test](#build-and-test) |
| Map a CPMD input line | Find the line in the [CPMD option mapping](docs/orgmode/reference/cpmd-options.org), as [Map a CPMD Input Line](#map-a-cpmd-input-line) describes |
| Link an OpenCPMD archive | Configure `-Dwith_cpmd=true` as shown in [OpenCPMD Archive Build](#opencpmd-archive-build) |
| Regenerate generated docs | Edit `docs/orgmode/`, then run `pixi run -e docs mkrst` and `pixi run -e docs sphinxbld` |
| Ask a question or report a bug | Open an issue as [the reporting section](#questions-and-bug-reports) describes |
| Send a fix or an improvement | Read the [contributing guide](docs/orgmode/contributing/index.org), then open a pull request |
| Cite cpmdc in a paper | Use the record in [`CITATION.cff`](CITATION.cff), as [Citation](#citation) describes |

New RPC and optimizer hosts should start with
`cpmdc_session_calculate_result()`. It keeps `CPMDParams` method setup in a
session, accepts one `ForceInput` per step, and writes the same
`PotentialResult` message used by the RPC layer.

## Python Quickstart

```bash
pip install "cpmdc[ase]"
```

Wheels cover Linux x86_64 and aarch64 and macOS arm64, and carry their own
`libcpmdc`: the default reference evaluator, without OpenCPMD or MPI. Run
[examples/python/ase_energy.py](examples/python/ase_energy.py) with
`CPMDC_PARAMS` set to a serialized `CPMDParams` file. To use OpenCPMD, point
`CPMDC_LIBRARY` at a `-Dwith_cpmd=true` build.
`cpmdc.ase.CPMDCCalculator` gives energy, free energy, forces, and, for a
fully periodic cell, stress in ASE units, and works with ASE optimizers, cell
filters, and trajectories;
[Use cpmdc as an ASE calculator](docs/orgmode/howto/ase.org) has the property
table and the rules for units, boundaries, caching, and the working directory.
The [install page](docs/orgmode/howto/install.org#install-the-python-package)
lists what the wheel needs and how other platforms build from the sdist.

## Build and Test

The default build does not require an OpenCPMD checkout. It builds the public C
ABI, Cap'n Proto readers, CPMD deck renderer, feature table, shared library, and
a deterministic reference evaluator used by the tests.

```bash
meson setup build -Dwith_tests=true
meson compile -C build
meson test -C build --print-errorlogs
```

With Pixi, the same default stub-library check is:

```bash
pixi run test-stub
```

A healthy default run passes the parser, ABI, feature inventory, session,
result-buffer, unit-conversion, and reference evaluator tests.

## Repository Map

| Path | Purpose |
| --- | --- |
| `include/cpmdc.h` | public C ABI, sessions, one-shot calls, result sizing |
| `include/cpmdc_restart.h` | C reader and writer for a CPMD `RESTART` file |
| `tools/cpmdc_restart.c` | `cpmdc-restart` command for coordinates, velocities, and cell |
| `tools/cpmdc_fd_check.c` | `cpmdc-fd-check` finite-difference check of session forces and stress |
| `tools/lammps/` | `cpmdcplugin.so`, the LAMMPS plugin with `pair_style cpmdc` |
| `python/cpmdc/` | the Python package: `Session`, the ASE calculator, the i-PI driver, `cpmdc.lammps`, and `cpmdc.gromacs` |
| `tools/gromacs/` | `libcpmdc_gromacs`, cpmdc behind the libcp2k API of GROMACS' QM/MM interface, and its `libcp2k.pc` |
| `bench/lammps_routes.py` | per-step overhead of the LAMMPS and i-PI routes |
| `bench/gromacs_routes.py` | per-step overhead of the GROMACS route beside the LAMMPS plugin |
| `include/cpmdc_features.h` | runtime feature discovery table |
| `include/cpmdc_export.h` | `CPMDC_API`, the mark on every exported function |
| `schema/Potentials.capnp` | `CPMDParams`, `ForceInput`, `PotentialResult`, `PotentialConfig` |
| `src/cpmdc_params.c` | Cap'n Proto decode and CPMD `INPUT` rendering |
| `src/cpmdc.c` | session lifecycle, topology checks, units, result writing |
| `src/cpmd_embed_c_api.F90` | Fortran `iso_c_binding` bridge and OpenCPMD archive path |
| `tests/` | encoded fixtures and C/cmocka contract tests |
| `docs/orgmode/` | authoritative documentation sources |

The docs build writes generated RST under `docs/source/`. Those files are not committed.

## Main API Choices

| Caller need | Entry point |
| --- | --- |
| Long-running driver with serialized input and output | `cpmdc_session_create`, then `cpmdc_session_calculate_result` |
| Long-running native C force loop | `cpmdc_session_energy_forces` |
| Native C, Fortran, or plugin loop with a periodic cell, in atomic units | `cpmdc_session_calculate_au` |
| LAMMPS integrating on cpmdc forces | `pair_style cpmdc` or `cpmdc.lammps.FixExternal`, both on `cpmdc_session_calculate_au` |
| GROMACS integrating on cpmdc forces, QM group or whole system | `libcpmdc_gromacs` linked as libcp2k (`-DGMX_CP2K=ON`), on `cpmdc_session_calculate_au` |
| Single serialized calculation | `cpmdc_calculate_result` |
| Native C compatibility call | `cpmdc_energy_forces` |
| Runtime capability checks | `cpmdc_feature_find` or `cpmdc_feature_table` |

## Map a CPMD Input Line

The [CPMD option mapping](docs/orgmode/reference/cpmd-options.org) gives the
writable `CPMDParams` carrier, the feature ID namespaces, and the discovery
row for each CPMD manual line. It also gives the order in which to choose a
carrier, from top-level fields to typed sections to raw deck text.

## Public ABI Surface

The [C ABI reference](docs/orgmode/reference/c-abi.org) lists every exported
function by group, with its units and session rules.
[Embedding cpmdc](docs/orgmode/howto/embedding.org) shows a session loop that
sizes the `PotentialResult` buffer and calls `cpmdc_session_calculate_result()`.

## Host Example

`examples/host_step.c` is the smallest complete host-side program in the tree.
It reads one serialized `CPMDParams` file and one serialized `ForceInput` file,
creates a session, sizes the `PotentialResult` output, evaluates one step, and
prints:

```text
energy_h=...
potential_result_size_bytes=...
message=...
```

Run it through Meson:

```bash
meson test -C build example-host-step --print-errorlogs
```

The `example-host-step` test uses generated fixture binaries from the same
Cap'n Proto text fixtures used by the E2E suites.

## RESTART Files

`cpmdc-restart` reads and writes the Fortran unformatted `RESTART.n` file that
OpenCPMD stores through `wv30` / `rv30`. Coordinates in the file are Bohr, in
species order. Replacing them rewrites those records and copies the
wavefunction records unchanged, so a host such as eOn can keep converged
orbitals and supply the next geometry.

```bash
cpmdc-restart info RESTART.1
cpmdc-restart positions RESTART.1 --angstrom > pos.txt
cpmdc-restart patch-positions RESTART.1 pos.txt RESTART.1.next --angstrom
```

Link `libcpmdc_restart` and include `cpmdc_restart.h` for the same operations
from C or C++. The library does not link OpenCPMD.

A build with `-Dwith_readcon=true` adds `export-con` and `import-con`. That
option is off by default, so the CPMD-only build does not link Rust.
readcon-core is pinned in `subprojects/readcon-core.wrap` to git main
`694e808080c49979ee3099afcfa18e5420c40433` until `v0.17.1` is tagged. The
bridge does not link readcon-db. `export-con` multiplies coordinates,
velocities, and the cell lengths by `CPMDC_RESTART_ANGSTROM_PER_BOHR` and
writes one RoundTrip CON with canonical left off and no decimal precision.
`import-con` divides by that same macro and calls
`cpmdc_restart_set_coordinates`, `cpmdc_restart_set_velocities`, and
`cpmdc_restart_set_cell`. The CON file is the geometry that is resumed. The
binary `RESTART` stays the wavefunction.

```bash
cpmdc-restart export-con RESTART.1 water.con --species O,H
cpmdc-restart import-con RESTART.1 water.con RESTART.1.next
```

## Inventory Verification

The default Meson suite includes checks that keep the schema, the feature
table, the renderer, and the reference pages aligned. The
[feature inventory reference](docs/orgmode/reference/feature-inventory.org)
lists each check, what makes it fail, and the command that runs them together.

## OpenCPMD Archive Build

The default build uses a reference evaluator. To link the same C ABI against a
completed OpenCPMD archive build:

```bash
export CPMD_ROOT=/path/to/OpenCPMD/CPMD
export CPMDC_PSEUDO_DIR=/path/to/Regtests/tests/PP_LIBRARY

meson setup build-cpmd \
  -Dwith_cpmd=true \
  -Dcpmd_root="$CPMD_ROOT" \
  -Dwith_tests=true
meson compile -C build-cpmd
meson test -C build-cpmd --print-errorlogs
```

`CPMD_ROOT` must contain `lib/libcpmd.a`. If runtime parameters use
library-style pseudopotential names such as `O_MT_BLYP.psp`, set
`CPMDC_PSEUDO_DIR` or provide absolute paths in `atoms.pseudopotentials`.

To see the deck cpmdc renders from `CPMDParams`, set
`CPMDC_DECK_OUT=/path/to/deck.inp`. Each configuration overwrites the file.
The OpenCPMD path then rebuilds `&ATOMS` from the step, using each listed
pseudopotential entry for the file, `LMAX`, `LOC`, and `KLEINMAN-BYLANDER`,
and adds `CELL` and `MAXITER` when the deck lacks them; see
[Debug an input with CPMDC_DECK_OUT](docs/orgmode/howto/debug-deck.org).

## Documentation

The published site is [cpmdc.rgoswami.me](https://cpmdc.rgoswami.me).

Agent and LLM index files live at the repository root:
[`AGENTS.md`](AGENTS.md) for layout, build, test, and docs workflow, and
[`llms.txt`](llms.txt) plus [`llms-full.txt`](llms-full.txt) generated by
`scripts/generate_llms_txt.py`.

Read the docs in this order:

- Tutorials: [first energy and forces](docs/orgmode/tutorials/first-energy.org),
  a C host from the default build to a real OpenCPMD energy;
  [an energy and forces from ASE](docs/orgmode/tutorials/ase-energy.org);
  [eOn through rgpot](docs/orgmode/tutorials/eon-rgpot.org), a minimisation
  driven by eOn with a `CPMDParams` file
- How-tos: [install and build](docs/orgmode/howto/install.org),
  [build and test cpmdc](docs/orgmode/howto/quickstart.org),
  [OpenCPMD archive](docs/orgmode/howto/opencpmd-archive.org),
  [write a CPMDParams message](docs/orgmode/howto/write-cpmdparams.org),
  [embedding cpmdc](docs/orgmode/howto/embedding.org),
  [call cpmdc from Python](docs/orgmode/howto/python.org),
  [use cpmdc as an ASE calculator](docs/orgmode/howto/ase.org),
  [relax a cell at a constant cutoff](docs/orgmode/howto/constant-cutoff.org),
  [drive cpmdc over the i-PI socket protocol](docs/orgmode/howto/socket-drivers.org),
  [run LAMMPS on cpmdc forces](docs/orgmode/howto/lammps.org),
  [run GROMACS on cpmdc forces](docs/orgmode/howto/gromacs.org),
  [run under mpirun](docs/orgmode/howto/mpi.org),
  [choose the wavefunction optimiser](docs/orgmode/howto/wavefunction-optimiser.org),
  [extrapolate the orbitals](docs/orgmode/howto/extrapolate-orbitals.org),
  [debug a deck](docs/orgmode/howto/debug-deck.org),
  [check forces and stress](docs/orgmode/howto/check-forces.org),
  [time the plane-wave scaffold](docs/orgmode/howto/pwscale-ranks.org),
  [troubleshooting](docs/orgmode/howto/troubleshooting.org)
- Reference: [C ABI](docs/orgmode/reference/c-abi.org),
  [plane-wave scaffold](docs/orgmode/reference/pwscale.org),
  [Python API](docs/orgmode/reference/python.org),
  [GROMACS interface](docs/orgmode/reference/gromacs.org),
  [LAMMPS interface](docs/orgmode/reference/lammps.org),
  [environment variables](docs/orgmode/reference/environment.org),
  [build options, artifacts, and Pixi tasks](docs/orgmode/reference/build-options.org),
  [CPMDParams fields](docs/orgmode/reference/cpmdparams-schema.org),
  [CPMD option mapping](docs/orgmode/reference/cpmd-options.org),
  [feature inventory](docs/orgmode/reference/feature-inventory.org),
  [public C surface map](docs/orgmode/reference/surface-map.org),
  [API contracts](docs/orgmode/reference/contracts.org),
  [glossary](docs/orgmode/reference/glossary.org)
- Explanation: [architecture](docs/orgmode/explanation/architecture.org),
  [file route against in-process route](docs/orgmode/explanation/routes.org),
  [wavefunction state](docs/orgmode/explanation/wavefunction-state.org),
  [MPI model](docs/orgmode/explanation/mpi-model.org),
  [atoms in an isolated cell](docs/orgmode/explanation/isolated-cells.org),
  [cpmdc inside MD engines](docs/orgmode/explanation/md-engines.org),
  [density of water at a constant cutoff](docs/orgmode/explanation/constant-cutoff-npt.org),
  [how the OpenCPMD patches work](docs/orgmode/explanation/opencpmd-patches.org),
  [SIMD paths and dispatch](docs/orgmode/explanation/simd-paths.org)
- [Contributing](docs/orgmode/contributing/index.org): the checks behind each
  contract, and the [release checklist](docs/orgmode/contributing/developer/release.org)

Build the documentation site with:

```bash
pixi run -e docs docbld
```

After editing Org sources, regenerate RST with:

```bash
pixi run -e docs mkrst
```

## Questions and bug reports

Ask questions and report bugs on the [GitHub issue tracker](https://github.com/OmniPotentRPC/cpmdc/issues).
Check [troubleshooting](docs/orgmode/howto/troubleshooting.org) first, because it lists known failures with their fixes.
Its last section lists what to put in a new report.

## Contributing

Bug reports, fixes, tests, documentation, and new typed keywords for the input deck are welcome.
Propose a change on the [GitHub issue tracker](https://github.com/OmniPotentRPC/cpmdc/issues).
Send the change itself as a pull request to this repository.
The [contributing guide](docs/orgmode/contributing/index.org) lists the checks to run first and the files that each kind of change touches.

## Version

<!-- [[[cog cog.outl(f'This README describes cpmdc {open("version.txt").read().strip()}.') ]]] -->
This README describes cpmdc 0.5.0.
<!-- [[[end]]] -->
`version.txt` holds the only hand-written copy of the version.
`pixi run version-sync` regenerates the other copies, and CI fails when
one differs.
[`CHANGELOG.md`](CHANGELOG.md) lists the changes in each release, and
`docs/newsfragments/` holds the entries for the next one.

## Citation

Cite cpmdc with the record in [`CITATION.cff`](CITATION.cff).
It gives the author, the version, the license, and the repository in the
[Citation File Format](https://citation-file-format.github.io).
GitHub reads that file and offers the record as "Cite this repository"
on the repository page, as the
[GitHub documentation on citation files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
describes.
The record carries no digital object identifier (DOI).

An OpenCPMD archive build computes its energies and forces with
OpenCPMD, so a paper that reports them cites OpenCPMD as well.
OpenCPMD ships no citation file; name it with its repository,
[github.com/OpenCPMD/CPMD](https://github.com/OpenCPMD/CPMD), and the
version or commit that the archive was built from.

## License

cpmdc is released under the `MIT License`.
The full terms are in [`LICENSE`](LICENSE).
