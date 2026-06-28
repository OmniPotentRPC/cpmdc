# cpmdc

`cpmdc` is the C ABI layer for embedding
[OpenCPMD](https://github.com/OpenCPMD/CPMD) in OmniPotentRPC tools. It gives
non-Fortran callers one stable boundary for CPMD setup, geometry steps, feature
discovery, and result messages.

Use this repository when you need to:

- pass CPMD method setup as a serialized Cap'n Proto `CPMDParams` message
- pass each geometry step as a serialized `ForceInput` message
- receive native C results or a serialized `PotentialResult` message
- load the backend from C, C++, Rust, Python, or an RPC host without exposing
  Fortran or C++ internals

`CPMDParams` owns method and backend setup. `ForceInput` owns per-step
coordinates, species, cell, and requested units. `PotentialResult` is the
portable result carrier used by RPC-facing callers.

## First Build

The default build does not require an OpenCPMD checkout. It builds the public C
ABI, Cap'n Proto readers, CPMD deck renderer, feature table, shared library, and
a deterministic reference evaluator used by the tests.

```bash
meson setup build -Dwith_tests=true
meson compile -C build
meson test -C build --print-errorlogs
```

With Pixi:

```bash
pixi run test-stub
```

A healthy default run passes the parser, ABI, feature inventory, session,
result-buffer, unit-conversion, and reference evaluator tests.

## Repository Map

| Path | Purpose |
| --- | --- |
| `include/cpmdc.h` | public C ABI, sessions, one-shot calls, result sizing |
| `include/cpmdc_features.h` | runtime feature discovery table |
| `schema/Potentials.capnp` | `CPMDParams`, `ForceInput`, `PotentialResult`, `PotentialConfig` |
| `src/cpmdc_params.c` | Cap'n Proto decode and CPMD `INPUT` rendering |
| `src/cpmdc.c` | session lifecycle, topology checks, units, result writing |
| `src/cpmd_embed_c_api.F90` | Fortran `iso_c_binding` bridge and OpenCPMD archive path |
| `tests/` | encoded fixtures and C/cmocka contract tests |
| `docs/orgmode/` | authoritative documentation sources |

Generated Sphinx sources live under `docs/source/`.

## Main API Choices

| Caller need | Entry point |
| --- | --- |
| Long-running driver with serialized input and output | `cpmdc_session_create`, then `cpmdc_session_calculate_result` |
| Long-running native C force loop | `cpmdc_session_energy_forces` |
| Single serialized calculation | `cpmdc_calculate_result` |
| Native C compatibility call | `cpmdc_energy_forces` |
| Runtime capability checks | `cpmdc_feature_find` or `cpmdc_feature_table` |

## Map a CPMD Input Line

When translating a CPMD manual line, separate the writable schema carrier from
the discovery row:

| Manual line | Writable carrier | Discovery row |
| --- | --- | --- |
| `OPTIMIZE GEOMETRY CLASSICAL XYZ SAMPLE` | `params.inputSections.cpmd.optimizeGeometryOptions` and `optimizeGeometrySample` | `catalog.cpmd.OPTIMIZE_GEOMETRY_CLASSICAL`, `catalog.cpmd.OPTIMIZE_GEOMETRY_XYZ`, `catalog.cpmd.OPTIMIZE_GEOMETRY_SAMPLE_OPTION` |
| `NOSE IONS ULTRA MASSIVE T0` | `params.inputSections.cpmd.noseIonsOptions` | `catalog.cpmd.NOSE_IONS_ULTRA`, `catalog.cpmd.NOSE_IONS_MASSIVE`, `catalog.cpmd.NOSE_IONS_T0` |
| `BOX WALLS` inside `&SYSTEM` | `params.inputSections.system.boxWalls` | `params.inputSections.system.boxWalls` |
| `&VDW ... EMPIRICAL CORRECTION ...` | `inputSections.vdw.directives` and `subsections` | `catalog.section.VDW` |

Use the `params.*` feature ID when the schema owns a typed field. Use
`catalog.cpmd.*` for recognized `&CPMD` keywords or inline option spellings
that render through an existing text field. Use `catalog.section.*` only to
check section support.

For a complete mapping table, see
[`docs/orgmode/reference/cpmd-options.org`](docs/orgmode/reference/cpmd-options.org).

## Public ABI Surface

| Group | Symbols |
| --- | --- |
| Library status | `cpmdc_version`, `cpmdc_available`, `cpmdc_finalize` |
| Feature discovery | `cpmdc_feature_count`, `cpmdc_feature_table`, `cpmdc_feature_find` |
| Global params and coordinate arrays | `cpmdc_set_params`, `cpmdc_energy`, `cpmdc_energy_gradient`, `cpmdc_energy_forces` |
| Session lifecycle | `cpmdc_session_create`, `cpmdc_session_set_params`, `cpmdc_session_destroy` |
| Session coordinate arrays | `cpmdc_session_energy`, `cpmdc_session_energy_gradient`, `cpmdc_session_energy_forces` |
| Session Cap'n Proto steps | `cpmdc_session_calculate_forces`, `cpmdc_session_calculate_result` |
| One-shot Cap'n Proto steps | `cpmdc_calculate_result`, `cpmdc_potential_result_size_for_force_input` |

New RPC-style callers should use the session result path:

```c
CPMDCSession *session =
    cpmdc_session_create(params_bytes, params_size);

size_t capacity = cpmdc_potential_result_size_for_force_input(
    force_input_bytes, force_input_size);

CPMDCResult result = cpmdc_session_calculate_result(
    session,
    force_input_bytes,
    force_input_size,
    potential_result_bytes,
    capacity,
    &potential_result_size);

cpmdc_session_destroy(session);
```

`cpmdc_session_create()` copies the serialized params buffer. The first
successful session evaluation fixes atom count and ordered atomic numbers.
Coordinates, units, and the 3x3 cell may change between later steps; atom count
or species changes require a new session.

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

## CPMD Input Model

Route configuration in this order:

| Need | Carrier |
| --- | --- |
| common method defaults | top-level `CPMDParams` fields |
| typed `&SYSTEM`, `&CPMD`, `&DFT`, and `&ATOMS` controls | `inputSections.system`, `inputSections.cpmd`, `inputSections.dft`, `inputSections.atoms` |
| catalog sections such as `&PIMD`, `&VDW`, `&LINRES`, `&TDDFT` | matching typed section arm with `directives` and `subsections` |
| one merge-only keyword inside a named section | `inputSections.set` |
| a non-catalog section alias | `inputSections.generic` |
| a deck fragment that must remain literal text | `inputBlocks` or `inputSections.raw` |

The feature table exposes stable IDs for the ABI, schema fields, CPMD catalog
keywords, and section arms. Examples:

- `abi.cpmdc_session_calculate_result`
- `params.inputSections.cpmd.maxIter`
- `params.inputSections.system.cell`
- `catalog.cpmd.OPTIMIZE_WAVEFUNCTION`
- `catalog.section.PIMD`

To map a CPMD manual line, first look for a typed
`params.inputSections.<section>.<field>` row. Use the matching
`catalog.cpmd.*` row when the line is an `&CPMD` keyword or an inline option
spelling that the renderer emits through an existing text field. Fall back to
`inputSections.set`, `inputSections.generic`, `inputSections.raw`, or
`inputBlocks` only when the typed schema and catalog rows do not represent the
deck form.

## Feature ID namespaces

Feature IDs separate writable schema fields from parser/catalog capability
rows:

| Namespace | Meaning |
| --- | --- |
| `abi.*` | exported C ABI symbols |
| `params.*` | top-level `CPMDParams` fields |
| `params.inputSections.*` | typed fields and directive carriers inside `CPMDInputSection` |
| `catalog.cpmd.*` | named `&CPMD` keyword capability rows |
| `catalog.dft.*` | named `&DFT` keyword capability rows |
| `catalog.section.*` | OpenCPMD section names discovered from parser `inscan('&SECTION')` calls |

Use `params.inputSections.system.*` for structured `&SYSTEM` controls such as
cell, cutoff, Poisson, KPOINTS, CDFT, pressure, and stress fields.
`catalog.section.SYSTEM` means the `&SYSTEM` deck section is recognized as a
section kind; individual `&SYSTEM` controls stay under
`params.inputSections.system.*` rather than a separate `catalog.system.*`
namespace.
Use `params.inputSections.atoms.*` for pseudopotential entries and fixed
non-coordinate `&ATOMS` keywords; geometry coordinates remain in each
`ForceInput` step.

Manual keyword names do not always map to a `catalog.cpmd.*` feature ID. Use
the field namespace when the renderer owns a typed section field:

| Manual name | Feature or field ID | Rendered deck form |
| --- | --- | --- |
| `BOX WALLS` | `params.inputSections.system.boxWalls` | `&SYSTEM` / `BOX WALLS` |
| `MODIFIED GOEDECKER` | `params.inputSections.system.modifiedGoedecker` | `&SYSTEM` / `MODIFIED GOEDECKER` |
| `MODIFIED GOEDECKER PARAMETERS` | `params.inputSections.system.modifiedGoedeckerParameters` | `&SYSTEM` / `MODIFIED GOEDECKER PARAMETERS` |
| `SOC` | `catalog.cpmd.SOC` and `params.inputSections.cpmd.spinOrbitCouplingStates` | `&CPMD` / `SPIN-ORBIT COUPLING` |

`FFTW WISDOM` is not a typed feature row because the OpenCPMD manual source
comments out that keyword and the parser tree used for inventory probing has no
matching branch. Use `inputSections.raw` or `inputBlocks` only for a downstream
tree that still accepts it.

## OpenCPMD Archive Build

The default build uses a reference evaluator. To link the same C ABI against a
completed OpenCPMD archive build:

```bash
export CPMD_ROOT=/path/to/OpenCPMD/CPMD
export CPMDC_PSEUDO_DIR=/path/to/CPMD-Regtests/tests/PP_LIBRARY

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

## Documentation

Start here:

- [Quickstart](docs/orgmode/tutorials/quickstart.org): build, tests, and common
  environment failures
- [Embedding cpmdc](docs/orgmode/howto/embedding.org): C ABI call flow,
  sessions, result buffers, units
- [CPMD option mapping](docs/orgmode/reference/cpmd-options.org): schema fields
  and feature IDs for rendered CPMD controls
- [Architecture](docs/orgmode/reference/architecture.org): layer boundaries and
  OpenCPMD archive link path

Build the documentation site with:

```bash
pixi run -e docs docbld
```

After editing Org sources, regenerate checked-in RST with:

```bash
pixi run -e docs mkrst
```

## Version

0.1.0 covers the Cap'n Proto contract, stable C ABI, Fortran bridge, shared
`libcpmdc`, deterministic reference evaluator, OpenCPMD archive link path,
parser tests, session tests, and rgpot-compatible result carriers.
