# cpmdc

`cpmdc` is the C ABI layer for embedding
[OpenCPMD](https://github.com/OpenCPMD/CPMD) in OmniPotentRPC tools. It gives
non-Fortran callers one stable boundary:

- method setup is an unpacked flat Cap'n Proto `CPMDParams` message
- per-step geometry is an unpacked flat Cap'n Proto `ForceInput` message
- per-step output can be returned as native C values or as a flat
  `PotentialResult` message

The ABI does not expose C++ types, Rust types, or a second JSON/TOML options
language. All portable configuration lives in `schema/Potentials.capnp`; the C
library decodes those bytes with generated `capnp-c` readers and renders the
CPMD `INPUT` deck internally.

## What Ships

| Component | Path | Purpose |
| --- | --- | --- |
| Public C ABI | `include/cpmdc.h` | `CPMDCResult`, opaque `CPMDCSession`, one-shot and session calls |
| Feature table | `include/cpmdc_features.h` | Runtime discovery for ABI, schema, and CPMD catalog support |
| Wire schema | `schema/Potentials.capnp` | `ForceInput`, `PotentialResult`, `CPMDParams`, `PotentialConfig` |
| Parser/deck renderer | `src/cpmdc_params.c` | Reads `CPMDParams` and renders structured CPMD input sections |
| Session runtime | `src/cpmdc.c` | Owns persistent params, topology checks, unit conversion, result writing |
| Fortran bridge | `src/cpmd_embed_c_api.F90` | `iso_c_binding` shell around the embedded CPMD evaluator |
| Tests and fixtures | `tests/` | cmocka tests plus encoded Cap'n Proto fixtures |

## Feature Discovery

`cpmdc_feature_table()` exposes ABI calls, renderable CPMD catalog entries, and
schema fields. Top-level fields use IDs such as `params.functional`;
structured section fields use `params.inputSections.<section>.<field>`, for
example `params.inputSections.cpmd.maxIter` and
`params.inputSections.dft.hfxScreening`.

## Which Interface To Use

| Caller need | API surface | Data contract |
| --- | --- | --- |
| Configure one long-running driver | `cpmdc_session_create` | one serialized `CPMDParams` |
| Evaluate repeated geometry steps | `cpmdc_session_calculate_result` | one serialized `ForceInput` per step, one serialized `PotentialResult` out |
| Get native C forces without result serialization | `cpmdc_session_energy_forces` | arrays in Angstrom, forces in Hartree/Bohr |
| Compatibility one-shot call | `cpmdc_calculate_result` | serialized `CPMDParams` + serialized `ForceInput` |
| Check support at runtime | `cpmdc_feature_find` / `cpmdc_feature_table` | stable `abi.*`, `params.*`, and `catalog.*` IDs |

## Public ABI Surface

The exported C surface is intentionally small. Use the Cap'n Proto entry points
for new callers; the coordinate-array calls are compatibility shims around the
same session state and native CPMD units.

| Group | Symbols |
| --- | --- |
| Library status | `cpmdc_version`, `cpmdc_available`, `cpmdc_finalize` |
| Feature discovery | `cpmdc_feature_count`, `cpmdc_feature_table`, `cpmdc_feature_find` |
| Global params and coordinate arrays | `cpmdc_set_params`, `cpmdc_energy`, `cpmdc_energy_gradient`, `cpmdc_energy_forces` |
| Session lifecycle | `cpmdc_session_create`, `cpmdc_session_set_params`, `cpmdc_session_destroy` |
| Session coordinate arrays | `cpmdc_session_energy`, `cpmdc_session_energy_gradient`, `cpmdc_session_energy_forces` |
| Session Cap'n Proto steps | `cpmdc_session_calculate_forces`, `cpmdc_session_calculate_result` |
| One-shot Cap'n Proto steps | `cpmdc_calculate_result`, `cpmdc_potential_result_size_for_force_input` |

## Data Flow

There are two wire messages in normal use:

| Message | Lifetime | Contents |
| --- | --- | --- |
| `CPMDParams` | session setup | method, CPMD sections, pseudopotentials, engine hints |
| `ForceInput` | one calculation step | positions, atomic numbers, optional cell, requested output units |

The host builds `CPMDParams` once, serializes it as an unpacked flat Cap'n Proto
message, and creates a `CPMDCSession`. Each geometry step is a separate
`ForceInput`; `cpmdc_session_calculate_result()` returns a flat
`PotentialResult` that uses the requested `ForceInput.energyUnit` and
`ForceInput.lengthUnit`.

Use typed `inputSections` when a field exists. Use `set` for one extra keyword
that should merge into a typed section, `generic` for a complete unsupported
section expressed as keyword/argument pairs, and `raw` only when preserving
existing deck text is more important than structure.

## Runtime Contracts

Session callers should rely on these invariants:

| Concern | Contract |
| --- | --- |
| Parameter lifetime | `cpmdc_session_create()` copies the serialized `CPMDParams`; the caller may release its input buffer after creation. |
| Topology | the first successful session evaluation fixes atom count and ordered atomic numbers; coordinate, cell, and unit changes are allowed after that. |
| Parameter replacement | `cpmdc_session_set_params()` is only valid before the session accepts a topology. |
| Native C units | `CPMDCResult.energy_h` is Hartree; gradient and force buffers are Hartree/Bohr. |
| Result-message units | `PotentialResult.energy` uses `ForceInput.energyUnit`; `PotentialResult.forces` use `energyUnit / lengthUnit`. |
| Result sizing | call `cpmdc_potential_result_size_for_force_input()` before writing a `PotentialResult`; a too-small buffer returns `ok == 0` and reports the required byte count. |

`ForceInput.pos` is a flat `natoms * 3` coordinate array and
`ForceInput.atmnrs` is a `natoms` atomic-number array. The optional
`ForceInput.box` is a row-major 3x3 cell. Unit strings default to `angstrom`
and `eV`.

## Build And Test Matrix

The default build does not need an OpenCPMD checkout. It builds the ABI,
parser, shared `libcpmdc`, and a deterministic reference evaluator used by
the test suite.

| Path | Commands | What it proves |
| --- | --- | --- |
| Default ABI/test build | `meson setup build -Dwith_tests=true`; `meson compile -C build`; `meson test -C build --print-errorlogs` | parser, C ABI, feature inventory, session contract, deterministic reference evaluator |
| cmocka render suites | `meson test -C build --suite cmocka --print-errorlogs` | Cap'n Proto decode and CPMD deck rendering |
| Session/result E2E suites | `meson test -C build --suite e2e --print-errorlogs` | `CPMDCSession`, topology checks, unit conversion, `PotentialResult` sizing |
| Documentation site | `pixi run -e docs docbld` | Org export, Doxygen/Doxyrest API pages, Sphinx HTML |
| OpenCPMD archive build | `meson setup build-cpmd -Dwith_cpmd=true -Dcpmd_root="$CPMD_ROOT" -Dwith_tests=true`; `meson compile -C build-cpmd`; `CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY meson test -C build-cpmd --print-errorlogs` | same C ABI linked to `libcpmd.a`, including live water parity when pseudopotentials are present |

```bash
meson setup build -Dwith_tests=true
meson compile -C build
meson test -C build --print-errorlogs
```

The reference evaluator exercises the same session, `ForceInput`, and
`PotentialResult` paths as the OpenCPMD build. The standalone
`libcpmdc_stub.a` target still reports `cpmdc_available() == 0`; it exists for
frontend links that need symbols but do not have an engine.

Build against a completed OpenCPMD tree with:

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

`cpmd_root` must contain `lib/libcpmd.a`. If runtime decks use library-style
pseudopotential names such as `O_MT_BLYP.psp`, put those files in the working
directory or set `CPMDC_PSEUDO_DIR`.

Common diagnostics:

| Symptom | Check |
| --- | --- |
| `with_cpmd requires ... lib/libcpmd.a` | build OpenCPMD first and pass the OpenCPMD tree as `-Dcpmd_root`, not the `cpmdc` tree. |
| pseudopotential open/read errors | set `CPMDC_PSEUDO_DIR` to the directory containing the named `.psp` files or use absolute paths in `atoms.pseudopotentials`. |
| topology-change errors | create a new `CPMDCSession` when atom count or ordered atomic numbers change. |
| `PotentialResult buffer too small` | resize to the byte count returned through `potential_result_capnp_size_bytes`. |
| unsupported unit errors | use the unit strings accepted by the schema helpers, such as `angstrom`, `bohr`, `eV`, and `hartree`. |

## Call Model

Long-running callers create one session from `CPMDParams`, then pass
`ForceInput` bytes for each geometry step:

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

`cpmdc_session_calculate_result()` converts energy and forces to
`ForceInput.energyUnit` and `energyUnit / lengthUnit` in the returned
`PotentialResult`. The lower-level force buffers and `CPMDCResult.energy_h`
remain in CPMD native units: Hartree and Hartree/Bohr.

The first accepted session evaluation fixes atom count and ordered atomic
numbers. Later steps may change coordinates, units, and the 3x3 cell. Topology
changes require a new session.

## CPMD Configuration

`CPMDParams` carries backend setup only. Geometry belongs in `ForceInput`.

Common scalar fields include:

| Field | Default | Meaning |
| --- | --- | --- |
| `functional` | `BLYP` | default DFT functional |
| `cutOffRy` | `70.0` | plane-wave cutoff in Rydberg |
| `charge` | `0` | system charge |
| `multiplicity` | `1` | spin multiplicity |
| `cpmdRoot` | empty | `CPMD_ROOT` hint for the embed layer |
| `enginePath` | empty | frontend hint used by loaders such as rgpot |

Long-tail CPMD input is represented with `inputSections`:

- `system` for `&SYSTEM` fields such as `CELL`, `CUTOFF`, `SCALE`, `CHARGE`
- `cpmd` for `&CPMD` controls such as `OPTIMIZE WAVEFUNCTION`, `MAXSTEP`,
  `OPTIMIZE GEOMETRY`, `MAXITER`, `CONVERGENCE GEOMETRY`, `EMASS`,
  `MOLECULAR DYNAMICS CP/BO/EH/PT/CLASSICAL`, `NOSE`, `BERENDSEN`,
  `LANGEVIN`, `ANNEALING`, `RATTLE`, `SHAKE`, `RESTART`,
  `PRINT`, `STORE`, `CENTER MOLECULE`, `DIIS`, `ODIIS`, `PCG`,
  `DIAGONALIZATION`, `FREE-ENERGY`, `QMMM`, `CDFT`, `PROPERTIES`,
  `RESTART WAVEFUNCTION`, `TRAJECTORY`
- `dft` for `&DFT` controls such as `FUNCTIONAL`, `LSD`, `GC-CUTOFF`,
  `XC_DRIVER`, `LIBXC`, `LR KERNEL`, `HFX`, `HFX-SCREENING`, `HUBBARD`,
  `ALPHA`, `BETA`, `OLDCODE`, `NEWCODE`, `CORRELATION`, `EXCHANGE`, and
  `BECKE88`
- `atoms` for pseudopotential entries and non-coordinate `&ATOMS` directives
- `generic`, `set`, and `raw` for CPMD sections that do not need a typed arm

Typed `atoms` sections group `ForceInput` coordinates into CPMD `&ATOMS`
entries by element symbol. Every atomic number in a geometry step needs a
matching pseudopotential entry. When `atoms` is omitted, built-in BLYP defaults
cover H and O only.

## rgpot Integration

`rgpot` uses the same schema. Its `PotentialConfig.cpmd` arm contains
`CPMDParams`; `calculate()` still receives `ForceInput`. The `CPMDPot`
frontend loads `libcpmdc` with `dlopen`, calls this C ABI, and returns the
same `PotentialResult` carrier used by other OmniPotentRPC potentials.

## Documentation

The Sphinx documentation is generated from Org sources under
`docs/orgmode/`. Start with:

- `docs/orgmode/tutorials/quickstart.org`
- `docs/orgmode/howto/embedding.org`
- `docs/orgmode/reference/cpmd-options.org`

Generated `.rst` files live under `docs/source/`.

## Version

0.1.0 - Cap'n Proto contract, stable C ABI, Fortran `iso_c_binding` shell,
shared `libcpmdc`, deterministic reference evaluator, OpenCPMD archive link
path, parser tests, session tests, and rgpot-compatible result carriers.
