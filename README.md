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

## Quick Build

The default build does not need an OpenCPMD checkout. It builds the ABI,
parser, shared `libcpmdc`, and a deterministic reference evaluator used by
the test suite.

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
meson setup build-cpmd \
  -Dwith_cpmd=true \
  -Dcpmd_root="$CPMD_ROOT" \
  -Dwith_tests=true
meson compile -C build-cpmd
```

`cpmd_root` must contain `lib/libcpmd.a`. If runtime decks use library-style
pseudopotential names such as `O_MT_BLYP.psp`, put those files in the working
directory or set `CPMDC_PSEUDO_DIR`.

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
  `RESTART WAVEFUNCTION`, `TRAJECTORY`
- `dft` for `&DFT` controls such as `FUNCTIONAL` and `LSD`
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
