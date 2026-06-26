# cpmdc

Stable C ABI for embedding [OpenCPMD](https://github.com/OpenCPMD/CPMD) from
language-neutral Cap'n Proto `CPMDParams` messages.

`cpmdc` mirrors the OmniPotentRPC `nwchemc` layout: a C ABI layer, modern
Fortran `iso_c_binding` / `iso_fortran_env` bridge code, and Cap'n Proto
`ForceInput` / `PotentialResult` carriers for the direct-call "socket" style
used by rgpot (session + per-step messages, no parallel user config language).

The public ABI does not expose C++ or Rust types:

```c
int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes);
CPMDCResult cpmdc_energy_gradient(
    int n_atoms, const double *positions_ang, const int *atomic_numbers,
    const void *params_capnp, size_t params_capnp_size_bytes,
    double *grad_h_bohr);
CPMDCSession *cpmdc_session_create(
    const void *params_capnp, size_t params_capnp_size_bytes);
size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes);
CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);
CPMDCResult cpmdc_calculate_result(
    const void *params_capnp, size_t params_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);
void cpmdc_session_destroy(CPMDCSession *session);
```

`params_capnp` is an unpacked flat Cap'n Proto message whose root is
`CPMDParams` from `schema/Potentials.capnp`. Geometry for each evaluation is a
serialized `ForceInput`. Long-running callers create one `CPMDCSession` from
`CPMDParams`, then pass a `ForceInput` per step and receive `PotentialResult`
bytes (energy/forces converted to `ForceInput.energyUnit` and
`energyUnit / lengthUnit`). Native evaluation units remain Hartree and
Hartree/Bohr on `CPMDCResult.energy_h` and force buffers.

The first accepted session evaluation fixes atom count and ordered atomic
numbers; later steps may change coordinates, units, and cell vectors. Topology
changes require a new session.

## Cap'n Proto (0.1 focus)

| Type | Role |
| --- | --- |
| `ForceInput` | Per-step geometry + units (shared with nwchemc / rgpot) |
| `PotentialResult` | Per-step energy + forces |
| `CPMDParams` | Method/backend setup (functional, cutoff Ry, sections, …) |
| `CPMDInputSection` | Structured `&CPMD` / `&SYSTEM` / `&DFT` / `&ATOMS` / generic / raw |
| `PotentialConfig` | Tagged union with `cpmd @2 :CPMDParams` arm for rgpot |

`cpmdc_params_render_input_deck()` turns `CPMDParams` into a CPMD-style
`&SECTION` … `&END` deck for embed configuration and debugging.

## Build

Frontend / stub build without OpenCPMD:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
meson test -C build --print-errorlogs
```

Optional embed shell shared library (ISO_C ABI present; evaluation still
requires a future OpenCPMD link-up):

```bash
meson setup build-cpmd -Dwith_cpmd=true -Dcpmd_root=/path/to/OpenCPMD/CPMD
meson compile -C build-cpmd
```

## Layout

```
include/cpmdc.h              Public C ABI
src/cpmdc.c                  Session + Cap'n Proto socket path
src/cpmdc_stub.c             Stub ABI when embed is not linked
src/cpmdc_params.c           Cap'n Proto parse / deck render / result write
src/cpmd_embed_c_api.f90     iso_c_binding / iso_fortran_env bridge
schema/Potentials.capnp      Wire schema
tests/                       cmocka + encoded Cap'n Proto fixtures
```

## Version

0.1.0 — Cap'n Proto contract, C ABI, Fortran ISO_C shell, stub + parser tests.
Full OpenCPMD energy/force evaluation is the next embed milestone.
