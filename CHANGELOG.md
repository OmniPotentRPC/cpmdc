# Changelog

<!-- towncrier release notes start -->

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
