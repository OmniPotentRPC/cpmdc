# Changelog

## 0.1.0

- OpenCPMD driver in the ISO_C embed shell for working E2E single-point and
  multi-step session tests without OpenCPMD archives.
- Docs tree with Antics analytics, pixi/prek/CI wiring.
- Initial OmniPotentRPC `cpmdc` package modeled on `nwchemc`.
- Cap'n Proto `ForceInput` / `PotentialResult` / `CPMDParams` / `PotentialConfig`.
- Stable C ABI with session direct-call socket entry points.
- Fortran `iso_c_binding` embed shell (`cpmd_embed_c_api.f90`).
- Stub ABI and cmocka tests for parser + result sizing without OpenCPMD.
