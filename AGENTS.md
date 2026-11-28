# Agent guide for cpmdc

This file orients automated agents and contributors to the repository layout,
build and test workflow, documentation pipeline, commit conventions, and how
cpmdc connects to rgpot.

## Repository layout

| Path | Role |
|------|------|
| `include/cpmdc.h` | Public C ABI: sessions, Cap'n Proto calculate paths, MPI bind and adopt |
| `include/cpmdc_features.h` | Runtime feature table (`cpmdc_feature_find`) |
| `include/cpmdc_restart.h` | RESTART file library (`libcpmdc_restart`) |
| `include/cpmdc_pwscale.h` | Opt-in plane-wave scaling scaffold, not part of `libcpmdc` |
| `include/cpmdc_csvr.h` | Opt-in CSVR thermostat for ionic velocities, not part of `libcpmdc` |
| `schema/Potentials.capnp` | Wire schema: `CPMDParams`, `ForceInput`, `PotentialResult` |
| `schema/inventory/cpmd_features.json` | Feature inventory mirrored in `src/cpmdc_features.c` |
| `src/cpmdc.c` | Sessions, topology guard, unit conversion, result encoding |
| `src/cpmdc_params.c` | Decode `CPMDParams` and render the CPMD INPUT deck |
| `src/cpmd_embed_c_api.F90` | Fortran bridge into OpenCPMD or the reference evaluator |
| `src/cpmdc_embed_host.c`, `src/cpmdc_memfd.c` | Host cwd, pseudopotentials, anonymous memfd deck |
| `python/cpmdc/` | Python `Session`, `Potential`, ASE calculator, drivers |
| `docs/orgmode/` | Authoritative documentation source (export to Sphinx) |
| `scripts/generate_llms_txt.py` | Generates `llms.txt` and `llms-full.txt` at the repo root |
| `examples/` | Runnable hosts, including `host_step.c` |

## Build

Default stub build (no OpenCPMD checkout):

```bash
meson setup build -Dwith_tests=true
meson compile -C build
```

With Pixi:

```bash
pixi run test-stub
```

OpenCPMD archive build:

```bash
meson setup build-cpmd -Dwith_tests=true -Dwith_cpmd=true -Dcpmd_root=/path/to/opencpmd-build
meson compile -C build-cpmd
```

Plane-wave scaffold (benchmarks, not installed):

```bash
meson setup build -Dwith_tests=true -Dwith_pwscale=true
```

## Test

Full default suite:

```bash
meson test -C build --print-errorlogs
```

Focused suites: `--suite docs`, `--suite inventory`, `--suite e2e`,
`--suite cmocka`.
Prose wrapping check: `pixi run snapper-check`.
It fails on files that would reformat and on fused lines, and reports the `long`
and `wrap` diagnostics; `SNAPPER_FAIL_ON=long,wrap` enforces those too.
Rewrap with `bash tools/check_snapper.sh --write`, which refuses symlinks and
paths outside the repository.
Meson build file check: `pixi run fmt-meson-check`; `pixi run fmt-meson` formats
in place with `meson.format`.

Documentation and API coverage checks include `reference-contract`,
`public-api-docs`, `cpmd-embedding-docs-coverage`, `llms-txt`, and
`docs-coverage-inventory`.

## Documentation build

Edit `docs/orgmode/` only.
Generated `docs/source/**/*.rst` is not committed.

```bash
pixi run -e docs docbld
```

Steps: `mkrst` (Org to RST), `doxybuild` (C API pages), `sphinxbld` (HTML to
`docs/build/`).

After changing Org sources or the public inventory, regenerate LLM index files:

```bash
python3 scripts/generate_llms_txt.py
```

`llms.txt` follows [llmstxt.org](https://llmstxt.org/); `llms-full.txt` holds
plain text from every Org page.

Coverage report (JSON on stdout):

```bash
python3 scripts/docs_coverage_inventory.py
```

## Commit messages

Use conventional prefixes: `feat`, `fix`, `chore`, `docs`, `refactor`, `test`.

Keep the subject at most 72 characters.
Wrap body lines at 72 characters.
Use backticks for code identifiers.
Do not use dashes as sentence punctuation in subjects or bodies.
Do not put pull request numbers, SHAs, or ticket ids in messages.
Do not use `!` or breaking-change framing unless the project explicitly requires
it.

Reference GitLab issues with `#number` when applicable.

## cpmdc and rgpot

[rgpot](https://github.com/OmniPotentRPC/rgpot) loads `libcpmdc.so` through
`CPMDPot`: it maps coordinates into `ForceInput`, keeps one `CPMDCSession` per
calculator group, and reads `PotentialResult` from
`cpmdc_session_calculate_result()`.

Environment variables tried in order:
`CPMDC_LIBRARY`, `RGPOT_CPMDC_ENGINE`, `RGPOT_CPMD_ENGINE`.

Python mirrors the contract with `cpmdc.Potential` for in-process rgpot-style
providers.
eOn uses rgpot via `[RgpotPot]`; see `docs/orgmode/tutorials/eon-rgpot.org`.

Cap'n Proto messages and units match the shared OmniPotentRPC `Potentials.capnp`
schema used across backends.
