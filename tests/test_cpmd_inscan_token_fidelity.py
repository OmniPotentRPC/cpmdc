#!/usr/bin/env python3
"""Assert shipped CPMD render tokens appear in OpenCPMD input-scan sources.

Epic cpmdc-3ba9: typed-section spellings must match OpenCPMD keyword scans.
Reads tokens from shipped cpmdc_params.c (render path) and checks each
required token string against OpenCPMD sources under CPMD_SRC (or
--cpmd-src).

Does not re-implement render; it only greps the shipped renderer source for
the token strings and greps OpenCPMD for the same literals (keyword tables /
comments / keyword_contains).
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

# Tokens the typed-section epic requires (from shipped render strings).
# Literals that must appear in shipped render AND OpenCPMD input scanners.
REQUIRED_TOKENS = [
    "CONSTRAINTS",
    "END CONSTRAINTS",
    "ISOTOPE",
    "VELOCITIES",
    "END VELOCITIES",
    "DUMMY ATOMS",
    "HUBBARD",  # OpenCPMD keyword_contains('HUBBARD')
    "HFX_BLOCK_SIZE",
    "HFX_DISTRIBUTION",
    "HFX",  # screening / hybrid path
    # Typed long-tail sections (3ba9 next slice)
    "EMPIRICAL CORRECTION",
    "GRIMME",
    "DIPOLE MOMENT",
    "LOCALIZE",
    "TROTTER DIMENSION",
    "REPLICA NUMBER",
    "TAMM-DANCOFF",
    "MAXSTEP",  # linres / general
    "HTHRS",
]


def collect_opencpmd_blob(src_root: Path) -> str:
    parts: list[str] = []
    for p in src_root.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in {".f", ".f90", ".inc", ".h", ".c"} and not p.name.endswith(
            ".mod.F90"
        ):
            # OpenCPMD uses .mod.F90 heavily
            if ".F90" not in p.name and ".F" not in p.name:
                continue
        try:
            parts.append(p.read_text(encoding="utf-8", errors="replace"))
        except OSError:
            pass
    return "\n".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cpmdc-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="cpmdc repo root (for src/cpmdc_params.c)",
    )
    ap.add_argument(
        "--cpmd-src",
        type=Path,
        default=None,
        help="OpenCPMD src tree (default: env CPMD_SRC or common terra path)",
    )
    args = ap.parse_args()
    params_c = args.cpmdc_root / "src" / "cpmdc_params.c"
    if not params_c.is_file():
        print(f"missing shipped renderer: {params_c}", file=sys.stderr)
        return 2
    render_src = params_c.read_text(encoding="utf-8", errors="replace")

    cpmd_src = args.cpmd_src
    if cpmd_src is None:
        env = os.environ.get("CPMD_SRC") or os.environ.get("OPENCPMD_SRC")
        if env:
            cpmd_src = Path(env)
        else:
            cand = Path.home() / "work/OpenCPMD/CPMD/src"
            if cand.is_dir():
                cpmd_src = cand
    if cpmd_src is None or not Path(cpmd_src).is_dir():
        print(
            "SKIP: OpenCPMD src not found (set CPMD_SRC or --cpmd-src). "
            "Token presence in shipped renderer still checked.",
            file=sys.stderr,
        )
        missing_render = [t for t in REQUIRED_TOKENS if t not in render_src]
        if missing_render:
            print("FAIL: tokens missing from shipped cpmdc_params.c:", missing_render)
            return 1
        print("OK: all required tokens present in shipped renderer (inscan skipped)")
        return 0

    opencpmd = collect_opencpmd_blob(Path(cpmd_src))
    failed: list[str] = []
    for tok in REQUIRED_TOKENS:
        in_render = tok in render_src
        # OpenCPMD often matches via keyword_contains fragments; allow token
        # without END- prefix as alternative for END lines.
        key = tok.replace("END ", "").strip()
        in_inscan = tok in opencpmd or key in opencpmd
        status = "OK" if (in_render and in_inscan) else "FAIL"
        print(f"{status}: {tok!r} render={in_render} opencpmd={in_inscan}")
        if not in_render or not in_inscan:
            failed.append(tok)
    if failed:
        print("FAILED tokens:", failed, file=sys.stderr)
        return 1
    print(f"OK: {len(REQUIRED_TOKENS)} tokens match shipped render + OpenCPMD src")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
