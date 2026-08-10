#!/usr/bin/env python3
"""Validate the vendored OpenCPMD patches as portable unified diffs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PATCHES = {
    "opencpmd_keep_fion.patch": {"src/rwfopt_utils.mod.F90"},
    "opencpmd_warm_orbitals.patch": {"src/rwfopt_utils.mod.F90"},
    "opencpmd_converged_state.patch": {"src/updwf_utils.mod.F90"},
    "opencpmd_gfortran14_array_sections.patch": {
        "src/rv30_utils.mod.F90",
        "src/wv30_utils.mod.F90",
    },
}


def patch_numstat(repo: Path, patch: Path) -> dict[str, tuple[int, int]]:
    result = subprocess.run(
        ["git", "apply", "--numstat", str(patch)],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [line.split("\t", 2) for line in result.stdout.splitlines() if line]
    if any(len(row) != 3 for row in rows):
        raise AssertionError(f"{patch.name}: malformed numstat rows {rows!r}")
    return {target: (int(added), int(removed)) for added, removed, target in rows}


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    for name, expected_targets in PATCHES.items():
        patch = repo / "tools" / name
        stats = patch_numstat(repo, patch)
        if set(stats) != expected_targets:
            raise AssertionError(
                f"{name}: targets {set(stats)!r}, expected {expected_targets!r}"
            )
        for target, (added, removed) in stats.items():
            if added <= 0:
                raise AssertionError(f"{name}: {target} must add integration code")
            if removed <= 0:
                raise AssertionError(f"{name}: {target} must replace upstream code")
    print(f"OK: {len(PATCHES)} portable OpenCPMD patches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
