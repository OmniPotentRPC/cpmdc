#!/usr/bin/env python3
"""Validate the vendored OpenCPMD patches as portable unified diffs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PATCHES = {
    "opencpmd_keep_fion.patch": "src/rwfopt_utils.mod.F90",
    "opencpmd_warm_orbitals.patch": "src/rwfopt_utils.mod.F90",
    "opencpmd_converged_state.patch": "src/updwf_utils.mod.F90",
    "opencpmd_stopgm_return.patch": "src/error_handling.mod.F90",
}

# The stopgm patch only inserts cpmd_stopgm_hook. It deletes no upstream line.
INSERT_ONLY = {"opencpmd_stopgm_return.patch"}


def patch_numstat(repo: Path, patch: Path) -> tuple[int, int, str]:
    result = subprocess.run(
        ["git", "apply", "--numstat", str(patch)],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [line.split("\t", 2) for line in result.stdout.splitlines() if line]
    if len(rows) != 1 or len(rows[0]) != 3:
        raise AssertionError(f"{patch.name}: expected one numstat row, got {rows!r}")
    added, removed, target = rows[0]
    return int(added), int(removed), target


def assert_stopgm_hook(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    if "NAME='cpmd_stopgm_hook'" not in text:
        raise AssertionError(f"{patch.name}: must publish cpmd_stopgm_hook")
    if "INTEGER(c_int), VALUE :: code" not in text:
        raise AssertionError(f"{patch.name}: hook must take the stop code by value")
    for symbol in ("cpmdc_embed_catch", "cpmdc_note_stop"):
        if symbol in text:
            raise AssertionError(
                f"{patch.name}: references {symbol}, so a standalone cpmd.x cannot link"
            )


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    for name, expected_target in PATCHES.items():
        patch = repo / "tools" / name
        added, removed, target = patch_numstat(repo, patch)
        if added <= 0:
            raise AssertionError(f"{name}: patch must add embed integration code")
        if name in INSERT_ONLY:
            if removed != 0:
                raise AssertionError(
                    f"{name}: hook patch must only add lines, removed {removed}"
                )
            assert_stopgm_hook(patch)
        elif removed <= 0:
            raise AssertionError(f"{name}: patch must replace upstream code")
        if target != expected_target:
            raise AssertionError(
                f"{name}: target {target!r}, expected {expected_target!r}"
            )
    print(f"OK: {len(PATCHES)} portable OpenCPMD patches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
