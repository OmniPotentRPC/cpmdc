#!/usr/bin/env python3
"""Validate the vendored OpenCPMD patches as portable unified diffs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PATCHES = {
    "opencpmd_embed_geometry.patch": (
        "src/SOURCES",
        "src/embed_ctrl.mod.F90",
        "src/initrun_driver.mod.F90",
    ),
    "opencpmd_embed_rwfopt.patch": "src/rwfopt_utils.mod.F90",
    "opencpmd_converged_state.patch": "src/updwf_utils.mod.F90",
    "opencpmd_kpoints_inputfile.patch": "src/rkpnt_utils.mod.F90",
    "opencpmd_stopgm_return.patch": "src/error_handling.mod.F90",
    "opencpmd_tistopgm.patch": "src/timer.mod.F90",
    "opencpmd_c_mem_addrs.patch": "src/c_mem_utils.c",
}

# The stopgm patch only inserts cpmd_stopgm_hook. It deletes no upstream line.
INSERT_ONLY = {"opencpmd_stopgm_return.patch"}


def patch_numstat(repo: Path, patch: Path) -> list[tuple[int, int, str]]:
    result = subprocess.run(
        ["git", "apply", "--numstat", str(patch)],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    rows = [line.split("\t", 2) for line in result.stdout.splitlines() if line]
    parsed: list[tuple[int, int, str]] = []
    for row in rows:
        if len(row) != 3:
            raise AssertionError(f"{patch.name}: bad numstat row {row!r}")
        added, removed, target = row
        parsed.append((int(added), int(removed), target))
    if not parsed:
        raise AssertionError(f"{patch.name}: empty numstat")
    return parsed


def assert_embed_rwfopt(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "PUBLIC :: embed_set_warm_orbitals",
        "PUBLIC :: embed_set_need_forces",
        "PUBLIC :: embed_reset_warm_orbitals",
        "PUBLIC :: embed_set_write_files",
        "CALL embed_reset_warm_orbitals",
        "IF (embed_warm_orbitals .AND. ropt_mod%convwf) THEN",
        "embed_have_partial = .TRUE.",
        "cntl%pcgmin .AND. .NOT. cntl%diis",
        "IF (.NOT.embed_need_forces) THEN",
        "IF (ALLOCATED(fion)) DEALLOCATE(fion)",
        "CALL stopgm('embed_save_orbitals', 'allocation problem'",
        "CALL stopgm('embed_restore_orbitals'",
        "saved orbitals do not match c0; call embed_reset_warm_orbitals",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    for symbol in (
        "cpmdc_set_warm_orbitals",
        "cpmdc_set_need_forces",
        "cpmdc_reset_warm_orbitals",
    ):
        if symbol in text:
            raise AssertionError(f"{patch.name}: still names {symbol}")


def assert_converged_state(patch: Path) -> None:
    lines = patch.read_text(encoding="utf-8").splitlines()
    gemax = [
        i
        for i, line in enumerate(lines)
        if line.startswith("+") and "gemax=2.0_real_8*cntr%tolog" in line
    ]
    check = [
        i
        for i, line in enumerate(lines)
        if line.startswith(" ") and "CALL tol_chk_cnvgrad" in line
    ]
    if len(gemax) != 1 or len(check) != 1 or gemax[0] > check[0]:
        raise AssertionError(
            f"{patch.name}: steepest-descent gemax reset must precede the convergence check"
        )
    if not any(
        line.startswith("+")
        and "Steepest descent with iproj <= 1 never converges on gemax." in line
        for line in lines
    ):
        raise AssertionError(f"{patch.name}: missing the steepest-descent gemax note")


def assert_tistopgm(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    if "trace_depth is HUGE(0) until tistart" not in text:
        raise AssertionError(f"{patch.name}: missing the pre-tistart note")
    needle = (
        "tname%trace_depth.GT.0 .AND. "
        "tname%trace_depth.LE.SIZE(tname%trace_names)"
    )
    if needle not in text:
        raise AssertionError(f"{patch.name}: missing the pre-tistart stack guard")


def assert_stopgm_hook(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    if "NAME='cpmd_stopgm_hook'" not in text:
        raise AssertionError(f"{patch.name}: must publish cpmd_stopgm_hook")
    if "INTEGER(c_int), VALUE :: code" not in text:
        raise AssertionError(f"{patch.name}: hook must take the stop code by value")
    if "its results and set CPMD up again" not in text:
        raise AssertionError(f"{patch.name}: missing the invalid-result contract")
    if "or abort the others itself" not in text:
        raise AssertionError(f"{patch.name}: missing the multi-rank contract")
    for symbol in ("cpmdc_embed_catch", "cpmdc_note_stop"):
        if symbol in text:
            raise AssertionError(
                f"{patch.name}: references {symbol}, so a standalone cpmd.x cannot link"
            )


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    for name, expected_target in PATCHES.items():
        patch = repo / "tools" / name
        rows = patch_numstat(repo, patch)
        added = sum(row[0] for row in rows)
        removed = sum(row[1] for row in rows)
        targets = tuple(row[2] for row in rows)
        if added <= 0:
            raise AssertionError(f"{name}: patch must add embed integration code")
        if name in INSERT_ONLY:
            if removed != 0:
                raise AssertionError(
                    f"{name}: hook patch must only add lines, removed {removed}"
                )
            assert_stopgm_hook(patch)
        elif name == "opencpmd_embed_rwfopt.patch":
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
            assert_embed_rwfopt(patch)
        elif name == "opencpmd_embed_geometry.patch":
            if "embed_write_files" not in patch.read_text(encoding="utf-8"):
                raise AssertionError(f"{name}: missing embed_write_files")
        elif name == "opencpmd_converged_state.patch":
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
            assert_converged_state(patch)
        elif name == "opencpmd_tistopgm.patch":
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
            assert_tistopgm(patch)
        elif name == "opencpmd_c_mem_addrs.patch":
            text = patch.read_text(encoding="utf-8")
            if "return (size_t) pp;" not in text:
                raise AssertionError(f"{name}: missing the address cast")
            if "return pp;" not in text:
                raise AssertionError(f"{name}: missing the stock return")
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
        elif removed <= 0:
            raise AssertionError(f"{name}: patch must replace upstream code")
        expected = expected_target if isinstance(expected_target, tuple) else (expected_target,)
        if targets != expected:
            raise AssertionError(
                f"{name}: targets {targets!r}, expected {expected!r}"
            )
    print(f"OK: {len(PATCHES)} portable OpenCPMD patches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
