#!/usr/bin/env python3
"""Validate the vendored OpenCPMD patches as portable unified diffs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PATCHES = {
    "opencpmd_embed_rinitwf.patch": "src/rinitwf_driver.mod.F90",
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
    "opencpmd_embed_teardown.patch": (
        "src/SOURCES",
        "src/broyden_utils.mod.F90",
        "src/calc_alm_utils.mod.F90",
        "src/chksym_utils.mod.F90",
        "src/detdof_utils.mod.F90",
        "src/drhov_utils.mod.F90",
        "src/ehpsi_utils.mod.F90",
        "src/embed_ctrl.mod.F90",
        "src/embed_teardown.mod.F90",
        "src/fftnew_utils.mod.F90",
        "src/forcedr_driver.mod.F90",
        "src/hfx_utils.mod.F90",
        "src/initclust_utils.mod.F90",
        "src/k_odiis_utils.mod.F90",
        "src/k_pcgrad_utils.mod.F90",
        "src/mixing_g_utils.mod.F90",
        "src/mixing_r_utils.mod.F90",
        "src/moverho_utils.mod.F90",
        "src/nlccset_utils.mod.F90",
        "src/numpw_utils.mod.F90",
        "src/phfac_utils.mod.F90",
        "src/pw_hfx.mod.F90",
        "src/qvan2_utils.mod.F90",
        "src/recpnew_utils.mod.F90",
        "src/rhodiis_utils.mod.F90",
        "src/rnlsmd_utils.mod.F90",
        "src/rpiiint_utils.mod.F90",
        "src/rwswap_utils.mod.F90",
        "src/setbasis_utils.mod.F90",
        "src/symtrz_utils.mod.F90",
        "src/tauofr_utils.mod.F90",
        "src/testex_utils.mod.F90",
        "src/updrho_utils.mod.F90",
        "src/vdw_utils.mod.F90",
        "src/vpsi_utils.mod.F90",
    ),
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


def assert_embed_rinitwf(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "LOGICAL, SAVE, PUBLIC :: embed_have_orbitals = .FALSE.",
        "CALL phfac(tau0)",
        "IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)",
        "CALL tihalt(procedureN,isub)",
        "RETURN",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    if "cntl%embed_warm_orbitals" in text:
        raise AssertionError(
            f"{patch.name}: must not use cntl%embed_warm_orbitals"
        )


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
        "USE rinitwf_driver,                  ONLY: embed_have_orbitals",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    if "LOGICAL, SAVE :: embed_have_orbitals" in text:
        raise AssertionError(
            f"{patch.name}: embed_have_orbitals must live in rinitwf_driver"
        )
    lines = text.splitlines()
    init_at = [
        i
        for i, line in enumerate(lines)
        if line.startswith(" ")
        and "CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)" in line
    ]
    restore_at = [
        i
        for i, line in enumerate(lines)
        if line.startswith("+") and "CALL embed_restore_orbitals(c0)" in line
    ]
    if len(init_at) != 1 or len(restore_at) != 2 or not (
        restore_at[0] < init_at[0] < restore_at[1]
    ):
        raise AssertionError(
            f"{patch.name}: stored c0 must be restored before and after initrun"
        )
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


def assert_embed_teardown(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "INTEGER, PUBLIC, SAVE :: embed_teardowns = 0",
        "PUBLIC :: embed_teardown_state",
        "CALL embed_reset_warm_orbitals()",
        "CALL free_cp_groups()",
        "embed_teardowns = embed_teardowns + 1",
        "embed_teardown.mod.F90",
        "IF (ALLOCATED(crge%f)) DEALLOCATE(crge%f)",
        "IF (ALLOCATED(eigr)) DEALLOCATE(eigr)",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    # Each guarded routine resets its first-call state once per teardown.
    guards = text.count("+    IF (embed_seen /= embed_teardowns) THEN") + text.count(
        "+  IF (embed_seen /= embed_teardowns) THEN"
    )
    if guards < 30:
        raise AssertionError(f"{patch.name}: {guards} first-call guards, expected 30 or more")
    # tistart executes STOP on a second call; the teardown keeps its state.
    if "USE timer" in text.split("+MODULE embed_teardown", 1)[1].split("END MODULE", 1)[0]:
        raise AssertionError(f"{patch.name}: embed_teardown must not touch the timers")


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
        elif name == "opencpmd_embed_rinitwf.patch":
            if removed != 0:
                raise AssertionError(
                    f"{name}: guess skip must only add lines, removed {removed}"
                )
            assert_embed_rinitwf(patch)
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
        elif name == "opencpmd_embed_teardown.patch":
            assert_embed_teardown(patch)
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
