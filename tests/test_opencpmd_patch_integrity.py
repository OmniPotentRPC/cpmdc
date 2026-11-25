#!/usr/bin/env python3
"""Validate the vendored OpenCPMD patches as portable unified diffs."""

from __future__ import annotations

import re
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
        "src/localize_utils.mod.F90",
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
        "src/stress_utils.mod.F90",
        "src/symtrz_utils.mod.F90",
        "src/tauofr_utils.mod.F90",
        "src/testex_utils.mod.F90",
        "src/updrho_utils.mod.F90",
        "src/vdw_utils.mod.F90",
        "src/vpsi_utils.mod.F90",
    ),
    "opencpmd_embed_extrapolate.patch": (
        "src/control_test_utils.mod.F90",
        "src/embed_ctrl.mod.F90",
        "src/rwfopt_utils.mod.F90",
    ),
    "opencpmd_embed_dipole.patch": (
        "src/ddipo_utils.mod.F90",
        "src/wfopts_utils.mod.F90",
    ),
    "opencpmd_constant_cutoff_precon.patch": "src/hesele_utils.mod.F90",
    "opencpmd_atomwf_kpts_ortho.patch": (
        "src/atomwf_utils.mod.F90",
        "src/gsortho_utils.mod.F90",
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
        "PUBLIC :: embed_abort_orbital_slot",
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
        "PUBLIC :: embed_set_orbital_slot_limit",
        "CPMDC_ORBITAL_SLOT_LIMIT",
        "PUBLIC :: embed_copy_orbital_slot",
        "PUBLIC :: embed_park_partial_orbital_slot",
        "SUBROUTINE embed_replace_key_store(c0, geo, key)",
        "LOGICAL, SAVE :: embed_store_converged = .FALSE.",
        "IF (.NOT. embed_slots(src)%converged) RETURN",
        "embed_slots(i)%converged = embed_store_converged",
        "embed_store_converged = embed_slots(i)%converged",
        "embed_store_converged = embed_slots(best)%converged",
        "ELSE IF (embed_warm_orbitals .AND. embed_have_orbitals) THEN",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    abort_begin = text.index("+  SUBROUTINE embed_abort_orbital_slot()")
    abort_end = text.index("+  END SUBROUTINE embed_abort_orbital_slot")
    abort = text[abort_begin:abort_end]
    for needle in (
        "DEALLOCATE(embed_pending_slots)",
        "DEALLOCATE(embed_pending_c0)",
        "DEALLOCATE(embed_pending_geo)",
        "embed_pending_action = .FALSE.",
        "embed_pending_index = 0",
        "embed_pending_grow = .FALSE.",
    ):
        if needle not in abort:
            raise AssertionError(f"{patch.name}: orbital abort missing {needle}")
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
    # localize2 keeps SAVE maps sized by the first basis; a second basis
    # overruns them in opeigr, so it resets like ddipo does.
    block = text.split("diff --git a/src/localize_utils.mod.F90", 1)
    if len(block) != 2:
        raise AssertionError(f"{patch.name}: localize_utils.mod.F90 is not patched")
    block = block[1].split("\ndiff --git ", 1)[0]
    for needle in (
        "+  USE embed_ctrl,                      ONLY: embed_teardowns",
        "+       IF (ALLOCATED(mapcol)) DEALLOCATE(mapcol)",
        "+       IF (ALLOCATED(mapful)) DEALLOCATE(mapful)",
        "+       IF (ALLOCATED(rotlsd)) DEALLOCATE(rotlsd)",
        "+       IF (ALLOCATED(rotmat)) DEALLOCATE(rotmat)",
        "+       IF (ALLOCATED(xyzmat)) DEALLOCATE(xyzmat)",
        "+       i_comp = 0",
        "+       ifirst = 0",
        "+       embed_seen = embed_teardowns",
    ):
        if needle not in block:
            raise AssertionError(f"{patch.name}: localize2 reset missing {needle}")
    # tistart executes STOP on a second call; the teardown keeps its state.
    if "USE timer" in text.split("+MODULE embed_teardown", 1)[1].split("END MODULE", 1)[0]:
        raise AssertionError(f"{patch.name}: embed_teardown must not touch the timers")


def assert_embed_extrapolate(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "LOGICAL, PUBLIC, SAVE :: embed_extrap_wfn = .FALSE.",
        "+    embed_extrap_wfn = cntl%textrap",
        "PUBLIC :: embed_extrapolation_used",
        "PUBLIC :: embed_clear_orbital_key",
        "SUBROUTINE embed_clear_orbital_key(drop_history)",
        "PUBLIC :: embed_extrapolation_count",
        "embed_have_key = .FALSE.",
        "embed_current_key = 0",
        "CALL extrapwf(0, cx, gam, embed_xhist, embed_xnow, embed_xcount, n, m)",
        "+    CALL embed_extrapolate_orbitals(c0, cx)",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    clear = text.split("+  SUBROUTINE embed_clear_orbital_key(drop_history)", 1)
    if len(clear) != 2 or (
        "+    IF (embed_have_key .OR. drop_history) CALL embed_break_history()"
    ) not in (
        clear[1].split("END SUBROUTINE embed_clear_orbital_key", 1)[0]
    ):
        raise AssertionError(
            f"{patch.name}: clearing a selected key must empty its "
            "extrapolation history"
        )
    stale = text.split(
        "ELSE IF (embed_warm_orbitals .AND. embed_have_orbitals) THEN", 1
    )
    if len(stale) != 2 or "CALL embed_break_history()" not in stale[1].split(
        "ELSE IF", 1
    )[0]:
        raise AssertionError(
            f"{patch.name}: an unconverged pass that already holds orbitals "
            "must break the extrapolation history"
        )
    # control_test still turns cntl%textrap off, so cpmd.x keeps its jobs.
    if "-       cntl%textrap=.FALSE." in text:
        raise AssertionError(f"{patch.name}: must not change cpmd.x's textrap")


def assert_embed_dipole(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "PUBLIC :: embed_get_dipole",
        "PUBLIC :: embed_set_dipole_snapshot",
        "embed_dipole_valid = .FALSE.",
        "CALL ddipo_embed(tau0,c0(:,:,1),dip_cm,c2,sc0(:,:,1),crge%n,center)",
        "PUBLIC :: ddipo_embed",
        "SUBROUTINE ddipo_embed(tau0,c0,cm,c2,sc0,nstate,center)",
        "embed_dipole = moment%dmom",
        ".NOT.isos1%tclust",
        "(cntl%tdipd .OR. prop1%dberry)",
        "USE embed_ctrl,                      ONLY: embed_teardowns",
        "IF (ALLOCATED(mapful)) DEALLOCATE(mapful)",
        "embed_seen = embed_teardowns",
        "wannier_recomput = wannier_recomput .OR. embed_force_recompute",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")
    # ddipo keeps the argument list that its other callers were compiled with.
    for needle in ("force_recompute)", "PRESENT(force_recompute)", "OPTIONAL"):
        if any(
            line.startswith("+") and needle in line and "embed_force_recompute" not in line
            for line in text.splitlines()
        ):
            raise AssertionError(f"{patch.name}: ddipo must not gain an argument ({needle})")


def assert_constant_cutoff_precon(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "USE prcp,                            ONLY: prcp_com",
        "IF (prcp_com%akin.GT.1.e-10_real_8) xskin=1._real_8/prcp_com%gskin",
        "(1._real_8+cp_erf((hg(j)-prcp_com%gckin)*xskin))",
        "-          vpp(j)=ftpi*hg(j)+vpp(j)",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")


def assert_atomwf_kpts_ortho(patch: Path) -> None:
    text = patch.read_text(encoding="utf-8")
    for needle in (
        "ALLOCATE(smat_c(nleft,MAX(atwp%nattot,nleft)),STAT=ierr)",
        "c0(:,atwp%nattot+1:,ikind),nleft,smat_c)",
        "-       CALL mp_sum(smat,2*np*n0,parai%allgrp)",
        "+       CALL mp_sum(smat,np*n0,parai%allgrp)",
    ):
        if needle not in text:
            raise AssertionError(f"{patch.name}: missing {needle}")


def assert_patch_application_order(repo: Path) -> None:
    script = (repo / "scripts" / "ci" / "build_opencpmd.sh").read_text(encoding="utf-8")
    match = re.search(r"for p in(.*?); do", script, flags=re.S)
    if not match:
        raise AssertionError("build_opencpmd.sh: missing patch application loop")
    tokens = [tok for tok in match.group(1).replace("\\", " ").split() if tok]
    try:
        rwfopt_at = tokens.index("embed_rwfopt")
        extrap_at = tokens.index("embed_extrapolate")
    except ValueError as exc:
        raise AssertionError(
            "build_opencpmd.sh: missing embed_rwfopt or embed_extrapolate in patch loop"
        ) from exc
    if rwfopt_at >= extrap_at:
        raise AssertionError(
            "build_opencpmd.sh: embed_rwfopt must apply before embed_extrapolate"
        )
    extrap = (repo / "tools" / "opencpmd_embed_extrapolate.patch").read_text(
        encoding="utf-8"
    )
    if "diff --git a/src/rwfopt_utils.mod.F90 b/src/rwfopt_utils.mod.F90" not in extrap:
        raise AssertionError(
            "opencpmd_embed_extrapolate.patch: must patch rwfopt_utils.mod.F90"
        )


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    assert_patch_application_order(repo)
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
        elif name == "opencpmd_embed_extrapolate.patch":
            assert_embed_extrapolate(patch)
        elif name == "opencpmd_embed_dipole.patch":
            assert_embed_dipole(patch)
        elif name == "opencpmd_constant_cutoff_precon.patch":
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
            assert_constant_cutoff_precon(patch)
        elif name == "opencpmd_atomwf_kpts_ortho.patch":
            if removed <= 0:
                raise AssertionError(f"{name}: patch must replace upstream code")
            assert_atomwf_kpts_ortho(patch)
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
