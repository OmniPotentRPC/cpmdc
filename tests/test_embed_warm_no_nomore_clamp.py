#!/usr/bin/env python3
"""Assert shipped embed warm path does not clamp nomore_iter to 0 or 1."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EMBED = ROOT / "src" / "cpmd_embed_c_api.F90"


def test_no_warm_nomore_iter_clamp() -> None:
    text = EMBED.read_text(encoding="utf-8", errors="replace")
    # Live warm eval must not assign nomore_iter to 0 or 1.
    # Allow only non-clamp assignments (e.g. 40 in dead helpers) — forbid 0/1 entirely
    # on the warm-eval subroutine body.
    m = re.search(
        r"SUBROUTINE embed_eval_energy_grad.*?END SUBROUTINE",
        text,
        re.S | re.I,
    )
    assert m, "embed_eval_energy_grad not found"
    body = m.group(0)
    bad = re.findall(r"nomore_iter\s*=\s*([01])\b", body)
    assert not bad, f"warm eval clamps nomore_iter to {bad} (not physical SCF)"
    assert "embed_set_tau0_from_pos" in body
    assert "phfac" in body
    assert "embed_set_warm_orbitals(.TRUE.)" in body
    assert "embed_set_warm_orbitals(.FALSE.)" not in body
    assert body.count("CALL wfopts") == 2
    assert body.count("CALL note_scf_steps(iteropt%nfi)") == 2
    assert "CPMDC_SCF_STEPS" in text
    assert "continuing with PCG MINIMIZE" in body
    assert "cntl%pcgmin = .TRUE." in body
    assert "cpmdc_stop_code() == 0_c_int" in body
    assert "embed_reset_warm_orbitals" not in body
    assert "embed_basis_latched" in body


def test_calculator_once_flag_and_stop_reset() -> None:
    embed = EMBED.read_text(encoding="utf-8", errors="replace")
    assert "mp_comm_set" not in embed
    assert "LOGICAL, SAVE :: embed_calculator_bound = .FALSE." in embed
    bind = re.search(
        r"FUNCTION cpmdc_embed_bind_calculator.*?END FUNCTION",
        embed,
        re.S,
    )
    assert bind, "cpmdc_embed_bind_calculator not found"
    body = bind.group(0)
    assert "USE mp_interface, ONLY: mp_comm_world" in body
    assert "IF (embed_calculator_bound) RETURN" in body
    assert "embed_calculator_bound = .TRUE." in body
    host = (ROOT / "src" / "cpmdc.c").read_text(encoding="utf-8", errors="replace")
    assert "if (stop != 0)" in host
    assert re.search(
        r"if \(stop != 0\) \{\s*.*?cpmdc_embed_reset_state\(image\);",
        host,
        re.S,
    ), "a stopgm return must reset the image before the next call"
    assert "cpmdc_embed_abort_other_ranks" not in host
    assert "MPI_Abort" not in host
    assert re.search(
        r'if \(stop != 0\) \{\s*.*?return fail_msg\("CPMD stopgm during embed SCF"\);',
        host,
        re.S,
    ), "a stopgm return must come back as an error the caller handles"
    scf = re.search(
        r"SUBROUTINE run_embed_scf.*?END SUBROUTINE",
        embed,
        re.S,
    )
    assert scf, "run_embed_scf not found"
    scf_body = scf.group(0)
    assert "embed_basis_changed" in scf_body
    assert "CALL embed_reset_warm_orbitals()" in scf_body
    assert "latch_embed_basis" in scf_body
    assert "embed_basis_latched" in embed
    # A failed warm SCF increments nothing and does not clear the counter.
    assert "IF (ok /= 0_c_int) THEN" in scf_body
    assert "image%cfg_warm_steps = 0" in scf_body


def test_new_session_keeps_matching_orbitals() -> None:
    embed = EMBED.read_text(encoding="utf-8", errors="replace")
    host = (ROOT / "src" / "cpmdc.c").read_text(encoding="utf-8", errors="replace")
    apply = re.search(
        r"static int embed_apply_from_wire\(.*?^\}",
        host,
        re.S | re.M,
    )
    assert apply, "embed_apply_from_wire not found"
    apply_body = apply.group(0)
    assert "cpmdc_embed_detach_image" in apply_body
    assert "cpmdc_embed_reset_state" not in apply_body
    basis = re.search(
        r"LOGICAL FUNCTION embed_basis_changed.*?END FUNCTION",
        embed,
        re.S,
    )
    assert basis, "embed_basis_changed not found"
    basis_body = basis.group(0)
    assert "same_element_counts" in basis_body
    assert "embed_saved_functional" in basis_body
    assert "embed_saved_deck" in basis_body
    assert "embed_saved_cell" in basis_body
    assert "embed_saved_z(i)" not in basis_body
    assert "FUNCTION cpmdc_embed_detach_image" in embed
    assert "CALL embed_reset_warm_orbitals" not in re.search(
        r"FUNCTION cpmdc_embed_detach_image.*?END FUNCTION",
        embed,
        re.S,
    ).group(0)


if __name__ == "__main__":
    test_no_warm_nomore_iter_clamp()
    test_calculator_once_flag_and_stop_reset()
    test_new_session_keeps_matching_orbitals()
    print("ok")
