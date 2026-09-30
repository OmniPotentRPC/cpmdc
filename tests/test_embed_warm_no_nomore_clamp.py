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
    assert "continuing with PCG MINIMIZE" in body
    assert "cntl%pcgmin = .TRUE." in body
    assert "cpmdc_stop_code() == 0_c_int" in body
    assert "embed_reset_warm_orbitals" not in body


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
    assert re.search(
        r"if \(stop != 0\) \{\s*.*?cpmdc_embed_abort_other_ranks\(\);",
        host,
        re.S,
    ), "a stopgm return on more than one rank must abort the others"
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


def test_cold_start_uses_pcg_minimize() -> None:
    text = EMBED.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"SUBROUTINE embed_eval_energy_grad.*?END SUBROUTINE",
        text,
        re.S | re.I,
    )
    assert match, "embed_eval_energy_grad not found"
    body = match.group(0)
    first = body.index("CALL wfopts")
    pre = body[:first]
    post = body[first:]
    assert "image%cfg_warm_steps <= 0_c_int" in pre
    assert "cntl%diis = .FALSE." in pre
    assert "cntl%pcg = .TRUE." in pre
    assert "cntl%pcgmin = .TRUE." in pre
    assert "used_opt = 'PCG MINIMIZE'" in pre
    assert "cold start; using PCG MINIMIZE" in pre
    assert "wf_optimiser_label" in pre
    assert "continuing with PCG MINIMIZE" in post
    assert "used_opt = 'PCG MINIMIZE'" in post
    assert "IF (cold_start) THEN" in post
    assert "ropt_mod%spcg = was_spcg" in post
    assert "CALL note_wf_optimiser" in post
    assert body.count("CALL wfopts") == 2
    assert "FUNCTION wf_optimiser_label" in text
    assert "label = 'ODIIS'" in text
    assert "label = 'PCG MINIMIZE'" in text
    assert "deck_contains_token(image, 'ODIIS')" in text
    assert "CHARACTER(KIND=c_char) :: name(32)" in text


if __name__ == "__main__":
    test_no_warm_nomore_iter_clamp()
    test_calculator_once_flag_and_stop_reset()
    test_cold_start_uses_pcg_minimize()
    print("ok")
