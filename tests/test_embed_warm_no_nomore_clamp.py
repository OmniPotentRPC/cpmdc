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
    assert "embed_set_warm_orbitals" in body
    assert "wfopts" in body


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


if __name__ == "__main__":
    test_no_warm_nomore_iter_clamp()
    test_calculator_once_flag_and_stop_reset()
    print("ok")
