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
    assert "cpmdc_set_warm_orbitals" in body
    assert "wfopts" in body


if __name__ == "__main__":
    test_no_warm_nomore_iter_clamp()
    print("ok")
