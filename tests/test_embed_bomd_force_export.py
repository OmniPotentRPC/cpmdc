#!/usr/bin/env python3
"""Guard: embed PEF/BOMD path requests nuclear forces and copies fion."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EMBED = ROOT / "src" / "cpmd_embed_c_api.F90"


def test_bomd_force_export_path() -> None:
    text = EMBED.read_text(encoding="utf-8", errors="replace")
    m = re.search(
        r"SUBROUTINE embed_eval_energy_grad.*?END SUBROUTINE",
        text,
        re.S | re.I,
    )
    assert m, "embed_eval_energy_grad not found"
    body = m.group(0)
    assert "embed_set_need_forces" in body, "must request OpenCPMD embed_need_forces"
    assert "embed_set_need_forces(.TRUE.)" in body or "embed_set_need_forces(.true.)" in body.lower()
    write_at = body.lower().find("embed_set_write_files(.false.)")
    wf_at = body.lower().find("call wfopts")
    assert write_at != -1 and wf_at != -1 and write_at < wf_at
    assert "iprint_force" in body
    assert "ALLOCATED(fion)" in body
    assert "-fion" in body or "-fion(" in body.replace(" ", "")
    assert "wfopts" in body
    # success must not ignore forces path machinery
    assert "grad(idx)" in body


if __name__ == "__main__":
    test_bomd_force_export_path()
    print("ok")
