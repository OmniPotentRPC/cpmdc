#!/usr/bin/env python3
"""Guard: CPMDC_TIPRI_EVERY prints CPMD's timers, and embed_scf times whole evaluations."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EMBED = ROOT / "src" / "cpmd_embed_c_api.F90"


def routine(text: str, name: str) -> str:
    m = re.search(
        rf"^\s*(?:INTEGER FUNCTION|SUBROUTINE) {name}\(.*?^\s*END (?:SUBROUTINE|FUNCTION)",
        text,
        re.S | re.M,
    )
    assert m, f"{name} not found"
    return m.group(0)


def test_the_timer_wraps_each_evaluation_and_prints_only_when_asked() -> None:
    text = EMBED.read_text(encoding="utf-8", errors="replace")
    assert text.count("CALL run_timed_embed_scf(") == 1
    assert text.count("CALL run_embed_scf(") == 1
    wrapper = routine(text, "run_timed_embed_scf")
    start = wrapper.find("IF (timed) CALL tiset('embed_scf', isub)")
    call = wrapper.find("CALL run_embed_scf(")
    stop = wrapper.find("IF (timed) CALL tihalt('embed_scf', isub)")
    assert -1 < start < call < stop, "embed_scf must time the whole evaluation, once CPMD's timer runs"
    assert wrapper.find("timed = timer_started") < start
    printed = wrapper.find("IF (MOD(evaluations, every) == 0) CALL tipri")
    assert printed > stop
    assert "IF (every > 0 .AND. timer_started) THEN" in wrapper[stop:printed]
    every = routine(text, "embed_tipri_every")
    assert "GET_ENVIRONMENT_VARIABLE('CPMDC_TIPRI_EVERY', v, STATUS=st)" in every
    assert "IF (ios == 0 .AND. n > 0) embed_tipri_every = n" in every


if __name__ == "__main__":
    test_the_timer_wraps_each_evaluation_and_prints_only_when_asked()
    print("ok")
