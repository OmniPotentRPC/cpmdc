#!/usr/bin/env python3
"""Guard: CPMDC_TIPRI_EVERY prints CPMD's timers, and embed_scf times whole evaluations."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EMBED = ROOT / "src" / "cpmd_embed_c_api.F90"
AGREE = ROOT / "src" / "cpmdc_orbital_limit.c"


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
    agree = wrapper.find("every = embed_tipri_every(failed_tipri)")
    failed = wrapper.find("IF (failed_tipri /= 0) THEN")
    start = wrapper.find("IF (timed) CALL tiset('embed_scf', isub)")
    call = wrapper.find("CALL run_embed_scf(")
    stop = wrapper.find("IF (timed) CALL tihalt('embed_scf', isub)")
    printed = wrapper.find("IF (MOD(evaluations, every) == 0) CALL tipri")
    assert -1 < agree < failed < start < call < stop < printed
    failed_body = wrapper[failed:start]
    assert "ok = 0_c_int" in failed_body
    assert "energy_h = 0.0_c_double" in failed_body
    assert "RETURN" in failed_body
    assert "CALL run_embed_scf(" not in failed_body
    assert wrapper.find("timed = timer_started") < start
    assert "IF (every > 0 .AND. timer_started) THEN" in wrapper[stop:printed]
    start = text.find("INTEGER FUNCTION embed_tipri_every")
    every = text[start:text.find("SUBROUTINE run_embed_scf", start)]
    assert "cpmdc_agree_tipri_every" in every
    assert "GET_ENVIRONMENT_VARIABLE" not in every
    assert "IF (rc /= 0) THEN" in every
    agree = AGREE.read_text(encoding="utf-8", errors="replace")
    assert "strlen(raw) > 16" in agree
    assert "MPI_Bcast(pair, 2, MPI_INT, 0, comm)" in agree
    assert "if (rank == 0)" in agree and "parse_tipri_every" in agree


if __name__ == "__main__":
    test_the_timer_wraps_each_evaluation_and_prints_only_when_asked()
    print("ok")
