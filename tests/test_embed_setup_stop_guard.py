#!/usr/bin/env python3
"""A stopgm in the cold setup list has to return before the next init.

The list used to run from init_fileopen through vdw_wf_alloc with
cpmd_set_up already true, and the stop test lived at the start of
embed_eval_energy_grad, after those calls. A fault in a later init
never reached the C handler. Each of those calls is now followed by
the stop test, and cpmd_set_up is set only after the list finishes.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

SETUP_CALLS = (
    "init_fileopen",
    "startpa",
    "New(bicanonicalCpmdInputConfig)",
    "init_pinf_pointers",
    "envir",
    "setcnst",
    "control",
    "dftin",
    "sysin",
    "setsc",
    "detsp",
    "mm_init",
    "ratom",
    "vdwin",
    "propin",
    "setsys",
    "New(bicanonicalCpmdConfig, bicanonicalCpmdInputConfig)",
    "genxc",
    "numpw",
    "rinit",
    "rinforce",
    "fft_init",
    "ortho_init",
    "initclust",
    "dg_init",
    "nosalloc",
    "exterp",
    "setbasis",
    "dqgalloc",
    "prnginit",
    "gle_alloc",
    "vdw_wf_alloc",
)


def routine_body(source: str, name: str) -> str:
    start = re.search(rf"^\s*SUBROUTINE\s+{name}\b", source, re.I | re.M)
    if not start:
        raise SystemExit(f"{name} not found in the Fortran bridge")
    end = re.compile(r"^\s*END\s+SUBROUTINE\b", re.I | re.M).search(source, start.end())
    if not end:
        raise SystemExit(f"{name} has no END SUBROUTINE")
    return source[start.start() : end.end()]


def code_lines(body: str) -> list[str]:
    lines = []
    for line in body.splitlines():
        code = line.split("!", 1)[0].strip()
        if code:
            lines.append(re.sub(r"\s+", " ", code))
    return lines


def call_indexes(lines: list[str]) -> list[int]:
    indexes = []
    cursor = 0
    for name in SETUP_CALLS:
        pattern = re.compile(rf"^CALL {re.escape(name)}(?:\s*\(.*\))?$", re.I)
        found = None
        for i in range(cursor, len(lines)):
            if pattern.match(lines[i]):
                found = i
                break
        if found is None:
            raise SystemExit(f"setup call {name} is missing from run_embed_scf")
        indexes.append(found)
        cursor = found + 1
    return indexes


def stop_checked(lines: list[str]) -> bool:
    return any(
        re.search(r"cold_setup_stopped|cpmdc_stop_code\s*\(", line, re.I)
        for line in lines
    )


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1]
    source = (root / "src" / "cpmd_embed_c_api.F90").read_text(encoding="utf-8")
    lines = code_lines(routine_body(source, "run_embed_scf"))
    indexes = call_indexes(lines)
    assignment = re.compile(r"cpmd_set_up\s*=\s*\.TRUE\.", re.I)
    assigned = [i for i, line in enumerate(lines) if assignment.search(line)]
    if not assigned:
        raise SystemExit("run_embed_scf never sets cpmd_set_up")
    if assigned[0] < indexes[0]:
        raise SystemExit(
            "cpmd_set_up is set before init_fileopen; a stop in the "
            "setup list is then a setup that already claims success"
        )
    if assigned[0] < indexes[-1]:
        raise SystemExit("cpmd_set_up is set before vdw_wf_alloc returns")
    for name, start, end in zip(SETUP_CALLS, indexes, indexes[1:]):
        if not stop_checked(lines[start + 1 : end]):
            raise SystemExit(
                f"no stop test between CALL {name} and the next setup call"
            )
    if not stop_checked(lines[indexes[-1] + 1 : assigned[0]]):
        raise SystemExit("no stop test after CALL vdw_wf_alloc before cpmd_set_up")
    stop_fn = re.search(
        r"LOGICAL FUNCTION cold_setup_stopped\(\).*?END FUNCTION cold_setup_stopped",
        source,
        re.I | re.S,
    )
    if not stop_fn:
        raise SystemExit("cold_setup_stopped is missing")
    stop_body = code_lines(stop_fn.group(0))
    try:
        agree = next(
            i
            for i, line in enumerate(stop_body)
            if "cpmdc_agree_cold_stop" in line.lower()
        )
    except StopIteration:
        raise SystemExit("cold setup does not agree the stop across ranks") from None
    leave = [
        i
        for i, line in enumerate(stop_body)
        if re.search(r"\bRETURN\b|cpmdc_embed_teardown\s*\(", line, re.I)
    ]
    if not leave or agree > min(leave):
        raise SystemExit(
            "cold setup leaves the list or frees the group before the stop is agreed"
        )
    print("cold setup reports stopgm before the next init")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
