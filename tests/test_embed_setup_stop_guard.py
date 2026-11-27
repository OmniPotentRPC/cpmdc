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


def find_line(lines: list[str], pattern: str, start: int = 0) -> int:
    rx = re.compile(pattern, re.I)
    for i in range(start, len(lines)):
        if rx.search(lines[i]):
            return i
    raise SystemExit(f"missing {pattern}")


def leaves_before_agree(lines: list[str]) -> bool:
    agree_at = None
    for i, line in enumerate(lines):
        lowered = line.lower()
        if "cold_setup_stopped" in lowered or "cpmdc_agree_cold_stop" in lowered:
            agree_at = i
            break
    if agree_at is None:
        return True
    for line in lines[:agree_at]:
        if "cpmdc_embed_teardown" in line.lower() or re.search(r"\bRETURN\b", line):
            return True
    return False


def c_function(source: str, signature: str) -> str:
    start = 0
    while True:
        start = source.find(signature, start)
        if start < 0:
            raise SystemExit(f"{signature} not found")
        brace = source.find("{", start)
        semi = source.find(";", start)
        if brace >= 0 and (semi < 0 or brace < semi):
            break
        start += len(signature)
    depth = 0
    for i in range(brace, len(source)):
        if source[i] == "{":
            depth += 1
        elif source[i] == "}":
            depth -= 1
            if depth == 0:
                return source[brace : i + 1]
    raise SystemExit(f"{signature} has no closing brace")


def descriptor_problems(source: str) -> list[str]:
    problems = []
    match = re.search(
        r"INTEGER FUNCTION close_cold_deck\(\).*?END FUNCTION close_cold_deck",
        source,
        re.I | re.S,
    )
    if not match:
        return ["close_cold_deck is missing"]
    body = code_lines(match.group(0))
    seen_result = False
    cleared = False
    for line in body:
        if re.search(r"\brc\b", line) and re.search(r"==|/=", line):
            seen_result = True
        if re.search(r"\bmfd\s*=\s*-1\b", line):
            cleared = True
            if not seen_result:
                problems.append(
                    "close_cold_deck clears the descriptor before the close result is tested"
                )
    if not cleared:
        problems.append("close_cold_deck never clears the descriptor")
    if "-2" not in " ".join(body):
        problems.append(
            "close_cold_deck keeps a descriptor number that close has already released"
        )
    return problems


def memfd_problems(source: str) -> list[str]:
    body = c_function(source, "int cpmdc_memfd_close")
    problems = []
    if "-2" not in body:
        problems.append(
            "cpmdc_memfd_close does not report a descriptor that close released"
        )
    if "return -1" not in body:
        problems.append(
            "cpmdc_memfd_close does not report a descriptor that is still open"
        )
    return problems


def agreement_problems(fortran: str, library: str, lines: list[str]) -> list[str]:
    problems = []
    ratom = find_line(lines, r"^CALL ratom$")
    vdwin = find_line(lines, r"^CALL vdwin$", ratom)
    output_dir = find_line(lines, r"embed_use_output_dir", ratom)
    if output_dir > vdwin or leaves_before_agree(lines[output_dir:vdwin]):
        problems.append(
            "output directory failure tears down or returns before every rank agrees"
        )
    vdw = find_line(lines, r"^CALL vdw_wf_alloc$")
    closed = find_line(lines, r"close_cold_deck", vdw)
    setup = find_line(lines, r"cpmd_set_up\s*=\s*\.TRUE\.", closed)
    if leaves_before_agree(lines[closed:setup]):
        problems.append(
            "cold deck close failure tears down or returns before every rank agrees"
        )
    agree = c_function(library, "int cpmdc_agree_cold_stop")
    if "g_embed_detail" not in agree or "MPI_Bcast" not in agree:
        problems.append(
            "cold stop agreement does not publish the pending failure to every rank"
        )
    if "cpmdc_agree_cold_stop" not in fortran:
        problems.append("cold setup does not agree the stop across ranks")
    return problems


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
    library = (root / "src" / "cpmdc.c").read_text(encoding="utf-8")
    memfd = (root / "src" / "cpmdc_memfd.c").read_text(encoding="utf-8")
    problems = (
        agreement_problems(source, library, lines)
        + descriptor_problems(source)
        + memfd_problems(memfd)
    )
    if problems:
        raise SystemExit("\n".join(problems))
    print("cold setup reports stopgm before the next init")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
