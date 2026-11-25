#!/usr/bin/env python3
"""Check that session charges alone own the external energy sum.

OpenCPMD already folds the deck's EXTERNAL POTENTIAL ADD energy into the
total. Adding ``ener_com%eext`` again when no session charge is active
counts it twice, so the sum in ``embed_eval_energy_grad`` has to sit inside
``IF (embed_external_count > 0)``. A deck that needs this energy also needs
a grid file on the host, so the guard is checked on the Fortran source.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def routine_body(source: str, name: str) -> str:
    start = re.search(rf"^\s*SUBROUTINE\s+{name}\b", source, re.I | re.M)
    if not start:
        raise SystemExit(f"{name} not found in the Fortran bridge")
    end = re.compile(r"^\s*END\s+SUBROUTINE\b", re.I | re.M).search(
        source, start.end()
    )
    if not end:
        raise SystemExit(f"{name} has no END SUBROUTINE")
    return source[start.start() : end.end()]


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1]
    source = (root / "src" / "cpmd_embed_c_api.F90").read_text(encoding="utf-8")
    body = routine_body(source, "embed_eval_energy_grad")
    lines = [line.split("!", 1)[0] for line in body.splitlines()]

    adds = [
        i
        for i, line in enumerate(lines)
        if re.search(r"ener_com%etot\s*=\s*ener_com%etot\s*\+", line, re.I)
        and re.search(r"ener_com%eext", line, re.I)
    ]
    if not adds:
        raise SystemExit("the eext addition to etot is missing")

    for i in adds:
        depth = 0
        guarded = False
        for line in lines[: i + 1]:
            text = line.strip()
            if re.match(r"IF\s*\(.*\)\s*THEN\s*$", text, re.I):
                depth += 1
                if re.search(r"embed_external_count\s*(>\s*0|/=\s*0|>=\s*1)", text, re.I):
                    guarded = True
            elif re.match(r"END\s*IF\b", text, re.I):
                depth -= 1
                if depth == 0:
                    guarded = False
            elif re.match(r"ELSE", text, re.I):
                guarded = False
        if not guarded:
            raise SystemExit(
                "ener_com%eext is added to etot outside "
                "IF (embed_external_count > 0): an EXTERNAL POTENTIAL ADD deck "
                "without session charges would count it twice"
            )
    print("embed_eval_energy_grad guards the eext sum")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
