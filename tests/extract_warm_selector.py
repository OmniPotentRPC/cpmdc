#!/usr/bin/env python3
"""Compile the exact orbital-history routines from the native source patch."""
from pathlib import Path
import sys

patch = Path(sys.argv[1]).read_text().splitlines()
start = next(i for i, line in enumerate(patch) if line == "   PUBLIC :: rwfopt")
end = next(i for i in range(start, len(patch)) if patch[i].startswith("@@ "))
body = []
for line in patch[start:end]:
    if line.startswith("+"):
        body.append(line[1:])
    elif line.strip() == "CONTAINS":
        body.append("CONTAINS")
Path(sys.argv[2]).write_text(
    "MODULE warm_selector\n"
    "IMPLICIT NONE\n"
    "INTEGER, PARAMETER :: real_8 = KIND(1.0D0)\n"
    "LOGICAL :: embed_have_orbitals = .FALSE., embed_write_files = .TRUE.\n"
    + "\n".join(body)
    + "\nSUBROUTINE stopgm(procedure, reason, line, file)\n"
    "CHARACTER(*), INTENT(IN) :: procedure, reason, file\n"
    "INTEGER, INTENT(IN) :: line\n"
    "PRINT *, procedure, reason, file, line\n"
    "ERROR STOP 99\nEND SUBROUTINE\nEND MODULE\n"
)
