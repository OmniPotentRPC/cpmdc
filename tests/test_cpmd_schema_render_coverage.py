#!/usr/bin/env python3
"""Check that CPMDCpmdSection fields have render coverage mappings."""

from __future__ import annotations

import re
import sys
from pathlib import Path


CPMD_SECTION = re.compile(r"struct CPMDCpmdSection \{(.*?)\n\}", re.S)
FIELD = re.compile(r"^\s*(\w+)\s+@\d+\s*:", re.M)
COVERAGE_FIELD = re.compile(r'"catalog\.cpmd\.[^"]+"\s*,\s*"(\w+)"')

EXEMPT_FIELDS = {
    "directives",
}


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    schema = (repo / "schema/Potentials.capnp").read_text()
    tests = (repo / "tests/test_params_cp_dft_render.c").read_text()

    match = CPMD_SECTION.search(schema)
    if not match:
        print("schema missing CPMDCpmdSection")
        return 1

    fields = set(FIELD.findall(match.group(1))) - EXEMPT_FIELDS
    covered = set(COVERAGE_FIELD.findall(tests))
    missing = sorted(fields - covered)

    if missing:
        print("CPMDCpmdSection fields without render coverage mapping:")
        print("\n".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
