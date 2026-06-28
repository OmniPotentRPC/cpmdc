#!/usr/bin/env python3
"""Check that typed CPMD render fixtures mention every structured field."""

from __future__ import annotations

import re
import sys
from pathlib import Path


STRUCTS = {
    "cpmd": "CPMDCpmdSection",
    "system": "CPMDSystemSection",
    "dft": "CPMDDftSection",
    "atoms": "CPMDAtomsSection",
}
EXEMPT_FIELDS = {
    "directives",
    "subsections",
}
FIELD = re.compile(r"^\s*(\w+)\s+@\d+\s*:", re.M)


def struct_body(schema: str, struct_name: str) -> str:
    match = re.search(rf"struct {struct_name} \{{(.*?)\n\}}", schema, re.S)
    if not match:
        raise RuntimeError(f"schema missing struct {struct_name}")
    return match.group(1)


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    schema = (repo / "schema/Potentials.capnp").read_text()
    haystack = (repo / "tests/test_params_cp_dft_render.c").read_text()
    haystack += "\n".join(
        fixture.read_text() for fixture in sorted((repo / "tests").glob("*.capnp.txt"))
    )

    missing: list[str] = []
    for section, struct_name in STRUCTS.items():
        fields = sorted(set(FIELD.findall(struct_body(schema, struct_name))))
        for field in fields:
            if field in EXEMPT_FIELDS:
                continue
            if not re.search(rf"\b{re.escape(field)}\b", haystack):
                missing.append(f"{section}.{field}")

    if missing:
        print("typed CPMD fields without render fixture coverage:")
        print("\n".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
