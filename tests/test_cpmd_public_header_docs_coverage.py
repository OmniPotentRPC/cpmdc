#!/usr/bin/env python3
"""Check that public API docs describe the serialized CPMD input surface."""

from __future__ import annotations

import sys
from pathlib import Path


REQUIRED_TERMS = [
    "do not need C setter functions",
    "inputSections",
    "inputBlocks",
    "params.inputSections.cpmd.maxIter",
    "params.inputSections.system.cell",
    "params.inputSections.dft.hfxScreening",
    "params.inputSections.atoms.pseudopotentials",
    "catalog.section.VDW",
    "params.inputSections.raw",
    "cpmdc_session_create",
    "cpmdc_calculate_result",
]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_cpmd_public_header_docs_coverage.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    targets = [
        root / "include/cpmdc.h",
        root / "docs/source/api/global.rst",
    ]
    failures: list[str] = []
    for path in targets:
        text = path.read_text(encoding="utf-8")
        for term in REQUIRED_TERMS:
            if term not in text:
                failures.append(f"{path}: missing {term!r}")

    if failures:
        print("public API docs do not describe the CPMD input surface:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
