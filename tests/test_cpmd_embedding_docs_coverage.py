#!/usr/bin/env python3
"""Check that the embedding how-to names the public CPMD input surface."""

from __future__ import annotations

import sys
from pathlib import Path


REQUIRED_GROUPS = {
    "public ABI calls": [
        "cpmdc_set_params",
        "cpmdc_session_create",
        "cpmdc_session_set_params",
        "cpmdc_session_calculate_result",
        "cpmdc_calculate_result",
    ],
    "wire messages": [
        "CPMDParams",
        "ForceInput",
        "PotentialResult",
    ],
    "ABI feature IDs": [
        "abi.cpmdc_set_params",
        "abi.cpmdc_session_create",
        "abi.cpmdc_session_set_params",
        "abi.cpmdc_session_calculate_result",
        "abi.cpmdc_calculate_result",
    ],
    "typed params fields": [
        "params.inputSections.cpmd.maxIter",
        "params.inputSections.system.cell",
        "params.inputSections.dft.hfxScreening",
        "params.inputSections.atoms.pseudopotentials",
    ],
    "long-tail section feature IDs": [
        "catalog.section.PIMD",
        "catalog.section.VDW",
        "params.inputSections.pimd.directives",
        "params.inputSections.vdw.subsections",
    ],
    "escape-hatch fields": [
        "params.inputBlocks",
        "params.inputSections.generic.directives",
        "params.inputSections.set.key",
        "params.inputSections.raw",
    ],
}


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_cpmd_embedding_docs_coverage.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    embedding = root / "docs/orgmode/howto/embedding.org"
    text = embedding.read_text(encoding="utf-8")

    failures: list[str] = []
    for group, terms in REQUIRED_GROUPS.items():
        for term in terms:
            if term not in text:
                failures.append(f"{embedding}: missing {group} term {term!r}")

    if failures:
        print("embedding docs coverage is incomplete:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
