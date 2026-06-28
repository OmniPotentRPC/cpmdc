#!/usr/bin/env python3
"""Check that CPMD escape-hatch carriers stay covered and documented."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REQUIRED_SCHEMA_TERMS = [
    "inputBlocks",
    "generic",
    "set",
    "raw",
    "CPMDGenericSection",
    "CPMDSetDirective",
]

FIXTURE_TERMS = {
    "tests/params_escape_generic.capnp.txt": [
        "inputBlocks",
        "generic",
        'name = "PROPERTIES"',
        "directives",
        "raw",
    ],
    "tests/params_set_directives.capnp.txt": [
        "set",
        'key = "CPMD.PRINT FORCES ON"',
        'key = "SYSTEM.POISSON SOLVER HOCKNEY"',
        'key = "DFT.GC-CUTOFF"',
    ],
}

RENDER_TEST_TERMS = {
    "tests/test_params_escape_render.c": [
        "&INFO",
        "DEMO KEYWORD",
        "&PROPERTIES",
        "DIPOLE",
        "&EXTE",
    ],
    "tests/cmocka/test_set_directives_render.c": [
        "PRINT FORCES ON",
        "MAXSTEP",
        "POISSON SOLVER HOCKNEY",
        "GC-CUTOFF",
        "assert_known_sections_are_unique",
    ],
}

DOC_TERMS = [
    "inputBlocks",
    "inputSections.generic",
    "inputSections.set",
    "inputSections.raw",
]


def missing_terms(path: Path, terms: list[str]) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [f"{path}: missing {term}" for term in terms if term not in text]


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    failures: list[str] = []
    failures.extend(
        missing_terms(repo / "schema/Potentials.capnp", REQUIRED_SCHEMA_TERMS)
    )
    for relpath, terms in FIXTURE_TERMS.items():
        failures.extend(missing_terms(repo / relpath, terms))
    for relpath, terms in RENDER_TEST_TERMS.items():
        failures.extend(missing_terms(repo / relpath, terms))

    reference = repo / "docs/orgmode/reference/cpmd-options.org"
    failures.extend(missing_terms(reference, DOC_TERMS))
    reference_text = reference.read_text(encoding="utf-8")
    for feature_id in (
        "params.inputBlocks",
        "params.inputSections.generic.directives",
        "params.inputSections.set.key",
        "params.inputSections.raw",
    ):
        if not re.search(rf"~{re.escape(feature_id)}~", reference_text):
            failures.append(f"{reference}: missing {feature_id}")

    if failures:
        print("escape-hatch coverage is incomplete:")
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
