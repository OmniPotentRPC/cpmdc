#!/usr/bin/env python3
"""Check that every directive section arm has long-tail render coverage."""

from __future__ import annotations

import re
import sys
from pathlib import Path


SECTION_UNION = re.compile(r"struct CPMDInputSection \{(.*?)\n\}", re.S)
UNION_ARM = re.compile(r"^\s*(\w+)\s+@\d+\s*:\s*([^;]+);", re.M)
CATALOG_SECTION_ROW = re.compile(r"\|\s*~catalog\.section\.([A-Z0-9_]+)~")

EXEMPT_ARMS = {
    "generic",
    "cpmd",
    "system",
    "dft",
    "atoms",
    "set",
    "raw",
}


def directive_section_arms(schema: str) -> dict[str, str]:
    match = SECTION_UNION.search(schema)
    if not match:
        raise RuntimeError("schema missing CPMDInputSection")
    arms: dict[str, str] = {}
    for arm, arm_type in UNION_ARM.findall(match.group(1)):
        if arm in EXEMPT_ARMS:
            continue
        if arm_type.strip() == "CPMDDirectiveSection":
            arms[arm] = arm.upper()
    return arms


def catalog_sections(reference: str) -> set[str]:
    return set(CATALOG_SECTION_ROW.findall(reference))


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    schema = (repo / "schema/Potentials.capnp").read_text(encoding="utf-8")
    fixture = (repo / "tests/params_long_tail_sections.capnp.txt").read_text(
        encoding="utf-8"
    )
    render_test = (repo / "tests/test_params_cp_dft_render.c").read_text(
        encoding="utf-8"
    )
    reference = (repo / "docs/orgmode/reference/cpmd-options.org").read_text(
        encoding="utf-8"
    )

    missing: list[str] = []
    for arm, section in sorted(directive_section_arms(schema).items()):
        if not re.search(rf"\(\s*{re.escape(arm)}\s*=", fixture):
            missing.append(f"fixture missing inputSections.{arm}")
        if f'"&{section}"' not in render_test:
            missing.append(f"render assertion missing &{section}")
        if section not in catalog_sections(reference):
            missing.append(f"reference docs missing catalog.section.{section}")

    if missing:
        print("long-tail directive section coverage is incomplete:")
        print("\n".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
