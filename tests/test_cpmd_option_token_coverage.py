#!/usr/bin/env python3
"""Check that fixture inline option tokens have catalog render coverage."""

from __future__ import annotations

import re
import sys
from pathlib import Path


COVERAGE_ENTRY = re.compile(
    r'\{\s*"([^"]+)",\s*"([^"]+)",\s*"((?:\\.|[^"])*)"\s*\}',
    re.S,
)
OPTION_ASSIGNMENT = re.compile(r'(\w*Options)\s*=\s*"([^"]*)"')
NUMBER_TOKEN = re.compile(r"[0-9.+-]+")


def decode_c_string(value: str) -> str:
    return bytes(value, "utf-8").decode("unicode_escape")


def coverage_by_field(repo: Path) -> dict[str, list[tuple[str, str]]]:
    source = (repo / "tests/test_params_cp_dft_render.c").read_text()
    by_field: dict[str, list[tuple[str, str]]] = {}
    for match in COVERAGE_ENTRY.finditer(source):
        feature_id, field, token = match.groups()
        by_field.setdefault(field, []).append((feature_id, decode_c_string(token)))
    return by_field


def option_tokens(payload: str) -> list[str]:
    tokens = []
    for raw in payload.replace("\\n", " ").split():
        token = raw.split("=", 1)[0]
        if token and not NUMBER_TOKEN.fullmatch(token):
            tokens.append(token)
    return tokens


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    coverage = coverage_by_field(repo)
    missing: list[str] = []
    for fixture in sorted((repo / "tests").glob("*.capnp.txt")):
        text = fixture.read_text()
        for match in OPTION_ASSIGNMENT.finditer(text):
            field, payload = match.groups()
            field_coverage = coverage.get(field, [])
            for token in option_tokens(payload):
                if not any(token in render_token for _, render_token in field_coverage):
                    missing.append(
                        f"{fixture.relative_to(repo)}:{field}: missing {token!r}"
                    )

    if missing:
        print("CPMD option tokens without catalog render coverage:")
        print("\n".join(missing))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
