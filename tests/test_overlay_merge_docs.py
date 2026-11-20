#!/usr/bin/env python3
"""Fail when the overlay merge docs disagree with the merge code.

The `cpmdc_configure()` comment in include/cpmdc.h, the C ABI reference, and
the overlay release notes must state the resolved precedence, the eV to Ry
cutoff conversion, the explicit singlet, the typed section directives that
block an overlay scalar, and the message an unreadable arm list stores, in
the format src/cpmdc.c writes.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

STALE = (
    "only while the arm still has the schema default",
    "only where the arm still has the schema default",
)


def configure_comment(header: str) -> str:
    end = header.index("CPMDC_API int cpmdc_configure(")
    start = header.rindex("/**", 0, end)
    lines = header[start:end].splitlines()
    return " ".join(re.sub(r"^\s*/?\*+/?\s?", "", line) for line in lines)


def invalid_list_pattern(source: str) -> re.Pattern[str]:
    match = re.search(
        r"static int cpmd_reject_invalid_list\([^)]*\)\s*\{(.*?)\n\}", source, re.S
    )
    if not match:
        raise LookupError("src/cpmdc.c has no cpmd_reject_invalid_list()")
    call = re.search(r"snprintf\((.*?)\);", match.group(1), re.S)
    if not call:
        raise LookupError("cpmd_reject_invalid_list() writes no message")
    fmt = "".join(re.findall(r'"((?:[^"\\]|\\.)*)"', call.group(1)))
    fmt = fmt.removeprefix("common overlay: ")
    parts = [re.escape(part) for part in fmt.split("%s")]
    return re.compile(r"\S+(?: directives)?".join(parts))


def flat(text: str) -> str:
    return " ".join(text.split())


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_overlay_merge_docs.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    failures = []
    pattern: re.Pattern[str] | None = None
    try:
        pattern = invalid_list_pattern((root / "src" / "cpmdc.c").read_text())
    except LookupError as exc:
        failures.append(str(exc))
    fragments = root / "docs" / "newsfragments"
    docs = {
        "include/cpmdc.h cpmdc_configure()": flat(
            configure_comment((root / "include" / "cpmdc.h").read_text())
        ),
        "docs/orgmode/reference/c-abi.org": flat(
            (root / "docs" / "orgmode" / "reference" / "c-abi.org")
            .read_text()
            .replace("~", "`")
        ),
        "overlay release notes": flat(
            (fragments / "+overlay-merge-fields.fixed.md").read_text()
            + " "
            + (fragments / "+overlay-typed-presence.fixed.md").read_text()
        ),
    }
    required = (
        "beats the overlay",
        "beats a schema default",
        "explicit singlet",
        "eV",
        "Ry",
        "`FUNCTIONAL`",
        "`MULTIPLICITY`",
    )
    for name, text in docs.items():
        for phrase in required:
            if phrase not in text:
                failures.append(f"{name}: does not say {phrase!r}")
        for phrase in STALE:
            if phrase in text:
                failures.append(f"{name}: still says {phrase!r}")
        if pattern and not pattern.search(text):
            failures.append(
                f"{name}: no unreadable list message matching {pattern.pattern!r}"
            )
    for failure in failures:
        print(failure, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
