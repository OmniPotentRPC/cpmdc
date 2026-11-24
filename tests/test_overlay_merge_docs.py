#!/usr/bin/env python3
"""Fail when the overlay merge docs disagree with the merge code.

The `cpmdc_configure()` comment in include/cpmdc.h, the C ABI reference, and
the overlay release notes must state the resolved precedence, the eV to Ry
cutoff conversion, the explicit singlet, the typed section directives that
block an overlay scalar, that the carrier check runs for every set overlay
field, and the message an unreadable arm list stores, in the format
src/cpmdc.c writes. No overlay doc may claim that the scalars are resolved
once or come from one result, and the merge in src/cpmdc.c must not skip the
carrier check when a typed field sets the value. The pages that describe the
rendered `&SYSTEM` scalars must say a `CUTOFF`, `CHARGE`, or `MULTIPLICITY`
directive is rejected for a default value and for any other value.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

STALE = (
    "only while the arm still has the schema default",
    "only where the arm still has the schema default",
    "resolved once",
    "one resolved result",
    "one resolved configuration",
    "this one resolved result",
    "would fill a schema default",
    "for a field no typed field sets",
    "a scalar the overlay would copy",
    "nothing else, so they cannot disagree",
)

CARRIER_CHECK = "the carrier check runs for every set overlay field"
SYSTEM_DIRECTIVE = (
    "a `CUTOFF`, `CHARGE`, or `MULTIPLICITY` directive in the typed "
    "`system` section"
)
SYSTEM_VALUES = "rejected for a default value and for any other value"


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


def typed_guarded_carrier_checks(source: str) -> list[str]:
    """Carrier checks in cpmd_merge_common_into_arm() under a typed test.

    Each cpmd_reject_opaque_scalar() call must sit under a condition on the
    overlay field alone; a `typed.` test in any enclosing `if` skips the
    check when a typed value and an opaque carrier are both present.
    """
    match = re.search(
        r"static int cpmd_merge_common_into_arm\([^)]*\)\s*\{(.*?)\n\}",
        source,
        re.S,
    )
    if not match:
        raise LookupError("src/cpmdc.c has no cpmd_merge_common_into_arm()")
    body = match.group(1)
    calls = list(re.finditer(r"cpmd_reject_opaque_scalar\(", body))
    if len(calls) != 4:
        raise LookupError(
            f"cpmd_merge_common_into_arm() has {len(calls)} carrier checks, want 4"
        )
    guarded = []
    for call in calls:
        open_conditions: list[str] = []
        prefix = body[: call.start()]
        for token in re.finditer(r"\bif\s*\(([^;{}]*)\)\s*\{|\{|\}", prefix, re.S):
            text = token.group(0)
            if text.startswith("if"):
                open_conditions.append(token.group(1))
            elif text == "{":
                open_conditions.append("")
            else:
                if open_conditions:
                    open_conditions.pop()
        line = body[: call.start()].count("\n") + 1
        for cond in open_conditions:
            if "typed." in cond:
                guarded.append(f"line {line} of the body, under if ({flat(cond)})")
    return guarded


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
    try:
        for where in typed_guarded_carrier_checks(
            (root / "src" / "cpmdc.c").read_text()
        ):
            failures.append(
                f"src/cpmdc.c: carrier check skipped for a typed value at {where}"
            )
    except LookupError as exc:
        failures.append(str(exc))
    fragments = root / "docs" / "newsfragments"
    reference = root / "docs" / "orgmode" / "reference"
    overlay_fragments = (
        fragments / "+overlay-merge-fields.fixed.md",
        fragments / "+overlay-typed-presence.fixed.md",
    )
    if all(path.is_file() for path in overlay_fragments):
        overlay_notes = " ".join(path.read_text() for path in overlay_fragments)
    else:
        changelog = (root / "CHANGELOG.md").read_text()
        # The two overlay fragments were consumed into the 0.6.0 release
        # notes, so that section is where their wording lives now. The
        # current version in version.txt would look in the wrong section.
        section = re.search(
            r"## \[0\.6\.0\].*?(?=\n## \[|\Z)",
            changelog,
            re.S,
        )
        overlay_notes = section.group(0) if section else changelog
    docs = {
        "include/cpmdc.h cpmdc_configure()": flat(
            configure_comment((root / "include" / "cpmdc.h").read_text())
        ),
        "docs/orgmode/reference/c-abi.org": flat(
            (root / "docs" / "orgmode" / "reference" / "c-abi.org")
            .read_text()
            .replace("~", "`")
        ),
        "overlay release notes": flat(overlay_notes),
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
    wording = dict(docs)
    for page in ("cpmd-options.org", "cpmdparams-schema.org", "contracts.org"):
        wording[f"docs/orgmode/reference/{page}"] = flat(
            (reference / page).read_text().replace("~", "`")
        )
    for name, text in docs.items():
        for phrase in required:
            if phrase not in text:
                failures.append(f"{name}: does not say {phrase!r}")
        if pattern and not pattern.search(text):
            failures.append(
                f"{name}: no unreadable list message matching {pattern.pattern!r}"
            )
    for name, text in wording.items():
        if name != "docs/orgmode/reference/contracts.org" and (
            CARRIER_CHECK not in text.lower()
        ):
            failures.append(f"{name}: does not say {CARRIER_CHECK!r}")
        for phrase in STALE:
            if phrase in text.lower():
                failures.append(f"{name}: still says {phrase!r}")
    system_pages = {
        "include/cpmdc.h cpmdc_configure()": docs[
            "include/cpmdc.h cpmdc_configure()"
        ],
        "docs/orgmode/reference/c-abi.org": docs[
            "docs/orgmode/reference/c-abi.org"
        ],
        "docs/orgmode/reference/cpmd-options.org": wording[
            "docs/orgmode/reference/cpmd-options.org"
        ],
        "docs/orgmode/reference/cpmdparams-schema.org": wording[
            "docs/orgmode/reference/cpmdparams-schema.org"
        ],
        "docs/orgmode/reference/contracts.org": wording[
            "docs/orgmode/reference/contracts.org"
        ],
    }
    fragment = fragments / "+system-scalar-directives.fixed.md"
    if fragment.is_file():
        system_notes = fragment.read_text()
        system_notes_name = "docs/newsfragments/+system-scalar-directives.fixed.md"
    else:
        changelog = (root / "CHANGELOG.md").read_text()
        section = re.search(
            r"## \[0\.7\.0\].*?(?=\n## \[|\Z)",
            changelog,
            re.S,
        )
        system_notes = section.group(0) if section else changelog
        system_notes_name = "CHANGELOG.md [0.7.0]"
    system_pages[system_notes_name] = flat(system_notes)
    for name, text in system_pages.items():
        for phrase in (SYSTEM_DIRECTIVE, SYSTEM_VALUES):
            if phrase not in text:
                failures.append(f"{name}: does not say {phrase!r}")
    for failure in failures:
        print(failure, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
