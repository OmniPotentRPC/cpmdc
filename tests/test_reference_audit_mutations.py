#!/usr/bin/env python3
"""The reference audit rejects a contracts page or a Python surface that drifts."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

PAGE = "docs/orgmode/reference/contracts.org"
INIT = "python/cpmdc/__init__.py"


def drop_section(heading: str) -> Callable[[str], str]:
    def mutate(text: str) -> str:
        pattern = re.compile(
            rf"^\*+ ~{re.escape(heading)}~\n.*?(?=^\*+ )", re.MULTILINE | re.DOTALL
        )
        result, count = pattern.subn("", text, count=1)
        if count != 1:
            raise LookupError(f"no section {heading}")
        return result

    return mutate


def replace(old: str, new: str) -> Callable[[str], str]:
    def mutate(text: str) -> str:
        if old not in text:
            raise LookupError(f"absent: {old}")
        return text.replace(old, new, 1)

    return mutate


def empty_row(heading: str, item: str) -> Callable[[str], str]:
    def mutate(text: str) -> str:
        pattern = re.compile(
            rf"(^\*+ ~{re.escape(heading)}~\n.*?^\| {re.escape(item)} \|)[^\n]*",
            re.MULTILINE | re.DOTALL,
        )
        result, count = pattern.subn(r"\1 |", text, count=1)
        if count != 1:
            raise LookupError(f"no row {item} in {heading}")
        return result

    return mutate


def append(extra: str) -> Callable[[str], str]:
    return lambda text: text + extra


STRAY = """
*** ~cpmdc_no_such_call~

| Item | Description |
|--------+--------|
| Units | none |
| Layout | none |
| Ownership | none |
| Threads and MPI | none |
| Returns and errors | none |
| Python | none |
| Defaults | none |
"""

MUTATIONS: list[tuple[str, str, Callable[[str], str], str]] = [
    ("C function entry removed", PAGE, drop_section("cpmdc_energy_forces"), "cpmdc_energy_forces"),
    (
        "C function entry renamed",
        PAGE,
        replace("*** ~cpmdc_energy_forces~", "*** ~cpmdc_energy_forcez~"),
        "cpmdc_energy_forcez",
    ),
    (
        "C function row emptied",
        PAGE,
        empty_row("cpmdc_energy_forces", "Units"),
        "cpmdc_energy_forces",
    ),
    ("entry for a function that does not exist", PAGE, append(STRAY), "cpmdc_no_such_call"),
    ("C type entry removed", PAGE, drop_section("CPMDCSession"), "CPMDCSession"),
    ("Python entry removed", PAGE, drop_section("cpmdc.Session.close"), "cpmdc.Session.close"),
    (
        "Python signature changed",
        PAGE,
        replace("| Signature | ~close()~ |", "| Signature | ~close(force)~ |"),
        "Signature row is 'close(force)'",
    ),
    (
        "Python entry row emptied",
        PAGE,
        empty_row("cpmdc.Session.close", "Raises"),
        "cpmdc.Session.close",
    ),
    (
        "Python name added without an entry",
        INIT,
        lambda text: replace("__all__ = [\n", '__all__ = [\n    "undocumented_name",\n')(text)
        + "\n\ndef undocumented_name() -> None:\n    return None\n",
        "undocumented_name",
    ),
    (
        "deprecated alias row renamed",
        PAGE,
        replace("| cpmdc.CpmdcError | cpmdc.CPMDCError |", "| cpmdc.OldError | cpmdc.CPMDCError |"),
        "cpmdc.OldError",
    ),
]


def copied_tree(source: Path, destination: Path) -> None:
    for relative in ("include", "docs", "python", "schema"):
        path = source / relative
        if path.is_dir():
            shutil.copytree(path, destination / relative)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_reference_audit_mutations.py SOURCE_ROOT BUILD_DIR", file=sys.stderr)
        return 2
    source = Path(sys.argv[1]).resolve()
    build = Path(sys.argv[2]).resolve()
    audit = source / "tests" / "test_reference_audit.py"
    failures: list[str] = []

    clean = subprocess.run(
        [sys.executable, audit, source, build], check=False, capture_output=True, text=True
    )
    if clean.returncode != 0:
        print(f"the unmutated tree fails the audit:\n{clean.stderr}", file=sys.stderr)
        return 1

    for description, relative, mutate, expected in MUTATIONS:
        with tempfile.TemporaryDirectory() as tmp:
            tree = Path(tmp) / "tree"
            tree.mkdir()
            copied_tree(source, tree)
            path = tree / relative
            try:
                path.write_text(mutate(path.read_text(encoding="utf-8")), encoding="utf-8")
            except LookupError as error:
                failures.append(f"{description}: mutation target is absent ({error})")
                continue
            checked = subprocess.run(
                [sys.executable, audit, tree, build],
                check=False,
                capture_output=True,
                text=True,
            )
            if checked.returncode == 0:
                failures.append(f"{description}: the audit accepted the mutation")
            elif expected not in checked.stderr:
                failures.append(
                    f"{description}: the audit failed without naming {expected}: "
                    f"{checked.stderr.strip()[:400]}"
                )
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"{len(MUTATIONS)} mutations rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
