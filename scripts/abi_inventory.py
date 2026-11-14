#!/usr/bin/env python3
"""Write or check schema/inventory/surface_map_abi.txt.

The inventory joins the cpmdc_ symbols exported by libcpmdc with the
function declarations of the public headers. Every public declaration must
be exported, and every exported cpmdc_ symbol must be declared in a public
header or listed in scripts/abi_internal_exports.txt.

The library is read with ``nm -D --defined-only`` (``nm -gU`` on macOS);
set NM to use another nm, such as llvm-nm.
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADERS = ("include/cpmdc.h", "include/cpmdc_features.h")
INTERNAL = "scripts/abi_internal_exports.txt"
INVENTORY = "schema/inventory/surface_map_abi.txt"
PREFIX = "cpmdc_"

DECLARATION = re.compile(r"(?m)^(?!typedef\b)[A-Za-z_][\w\s*]*?\b(cpmdc_\w+)\s*\(")
DOC_BEFORE = re.compile(r"/\*\*((?:(?!\*/).)*)\*/\s*$", re.S)
SINCE = re.compile(r"[@\\]since\s+(?:ABI\s+)?(\d+)", re.I)
ABI_GENERATION = re.compile(r"(?m)^\s*#\s*define\s+CPMDC_ABI_VERSION\s+(\d+)\b")


@dataclass(frozen=True)
class Declaration:
    symbol: str
    header: str
    since: str | None


def blank_comments(text: str) -> str:
    """Replace comments and preprocessor lines with spaces, keeping offsets."""

    def spaces(match: re.Match[str]) -> str:
        return re.sub(r"[^\n]", " ", match.group(0))

    text = re.sub(r"/\*.*?\*/|//[^\n]*", spaces, text, flags=re.S)
    return re.sub(r"(?m)^[ \t]*#(?:[^\n]*\\\n)*[^\n]*", spaces, text)


def header_declarations(root: Path) -> tuple[list[Declaration], str | None]:
    declarations: list[Declaration] = []
    generation = None
    for relative in HEADERS:
        raw = (root / relative).read_text(encoding="utf-8")
        found = ABI_GENERATION.search(raw)
        if found:
            generation = found.group(1)
        code = blank_comments(raw)
        for match in DECLARATION.finditer(code):
            doc = DOC_BEFORE.search(raw, 0, match.start())
            since = SINCE.search(doc.group(1)) if doc else None
            declarations.append(
                Declaration(
                    match.group(1),
                    Path(relative).name,
                    since.group(1) if since else None,
                )
            )
    return declarations, generation


def nm_command(library: Path) -> list[str]:
    nm = shlex.split(os.environ.get("NM", "nm"))
    if sys.platform == "darwin":
        return [*nm, "-gU", str(library)]
    return [*nm, "-D", "--defined-only", str(library)]


def exported_symbols(library: Path) -> set[str]:
    command = nm_command(library)
    try:
        done = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        raise SystemExit(f"{command[0]} not found; set NM to an nm for {library}")
    if done.returncode != 0:
        raise SystemExit(f"{shlex.join(command)} failed:\n{done.stderr}")
    symbols = set()
    for line in done.stdout.splitlines():
        fields = line.split()
        if len(fields) < 3 or fields[-2] == "U" or not fields[-2].isupper():
            continue
        name = fields[-1].split("@", 1)[0]
        if sys.platform == "darwin" and name.startswith("_"):
            name = name[1:]
        if name.startswith(PREFIX):
            symbols.add(name)
    return symbols


def internal_exports(root: Path) -> list[str]:
    lines = (root / INTERNAL).read_text(encoding="utf-8").splitlines()
    return [line.strip() for line in lines if line.strip() and not line.startswith("#")]


def table(header: list[str], rows: list[list[str]]) -> list[str]:
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(len(header))]

    def line(cells: list[str]) -> str:
        return "| " + " | ".join(c.ljust(w) for c, w in zip(cells, widths)) + " |"

    rule = "|" + "+".join("-" * (w + 2) for w in widths) + "|"
    return [line(header), rule, *(line(row) for row in rows)]


def build(root: Path, exported: set[str]) -> tuple[str, list[str]]:
    declarations, generation = header_declarations(root)
    internal = internal_exports(root)
    public = {d.symbol for d in declarations}
    problems: list[str] = []

    seen: set[str] = set()
    for d in declarations:
        if d.symbol in seen:
            problems.append(f"declared more than once in the public headers: {d.symbol}")
        seen.add(d.symbol)
    seen = set()
    for symbol in internal:
        if symbol in seen:
            problems.append(f"listed more than once in {INTERNAL}: {symbol}")
        seen.add(symbol)

    for symbol in sorted(public - exported):
        problems.append(f"declared in a public header but not exported: {symbol}")
    for symbol in sorted(exported - public - set(internal)):
        problems.append(
            f"exported but declared in no public header: {symbol}"
            f" (declare it in include/ or list it in {INTERNAL})"
        )
    for symbol in sorted(set(internal) & public):
        problems.append(f"listed in {INTERNAL} but declared in a public header: {symbol}")
    for symbol in sorted(set(internal) - exported):
        problems.append(f"listed in {INTERNAL} but not exported: {symbol}")

    mark = {True: "Y", False: "N"}
    public_rows = [
        [f"={d.symbol}=", d.header, d.since or "-", mark[d.symbol in exported]]
        for d in sorted({d.symbol: d for d in declarations}.values(), key=lambda d: d.symbol)
    ]
    internal_rows = [[f"={s}=", mark[s in exported]] for s in sorted(set(internal))]

    out = [
        "#+TITLE: Exported C ABI inventory",
        "# Generated by scripts/abi_inventory.py from the exported symbols of",
        "# libcpmdc and the public headers. Do not edit; regenerate with",
        "# meson compile -C build abi-inventory",
        "",
        f"ABI generation (CPMDC_ABI_VERSION): {generation if generation is not None else '-'}",
        "",
        "* Public symbols",
        "Declared in " + " or ".join(HEADERS) + ".",
        "The abi column is the ABI generation a declaration records with",
        "@since, or - when it records none.",
        "",
        *table(["symbol", "header", "abi", "exported"], public_rows),
        "",
        "* Internal exports",
        f"Exported without a public declaration and listed in {INTERNAL}.",
        "",
        *table(["symbol", "exported"], internal_rows),
        "",
        "* Summary",
        f"Public declarations: {len(public)}",
        f"Exported {PREFIX} symbols: {len(exported)}",
        f"Internal exports: {len(set(internal))}",
        "Parity: " + ("OK" if not problems else "FAIL"),
        *(f"- {p}" for p in problems),
    ]
    return "\n".join(out) + "\n", problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="rewrite the inventory")
    mode.add_argument("--check", action="store_true", help="fail when it differs")
    parser.add_argument("library", type=Path, help="built libcpmdc shared library")
    parser.add_argument("--root", type=Path, default=ROOT, help="source tree")
    args = parser.parse_args()

    text, problems = build(args.root, exported_symbols(args.library))
    target = args.root / INVENTORY
    regenerate = (
        f"meson compile -C {args.library.resolve().parent} abi-inventory"
        f" (or python3 scripts/abi_inventory.py --write {args.library})"
    )
    status = 0
    if args.write:
        target.write_text(text, encoding="utf-8")
        print(f"wrote {INVENTORY}")
    else:
        current = target.read_text(encoding="utf-8") if target.is_file() else ""
        if current != text:
            sys.stdout.writelines(
                difflib.unified_diff(
                    current.splitlines(keepends=True),
                    text.splitlines(keepends=True),
                    f"a/{INVENTORY}",
                    f"b/{INVENTORY}",
                )
            )
            print(f"{INVENTORY} is stale; run {regenerate}", file=sys.stderr)
            status = 1
    for problem in problems:
        print(f"C ABI parity: {problem}", file=sys.stderr)
    if problems:
        print(
            "fix the header, the library, or the internal list, then run " + regenerate,
            file=sys.stderr,
        )
        status = 1
    if status == 0 and args.check:
        print(f"{INVENTORY} matches {args.library.name} and the public headers")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
