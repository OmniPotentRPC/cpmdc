#!/usr/bin/env python3
"""Every public C function and Python name is documented."""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path


def c_functions(text: str) -> set[str]:
    stripped = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    stripped = re.sub(r"//.*", " ", stripped)
    return set(re.findall(r"\b(cpmdc_[A-Za-z0-9_]+)\s*\(", stripped))


def documented_functions(text: str) -> set[str]:
    found = set()
    pattern = re.compile(
        r"/\*\*(?:[^*]|\*(?!/))*\*/\s*"
        r"(?:[\w][\w\s\*]*?)\b(cpmdc_[A-Za-z0-9_]+)\s*\(",
        re.S,
    )
    for match in pattern.finditer(text):
        found.add(match.group(1))
    return found


def python_public(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    missing: list[str] = []
    names: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
            names.append(node.name)
            if not ast.get_docstring(node):
                missing.append(f"{path.name}:{node.name} has no docstring")
        if isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            names.append(node.name)
            if not ast.get_docstring(node):
                missing.append(f"{path.name}:{node.name} has no docstring")
            for item in node.body:
                if not isinstance(item, ast.FunctionDef):
                    continue
                if item.name.startswith("_"):
                    continue
                names.append(f"{node.name}.{item.name}")
                if not ast.get_docstring(item):
                    missing.append(f"{path.name}:{node.name}.{item.name} has no docstring")
    return missing + [f"name:{item}" for item in names]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_public_api_docs.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    failures: list[str] = []
    for header in sorted((root / "include").glob("*.h")):
        text = header.read_text(encoding="utf-8")
        missing = sorted(c_functions(text) - documented_functions(text))
        for name in missing:
            failures.append(f"{header.relative_to(root)}: {name} has no doc comment")
    reference = (root / "docs/orgmode/reference/python.org").read_text(encoding="utf-8")
    howto = (root / "docs/orgmode/howto/python.org").read_text(encoding="utf-8")
    docs = reference + "\n" + howto
    init = (root / "python/cpmdc/__init__.py").read_text(encoding="utf-8")
    exported = re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"', init.split("__all__", 1)[1].split("]", 1)[0])
    for name in exported:
        if name not in docs:
            failures.append(f"{name} is exported but not named in the Python reference")
    if "CPMDCCalculator" not in docs:
        failures.append("CPMDCCalculator is not named in the Python reference")
    for path in sorted((root / "python/cpmdc").glob("*.py")):
        for item in python_public(path):
            if item.startswith("name:"):
                name = item.removeprefix("name:")
                if name not in docs and name.split(".")[-1] not in docs:
                    failures.append(f"{name} is not named in the Python reference")
            else:
                failures.append(item)
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
