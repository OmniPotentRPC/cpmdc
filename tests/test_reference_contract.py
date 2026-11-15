#!/usr/bin/env python3
"""The reference documents every public name, and the docs name no other.

Every function, type, enum constant, and macro of the installed headers
cpmdc.h, cpmdc_features.h, and cpmdc_restart.h must appear in
docs/orgmode/reference/c-abi.org; a function must appear with its
parameter list in a table row, which is where that page states a contract.
Every module of python/cpmdc must list its public functions and classes in
__all__. Every name in __all__, and every public member that an exported
class defines, must appear in docs/orgmode/reference/python.org; a member
appears as ~Class.member~ or ~module.Class.member~.

The README and the Org pages may not present any other cpmdc_ identifier
as a C symbol: a name that is hidden in libcpmdc or gone from the headers.
File names, library names, and the log tokens that CPMDC_SCF_STEPS writes
are not symbols. The changelog keeps the history and is not checked.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

HEADERS = ("include/cpmdc.h", "include/cpmdc_features.h", "include/cpmdc_restart.h")
# Headers whose names the docs may mention without the reference contract:
# the export macro header, and the optional plane-wave scaffold.
OTHER_HEADERS = ("include/cpmdc_export.h", "include/cpmdc_pwscale.h")
# Written to stderr under CPMDC_SCF_STEPS=1; see reference/environment.org.
LOG_TOKENS = {"cpmdc_scf_steps", "cpmdc_warm_start", "cpmdc_wfn_extrapolated"}


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def header_names(text: str) -> dict[str, set[str]]:
    code = strip_comments(text)
    return {
        "function": set(re.findall(r"\b(cpmdc_\w+)\s*\(", code)),
        "type": set(re.findall(r"\btypedef\b[^;]*?\b(\w+)\s*;", code))
        | set(re.findall(r"}\s*(\w+)\s*;", code)),
        "enum constant": set(re.findall(r"\b(CPMDC_[A-Z0-9_]+)\s*=", code)),
        "macro": {
            name
            for name in re.findall(r"(?m)^\s*#\s*define\s+(CPMDC\w*)", code)
            if not name.endswith("_H")
        },
    }


def table_rows(text: str) -> list[str]:
    return [line for line in text.splitlines() if line.lstrip().startswith("|")]


def check_reference(root: Path, failures: list[str]) -> set[str]:
    reference = (root / "docs/orgmode/reference/c-abi.org").read_text(encoding="utf-8")
    rows = table_rows(reference)
    public: set[str] = set()
    for rel in HEADERS:
        names = header_names((root / rel).read_text(encoding="utf-8"))
        for kind, found in names.items():
            public |= found
            for name in sorted(found):
                if kind == "function":
                    if not any(re.search(rf"\b{name}\(\w", row) for row in rows):
                        failures.append(
                            f"reference/c-abi.org: {name} ({rel}) has no table row with its "
                            "parameter list"
                        )
                elif not re.search(rf"~[^~\n]*\b{name}\b[^~\n]*~", reference):
                    failures.append(f"reference/c-abi.org: {kind} {name} ({rel}) is not named")
    for rel in OTHER_HEADERS:
        for found in header_names((root / rel).read_text(encoding="utf-8")).values():
            public |= found
    return public


def mentioned_symbols(text: str) -> set[str]:
    found = set()
    for match in re.finditer(r"(?<![\w/.-])(cpmdc_[a-z0-9_]*)(\*)?(?![\w.])", text):
        name, star = match.group(1), match.group(2)
        before = text[max(0, match.start() - 3) : match.start()]
        if before.endswith("lib") or (name.endswith("_") and not star):
            continue
        found.add(name + ("*" if star else ""))
    return found


def check_mentions(root: Path, public: set[str], failures: list[str]) -> None:
    pages = [root / "README.md"] + sorted((root / "docs/orgmode").rglob("*.org"))
    for page in pages:
        if page.name == "changelog.org":
            continue
        rel = page.relative_to(root).as_posix()
        text = page.read_text(encoding="utf-8")
        for name in sorted(mentioned_symbols(text)):
            if name.endswith("*"):
                prefix = name[:-1]
                if any(symbol.startswith(prefix) for symbol in public):
                    continue
            elif name in public or name in LOG_TOKENS:
                continue
            failures.append(f"{rel}: names {name}, which no public header declares")


def class_members(node: ast.ClassDef) -> set[str]:
    names = set()
    for item in node.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(item.name)
        elif isinstance(item, ast.Assign):
            names |= {target.id for target in item.targets if isinstance(target, ast.Name)}
        elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
            names.add(item.target.id)
    return {name for name in names if not name.startswith("_")}


def python_public(path: Path) -> tuple[set[str] | None, dict[str, set[str]]]:
    """__all__ of a module, None without one, and the public members of each class."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    exported = None
    defined: dict[str, set[str]] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets
        ):
            exported = {elt.value for elt in node.value.elts if isinstance(elt, ast.Constant)}
        elif isinstance(node, ast.ClassDef):
            defined[node.name] = class_members(node)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            defined[node.name] = set()
    return exported, defined


def check_python(root: Path, failures: list[str]) -> None:
    reference = (root / "docs/orgmode/reference/python.org").read_text(encoding="utf-8")
    for path in sorted((root / "python/cpmdc").glob("*.py")):
        rel = path.relative_to(root).as_posix()
        module = "cpmdc" if path.name == "__init__.py" else f"cpmdc.{path.stem}"
        exported, defined = python_public(path)
        if exported is None:
            failures.append(f"{rel}: no __all__ lists the public names")
            continue
        for name in sorted(set(defined) - exported):
            if not name.startswith("_"):
                failures.append(f"{rel}: public {name} is missing from __all__")
        for name in sorted(exported):
            named = rf"~(?:{re.escape(module)}\.)?{re.escape(name)}"
            if not re.search(rf"{named}\b", reference):
                failures.append(f"reference/python.org: {module}.{name} is not named")
            for member in sorted(defined.get(name, ())):
                if not re.search(rf"{named}\.{re.escape(member)}\b", reference):
                    failures.append(
                        f"reference/python.org: {module}.{name}.{member} is not named"
                    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_reference_contract.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    failures: list[str] = []
    public = check_reference(root, failures)
    check_mentions(root, public, failures)
    check_python(root, failures)
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"{len(public)} public C names documented")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
