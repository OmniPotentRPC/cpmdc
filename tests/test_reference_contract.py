#!/usr/bin/env python3
"""The reference documents every public name, and the docs name no other.

Every header of include/ must be assigned to a reference page in PAGES. Every
function, type, enum constant, macro, and struct field of a header must appear
in its page. A function must appear with its parameter list in a table row,
which is where the pages state a contract. A struct field, enum constant, or
macro must start a table row, so that a mention in prose does not count. The
page must be linked from the toctree of docs/orgmode/index.org.
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
import tempfile
from pathlib import Path

PAGES = {
    "include/cpmdc.h": "c-abi",
    "include/cpmdc_export.h": "c-abi",
    "include/cpmdc_features.h": "c-abi",
    "include/cpmdc_restart.h": "c-abi",
    "include/cpmdc_pwscale.h": "pwscale",
    "include/cpmdc_csvr.h": "csvr",
    "include/cpmdc_damp.h": "damp",
}
# Written to stderr under CPMDC_SCF_STEPS=1; see reference/environment.org.
LOG_TOKENS = {"cpmdc_scf_steps", "cpmdc_warm_start", "cpmdc_wfn_extrapolated"}


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def struct_fields(code: str) -> set[tuple[str, str]]:
    fields = set()
    structs = re.findall(
        r"\btypedef\s+struct(?:\s+\w+)?\s*\{([^{}]*)\}\s*(\w+)\s*;",
        code,
    )
    for body, struct in structs:
        for declaration in body.split(";"):
            declaration = declaration.strip()
            if not declaration:
                continue
            pointer = re.search(r"\(\s*\*\s*(\w+)\s*\)", declaration)
            if pointer:
                fields.add((struct, pointer.group(1)))
                continue
            for declarator in re.sub(r"\[[^\]]*\]", "", declaration).split(","):
                found = re.findall(r"\w+", declarator)
                if found:
                    fields.add((struct, found[-1]))
    return fields


def header_names(text: str) -> dict[str, set[str] | set[tuple[str, str]]]:
    code = strip_comments(text)
    return {
        "struct field": struct_fields(code),
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


ROW_KINDS = {"struct field", "enum constant", "macro"}


def row_names(rows: list[str]) -> set[str]:
    """Names in the ~code~ spans of the first cell of each table row."""
    names: set[str] = set()
    for row in rows:
        cells = row.strip().strip("|").split("|")
        for span in re.findall(r"~([^~\n]+)~", cells[0]):
            names |= set(re.findall(r"\w+", span))
    return names


def _table_cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def _is_separator_row(cells: list[str]) -> bool:
    if not cells:
        return False
    stripped = [cell.replace(":", "").strip() for cell in cells if cell.strip()]
    return bool(stripped) and all(
        re.fullmatch(r"-+(?:\+-+)*", cell) for cell in stripped
    )


def _is_field_table_header(cells: list[str]) -> bool:
    return bool(cells) and cells[0] == "Field"


def _field_names_from_cell(cell: str) -> set[str]:
    fields: set[str] = set()
    for span in re.findall(r"~([^~\n]+)~", cell):
        found = re.search(r"\w+", span)
        if found:
            fields.add(found.group())
    return fields


def _structs_in_text(fragment: str, structs: set[str]) -> set[str]:
    return set(re.findall(r"~(\w+)~", fragment)) & structs


def _singular_struct(candidates: set[str]) -> str | None:
    if len(candidates) == 1:
        return next(iter(candidates))
    return None


def _intro_struct(lines: list[str], structs: set[str]) -> str | None:
    for line in reversed(lines):
        if line.lstrip().startswith("Fields of "):
            found = _singular_struct(_structs_in_text(line, structs))
            if found is not None:
                return found
    return None


_ROW_STRUCT_REF = re.compile(r"(?:see|refers\s+to)\s+~(\w+)~", re.I)


def _row_struct_references(cells: list[str], structs: set[str]) -> set[str]:
    owners: set[str] = set()
    for cell in cells:
        match = _ROW_STRUCT_REF.fullmatch(cell.strip())
        if match and match.group(1) in structs:
            owners.add(match.group(1))
    return owners


def documented_struct_fields(text: str, structs: set[str]) -> set[tuple[str, str]]:
    """Struct-field pairs represented by field tables (header row starts with Field)."""
    pairs: set[tuple[str, str]] = set()
    lines = text.splitlines()
    last_heading_struct: str | None = None
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("*"):
            heading_names = set(re.findall(r"~(\w+)~", line)) & structs
            last_heading_struct = next(iter(heading_names), None) if heading_names else None
        if not line.lstrip().startswith("|"):
            index += 1
            continue
        header_cells = _table_cells(line)
        if not _is_field_table_header(header_cells):
            index += 1
            continue
        table_start = index
        has_struct_col = any(cell == "Struct" for cell in header_cells)
        struct_col = header_cells.index("Struct") if has_struct_col else None
        table_lines = [line]
        index += 1
        while index < len(lines) and lines[index].lstrip().startswith("|"):
            row_cells = _table_cells(lines[index])
            if _is_field_table_header(row_cells):
                break
            table_lines.append(lines[index])
            index += 1
        context_lines: list[str] = []
        scan = table_start - 1
        while scan >= 0:
            prior = lines[scan]
            if prior.lstrip().startswith("|") or prior.startswith("*"):
                break
            context_lines.insert(0, prior)
            scan -= 1
        context_block = "\n".join(context_lines)
        context_structs = _structs_in_text(context_block, structs)
        default_owner = (
            _intro_struct(context_lines, structs)
            or last_heading_struct
            or _singular_struct(context_structs)
        )
        for table_line in table_lines[1:]:
            row_cells = _table_cells(table_line)
            if _is_separator_row(row_cells) or _is_field_table_header(row_cells):
                continue
            if not row_cells or not row_cells[0].startswith("~"):
                continue
            fields = _field_names_from_cell(row_cells[0])
            if not fields:
                continue
            owners: set[str] = set()
            if struct_col is not None and struct_col < len(row_cells):
                owners = _structs_in_text(row_cells[struct_col], structs)
            elif struct_col is None:
                owners = _row_struct_references(row_cells[1:], structs)
            if not owners and default_owner is not None:
                owners = {default_owner}
            pairs |= {(owner, field) for owner in owners for field in fields}
    return pairs


def check_reference(root: Path, failures: list[str]) -> set[str]:
    pages: dict[str, tuple[str, list[str]] | None] = {}
    first_cells: dict[str, set[str]] = {}
    index = (root / "docs/orgmode/index.org").read_text(encoding="utf-8")
    public: set[str] = set()
    headers = sorted(path.relative_to(root).as_posix() for path in (root / "include").glob("*.h"))
    names_by_header = {
        rel: header_names((root / rel).read_text(encoding="utf-8")) for rel in headers
    }
    page_struct_fields: dict[str, set[tuple[str, str]]] = {}
    for rel, names in names_by_header.items():
        if rel in PAGES:
            page = f"reference/{PAGES[rel]}.org"
            page_struct_fields.setdefault(page, set()).update(names["struct field"])
    for rel in headers:
        if rel not in PAGES:
            failures.append(f"{rel}: no reference page is assigned in PAGES")
    for rel in sorted(set(PAGES) - set(headers)):
        failures.append(f"{rel}: listed in PAGES but not in include/")
    for rel in headers:
        stem = PAGES.get(rel)
        if stem is None:
            continue
        page = f"reference/{stem}.org"
        names = names_by_header[rel]
        for kind, found in names.items():
            if kind != "struct field":
                public |= found
        if page not in pages:
            source = root / "docs/orgmode" / page
            if not source.is_file():
                failures.append(f"{rel}: its page docs/orgmode/{page} does not exist")
                pages[page] = None
                continue
            text = source.read_text(encoding="utf-8")
            pages[page] = (text, table_rows(text))
            first_cells[page] = row_names(pages[page][1])
            if not re.search(rf"^\s+reference/{re.escape(stem)}\s*$", index, re.M):
                failures.append(f"index.org: toctree does not list reference/{stem}")
        if pages[page] is None:
            continue
        text, rows = pages[page]
        for kind, found in names.items():
            for name in sorted(found):
                if kind == "function":
                    if not any(re.search(rf"\b{name}\(\w", row) for row in rows):
                        failures.append(
                            f"{page}: {name} ({rel}) has no table row with its parameter list"
                        )
                elif kind == "struct field":
                    documented = documented_struct_fields(
                        text,
                        {struct for struct, _field in page_struct_fields[page]},
                    )
                    if name not in documented:
                        struct, field = name
                        failures.append(
                            f"{page}: struct field {struct}.{field} ({rel}) "
                            "has no table row for that struct"
                        )
                elif kind in ROW_KINDS:
                    if name not in first_cells[page]:
                        failures.append(
                            f"{page}: {kind} {name} ({rel}) has no table row that starts with it"
                        )
                elif not re.search(rf"~[^~\n]*\b{name}\b[^~\n]*~", text):
                    failures.append(f"{page}: {kind} {name} ({rel}) is not named")
    for page, expected in page_struct_fields.items():
        cached = pages.get(page)
        if cached is None:
            continue
        documented = documented_struct_fields(
            cached[0],
            {struct for struct, _field in expected},
        )
        for struct, field in sorted(documented - expected):
            failures.append(
                f"{page}: table documents nonexistent struct field {struct}.{field}"
            )
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


def test_org_separator_rows() -> None:
    for row in ("|---|", "|---+---|", "|-----+-----+----|"):
        assert _is_separator_row(_table_cells(row))
    assert not _is_separator_row(_table_cells("| Field | Type | Meaning |"))
    assert not _is_separator_row(_table_cells("| ~alpha~ | ~int~ | first member |"))
    text = (
        "Fields of ~DemoStruct~:\n\n"
        "| Field | Type | Meaning |\n"
        "|---|\n"
        "|---+---|\n"
        "|-----+-----+----|\n"
        "| ~alpha~ | ~int~ | first member |\n"
    )
    assert documented_struct_fields(text, {"DemoStruct"}) == {("DemoStruct", "alpha")}


def test_field_table_with_nonstandard_intro() -> None:
    text = (
        "Column overview for ~DemoStruct~ without the usual opening line:\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | ~int~ | first member |\n"
    )
    documented = documented_struct_fields(text, {"DemoStruct"})
    assert documented == {("DemoStruct", "alpha")}


def test_bogus_field_in_nonstandard_intro_table() -> None:
    valid_page = (
        "Column overview for ~DemoStruct~ without the usual opening line:\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | ~int~ | first member |\n"
    )
    bogus_page = valid_page + "| ~bogus~ | ~int~ | not in the header |\n"
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "include").mkdir(parents=True)
        (root / "docs/orgmode/reference").mkdir(parents=True)
        for rel in PAGES:
            (root / rel).write_text("", encoding="utf-8")
        (root / "include/cpmdc.h").write_text(
            "typedef struct DemoStruct { int alpha; } DemoStruct;\n",
            encoding="utf-8",
        )
        (root / "docs/orgmode/index.org").write_text(
            "\n".join(f"   reference/{stem}" for stem in sorted(set(PAGES.values())))
            + "\n",
            encoding="utf-8",
        )
        for stem in sorted(set(PAGES.values())):
            (root / "docs/orgmode/reference" / f"{stem}.org").write_text("", encoding="utf-8")
        c_abi = root / "docs/orgmode/reference/c-abi.org"
        c_abi.write_text(valid_page, encoding="utf-8")
        clean_failures: list[str] = []
        check_reference(root, clean_failures)
        assert not clean_failures, clean_failures
        c_abi.write_text(bogus_page, encoding="utf-8")
        failures: list[str] = []
        check_reference(root, failures)
        assert failures == [
            "reference/c-abi.org: table documents nonexistent struct field DemoStruct.bogus"
        ]


def test_meaning_cell_cross_reference_overrides_default_owner() -> None:
    text = (
        "* ~HeadingStruct~\n\n"
        "Fields of ~IntroStruct~:\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | ~int~ | first member |\n"
        "| ~beta~ | ~int~ | refers to ~OtherStruct~ |\n"
    )
    documented = documented_struct_fields(
        text,
        {"HeadingStruct", "IntroStruct", "OtherStruct"},
    )
    assert documented == {("IntroStruct", "alpha"), ("OtherStruct", "beta")}
    heading_only = (
        "* ~HeadingStruct~\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | ~int~ | first member |\n"
        "| ~beta~ | ~int~ | refers to ~OtherStruct~ |\n"
    )
    assert documented_struct_fields(
        heading_only,
        {"HeadingStruct", "OtherStruct"},
    ) == {("HeadingStruct", "alpha"), ("OtherStruct", "beta")}
    explicit_col = (
        "* ~HeadingStruct~\n\n"
        "Fields of ~IntroStruct~:\n\n"
        "| Field | Struct | Meaning |\n"
        "|-------+--------+---------|\n"
        "| ~alpha~ | ~ColumnStruct~ | refers to ~OtherStruct~ |\n"
    )
    assert documented_struct_fields(
        explicit_col,
        {"HeadingStruct", "IntroStruct", "ColumnStruct", "OtherStruct"},
    ) == {("ColumnStruct", "alpha")}


def test_row_type_or_meaning_struct_reference_without_struct_column() -> None:
    text = (
        "Fields of ~A~:\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | ~int~ | first member |\n"
        "| ~beta~ | see ~B~ | delegated in Type cell |\n"
        "| ~gamma~ | ~int~ | see ~B~ |\n"
    )
    assert documented_struct_fields(text, {"A", "B"}) == {
        ("A", "alpha"),
        ("B", "beta"),
        ("B", "gamma"),
    }


def test_prose_struct_mention_keeps_intro_owner() -> None:
    text = (
        "Fields of ~IntroStruct~:\n\n"
        "| Field | Type | Meaning |\n"
        "|-------+------+---------|\n"
        "| ~alpha~ | pointer to ~OtherStruct~ | first member |\n"
        "| ~beta~ | ~int~ | copied from ~OtherStruct~ on reset |\n"
    )
    documented = documented_struct_fields(
        text,
        {"IntroStruct", "OtherStruct"},
    )
    assert documented == {("IntroStruct", "alpha"), ("IntroStruct", "beta")}


def run_self_tests() -> None:
    tests = [
        obj
        for name, obj in globals().items()
        if name.startswith("test_") and callable(obj)
    ]
    assert tests, "no test_* functions discovered"
    for test in tests:
        test()


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_reference_contract.py SOURCE_ROOT", file=sys.stderr)
        return 2
    run_self_tests()
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
