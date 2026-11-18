#!/usr/bin/env python3
"""Every public name has a contract entry, and every entry names a real symbol.

``reference/contracts.org`` holds one entry per public name:

* each function of the installed C headers, with the rows ``Units``,
  ``Layout``, ``Ownership``, ``Threads and MPI``, ``Returns and errors``,
  ``Python`` and ``Defaults``;
* each type of those headers, with ``Layout``, ``Units``, ``Ownership`` and
  ``Python``;
* each macro and enum constant, in the constants table, with its value;
* each name of ``__all__`` in every module of ``python/cpmdc``, and each
  public member of an exported class, with ``Signature``, ``Units``,
  ``Shape and dtype``, ``Ownership and threads``, ``Raises`` and
  ``Defaults``;
* each property that ``CPMDCCalculator`` implements and each parameter it
  takes, in the ASE tables.

The C side is read from the headers that ``meson install`` places (the
install plan of the build directory), the Python side from the imported
package. The signature row of a Python entry must equal the signature
``inspect`` reports, so a changed default or a renamed parameter fails here
until the page follows. The ``Raises`` row must name every exception that
the function raises by name. An entry whose symbol is gone fails too.

usage: test_reference_audit.py SOURCE_ROOT [BUILD_DIR]
"""

from __future__ import annotations

import ast
import builtins
import importlib
import inspect
import json
import re
import subprocess
import sys
import warnings
from pathlib import Path

PAGE = "docs/orgmode/reference/contracts.org"
C_PAGE = "docs/orgmode/reference/c-abi.org"

C_PREFIXES = ("cpmdc_", "cp2k_")
FALLBACK_HEADERS = ("include/*.h", "tools/gromacs/libcp2k.h")
SKIPPED_HEADERS = {"cpmdc_pwscale.h"}

C_FUNCTION_ROWS = (
    "Units",
    "Layout",
    "Ownership",
    "Threads and MPI",
    "Returns and errors",
    "Python",
    "Defaults",
)
C_TYPE_ROWS = ("Layout", "Units", "Ownership", "Python")
PY_ROWS = (
    "Signature",
    "Units",
    "Shape and dtype",
    "Ownership and threads",
    "Raises",
    "Defaults",
)
CONSTANT_HEAD = ("Name", "Value", "Meaning", "Python")
PROPERTY_HEAD = ("Property", "Units", "Shape and dtype", "Sign and convention", "Raises")
PARAMETER_HEAD = ("Parameter", "Units", "Type", "Default", "Meaning")

UNIT_WORDS = re.compile(
    r"hartree|bohr|angstrom|rydberg|\bry\b|\bev\b|atomic unit|electrons?|"
    r"\bnm\b|kj/mol|kcal/mol|dimensionless|seconds?|bytes?|\bps\b|"
    r"lammps unit|per atom type",
    re.I,
)
OPAQUE_TYPES = {"CPMDCSession", "CPMDCRestart", "cpmdc_restart"}
NON_QUANTITY_STRUCTS = OPAQUE_TYPES | {"CPMDCFeatureEntry"}
EXCEPTION_NAME = re.compile(r"\b[A-Z]\w*(?:Error|Failed|Warning|Exception|Exit)\b")
WORD_NONE = re.compile(r"^\s*none\b", re.I)


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


# --- C headers ----------------------------------------------------------------


class CFunction:
    def __init__(self, name: str, ret: str, params: list[tuple[str, str]], header: str) -> None:
        self.name = name
        self.ret = ret
        self.params = params
        self.header = header

    @property
    def return_kind(self) -> str:
        ret = self.ret
        if "*" in ret:
            return "pointer"
        if ret == "void":
            return "void"
        if ret == "int":
            return "int"
        if ret == "size_t":
            return "size_t"
        if ret == "CPMDCResult":
            return "CPMDCResult"
        return "other"

    @property
    def pointer_params(self) -> list[str]:
        return [name for name, ctype in self.params if "*" in ctype or "[" in ctype]

    @property
    def layout_params(self) -> list[str]:
        out = []
        for name, ctype in self.params:
            base = re.sub(r"\b(const|struct)\b|\*|\[.*\]", " ", ctype).split()
            if ("*" in ctype or "[" in ctype) and not (base and base[0] in OPAQUE_TYPES):
                out.append(name)
        return out

    @property
    def has_quantity(self) -> bool:
        if self.return_kind == "CPMDCResult":
            return True
        for _name, ctype in self.params:
            words = re.sub(r"\b(const|struct)\b|\*|\[.*\]", " ", ctype).split()
            if not words:
                continue
            if words[0] == "double":
                return True
            if words[0].startswith("CPMDC") and words[0] not in NON_QUANTITY_STRUCTS and "*" in ctype:
                return True
        return False


class CHeaders:
    def __init__(self) -> None:
        self.functions: dict[str, CFunction] = {}
        self.types: dict[str, str] = {}
        self.constants: dict[str, str] = {}
        self.fields: set[tuple[str, str, str, str]] = set()
        self.files: list[str] = []


def split_params(args: str) -> list[tuple[str, str]]:
    args = args.strip()
    if args in ("", "void"):
        return []
    out = []
    for chunk in args.split(","):
        chunk = " ".join(chunk.split())
        dims = re.search(r"\[[^\]]*\]\s*$", chunk)
        stem = chunk[: dims.start()] if dims else chunk
        found = re.search(r"([A-Za-z_]\w*)\s*$", stem)
        if not found:
            raise ValueError(f"cannot read the parameter {chunk!r}")
        name = found.group(1)
        ctype = stem[: found.start()].strip() + (" []" if dims else "")
        out.append((name, ctype))
    return out


def parse_header(text: str, rel: str, into: CHeaders) -> None:
    code = strip_comments(text)
    macros = re.findall(r"(?m)^[ \t]*#[ \t]*define[ \t]+(\w+)(?:[ \t]+([^\n]*))?$", code)
    for name, value in macros:
        if name.endswith("_H") or "(" in name:
            continue
        if name.startswith(("CPMDC", "cp2k")):
            into.constants.setdefault(name, value.strip())
    for name, value in re.findall(r"\b(CPMDC_[A-Z0-9_]+)\s*=\s*([^,}\n]+)", code):
        into.constants[name] = value.strip()
    flat = re.sub(r"(?m)^[ \t]*#[^\n]*$", " ", code)
    flat = re.sub(r'extern\s+"C"\s*\{', " ", flat)
    for keyword, body, name in re.findall(
        r"\btypedef\s+(struct|enum)(?:\s+\w+)?\s*\{([^{}]*)\}\s*(\w+)\s*;", flat
    ):
        into.types[name] = keyword
        if keyword != "struct":
            continue
        for declaration in body.split(";"):
            declaration = " ".join(declaration.split())
            if not declaration:
                continue
            dims = re.search(r"\[([^\]]*)\]\s*$", declaration)
            stem = re.sub(r"\[[^\]]*\]\s*$", "", declaration)
            found = re.search(r"([A-Za-z_]\w*)\s*$", stem)
            if found:
                ctype = stem[: found.start()].strip()
                into.fields.add((name, found.group(1), ctype, dims.group(1) if dims else ""))
    for old, name in re.findall(r"\btypedef\s+struct\s+(\w+)\s+(\w+)\s*;", flat):
        into.types[name] = "opaque"
    for old, name in re.findall(r"\btypedef\s+(int|long|unsigned)\s+(\w+)\s*;", flat):
        into.types[name] = "handle"
    for ret, name, args in re.findall(
        r"((?:[A-Za-z_][\w \t\*]*?))\b((?:cpmdc|cp2k)_\w+)\s*\(([^()]*)\)\s*;", flat
    ):
        ret = re.sub(r"\bCPMDC_API\b", " ", ret)
        ret = " ".join(ret.replace("*", " * ").split()).replace(" * ", " *").replace("* ", "*")
        if "typedef" in ret:
            continue
        ret = re.sub(r"^(?:extern|static)\s+", "", ret)
        into.functions[name] = CFunction(name, ret.replace("const char*", "const char *"), split_params(args), rel)


def header_paths(root: Path, build: Path | None) -> list[Path]:
    plan = build / "meson-info" / "intro-install_plan.json" if build else None
    if plan is not None and plan.is_file():
        data = json.loads(plan.read_text(encoding="utf-8"))
        found = []
        for source, entry in sorted(data.get("headers", {}).items()):
            dest = entry.get("destination", "")
            if "/cpmdc/" in dest and dest.endswith(".h"):
                found.append(Path(source))
        if found:
            return found
    out: list[Path] = []
    for pattern in FALLBACK_HEADERS:
        out.extend(sorted(root.glob(pattern)))
    return out


def read_headers(root: Path, build: Path | None) -> CHeaders:
    into = CHeaders()
    for path in header_paths(root, build):
        if path.name in SKIPPED_HEADERS:
            continue
        rel = path.relative_to(root).as_posix() if path.is_relative_to(root) else path.name
        into.files.append(rel)
        parse_header(path.read_text(encoding="utf-8"), rel, into)
    return into


# --- the reference page -------------------------------------------------------


class Page:
    def __init__(self, text: str) -> None:
        self.sections: dict[str, dict[str, str]] = {}
        self.tables: dict[str, list[list[str]]] = {}
        title = ""
        heading = re.compile(r"^\*+\s+(.*?)\s*$")
        rows: list[str] = []

        def flush() -> None:
            if not title:
                return
            table = [parse_row(line) for line in rows]
            table = [row for row in table if not set("".join(row)) <= set("-+ ")]
            if not table:
                return
            self.tables[title] = table
            self.sections[title] = {
                row[0]: row[1] for row in table[1:] if len(row) >= 2
            } if len(table[0]) == 2 and table[0][0] == "Item" else {}

        for line in text.splitlines():
            match = heading.match(line)
            if match:
                flush()
                rows = []
                title = match.group(1).replace("~", "").strip()
                continue
            if line.lstrip().startswith("|"):
                rows.append(line)
            elif rows and line.strip() == "":
                continue
        flush()


def parse_row(line: str) -> list[str]:
    if re.match(r"^\s*\|[-+]+\|?\s*$", line):
        return ["---"]
    cells = line.strip().strip("|").split("|")
    return [cell.strip() for cell in cells]


def clean(cell: str) -> str:
    return re.sub(r"~([^~\n]*)~", r"\1", cell).strip()


# --- Python surface -----------------------------------------------------------


def module_names(root: Path) -> list[str]:
    names = []
    for path in sorted((root / "python/cpmdc").glob("*.py")):
        if path.name == "__init__.py":
            names.append("cpmdc")
        elif not path.name.startswith("_"):
            names.append(f"cpmdc.{path.stem}")
    return names


def init_attributes(cls: type) -> set[str]:
    try:
        tree = ast.parse(_dedent(inspect.getsource(cls)))
    except (OSError, TypeError, SyntaxError):
        return set()
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "__init__":
            for sub in ast.walk(node):
                targets: list[ast.expr] = []
                if isinstance(sub, ast.Assign):
                    targets = list(sub.targets)
                elif isinstance(sub, (ast.AnnAssign, ast.AugAssign)):
                    targets = [sub.target]
                for target in targets:
                    if (
                        isinstance(target, ast.Attribute)
                        and isinstance(target.value, ast.Name)
                        and target.value.id == "self"
                        and not target.attr.startswith("_")
                    ):
                        found.add(target.attr)
    return found


def _dedent(source: str) -> str:
    lines = source.splitlines()
    indent = min((len(line) - len(line.lstrip()) for line in lines if line.strip()), default=0)
    return "\n".join(line[indent:] for line in lines)


def class_members(cls: type) -> list[str]:
    own = {k for k in vars(cls) if not k.startswith("_")}
    own |= {k for k in getattr(cls, "__annotations__", {}) if not k.startswith("_")}
    own |= init_attributes(cls)
    return sorted(own)


def render_signature(obj: object, name: str, drop_first: bool) -> str:
    sig = inspect.signature(obj)  # type: ignore[arg-type]
    parts: list[str] = []
    seen_star = False
    items = list(sig.parameters.values())
    if drop_first and items and items[0].name in ("self", "cls"):
        items = items[1:]
    for index, param in enumerate(items):
        if param.kind is inspect.Parameter.KEYWORD_ONLY and not seen_star:
            parts.append("*")
            seen_star = True
        if param.kind is inspect.Parameter.VAR_POSITIONAL:
            seen_star = True
            text = f"*{param.name}"
        elif param.kind is inspect.Parameter.VAR_KEYWORD:
            text = f"**{param.name}"
        elif param.default is inspect.Parameter.empty:
            text = param.name
        else:
            text = f"{param.name}={param.default!r}"
        parts.append(text)
        del index
    return f"{name}({', '.join(parts)})"


def qualname_of(base: type) -> str:
    if base.__module__ == "builtins":
        return base.__name__
    return f"{base.__module__}.{base.__qualname__}"


def raised_names(func: object) -> set[str]:
    try:
        tree = ast.parse(_dedent(inspect.getsource(func)))  # type: ignore[arg-type]
    except (OSError, TypeError, SyntaxError):
        return set()
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Raise) and node.exc is not None:
            exc = node.exc
            if isinstance(exc, ast.Call):
                exc = exc.func
            if isinstance(exc, ast.Name) and EXCEPTION_NAME.fullmatch(exc.id):
                names.add(exc.id)
            elif isinstance(exc, ast.Attribute) and EXCEPTION_NAME.fullmatch(exc.attr):
                names.add(exc.attr)
    return names


def defined_public(root: Path, module: str) -> set[str]:
    """Names a module defines at top level without a leading underscore."""
    stem = "__init__" if module == "cpmdc" else module.split(".", 1)[1]
    tree = ast.parse((root / "python/cpmdc" / f"{stem}.py").read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return {name for name in names if not name.startswith("_")}


class PythonSurface:
    """Public names of the package with their kind, object, and members."""

    def __init__(self, root: Path) -> None:
        sys.path.insert(0, str(root / "python"))
        for stale in [m for m in sys.modules if m == "cpmdc" or m.startswith("cpmdc.")]:
            del sys.modules[stale]
        self.entries: dict[str, tuple[str, object, object | None]] = {}
        self.namespaces: dict[str, object] = {}
        self.problems: list[str] = []
        self.skipped: list[str] = []
        for module in module_names(root):
            try:
                mod = importlib.import_module(module)
            except ImportError as exc:
                self.skipped.append(f"{module}: {exc}")
                continue
            exported = getattr(mod, "__all__", None)
            if exported is None:
                self.problems.append(f"{module}: no __all__ lists the public names")
                continue
            for name in sorted(defined_public(root, module) - set(exported)):
                self.problems.append(f"{module}: public {name} is missing from __all__")
            self.namespaces[module] = mod
            for name in exported:
                obj = getattr(mod, name)
                key = f"{module}.{name}"
                if inspect.isclass(obj):
                    kind = "exception" if issubclass(obj, BaseException) else "class"
                    self.entries[key] = (kind, obj, None)
                    for member in class_members(obj):
                        self.add_member(key, kind, obj, member)
                elif callable(obj):
                    self.entries[key] = ("function", obj, None)
                else:
                    self.entries[key] = ("constant", obj, None)

    def add_member(self, key: str, kind: str, cls: type, member: str) -> None:
        raw = inspect.getattr_static(cls, member, None)
        path = f"{key}.{member}"
        if isinstance(raw, property):
            self.entries[path] = ("property", raw.fget, cls)
        elif isinstance(raw, (classmethod, staticmethod)):
            self.entries[path] = ("classmethod" if isinstance(raw, classmethod) else "staticmethod", raw.__func__, cls)
        elif inspect.isfunction(raw):
            self.entries[path] = ("method", raw, cls)
        else:
            self.entries[path] = ("attribute", raw, cls)

    def exception_known(self, name: str) -> bool:
        short = name.rsplit(".", 1)[-1]
        if isinstance(getattr(builtins, short, None), type):
            return True
        for mod in self.namespaces.values():
            if isinstance(getattr(mod, short, None), type):
                return True
        for module in ("ase.calculators.calculator", "ase.utils"):
            try:
                if isinstance(getattr(importlib.import_module(module), short, None), type):
                    return True
            except ImportError:
                continue
        return False

    def resolves(self, dotted: str) -> bool:
        return dotted in self.entries or dotted in self.namespaces


# --- the audit ----------------------------------------------------------------


def check_functions(headers: CHeaders, page: Page, surface: PythonSurface, fails: list[str]) -> None:
    for name, func in sorted(headers.functions.items()):
        rows = page.sections.get(name)
        if rows is None:
            fails.append(f"{PAGE}: {name} ({func.header}) has no contract entry")
            continue
        for key in C_FUNCTION_ROWS:
            if not rows.get(key, "").strip():
                fails.append(f"{PAGE}: {name} has no {key!r} row")
        if any(not rows.get(key, "").strip() for key in C_FUNCTION_ROWS):
            continue
        units = clean(rows["Units"])
        if func.has_quantity:
            if not UNIT_WORDS.search(units):
                fails.append(f"{PAGE}: {name} takes or returns quantities but its Units row names no unit")
        elif not WORD_NONE.match(units):
            fails.append(f"{PAGE}: {name} has no quantities, so its Units row must start with 'none'")
        layout = clean(rows["Layout"])
        if func.layout_params:
            for param in func.layout_params:
                if not re.search(rf"\b{re.escape(param)}\b", layout):
                    fails.append(f"{PAGE}: {name} Layout row does not describe {param}")
        elif not WORD_NONE.match(layout):
            fails.append(f"{PAGE}: {name} has no array or string argument, so its Layout row must start with 'none'")
        ownership = clean(rows["Ownership"])
        if func.pointer_params or func.return_kind == "pointer":
            blanket = re.search(
                r"\b(?:every|all|each)\b[^.;]{0,24}\b(?:buffers?|arrays?|pointers?|arguments?|inputs?)\b",
                ownership,
                re.I,
            )
            for param in func.pointer_params:
                if not blanket and not re.search(rf"\b{re.escape(param)}\b", ownership):
                    fails.append(f"{PAGE}: {name} Ownership row does not name {param}")
            if func.return_kind == "pointer" and not re.search(r"return", ownership, re.I):
                fails.append(f"{PAGE}: {name} returns a pointer, so its Ownership row must say who owns it")
        elif not WORD_NONE.match(ownership):
            fails.append(f"{PAGE}: {name} has no pointer, so its Ownership row must start with 'none'")
        threads = clean(rows["Threads and MPI"])
        if not re.search(r"thread", threads, re.I) or not re.search(r"MPI|rank", threads):
            fails.append(f"{PAGE}: {name} Threads and MPI row must say both what threads and what MPI ranks may do")
        errors = clean(rows["Returns and errors"])
        needle = {"int": "-1", "pointer": "NULL", "size_t": "0", "CPMDCResult": "ok"}.get(func.return_kind)
        if needle is not None and not WORD_NONE.match(errors) and needle not in errors:
            fails.append(f"{PAGE}: {name} returns {func.ret!r} but its Returns and errors row never says {needle!r}")
        if not WORD_NONE.match(errors) and func.header.endswith(("cpmdc.h", "cpmdc_features.h")):
            if not re.search(r"last_error", errors):
                fails.append(f"{PAGE}: {name} Returns and errors row must say whether cpmdc_last_error() is written")
        python = clean(rows["Python"])
        for span in re.findall(r"[A-Za-z_]\w*(?:\.\w+)*", python):
            if EXCEPTION_NAME.fullmatch(span.rsplit(".", 1)[-1]) and surface.namespaces:
                if not surface.exception_known(span):
                    fails.append(f"{PAGE}: {name} Python row names {span}, which is no exception of the package")
            elif span.startswith("cpmdc.") and surface.namespaces and not surface.resolves(span):
                fails.append(f"{PAGE}: {name} Python row names {span}, which does not exist")
    for name in sorted(set(page.sections) - set(headers.functions)):
        if name.startswith(C_PREFIXES) and set(page.sections[name]) >= {"Returns and errors"}:
            fails.append(f"{PAGE}: entry {name} names a C function that no installed header declares")


def check_types(headers: CHeaders, page: Page, fails: list[str]) -> None:
    for name in sorted(headers.types):
        rows = page.sections.get(name)
        if rows is None:
            fails.append(f"{PAGE}: type {name} has no contract entry")
            continue
        for key in C_TYPE_ROWS:
            if not rows.get(key, "").strip():
                fails.append(f"{PAGE}: type {name} has no {key!r} row")
    for name, rows in sorted(page.sections.items()):
        if set(rows) >= set(C_TYPE_ROWS) and "Signature" not in rows and "Returns and errors" not in rows:
            if name not in headers.types:
                fails.append(f"{PAGE}: entry {name} names a C type that no installed header declares")


def check_constants(headers: CHeaders, page: Page, fails: list[str]) -> None:
    table = None
    for rows in page.tables.values():
        if tuple(rows[0][: len(CONSTANT_HEAD)]) == CONSTANT_HEAD:
            table = rows
    if table is None:
        fails.append(f"{PAGE}: no constants table with the columns {' | '.join(CONSTANT_HEAD)}")
        return
    documented = {clean(row[0]): row for row in table[1:] if len(row) >= 4}
    for name, value in sorted(headers.constants.items()):
        row = documented.get(name)
        if row is None:
            fails.append(f"{PAGE}: constant {name} has no row in the constants table")
            continue
        if not all(cell.strip() for cell in row):
            fails.append(f"{PAGE}: constants row {name} has an empty cell")
        numeric = re.fullmatch(r"-?\d[\d.eE+-]*", value.strip())
        if numeric and name != "CPMDC_ABI_VERSION" and value.strip() not in clean(row[1]).split():
            fails.append(f"{PAGE}: constant {name} is {value.strip()} in the header, but its Value cell says {clean(row[1])!r}")
    for name in sorted(set(documented) - set(headers.constants)):
        fails.append(f"{PAGE}: constants row {name} names a macro or enum constant that no installed header defines")


def check_fields(headers: CHeaders, fails: list[str], root: Path) -> None:
    text = (root / C_PAGE).read_text(encoding="utf-8")
    rows = [parse_row(line) for line in text.splitlines() if line.lstrip().startswith("|")]
    for struct, field, ctype, dims in sorted(headers.fields):
        if ctype.split()[-1] != "double":
            continue
        label = f"{field}[{dims}]" if dims else field
        matches = [
            " ".join(row)
            for row in rows
            if row and re.search(rf"~{re.escape(label)}~", row[0])
        ]
        if not matches:
            fails.append(f"{C_PAGE}: {struct}.{label} has no row that starts with it")
        elif ctype.split()[-1] == "double" and not any(UNIT_WORDS.search(row) for row in matches):
            fails.append(f"{C_PAGE}: {struct}.{label} is a double, but its row names no unit")


def py_unit_hint(name: str) -> bool:
    return bool(re.search(r"_(h|bohr|ry|ang|ev|nm|h_bohr|bohr3|ev_per_angstrom)\b|_to_", name))


def check_python(page: Page, surface: PythonSurface, fails: list[str]) -> None:
    fails.extend(surface.problems)
    for key, (kind, obj, owner) in sorted(surface.entries.items()):
        rows = page.sections.get(key)
        if rows is None:
            fails.append(f"{PAGE}: {key} ({kind}) has no contract entry")
            continue
        for row_name in PY_ROWS:
            if not rows.get(row_name, "").strip():
                fails.append(f"{PAGE}: {key} has no {row_name!r} row")
        if any(not rows.get(row_name, "").strip() for row_name in PY_ROWS):
            continue
        signature = clean(rows["Signature"])
        short = key.rsplit(".", 1)[-1]
        if kind == "exception":
            bases = ", ".join(qualname_of(base) for base in obj.__bases__)  # type: ignore[attr-defined]
            expect = f"{short}({bases})"
        elif kind == "class":
            try:
                expect = render_signature(obj.__init__, short, True)  # type: ignore[attr-defined]
            except (TypeError, ValueError):
                expect = f"{short}()"
        elif kind in ("function",):
            expect = render_signature(obj, short, False)
        elif kind in ("method", "classmethod", "staticmethod"):
            expect = render_signature(obj, short, kind != "staticmethod")
        elif kind == "constant":
            if isinstance(obj, (int, float)):
                expect = f"{short} = {obj!r}"
            else:
                expect = f"{short} ({type(obj).__name__})"
        elif kind == "property":
            expect = f"{short} (property)"
        else:
            expect = f"{short} (attribute)"
        if signature != expect:
            fails.append(f"{PAGE}: {key} Signature row is {signature!r}, the code has {expect!r}")
        units = clean(rows["Units"])
        if not WORD_NONE.match(units) and not UNIT_WORDS.search(units):
            fails.append(f"{PAGE}: {key} Units row names no unit; say 'none' when there is none")
        if kind in ("function", "method", "classmethod", "staticmethod") and py_unit_hint(short):
            if WORD_NONE.match(units):
                fails.append(f"{PAGE}: {key} has a unit in its name but its Units row says none")
        if kind in ("class", "function", "method", "exception"):
            raises = clean(rows["Raises"])
            if kind in ("function", "method", "class"):
                target = obj.__init__ if kind == "class" else obj  # type: ignore[attr-defined]
                for name in sorted(raised_names(target)):
                    if not re.search(rf"\b{re.escape(name)}\b", raises):
                        fails.append(f"{PAGE}: {key} raises {name}, which its Raises row does not name")
            for span in set(EXCEPTION_NAME.findall(raises)):
                if not surface.exception_known(span):
                    fails.append(f"{PAGE}: {key} Raises row names {span}, which is no exception")
        if kind in ("class", "function", "method"):
            if not re.search(r"thread|stateless|pure function", clean(rows["Ownership and threads"]), re.I):
                fails.append(f"{PAGE}: {key} Ownership and threads row says nothing about threads")
        if kind == "class" and "__enter__" in vars(obj):  # type: ignore[arg-type]
            if not any(re.search(r"context manager", cell, re.I) for cell in rows.values()):
                fails.append(f"{PAGE}: {key} supports 'with' but no row says it is a context manager")
    for name in sorted(page.sections):
        if name.startswith("cpmdc.") and set(page.sections[name]) >= set(PY_ROWS):
            if name not in surface.entries and surface.namespaces:
                module = name.rsplit(".", 1)[0]
                if module in surface.namespaces or any(m.startswith(module) for m in surface.namespaces):
                    fails.append(f"{PAGE}: entry {name} names a Python symbol that does not exist")


def check_ase(page: Page, surface: PythonSurface, fails: list[str]) -> None:
    entry = surface.entries.get("cpmdc.ase.CPMDCCalculator")
    if entry is None:
        return
    cls = entry[1]
    for title, head, wanted, label in (
        ("ASE properties", PROPERTY_HEAD, list(cls.implemented_properties), "implemented property"),  # type: ignore[attr-defined]
        ("ASE parameters", PARAMETER_HEAD, list(cls.default_parameters), "parameter"),  # type: ignore[attr-defined]
    ):
        table = page.tables.get(title)
        if table is None or tuple(table[0][: len(head)]) != head:
            fails.append(f"{PAGE}: section {title!r} needs a table with the columns {' | '.join(head)}")
            continue
        rows = {clean(row[0]): row for row in table[1:] if row and len(row) >= len(head)}
        for name in wanted:
            row = rows.get(name)
            if row is None:
                fails.append(f"{PAGE}: ASE {label} {name} has no row in {title!r}")
            elif not all(cell.strip() for cell in row):
                fails.append(f"{PAGE}: ASE {label} {name} has an empty cell")
        for name in sorted(set(rows) - set(wanted)):
            fails.append(f"{PAGE}: {title!r} lists {name}, which CPMDCCalculator does not have")
    table = page.tables.get("ASE parameters")
    if table is not None:
        for row in table[1:]:
            name = clean(row[0])
            if name in cls.default_parameters and len(row) >= 4:  # type: ignore[attr-defined]
                default = repr(cls.default_parameters[name])  # type: ignore[attr-defined]
                if clean(row[3]) != default:
                    fails.append(f"{PAGE}: ASE parameter {name} defaults to {default}, the page says {clean(row[3])!r}")
    props = page.tables.get("ASE properties")
    if props is not None:
        for row in props[1:]:
            if clean(row[0]) in cls.implemented_properties and not UNIT_WORDS.search(clean(row[1])):  # type: ignore[attr-defined]
                fails.append(f"{PAGE}: ASE property {clean(row[0])} names no unit")


def check_aliases(page: Page, surface: PythonSurface, fails: list[str]) -> None:
    table = page.tables.get("Deprecated aliases")
    if table is None:
        fails.append(f"{PAGE}: no 'Deprecated aliases' table")
        return
    for row in table[1:]:
        if len(row) < 3 or "---" in row:
            continue
        old, new = clean(row[0]), clean(row[1])
        if old.startswith("cpmdc_") or "." not in old:
            continue
        module, _, attr = old.rpartition(".")
        mod = surface.namespaces.get(module)
        if mod is None:
            continue
        target = surface.entries.get(new)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                got = getattr(mod, attr)
            except AttributeError:
                fails.append(f"{PAGE}: deprecated alias {old} does not exist")
                continue
        if not any(issubclass(w.category, DeprecationWarning) for w in caught):
            fails.append(f"{PAGE}: deprecated alias {old} does not warn with DeprecationWarning")
        if target is None or got is not target[1]:
            fails.append(f"{PAGE}: deprecated alias {old} is not the same object as {new}")


def check_exports(headers: CHeaders, build: Path | None, fails: list[str]) -> None:
    if build is None:
        return
    library = build / "libcpmdc.so"
    if not library.is_file():
        return
    try:
        out = subprocess.run(
            ["nm", "-D", "--defined-only", str(library)], capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return
    exported = {line.split()[-1] for line in out.splitlines() if line.split() and line.split()[-1].startswith("cpmdc_")}
    declared = {name for name, func in headers.functions.items() if func.header.endswith(("cpmdc.h", "cpmdc_features.h"))}
    for name in sorted(exported - declared):
        fails.append(f"libcpmdc.so exports {name}, which no installed header declares")
    for name in sorted(declared - exported):
        fails.append(f"{name} is declared in an installed header but libcpmdc.so does not export it")


def audit(root: Path, build: Path | None = None) -> list[str]:
    fails: list[str] = []
    page_path = root / PAGE
    if page_path.is_file():
        page = Page(page_path.read_text(encoding="utf-8"))
    else:
        fails.append(f"{PAGE} does not exist")
        page = Page("")
    headers = read_headers(root, build)
    if not headers.functions:
        return ["no public C function was read from the installed headers"]
    surface = PythonSurface(root)
    check_functions(headers, page, surface, fails)
    check_types(headers, page, fails)
    check_constants(headers, page, fails)
    check_fields(headers, fails, root)
    check_python(page, surface, fails)
    check_ase(page, surface, fails)
    check_aliases(page, surface, fails)
    check_exports(headers, build, fails)
    for note in surface.skipped:
        print(f"skipped {note}", file=sys.stderr)
    return fails


def main() -> int:
    if len(sys.argv) not in (2, 3):
        print(__doc__.splitlines()[-1], file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    build = Path(sys.argv[2]).resolve() if len(sys.argv) == 3 else None
    fails = audit(root, build)
    if fails:
        print("\n".join(fails), file=sys.stderr)
        print(f"{len(fails)} contract failures", file=sys.stderr)
        return 1
    print("every public name has a contract entry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
