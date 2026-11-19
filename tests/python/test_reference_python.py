"""docs/orgmode/reference/python.org matches the public API of cpmdc."""

from __future__ import annotations

import ast
import importlib
import inspect
import re
import warnings
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "orgmode" / "reference" / "python.org"
PACKAGE_DIR = ROOT / "python" / "cpmdc"
CLASSES = ("Session", "Potential", "Evaluation", "ConstantCutoff")
CPMDC_SYMBOL = re.compile(r"(?<![\w/.-])cpmdc(?:\.[A-Za-z_]\w*)+")


def _import_cpmdc() -> Any:
    try:
        return importlib.import_module("cpmdc")
    except ImportError as exc:
        pytest.skip(f"cannot import cpmdc: {exc}")


def _public_members(cls: type) -> dict[str, Any]:
    members: dict[str, Any] = {}
    for name in sorted(set(dir(cls)) | set(getattr(cls, "__annotations__", {}))):
        if name.startswith("_"):
            continue
        members[name] = inspect.getattr_static(cls, name, None)
    return members


def _ast_module_symbols(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    symbols: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_"):
            symbols.add(node.name)
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            symbols.add(node.name)
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and not item.name.startswith("_"):
                    symbols.add(f"{node.name}.{item.name}")
                    if item.name == "__init__":
                        for sub in ast.walk(item):
                            if isinstance(sub, ast.Assign):
                                for target in sub.targets:
                                    if (
                                        isinstance(target, ast.Attribute)
                                        and isinstance(target.value, ast.Name)
                                        and target.value.id == "self"
                                        and not target.attr.startswith("_")
                                    ):
                                        symbols.add(f"{node.name}.{target.attr}")
                elif isinstance(item, ast.Assign):
                    for target in item.targets:
                        if isinstance(target, ast.Name) and not target.id.startswith("_"):
                            symbols.add(f"{node.name}.{target.id}")
                elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
                    if not item.target.id.startswith("_"):
                        symbols.add(f"{node.name}.{item.target.id}")
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and not target.id.startswith("_") and target.id != "__all__":
                    symbols.add(target.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if not node.target.id.startswith("_"):
                symbols.add(node.target.id)
    return symbols


def _known_cpmdc_symbols(cpmdc_mod: Any) -> set[str]:
    known: set[str] = {"cpmdc", "cpmdc.CpmdcError"}
    for name in cpmdc_mod.__all__:
        known.add(f"cpmdc.{name}")
        obj = getattr(cpmdc_mod, name)
        if inspect.isclass(obj):
            for member in _public_members(obj):
                known.add(f"cpmdc.{name}.{member}")
    for path in sorted(PACKAGE_DIR.glob("*.py")):
        if path.name.startswith("_"):
            continue
        submod_name = f"cpmdc.{path.stem}"
        known.add(submod_name)
        for item in _ast_module_symbols(path):
            known.add(f"{submod_name}.{item}")
        try:
            submod = importlib.import_module(submod_name)
        except ImportError:
            continue
        for name in getattr(submod, "__all__", ()):
            known.add(f"{submod_name}.{name}")
            obj = getattr(submod, name, None)
            if inspect.isclass(obj):
                for member in _public_members(obj):
                    known.add(f"{submod_name}.{name}.{member}")
    return known


def _has_token(text: str, token: str) -> bool:
    return bool(re.search(rf"(?<![\w.]){re.escape(token)}(?![\w])", text))


def _rows_by_symbol(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line in text.splitlines():
        if not line.lstrip().startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 2 or not cells[0].startswith("~"):
            continue
        for span in re.findall(r"~([^~\n]+)~", cells[0]):
            key = span.split("(", 1)[0].strip().removeprefix("cpmdc.")
            if key:
                rows[key] = rows.get(key, "") + "\n" + line
    return rows


def check_reference_text(text: str, cpmdc_mod: Any) -> list[str]:
    failures: list[str] = []
    rows = _rows_by_symbol(text)

    for name in cpmdc_mod.__all__:
        if not _has_token(text, name):
            failures.append(f"cpmdc.__all__ name {name} is not mentioned in reference/python.org")
        obj = getattr(cpmdc_mod, name)
        if inspect.isfunction(obj):
            scope = rows.get(name, text)
            for param in inspect.signature(obj).parameters:
                if not _has_token(scope, param):
                    failures.append(
                        f"parameter {param} of cpmdc.{name} is not mentioned in reference/python.org"
                    )

    for cls_name in CLASSES:
        cls = getattr(cpmdc_mod, cls_name)
        init_sig = inspect.signature(cls.__init__)
        cls_scope = rows.get(cls_name, text)
        for param_name in init_sig.parameters:
            if param_name == "self":
                continue
            if not _has_token(cls_scope, param_name):
                failures.append(
                    f"constructor parameter {param_name} of {cls_name} is not mentioned in reference/python.org"
                )
        for member_name, raw in _public_members(cls).items():
            qualified = f"{cls_name}.{member_name}"
            if not _has_token(text, qualified):
                failures.append(f"{qualified} is not mentioned in reference/python.org")
            func: Any = None
            if isinstance(raw, (classmethod, staticmethod)):
                func = raw.__func__
            elif inspect.isfunction(raw):
                func = raw
            if func is not None:
                sig = inspect.signature(func)
                method_scope = rows.get(qualified, text)
                for param_name in sig.parameters:
                    if param_name in ("self", "cls"):
                        continue
                    if not _has_token(method_scope, param_name):
                        failures.append(
                            f"parameter {param_name} of {qualified} is not mentioned in reference/python.org"
                        )

    known = _known_cpmdc_symbols(cpmdc_mod)
    for match in sorted(set(CPMDC_SYMBOL.findall(text))):
        if match not in known:
            failures.append(f"reference/python.org names {match}, which does not exist")

    return failures


def test_python_reference_matches_cpmdc() -> None:
    cpmdc_mod = _import_cpmdc()
    if not PAGE.is_file():
        pytest.skip("docs/orgmode/reference/python.org is not present")
    text = PAGE.read_text(encoding="utf-8")
    failures = check_reference_text(text, cpmdc_mod)
    assert failures == [], "\n".join(failures)
    assert len(re.findall(r"\bCpmdcError\b", text)) == 1
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        alias = cpmdc_mod.CpmdcError
    assert alias is cpmdc_mod.CPMDCError
    assert any(issubclass(item.category, DeprecationWarning) for item in caught)


@pytest.mark.parametrize(
    ("removed", "expected"),
    [
        ("HARTREE_TO_KCAL_PER_MOL", "cpmdc.__all__ name HARTREE_TO_KCAL_PER_MOL"),
        ("Session.copy_orbitals", "Session.copy_orbitals is not mentioned"),
        ("positions_bohr", "parameter positions_bohr of Session.calculate_au"),
        ("from_key", "parameter from_key of Session.copy_orbitals"),
        ("cutoff_ry", "constructor parameter cutoff_ry of ConstantCutoff"),
        ("forces_h_bohr", "parameter forces_h_bohr of cpmdc.forces_to_ev_per_angstrom"),
    ],
)
def test_python_reference_check_bites_when_name_is_removed(removed: str, expected: str) -> None:
    cpmdc_mod = _import_cpmdc()
    if not PAGE.is_file():
        pytest.skip("docs/orgmode/reference/python.org is not present")
    text = PAGE.read_text(encoding="utf-8")
    mutated = re.sub(rf"(?<![\w.]){re.escape(removed)}(?![\w])", "removed_token", text)
    failures = check_reference_text(mutated, cpmdc_mod)
    assert any(expected in item for item in failures), failures


@pytest.mark.parametrize(
    "bogus",
    [
        "cpmdc.no_such_symbol",
        "cpmdc.Session.no_such_method",
        "cpmdc.ipi.no_such_helper",
    ],
)
def test_python_reference_check_rejects_nonexistent_cpmdc_symbol(bogus: str) -> None:
    cpmdc_mod = _import_cpmdc()
    if not PAGE.is_file():
        pytest.skip("docs/orgmode/reference/python.org is not present")
    text = PAGE.read_text(encoding="utf-8") + f"\n| ~{bogus}~ | bogus |\n"
    failures = check_reference_text(text, cpmdc_mod)
    assert f"reference/python.org names {bogus}, which does not exist" in failures
