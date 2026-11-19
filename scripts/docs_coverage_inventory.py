#!/usr/bin/env python3
"""Report documentation coverage for public C symbols, params, and Python exports."""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DOC_PAGES = {
    "c_abi": ROOT / "docs/orgmode/reference/c-abi.org",
    "contracts": ROOT / "docs/orgmode/reference/contracts.org",
    "python": ROOT / "docs/orgmode/reference/python.org",
    "cpmd_options": ROOT / "docs/orgmode/reference/cpmd-options.org",
    "embedding": ROOT / "docs/orgmode/howto/embedding.org",
    "pwscale": ROOT / "docs/orgmode/reference/pwscale.org",
}


def load_headers() -> set[str]:
    names: set[str] = set()
    for header in ("cpmdc.h", "cpmdc_features.h", "cpmdc_restart.h"):
        text = (ROOT / "include" / header).read_text(encoding="utf-8")
        text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
        text = re.sub(r"//[^\n]*", " ", text)
        names |= set(re.findall(r"\b(cpmdc_\w+)\s*\(", text))
    pw = ROOT / "include/cpmdc_pwscale.h"
    if pw.is_file():
        text = re.sub(r"/\*.*?\*/", " ", pw.read_text(encoding="utf-8"), flags=re.S)
        names |= set(re.findall(r"\b(cpmdc_pwscale_\w+)\s*\(", text))
    return names


def python_exported() -> set[str]:
    names: set[str] = set()
    for path in (ROOT / "python" / "cpmdc").glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if path.name == "__init__.py":
            for node in tree.body:
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id == "__all__":
                            if isinstance(node.value, (ast.List, ast.Tuple)):
                                for elt in node.value.elts:
                                    if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                                        names.add(elt.value)
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                names.add(node.name)
            if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
                names.add(f"{path.stem}.{node.name}" if path.name != "__init__.py" else node.name)
    names.add("CPMDCCalculator")
    return names


def param_feature_ids() -> list[str]:
    data = json.loads((ROOT / "schema/inventory/cpmd_features.json").read_text())
    return sorted(
        f["feature_id"]
        for f in data.get("features", [])
        if f.get("kind") == "params_field" or str(f.get("feature_id", "")).startswith("params.")
    )


def mention(haystack: str, needle: str) -> bool:
    return needle in haystack


def report() -> dict[str, object]:
    c_docs = (
        DOC_PAGES["c_abi"].read_text(encoding="utf-8")
        + DOC_PAGES["contracts"].read_text(encoding="utf-8")
        + DOC_PAGES["pwscale"].read_text(encoding="utf-8")
    )
    py_docs = (
        DOC_PAGES["python"].read_text(encoding="utf-8")
        + (ROOT / "docs/orgmode/howto/python.org").read_text(encoding="utf-8")
    )
    opts = DOC_PAGES["cpmd_options"].read_text(encoding="utf-8")
    embed = DOC_PAGES["embedding"].read_text(encoding="utf-8")
    architecture = (ROOT / "docs/orgmode/explanation/architecture.org").read_text(
        encoding="utf-8"
    )

    c_symbols = sorted(load_headers())
    c_missing = [s for s in c_symbols if not mention(c_docs, s)]

    py_names = sorted(python_exported())
    py_missing = [n for n in py_names if not mention(py_docs, n.split(".")[-1])]

    params = param_feature_ids()
    params_missing = [p for p in params if not mention(opts, p)]

    embed_terms = [
        "cpmd_embed_c_api.F90",
        "cpmdc_embed_host",
        "memfd",
        "cpmdc_adopt_calculator_comm",
        "cpmdc_restart_read_mem",
        "cpmdc_pwscale_eval_cuda",
    ]
    embed_corpus = embed + architecture + c_docs
    embed_missing = [t for t in embed_terms if not mention(embed_corpus, t)]

    return {
        "c_symbols_total": len(c_symbols),
        "c_symbols_missing_docs": c_missing,
        "python_names_total": len(py_names),
        "python_names_missing_docs": py_missing,
        "params_total": len(params),
        "params_missing_cpmd_options": params_missing,
        "embedding_terms_missing": embed_missing,
    }


def main() -> int:
    data = report()
    gaps = (
        data["c_symbols_missing_docs"]
        or data["python_names_missing_docs"]
        or data["params_missing_cpmd_options"]
        or data["embedding_terms_missing"]
    )
    print(json.dumps(data, indent=2))
    return 1 if gaps else 0


if __name__ == "__main__":
    raise SystemExit(main())