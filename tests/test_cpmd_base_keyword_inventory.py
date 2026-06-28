#!/usr/bin/env python3
"""Check that base &CPMD keyword rows are required by the inventory audit."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
import sys
import tempfile
from pathlib import Path


REMOVED_FEATURE = "catalog.cpmd.OPTIMIZE_WAVEFUNCTION"


def load_checker(repo: Path):
    path = repo / "tools/check_feature_inventory.py"
    spec = importlib.util.spec_from_file_location("check_feature_inventory", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inventory_without_feature(repo: Path, tmpdir: Path, feature_id: str) -> Path:
    data = json.loads((repo / "schema/inventory/cpmd_features.json").read_text())
    kept = [feature for feature in data["features"] if feature["feature_id"] != feature_id]
    if len(kept) == len(data["features"]):
        raise RuntimeError(f"fixture feature not found: {feature_id}")
    data["features"] = kept
    path = tmpdir / "cpmd_features.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def run_checker_with_inventory(checker, inventory: Path) -> tuple[int, str]:
    old_inventory = checker.INVENTORY
    old_cpmd_root = os.environ.pop("CPMD_ROOT", None)
    old_cpmdc_root = os.environ.pop("CPMDC_CPMD_ROOT", None)
    stream = io.StringIO()
    try:
        checker.INVENTORY = inventory
        with contextlib.redirect_stdout(stream):
            code = checker.main()
    finally:
        checker.INVENTORY = old_inventory
        if old_cpmd_root is not None:
            os.environ["CPMD_ROOT"] = old_cpmd_root
        if old_cpmdc_root is not None:
            os.environ["CPMDC_CPMD_ROOT"] = old_cpmdc_root
    return code, stream.getvalue()


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    checker = load_checker(repo)
    with tempfile.TemporaryDirectory(prefix="cpmd-keyword-inventory-") as raw:
        inventory = write_inventory_without_feature(repo, Path(raw), REMOVED_FEATURE)
        code, output = run_checker_with_inventory(checker, inventory)

    expected = (
        "inventory missing CPMD base keyword from cpmd_cp_keywords.txt: "
        "OPTIMIZE WAVEFUNCTION"
    )
    if code == 0 or expected not in output:
        print("expected base keyword inventory failure")
        print(output)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
