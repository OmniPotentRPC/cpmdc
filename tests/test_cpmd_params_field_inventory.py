#!/usr/bin/env python3
"""Check that top-level CPMDParams fields are required as feature rows."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
import sys
import tempfile
from pathlib import Path


REMOVED_FEATURE = "params.functional"
DUPLICATED_FEATURE = "params.functional"


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


def write_inventory_with_duplicate_feature(
    repo: Path, tmpdir: Path, feature_id: str
) -> Path:
    data = json.loads((repo / "schema/inventory/cpmd_features.json").read_text())
    duplicate = next(
        (
            feature.copy()
            for feature in data["features"]
            if feature["feature_id"] == feature_id
        ),
        None,
    )
    if duplicate is None:
        raise RuntimeError(f"fixture feature not found: {feature_id}")
    data["features"].append(duplicate)
    path = tmpdir / "cpmd_features_duplicate.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def write_inventory_with_duplicate_list_entry(
    repo: Path, tmpdir: Path, key: str
) -> Path:
    data = json.loads((repo / "schema/inventory/cpmd_features.json").read_text())
    duplicate = data[key][0]
    if isinstance(duplicate, dict):
        duplicate = duplicate.copy()
    data[key].append(duplicate)
    path = tmpdir / f"cpmd_features_duplicate_{key}.json"
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
    with tempfile.TemporaryDirectory(prefix="cpmd-params-field-inventory-") as raw:
        tmpdir = Path(raw)
        missing_inventory = write_inventory_without_feature(
            repo, tmpdir, REMOVED_FEATURE
        )
        missing_code, missing_output = run_checker_with_inventory(
            checker, missing_inventory
        )
        duplicate_inventory = write_inventory_with_duplicate_feature(
            repo, tmpdir, DUPLICATED_FEATURE
        )
        duplicate_code, duplicate_output = run_checker_with_inventory(
            checker, duplicate_inventory
        )
        duplicate_lists = {
            "section_kinds": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "section_kinds"
                ),
            ),
            "params_fields": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "params_fields"
                ),
            ),
            "cpmd_sections": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "cpmd_sections"
                ),
            ),
            "abi_symbols": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "abi_symbols"
                ),
            ),
        }

    expected = "inventory missing top-level params feature from CPMDParams: functional"
    if missing_code == 0 or expected not in missing_output:
        print("expected top-level params feature inventory failure")
        print(missing_output)
        return 1

    expected_duplicate = "inventory feature_id duplicated: params.functional"
    if duplicate_code == 0 or expected_duplicate not in duplicate_output:
        print("expected duplicate feature_id inventory failure")
        print(duplicate_output)
        return 1

    expected_list_failures = {
        "section_kinds": "inventory section_kinds duplicated: generic",
        "params_fields": "inventory params_fields duplicated: functional",
        "cpmd_sections": "inventory cpmd_sections duplicated: ATOM",
        "abi_symbols": "inventory abi_symbols duplicated: cpmdc_set_params",
    }
    for key, expected_message in expected_list_failures.items():
        code, output = duplicate_lists[key]
        if code == 0 or expected_message not in output:
            print(f"expected duplicate {key} inventory failure")
            print(output)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
