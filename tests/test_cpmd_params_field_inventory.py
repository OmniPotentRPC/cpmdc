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
EXTRA_SECTION_FEATURE = "catalog.section.EXTRA"
DUPLICATE_C_FEATURE = "catalog.section.ATOM"


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


def write_inventory_with_removed_list_entry(
    repo: Path, tmpdir: Path, key: str
) -> Path:
    data = json.loads((repo / "schema/inventory/cpmd_features.json").read_text())
    data[key] = data[key][1:]
    path = tmpdir / f"cpmd_features_removed_{key}.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def write_inventory_with_extra_section_feature(repo: Path, tmpdir: Path) -> Path:
    data = json.loads((repo / "schema/inventory/cpmd_features.json").read_text())
    data["features"].append(
        {
            "feature_id": EXTRA_SECTION_FEATURE,
            "kind": "cpmd_section",
            "name": "EXTRA",
            "stub_applicable": True,
            "embed_applicable": True,
            "render_via": "set_or_generic_or_raw",
            "source": "inscan",
        }
    )
    path = tmpdir / "cpmd_features_extra_section.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def write_features_c_with_extra_section(repo: Path, tmpdir: Path) -> Path:
    text = (repo / "src/cpmdc_features.c").read_text(encoding="utf-8")
    row = f'    {{"{EXTRA_SECTION_FEATURE}", CPMDC_FEATURE_SECTION, 1, 1}},\n'
    path = tmpdir / "cpmdc_features_extra_section.c"
    path.write_text(
        text.replace(
            "\n};\n\n#undef CPMDC_PARAM_FEATURE",
            f"\n{row}}};\n\n#undef CPMDC_PARAM_FEATURE",
        ),
        encoding="utf-8",
    )
    return path


def write_features_c_with_duplicate_row(repo: Path, tmpdir: Path) -> Path:
    text = (repo / "src/cpmdc_features.c").read_text(encoding="utf-8")
    row = f'    {{"{DUPLICATE_C_FEATURE}", CPMDC_FEATURE_SECTION, 1, 1}},\n'
    path = tmpdir / "cpmdc_features_duplicate_row.c"
    path.write_text(
        text.replace(
            "\n};\n\n#undef CPMDC_PARAM_FEATURE",
            f"\n{row}}};\n\n#undef CPMDC_PARAM_FEATURE",
        ),
        encoding="utf-8",
    )
    return path


def write_allowlist_with_duplicate(repo: Path, tmpdir: Path) -> Path:
    text = (repo / "schema/inventory/opencpmd_sections.txt").read_text(
        encoding="utf-8"
    )
    first = next(line for line in text.splitlines() if line.strip())
    path = tmpdir / "opencpmd_sections_duplicate.txt"
    path.write_text(f"{text.rstrip()}\n{first}\n", encoding="utf-8")
    return path


def run_checker_with_inventory(checker, inventory: Path) -> tuple[int, str]:
    return run_checker_with_paths(checker, INVENTORY=inventory)


def run_checker_with_paths(checker, **paths: Path) -> tuple[int, str]:
    old_paths = {name: getattr(checker, name) for name in paths}
    old_cpmd_root = os.environ.pop("CPMD_ROOT", None)
    old_cpmdc_root = os.environ.pop("CPMDC_CPMD_ROOT", None)
    stream = io.StringIO()
    try:
        for name, path in paths.items():
            setattr(checker, name, path)
        with contextlib.redirect_stdout(stream):
            code = checker.main()
    finally:
        for name, path in old_paths.items():
            setattr(checker, name, path)
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
            "opencpmd_inscan_sections": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "opencpmd_inscan_sections"
                ),
            ),
            "abi_symbols": run_checker_with_inventory(
                checker,
                write_inventory_with_duplicate_list_entry(
                    repo, tmpdir, "abi_symbols"
                ),
            ),
        }
        mismatched_inscan_code, mismatched_inscan_output = run_checker_with_inventory(
            checker,
            write_inventory_with_removed_list_entry(
                repo, tmpdir, "opencpmd_inscan_sections"
            ),
        )
        extra_section_code, extra_section_output = run_checker_with_paths(
            checker,
            INVENTORY=write_inventory_with_extra_section_feature(repo, tmpdir),
            FEATURES_C=write_features_c_with_extra_section(repo, tmpdir),
        )
        duplicate_allow_code, duplicate_allow_output = run_checker_with_paths(
            checker,
            SEC_ALLOW=write_allowlist_with_duplicate(repo, tmpdir),
        )
        duplicate_c_row_code, duplicate_c_row_output = run_checker_with_paths(
            checker,
            FEATURES_C=write_features_c_with_duplicate_row(repo, tmpdir),
        )

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
        "opencpmd_inscan_sections": (
            "inventory opencpmd_inscan_sections duplicated: ATOM"
        ),
        "abi_symbols": "inventory abi_symbols duplicated: cpmdc_set_params",
    }
    for key, expected_message in expected_list_failures.items():
        code, output = duplicate_lists[key]
        if code == 0 or expected_message not in output:
            print(f"expected duplicate {key} inventory failure")
            print(output)
            return 1
    expected_mismatch = "opencpmd_inscan_sections != cpmd_sections"
    if (
        mismatched_inscan_code == 0
        or expected_mismatch not in mismatched_inscan_output
    ):
        print("expected opencpmd_inscan_sections mismatch failure")
        print(mismatched_inscan_output)
        return 1
    expected_extra_section = (
        f"inventory cpmd_section feature outside cpmd_sections: {EXTRA_SECTION_FEATURE}"
    )
    if (
        extra_section_code == 0
        or expected_extra_section not in extra_section_output
    ):
        print("expected extra cpmd section feature failure")
        print(extra_section_output)
        return 1
    expected_allow_duplicate = "inventory opencpmd_sections.txt duplicated: ATOM"
    if (
        duplicate_allow_code == 0
        or expected_allow_duplicate not in duplicate_allow_output
    ):
        print("expected duplicate opencpmd_sections.txt failure")
        print(duplicate_allow_output)
        return 1
    expected_c_duplicate = f"C feature table duplicated: {DUPLICATE_C_FEATURE}"
    if (
        duplicate_c_row_code == 0
        or expected_c_duplicate not in duplicate_c_row_output
    ):
        print("expected duplicate C feature table failure")
        print(duplicate_c_row_output)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
