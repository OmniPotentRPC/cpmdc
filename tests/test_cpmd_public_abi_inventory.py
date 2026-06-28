#!/usr/bin/env python3
"""Check that public cpmdc_* declarations require ABI inventory rows."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
import sys
import tempfile
from pathlib import Path


ADDED_SYMBOL = "cpmdc_uninventoried_demo"
UNIMPLEMENTED_SYMBOL = "cpmdc_unimplemented_demo"


def load_checker(repo: Path):
    path = repo / "tools/check_feature_inventory.py"
    spec = importlib.util.spec_from_file_location("check_feature_inventory", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_header_with_extra_symbol(repo: Path, tmpdir: Path) -> Path:
    header = (repo / "include/cpmdc.h").read_text(encoding="utf-8")
    path = tmpdir / "cpmdc.h"
    path.write_text(
        f"{header}\nint {ADDED_SYMBOL}(void);\n",
        encoding="utf-8",
    )
    return path


def write_inventory_with_extra_symbol(repo: Path, tmpdir: Path, symbol: str) -> Path:
    inventory = json.loads(
        (repo / "schema/inventory/cpmd_features.json").read_text(encoding="utf-8")
    )
    inventory["abi_symbols"].append(symbol)
    inventory["features"].append(
        {
            "feature_id": f"abi.{symbol}",
            "kind": "abi",
            "stub_applicable": True,
            "embed_applicable": True,
            "name": symbol,
            "source": "public C ABI",
        }
    )
    path = tmpdir / "cpmd_features.json"
    path.write_text(json.dumps(inventory), encoding="utf-8")
    return path


def write_features_c_with_extra_symbol(repo: Path, tmpdir: Path, symbol: str) -> Path:
    text = (repo / "src/cpmdc_features.c").read_text(encoding="utf-8")
    row = f'    {{"abi.{symbol}", CPMDC_FEATURE_ABI, 1, 1}},\n'
    path = tmpdir / "cpmdc_features.c"
    path.write_text(
        text.replace(
            "\n};\n\n#undef CPMDC_PARAM_FEATURE",
            f"\n{row}}};\n\n#undef CPMDC_PARAM_FEATURE",
        ),
        encoding="utf-8",
    )
    return path


def write_readme_with_extra_symbol(repo: Path, tmpdir: Path, symbol: str) -> Path:
    path = tmpdir / "README.md"
    path.write_text(
        (repo / "README.md").read_text(encoding="utf-8") + f"\n`{symbol}`\n",
        encoding="utf-8",
    )
    return path


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
    with tempfile.TemporaryDirectory(prefix="cpmd-public-abi-inventory-") as raw:
        tmpdir = Path(raw)
        header = write_header_with_extra_symbol(repo, tmpdir)
        code, output = run_checker_with_paths(checker, HEADER=header)

    expected = f"inventory abi_symbols missing public header function: {ADDED_SYMBOL}"
    if code == 0 or expected not in output:
        print("expected public ABI inventory failure")
        print(output)
        return 1

    with tempfile.TemporaryDirectory(prefix="cpmd-public-abi-implementation-") as raw:
        tmpdir = Path(raw)
        header = write_header_with_extra_symbol(repo, tmpdir)
        header.write_text(
            header.read_text(encoding="utf-8").replace(
                ADDED_SYMBOL, UNIMPLEMENTED_SYMBOL
            ),
            encoding="utf-8",
        )
        code, output = run_checker_with_paths(
            checker,
            HEADER=header,
            INVENTORY=write_inventory_with_extra_symbol(
                repo, tmpdir, UNIMPLEMENTED_SYMBOL
            ),
            FEATURES_C=write_features_c_with_extra_symbol(
                repo, tmpdir, UNIMPLEMENTED_SYMBOL
            ),
            README=write_readme_with_extra_symbol(
                repo, tmpdir, UNIMPLEMENTED_SYMBOL
            ),
        )

    expected = (
        f"src/cpmdc.c missing public ABI implementation: {UNIMPLEMENTED_SYMBOL}"
    )
    if code == 0 or expected not in output:
        print("expected public ABI implementation failure")
        print(output)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
