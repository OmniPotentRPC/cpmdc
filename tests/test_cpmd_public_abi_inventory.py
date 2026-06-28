#!/usr/bin/env python3
"""Check that public cpmdc_* declarations require ABI inventory rows."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import sys
import tempfile
from pathlib import Path


ADDED_SYMBOL = "cpmdc_uninventoried_demo"


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


def run_checker_with_header(checker, header: Path) -> tuple[int, str]:
    old_header = checker.HEADER
    old_cpmd_root = os.environ.pop("CPMD_ROOT", None)
    old_cpmdc_root = os.environ.pop("CPMDC_CPMD_ROOT", None)
    stream = io.StringIO()
    try:
        checker.HEADER = header
        with contextlib.redirect_stdout(stream):
            code = checker.main()
    finally:
        checker.HEADER = old_header
        if old_cpmd_root is not None:
            os.environ["CPMD_ROOT"] = old_cpmd_root
        if old_cpmdc_root is not None:
            os.environ["CPMDC_CPMD_ROOT"] = old_cpmdc_root
    return code, stream.getvalue()


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    checker = load_checker(repo)
    with tempfile.TemporaryDirectory(prefix="cpmd-public-abi-inventory-") as raw:
        header = write_header_with_extra_symbol(repo, Path(raw))
        code, output = run_checker_with_header(checker, header)

    expected = f"inventory abi_symbols missing public header function: {ADDED_SYMBOL}"
    if code == 0 or expected not in output:
        print("expected public ABI inventory failure")
        print(output)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
