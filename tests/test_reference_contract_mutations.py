#!/usr/bin/env python3
"""The reference checker rejects struct fields assigned to the wrong struct."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


MUTATIONS = {
    "missing CPMDCRankSpread.energy_h": (
        "docs/orgmode/reference/c-abi.org",
        "| ~energy_h~ | ~CPMDCRankSpread~ | largest absolute energy difference "
        "from rank 0 over the ranks, Hartree |\n",
        "",
        "CPMDCRankSpread.energy_h",
    ),
    "bogus cpmdc_pwscale_stats.natoms": (
        "docs/orgmode/reference/pwscale.org",
        "| ~natoms~ | ~int~ | atoms | 1 to 1000000 |\n",
        "| ~natoms~ | ~int~ | atoms | 1 to 1000000 |\n"
        "| ~natoms~ | ~cpmdc_pwscale_stats~ | bogus field |\n",
        "cpmdc_pwscale_stats.natoms",
    ),
}


def copied_tree(source: Path, destination: Path) -> None:
    for relative in ("include", "docs", "python", "README.md"):
        source_path = source / relative
        destination_path = destination / relative
        if source_path.is_dir():
            shutil.copytree(source_path, destination_path)
        else:
            destination_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_path, destination_path)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_reference_contract_mutations.py SOURCE_ROOT", file=sys.stderr)
        return 2
    source = Path(sys.argv[1]).resolve()
    failures = []
    for description, (relative, original, replacement, expected) in MUTATIONS.items():
        with tempfile.TemporaryDirectory() as tmp:
            tree = Path(tmp) / "tree"
            copied_tree(source, tree)
            path = tree / relative
            text = path.read_text(encoding="utf-8")
            if original not in text:
                failures.append(f"{description}: mutation target is absent")
                continue
            path.write_text(text.replace(original, replacement, 1), encoding="utf-8")
            checked = subprocess.run(
                [sys.executable, source / "tests/test_reference_contract.py", tree],
                check=False,
                capture_output=True,
                text=True,
            )
            if checked.returncode == 0:
                failures.append(f"{description}: checker accepted the mutation")
            elif expected not in checked.stderr:
                failures.append(
                    f"{description}: checker failed without naming {expected}: "
                    f"{checked.stderr.strip()}"
                )
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"{len(MUTATIONS)} struct-field mutations rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
