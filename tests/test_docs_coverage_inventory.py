#!/usr/bin/env python3
"""Fail when the coverage inventory reports documentation gaps."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_docs_coverage_inventory.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    script = root / "scripts" / "docs_coverage_inventory.py"
    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        print("documentation coverage gaps remain; see JSON above", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())