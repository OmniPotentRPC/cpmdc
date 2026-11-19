#!/usr/bin/env python3
"""Fail when llms.txt or llms-full.txt are stale."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_llms_txt.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    script = root / "scripts" / "generate_llms_txt.py"
    result = subprocess.run(
        [sys.executable, str(script), "--check"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())