#!/usr/bin/env python3
"""python/cpmdc passes ``mypy --strict`` with the configuration of pyproject.toml.

Exit status 77 (a Meson skip) when mypy is not installed.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

SKIP = 77


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_mypy_strict.py SOURCE_ROOT", file=sys.stderr)
        return 2
    if importlib.util.find_spec("mypy") is None:
        print("mypy is not installed", file=sys.stderr)
        return SKIP
    root = Path(sys.argv[1])
    done = subprocess.run(
        [sys.executable, "-m", "mypy", "--strict", "--config-file", str(root / "pyproject.toml")],
        cwd=root,
        check=False,
    )
    return done.returncode


if __name__ == "__main__":
    raise SystemExit(main())
