#!/usr/bin/env python3
"""Validate CPMD_ROOT / CPMDC_CPMD_EXE for embed builds."""
import os
import pathlib
import sys


def main() -> int:
    exe = os.environ.get("CPMDC_CPMD_EXE", "")
    root = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CPMD_ROOT", "")
    candidates = []
    if exe:
        candidates.append(pathlib.Path(exe))
    if root:
        r = pathlib.Path(root)
        candidates += [r / "bin" / "cpmd.x", r / "cpmd.x"]
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            print(c)
            return 0
    print("cpmd.x not found; set CPMDC_CPMD_EXE or pass cpmd_root with bin/cpmd.x",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
