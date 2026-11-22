#!/usr/bin/env python3
"""Verify LeakSanitizer reports allocations leaked inside a dlopen constructor."""

from __future__ import annotations

import os
import subprocess
import sys


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "usage: test_dlopen_leak.py DRIVER MODULE",
            file=sys.stderr,
        )
        return 2
    driver, module = sys.argv[1], sys.argv[2]
    asan_opts = os.environ.get("ASAN_OPTIONS", "")
    if "detect_leaks=0" in asan_opts.split(":"):
        print("skipped: leak detection disabled via ASAN_OPTIONS")
        return 77
    proc = subprocess.run(
        [driver, module],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    output = (proc.stderr or "") + (proc.stdout or "")
    if proc.returncode == 0:
        print(
            "expected LeakSanitizer failure for 77-byte dlopen constructor leak, "
            "but driver exited 0\n" + output,
            file=sys.stderr,
        )
        return 1
    if "LeakSanitizer" not in output or "77 byte(s)" not in output:
        print(
            "expected 77-byte LeakSanitizer report, got:\n" + output,
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
