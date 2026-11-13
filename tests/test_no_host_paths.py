#!/usr/bin/env python3
"""Tracked files must not name host-local directories."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def needles() -> list[bytes]:
    slash = bytes([47])
    return [slash + name for name in (b"tmp/", b"home/", b"opt/")]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_no_host_paths.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    listed = subprocess.check_output(["git", "-C", str(root), "ls-files", "-z"])
    bad: list[str] = []
    wanted = needles()
    for raw in listed.split(b"\0"):
        if not raw:
            continue
        path = root / raw.decode("utf-8")
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data[:4096]:
            continue
        for needle in wanted:
            if needle in data:
                bad.append(f"{path.relative_to(root)} contains {needle.decode()}")
                break
    if bad:
        print("\n".join(bad), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
