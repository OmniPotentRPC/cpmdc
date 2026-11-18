#!/usr/bin/env python3
"""Check that .lychee.toml exclude patterns match doi.org URLs."""

from __future__ import annotations

import re
import sys
from pathlib import Path

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_lychee_toml.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    config_path = root / ".lychee.toml"
    data = tomllib.loads(config_path.read_text(encoding="utf-8"))
    patterns = data.get("exclude", [])
    if not patterns:
        print(f"{config_path}: exclude list is empty", file=sys.stderr)
        return 1

    compiled = [re.compile(pattern) for pattern in patterns]
    doi_url = "https://doi.org/10.1016/0022-3697(94)00228-2"
    github_url = "https://github.com/OmniPotentRPC/cpmdc"

    if not any(regex.search(doi_url) for regex in compiled):
        print(
            f"{config_path}: no exclude pattern matches a doi.org URL ({doi_url!r})",
            file=sys.stderr,
        )
        for pattern in patterns:
            print(f"  pattern {pattern!r}", file=sys.stderr)
        return 1

    if any(regex.search(github_url) for regex in compiled):
        print(
            f"{config_path}: exclude pattern wrongly matches {github_url!r}",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())