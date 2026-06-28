#!/usr/bin/env python3
"""Check that the README gives readers a clear entry path."""

from __future__ import annotations

import sys
from pathlib import Path


def require_text(path: Path, needles: list[str]) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [
        f"{path} does not mention {needle!r}"
        for needle in needles
        if needle not in text
    ]


def reject_text(path: Path, needles: list[str]) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [
        f"{path} still mentions {needle!r}"
        for needle in needles
        if needle in text
    ]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_readme_navigation.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    readme = root / "README.md"
    failures = require_text(
        readme,
        [
            "Start Here",
            "Choose the API path first",
            "Build and test the stub library",
            "Map a CPMD input line",
            "Link an OpenCPMD archive",
            "Regenerate generated docs",
        ],
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/tutorials/quickstart.org",
            [
                "Work Loop",
                "Build once, then iterate with focused tests",
                "Reconfigure only when build definitions change",
                "Run the full suite before publishing",
            ],
        )
    )
    failures.extend(
        require_text(root / "docs/orgmode/index.org", ["Build and Test"])
    )
    failures.extend(
        reject_text(root / "docs/orgmode/index.org", ["First Command"])
    )
    failures.extend(
        reject_text(root / "docs/orgmode/tutorials/quickstart.org", ["green Meson"])
    )
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
