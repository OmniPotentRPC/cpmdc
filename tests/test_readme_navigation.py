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
            root / "docs/orgmode/howto/quickstart.org",
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
        reject_text(root / "docs/orgmode/howto/quickstart.org", ["green Meson"])
    )
    _, found, rest = readme.read_text(encoding="utf-8").partition("\n## Documentation\n")
    if not found:
        failures.append("README.md has no Documentation section")
    documentation = rest.split("\n## ", 1)[0]
    base = root / "docs" / "orgmode"
    for page in sorted(base.rglob("*.org")):
        rel = page.relative_to(root).as_posix()
        # changelog.org wraps CHANGELOG.md, which the README links itself.
        if page.name in ("index.org", "changelog.org") and page.parent == base:
            continue
        # An orphan page keeps a moved page's address and links to the new one.
        if "\n:orphan:\n" in page.read_text(encoding="utf-8"):
            continue
        if f"]({rel})" not in documentation:
            failures.append(f"README.md Documentation does not link {rel}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
