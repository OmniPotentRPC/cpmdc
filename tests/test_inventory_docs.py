#!/usr/bin/env python3
"""Check that inventory coverage docs name the runnable guards."""

from __future__ import annotations

import sys
from pathlib import Path


def require_text(path: Path, needles: list[str]) -> list[str]:
    if not path.is_file():
        return [f"missing {path}"]
    text = path.read_text(encoding="utf-8")
    return [
        f"{path} does not mention {needle!r}"
        for needle in needles
        if needle not in text
    ]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_inventory_docs.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    inventory_guards = [
        "feature-inventory",
        "cpmd-base-keyword-inventory",
        "cpmd-params-field-inventory",
        "cpmd-public-abi-inventory",
        "stub-abi-symbol-coverage",
        "shared-dlopen-symbol-coverage",
        "cpmd-schema-render-coverage",
        "cpmd-option-token-coverage",
        "cpmd-typed-render-field-coverage",
        "cpmd-long-tail-section-coverage",
        "cpmd-escape-hatch-coverage",
    ]
    docs_guards = [
        "examples-documented",
        "readme-navigation",
    ]
    failures: list[str] = []
    failures.extend(require_text(root / "README.md", inventory_guards))
    failures.extend(require_text(root / "README.md", docs_guards))
    failures.extend(
        require_text(
            root / "README.md",
            [
                "src/cpmdc.c",
                "src/cpmdc_stub.c",
                "src/cpmdc_features.c",
            ],
        )
    )
    failures.extend(
        require_text(root / "docs/orgmode/contributing/index.org", inventory_guards)
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/contributing/index.org",
            ["C feature table duplicates"],
        )
    )
    failures.extend(
        require_text(root / "docs/orgmode/contributing/index.org", docs_guards)
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/contributing/index.org",
            [
                "src/cpmdc.c",
                "src/cpmdc_stub.c",
                "src/cpmdc_features.c",
            ],
        )
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/reference/cpmd-options.org",
            [
                "schema/inventory/cpmd_cp_keywords.txt",
                "cpmd-base-keyword-inventory",
                "cpmd-params-field-inventory",
                "cpmd-public-abi-inventory",
                "stub-abi-symbol-coverage",
                "shared-dlopen-symbol-coverage",
                "cpmd-schema-render-coverage",
                "cpmd-option-token-coverage",
                "cpmd-typed-render-field-coverage",
                "cpmd-long-tail-section-coverage",
                "cpmd-escape-hatch-coverage",
                "src/cpmdc.c",
                "src/cpmdc_stub.c",
                "src/cpmdc_features.c",
                "opencpmd_inscan_sections",
                "cpmd_sections",
                "C feature table duplicates",
            ],
        )
    )

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
