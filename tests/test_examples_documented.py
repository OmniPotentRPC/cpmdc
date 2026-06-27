#!/usr/bin/env python3
"""Check that public docs point to a runnable host-side example."""

from __future__ import annotations

import sys
from pathlib import Path


def require_text(path: Path, needles: list[str]) -> list[str]:
    if not path.is_file():
        return [f"missing {path}"]
    text = path.read_text(encoding="utf-8")
    return [f"{path} does not mention {needle!r}" for needle in needles if needle not in text]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_examples_documented.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    failures: list[str] = []
    failures.extend(
        require_text(
            root / "examples/host_step.c",
            [
                "cpmdc_session_create",
                "cpmdc_potential_result_size_for_force_input",
                "cpmdc_session_calculate_result",
                "energy_h=",
                "potential_result_size_bytes=",
            ],
        )
    )
    failures.extend(
        require_text(
            root / "README.md",
            ["examples/host_step.c", "example-host-step", "potential_result_size_bytes="],
        )
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/tutorials/quickstart.org",
            ["examples/host_step.c", "example-host-step", "potential_result_size_bytes="],
        )
    )
    failures.extend(
        require_text(
            root / "docs/orgmode/howto/embedding.org",
            ["examples/host_step.c", "potential_result_size_bytes="],
        )
    )

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
