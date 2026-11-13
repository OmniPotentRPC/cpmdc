#!/usr/bin/env python3
"""The reference-well figure is drawn from the evaluator formula."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_reference_well.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    script = root / "docs" / "render_reference_well.py"
    with tempfile.TemporaryDirectory(prefix="cpmdc-well-") as tmp:
        out = Path(tmp) / "reference-well.svg"
        proc = subprocess.run(
            [sys.executable, str(script), "--out", str(out), "--check"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr)
            return proc.returncode
        text = out.read_text(encoding="utf-8")
    if "placeholder" in text.lower() or "<polyline" not in text:
        print("reference well is not a drawn curve", file=sys.stderr)
        return 1
    if text.count(",") < 16:
        print("reference well has too few points", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
