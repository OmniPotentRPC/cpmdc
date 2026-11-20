#!/usr/bin/env python3
"""Keep orbital storage bookkeeping non-failing."""

from __future__ import annotations

import re
import sys
from pathlib import Path


def main() -> int:
    repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    frontend = (repo / "src/cpmdc.c").read_text(encoding="utf-8")
    ledger = (repo / "src/cpmdc_orbital_limit.c").read_text(encoding="utf-8")
    if "void cpmdc_orbital_note_stored(void)" not in frontend:
        raise AssertionError("the frontend declaration must be non-failing")
    if "void cpmdc_orbital_note_stored(void)" not in ledger:
        raise AssertionError("the ledger implementation must be non-failing")
    if re.search(r"if\s*\(\s*cpmdc_orbital_note_stored\s*\(", frontend):
        raise AssertionError("the storage note must not have an error branch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
