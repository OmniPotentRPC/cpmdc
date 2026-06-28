#!/usr/bin/env python3
"""Check that the stub ABI test covers every public ABI feature row."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    repo = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    inventory = json.loads(
        (repo / "schema/inventory/cpmd_features.json").read_text(encoding="utf-8")
    )
    source = (repo / "tests/test_stub_abi.c").read_text(encoding="utf-8")
    missing = [
        symbol
        for symbol in inventory.get("abi_symbols", [])
        if f'"abi.{symbol}"' not in source
    ]
    if missing:
        print("stub ABI test does not assert ABI feature rows:")
        for symbol in missing:
            print(f"abi.{symbol}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
