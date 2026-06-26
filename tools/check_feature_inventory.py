#!/usr/bin/env python3
"""Cross-check cpmdc feature inventory vs schema and public ABI header."""
from __future__ import annotations
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schema" / "Potentials.capnp"
INVENTORY = ROOT / "schema" / "inventory" / "cpmd_features.json"
HEADER = ROOT / "include" / "cpmdc.h"

def main() -> int:
    schema = SCHEMA.read_text(encoding="utf-8")
    inv = json.loads(INVENTORY.read_text(encoding="utf-8"))
    header = HEADER.read_text(encoding="utf-8")
    errors: list[str] = []

    kinds = set(re.findall(r"enum CPMDSectionKind \{(.*?)\}", schema, re.S)[0]
                and re.findall(r"(\w+)\s+@\d+;", re.search(r"enum CPMDSectionKind \{(.*?)\}", schema, re.S).group(1)))
    inv_kinds = {s["kind"] for s in inv["section_kinds"]}
    if kinds - inv_kinds:
        errors.append(f"inventory missing section kinds: {sorted(kinds - inv_kinds)}")
    if inv_kinds - kinds:
        errors.append(f"inventory extra section kinds: {sorted(inv_kinds - kinds)}")

    fields = set(re.findall(r"(\w+)\s+@\d+\s*:", re.search(r"struct CPMDParams \{(.*?)\}", schema, re.S).group(1)))
    inv_fields = {f["name"] for f in inv["params_fields"]}
    if fields - inv_fields:
        errors.append(f"inventory missing params fields: {sorted(fields - inv_fields)}")
    if inv_fields - fields:
        errors.append(f"inventory extra params fields: {sorted(inv_fields - fields)}")

    for sym in inv["abi_symbols"]:
        if sym not in header:
            errors.append(f"header missing ABI symbol mention: {sym}")

    if errors:
        print("FAIL")
        for e in errors:
            print(e)
        return 1
    print("OK inventory consistent with schema and cpmdc.h")
    return 0

if __name__ == "__main__":
    sys.exit(main())
