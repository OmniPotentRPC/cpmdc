#!/usr/bin/env python3
"""Cross-check inventory JSON, C table, schema, headers, OpenCPMD allowlists."""
from __future__ import annotations
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schema" / "Potentials.capnp"
INVENTORY = ROOT / "schema" / "inventory" / "cpmd_features.json"
HEADER = ROOT / "include" / "cpmdc.h"
FEATURES_H = ROOT / "include" / "cpmdc_features.h"
FEATURES_C = ROOT / "src" / "cpmdc_features.c"
SEC_ALLOW = ROOT / "schema" / "inventory" / "opencpmd_sections.txt"
CP_ALLOW = ROOT / "schema" / "inventory" / "cpmd_cp_keywords.txt"
DFT_ALLOW = ROOT / "schema" / "inventory" / "dft_functionals.txt"

def main() -> int:
    inv = json.loads(INVENTORY.read_text(encoding="utf-8"))
    schema = SCHEMA.read_text(encoding="utf-8")
    header = HEADER.read_text(encoding="utf-8")
    features_c = FEATURES_C.read_text(encoding="utf-8")
    features_h = FEATURES_H.read_text(encoding="utf-8")
    errors: list[str] = []

    kinds = set(re.findall(r"(\w+)\s+@\d+;", re.search(r"enum CPMDSectionKind \{(.*?)\}", schema, re.S).group(1)))
    inv_kinds = {s["kind"] for s in inv["section_kinds"]}
    if kinds != inv_kinds:
        errors.append(f"section_kinds mismatch schema={sorted(kinds)} inv={sorted(inv_kinds)}")

    fields = set(re.findall(r"(\w+)\s+@\d+\s*:", re.search(r"struct CPMDParams \{(.*?)\}", schema, re.S).group(1)))
    inv_fields = {f["name"] for f in inv["params_fields"]}
    if fields != inv_fields:
        errors.append(f"params_fields mismatch schema={sorted(fields)} inv={sorted(inv_fields)}")

    fids = {f["feature_id"] for f in inv["features"]}
    for sec in inv.get("cpmd_sections", []):
        fid = f"catalog.section.{sec}"
        if fid not in fids:
            errors.append(f"missing catalog section {fid}")
        if f'"{fid}"' not in features_c:
            errors.append(f"C table missing {fid}")

    for kw in inv.get("cpmd_keywords", []):
        if not any(f.get("name") == kw for f in inv["features"] if f["kind"] == "cpmd_keyword"):
            errors.append(f"missing cpmd keyword {kw}")

    for fn in inv.get("dft_functionals", []):
        if not any(f.get("name") == fn for f in inv["features"] if f["kind"] == "dft_functional"):
            errors.append(f"missing dft functional {fn}")

    for sym in inv.get("abi_symbols", []):
        if sym.startswith("cpmdc_feature"):
            if "cpmdc_feature" not in features_h:
                errors.append(f"features header missing {sym}")
        elif sym not in header and sym not in features_h:
            errors.append(f"header missing {sym}")
        if f'"abi.{sym}"' not in features_c:
            errors.append(f"C table missing abi.{sym}")

    if SEC_ALLOW.exists():
        allow = {l.strip() for l in SEC_ALLOW.read_text().splitlines() if l.strip()}
        if allow != set(inv["cpmd_sections"]):
            errors.append("opencpmd_sections.txt diverges from inventory cpmd_sections")
    if CP_ALLOW.exists():
        allow = [l.strip() for l in CP_ALLOW.read_text().splitlines() if l.strip()]
        if allow != inv["cpmd_keywords"]:
            errors.append("cpmd_cp_keywords.txt diverges from inventory")
    if DFT_ALLOW.exists():
        allow = [l.strip() for l in DFT_ALLOW.read_text().splitlines() if l.strip()]
        if allow != inv["dft_functionals"]:
            errors.append("dft_functionals.txt diverges from inventory")

    if "inputBlocks" not in schema:
        errors.append("schema missing inputBlocks")
    if "raw" not in schema:
        errors.append("schema missing raw")

    if errors:
        print("FAIL")
        for e in errors[:50]:
            print(e)
        return 1
    print(
        f"OK inventory: {len(inv['features'])} features, "
        f"{len(inv['cpmd_sections'])} sections, "
        f"{len(inv['cpmd_keywords'])} CP keywords, "
        f"{len(inv['dft_functionals'])} DFT functionals"
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
