#!/usr/bin/env python3
"""Cross-check inventory vs schema/C table; probe OpenCPMD inscan when tree present."""
from __future__ import annotations
import json, os, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "schema" / "Potentials.capnp"
INVENTORY = ROOT / "schema" / "inventory" / "cpmd_features.json"
HEADER = ROOT / "include" / "cpmdc.h"
FEATURES_H = ROOT / "include" / "cpmdc_features.h"
FEATURES_C = ROOT / "src" / "cpmdc_features.c"
SEC_ALLOW = ROOT / "schema" / "inventory" / "opencpmd_sections.txt"
CPMD_OPTIONS_DOC = ROOT / "docs" / "orgmode" / "reference" / "cpmd-options.org"
README = ROOT / "README.md"

def probe_inscan_sections(cpmd_root: Path) -> set[str]:
    secs: set[str] = set()
    src = cpmd_root / "src"
    if not src.is_dir():
        return secs
    for p in src.rglob("*.F90"):
        t = p.read_text(errors="replace")
        for m in re.finditer(r"inscan\s*\([^)]*'&([A-Z][A-Z0-9_]*)'", t, re.I):
                secs.add(m.group(1).upper())
    return secs

def struct_fields(schema: str, struct_name: str) -> list[str]:
    match = re.search(rf"struct {struct_name} \{{(.*?)\n\}}", schema, re.S)
    if not match:
        raise ValueError(f"schema missing struct {struct_name}")
    return re.findall(r"(\w+)\s+@\d+\s*:", match.group(1))

def structured_param_feature_ids(schema: str) -> list[str]:
    section_structs = {
        "generic": "CPMDGenericSection",
        "atom": "CPMDDirectiveSection",
        "system": "CPMDSystemSection",
        "cpmd": "CPMDCpmdSection",
        "dft": "CPMDDftSection",
        "atoms": "CPMDAtomsSection",
        "basis": "CPMDDirectiveSection",
        "clas": "CPMDDirectiveSection",
        "eam": "CPMDDirectiveSection",
        "exte": "CPMDDirectiveSection",
        "hardness": "CPMDDirectiveSection",
        "info": "CPMDDirectiveSection",
        "linres": "CPMDDirectiveSection",
        "molstates": "CPMDDirectiveSection",
        "mts": "CPMDDirectiveSection",
        "nlcc": "CPMDDirectiveSection",
        "path": "CPMDDirectiveSection",
        "pimd": "CPMDDirectiveSection",
        "potential": "CPMDDirectiveSection",
        "prop": "CPMDDirectiveSection",
        "ptddft": "CPMDDirectiveSection",
        "resp": "CPMDDirectiveSection",
        "tddft": "CPMDDirectiveSection",
        "vdw": "CPMDDirectiveSection",
        "vectors": "CPMDDirectiveSection",
        "wavefunction": "CPMDDirectiveSection",
        "set": "CPMDSetDirective",
    }
    ids: list[str] = []
    for section, struct_name in section_structs.items():
        for field in struct_fields(schema, struct_name):
            ids.append(f"params.inputSections.{section}.{field}")
    ids.append("params.inputSections.raw")
    return ids

def main() -> int:
    inv = json.loads(INVENTORY.read_text(encoding="utf-8"))
    schema = SCHEMA.read_text(encoding="utf-8")
    header = HEADER.read_text(encoding="utf-8")
    features_c = FEATURES_C.read_text(encoding="utf-8")
    features_h = FEATURES_H.read_text(encoding="utf-8")
    cpmd_options_doc = CPMD_OPTIONS_DOC.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")
    errors: list[str] = []

    reference = str(inv.get("opencpmd_reference") or "")
    if reference.startswith(("/", "~")):
        errors.append("opencpmd_reference must describe the reference source, not a host-local path")

    kinds = set(re.findall(r"(\w+)\s+@\d+;", re.search(r"enum CPMDSectionKind \{(.*?)\}", schema, re.S).group(1)))
    inv_kinds = {s["kind"] for s in inv["section_kinds"]}
    if kinds != inv_kinds:
        errors.append(f"section_kinds mismatch schema={sorted(kinds)} inv={sorted(inv_kinds)}")

    fields = set(re.findall(r"(\w+)\s+@\d+\s*:", re.search(r"struct CPMDParams \{(.*?)\}", schema, re.S).group(1)))
    inv_fields = {f["name"] for f in inv["params_fields"]}
    if fields != inv_fields:
        errors.append(f"params_fields mismatch")

    inv_secs = set(inv.get("cpmd_sections", []))
    allow = {l.strip() for l in SEC_ALLOW.read_text().splitlines() if l.strip()} if SEC_ALLOW.exists() else set()
    if allow and allow != inv_secs:
        errors.append(f"opencpmd_sections.txt != inventory cpmd_sections: {sorted(allow ^ inv_secs)}")

    # Live OpenCPMD probe (non-circular completeness)
    cpmd_root = os.environ.get("CPMD_ROOT") or os.environ.get("CPMDC_CPMD_ROOT")
    live_probe_status = "live probe skipped"
    if not cpmd_root:
        # meson often uses /tmp/OpenCPMD-gf14
        for cand in (Path("/tmp/OpenCPMD-gf14"), Path(inv.get("opencpmd_reference") or "")):
            if cand.is_dir() and (cand / "src").is_dir():
                cpmd_root = str(cand)
                break
    if cpmd_root and Path(cpmd_root).joinpath("src").is_dir():
        live = probe_inscan_sections(Path(cpmd_root))
        missing = live - inv_secs
        extra = inv_secs - live
        if missing:
            errors.append(f"inventory missing OpenCPMD inscan sections: {sorted(missing)}")
        if extra:
            errors.append(f"inventory has non-inscan names as sections: {sorted(extra)}")
        if allow and live - allow:
            errors.append(f"allowlist missing live inscan: {sorted(live - allow)}")
        print(f"probed {cpmd_root}: {len(live)} inscan sections")
        live_probe_status = "live probe passed"
    else:
        print("WARN: no OpenCPMD tree to probe (set CPMD_ROOT); skipping live completeness")

    fids = {f["feature_id"] for f in inv["features"]}
    for sec in inv_secs:
        fid = f"catalog.section.{sec}"
        if fid not in fids:
            errors.append(f"missing catalog section {fid}")
        if f'"{fid}"' not in features_c:
            errors.append(f"C table missing {fid}")
        if fid not in cpmd_options_doc:
            errors.append(f"cpmd-options docs missing {fid}")

    for feature in inv["features"]:
        fid = feature["feature_id"]
        if fid.startswith(("catalog.cpmd.", "catalog.dft.")):
            if fid not in cpmd_options_doc:
                errors.append(f"cpmd-options docs missing {fid}")
        if fid.startswith("params."):
            if fid not in cpmd_options_doc:
                errors.append(f"cpmd-options docs missing {fid}")

    c_flags = {
        fid: (bool(int(stub)), bool(int(embed)))
        for fid, stub, embed in re.findall(
            r'\{"([^"]+)",\s*CPMDC_FEATURE_\w+,\s*([01]),\s*([01])\}',
            features_c,
        )
    }
    for fid in re.findall(r'CPMDC_PARAM_FEATURE\(\s*"([^"]+)"\s*\)', features_c):
        c_flags[fid] = (True, True)

    for fid in structured_param_feature_ids(schema):
        if fid not in fids:
            errors.append(f"inventory missing structured params feature {fid}")
        if fid not in cpmd_options_doc:
            errors.append(f"cpmd-options docs missing {fid}")
        if fid not in c_flags:
            errors.append(f"C table missing structured params feature {fid}")
            continue
        stub, embed = c_flags[fid]
        if not stub or not embed:
            errors.append(f"{fid} structured params applicability should be stub+embed")

    for feature in inv["features"]:
        fid = feature["feature_id"]
        if fid not in c_flags:
            errors.append(f"C table missing inventory feature {fid}")
            continue
        stub, embed = c_flags[fid]
        if bool(feature["stub_applicable"]) != stub:
            errors.append(f"{fid} stub_applicable inventory={feature['stub_applicable']} C={stub}")
        if bool(feature["embed_applicable"]) != embed:
            errors.append(f"{fid} embed_applicable inventory={feature['embed_applicable']} C={embed}")

    for section in inv["section_kinds"]:
        fid = section["feature_id"]
        if fid not in c_flags:
            errors.append(f"section kind C table missing {fid}")
            continue
        stub, embed = c_flags[fid]
        if bool(section["stub_applicable"]) != stub:
            errors.append(f"{fid} section stub_applicable inventory={section['stub_applicable']} C={stub}")
        if bool(section["embed_applicable"]) != embed:
            errors.append(f"{fid} section embed_applicable inventory={section['embed_applicable']} C={embed}")

    # Must include skeptic-required
    for req in ("EAM", "MOLSTATES", "NLCC", "VECTORS"):
        if req not in inv_secs:
            errors.append(f"required inscan section absent from inventory: {req}")

    typed_catalog_sections = inv_secs
    for feature in inv["features"]:
        if feature.get("kind") != "cpmd_section":
            continue
        name = feature["name"]
        expected_render = (
            "typed_or_set_or_generic_or_raw"
            if name in typed_catalog_sections
            else "set_or_generic_or_raw"
        )
        if feature.get("render_via") != expected_render:
            errors.append(
                f"{feature['feature_id']} render_via={feature.get('render_via')} "
                f"expected {expected_render}"
            )

    # Must NOT list keyword-only as sections
    for bad in ("QMMM", "PROPERTIES", "BICANONICAL", "CDFT", "CLASSIC", "F_B"):
        if bad in inv_secs:
            errors.append(f"non-deck name listed as section: {bad}")

    for sym in inv.get("abi_symbols", []):
        if sym.startswith("cpmdc_feature"):
            if "cpmdc_feature" not in features_h:
                errors.append(f"features header missing {sym}")
        elif sym not in header and sym not in features_h:
            errors.append(f"header missing {sym}")
        if f'"abi.{sym}"' not in features_c:
            errors.append(f"C table missing abi.{sym}")
        if sym not in readme:
            errors.append(f"README missing ABI symbol {sym}")

    if "inputBlocks" not in schema or "raw" not in schema:
        errors.append("schema missing escape hatches")

    if errors:
        print("FAIL")
        for e in errors[:60]:
            print(e)
        return 1
    print(
        f"OK inventory: {len(inv['features'])} features, "
        f"{len(inv_secs)} inscan sections (incl EAM/MOLSTATES/NLCC/VECTORS), "
        f"{live_probe_status}"
    )
    return 0

if __name__ == "__main__":
    sys.exit(main())
