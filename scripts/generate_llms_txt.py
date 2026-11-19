#!/usr/bin/env python3
"""Generate llms.txt and llms-full.txt from docs/orgmode and the public inventory.

Follows https://llmstxt.org/ for the summary file. The full file concatenates plain
text derived from every Org page under docs/orgmode/.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORG_ROOT = ROOT / "docs" / "orgmode"
INVENTORY = ROOT / "schema" / "inventory" / "cpmd_features.json"
SITE = "https://cpmdc.rgoswami.me"
SUMMARY_PATH = ROOT / "llms.txt"
FULL_PATH = ROOT / "llms-full.txt"

TOCTREE_CAPTION = re.compile(r":caption:\s*(.+)")
TOCTREE_ENTRY = re.compile(r"^\s+(\S+)\s*$")
TITLE_LINE = re.compile(r"^#\+TITLE:\s*(.+)", re.I)


def repo_rel_doc_url(org_path: Path) -> str:
    rel = org_path.relative_to(ORG_ROOT).with_suffix("").as_posix()
    return f"{SITE}/{rel}"


def page_title(org_path: Path) -> str:
    for line in org_path.read_text(encoding="utf-8").splitlines():
        match = TITLE_LINE.match(line.strip())
        if match:
            return match.group(1).strip()
    return org_path.stem.replace("-", " ").title()


def org_plaintext(text: str) -> str:
    """Rough Org to plain text for llms-full.txt (not a full exporter)."""
    lines: list[str] = []
    in_export = False
    for raw in text.splitlines():
        line = raw.rstrip()
        if line.startswith("#+begin_export") or line.startswith("#+BEGIN_EXPORT"):
            in_export = True
            continue
        if line.startswith("#+end_export") or line.startswith("#+END_EXPORT"):
            in_export = False
            continue
        if in_export:
            lines.append(line)
            continue
        if line.startswith("#+"):
            continue
        if line.startswith("*"):
            level = len(line) - len(line.lstrip("*"))
            heading = line.lstrip("*").strip()
            lines.append("#" * min(level, 6) + " " + heading)
            continue
        line = re.sub(r"~\([^)]+\)~", "", line)
        line = re.sub(r"~([^~]+)~", r"`\1`", line)
        line = re.sub(r"\[\[([^]]+)\]\[([^]]*)\]\]", r"\2 (\1)", line)
        line = re.sub(r"\[\[([^]]+)\]\]", r"\1", line)
        lines.append(line)
    return "\n".join(lines).strip() + "\n"


def discover_sections() -> list[tuple[str, list[Path]]]:
    """Sections and pages from index.org toctree blocks."""
    index = (ORG_ROOT / "index.org").read_text(encoding="utf-8")
    sections: list[tuple[str, list[Path]]] = []
    current_caption = "Documentation"
    current_pages: list[Path] = []
    in_toctree = False
    for line in index.splitlines():
        if line.strip().startswith(".. toctree::"):
            in_toctree = True
            current_pages = []
            continue
        if in_toctree:
            cap = TOCTREE_CAPTION.search(line)
            if cap:
                current_caption = cap.group(1).strip()
                continue
            entry = TOCTREE_ENTRY.match(line)
            if entry:
                stem = entry.group(1)
                org = ORG_ROOT / f"{stem}.org"
                if not org.is_file() and "/" in stem:
                    org = ORG_ROOT / f"{stem}.org"
                if org.is_file():
                    current_pages.append(org)
                continue
            if line.strip() == "" and current_pages:
                sections.append((current_caption, list(current_pages)))
                current_pages = []
                in_toctree = False
                continue
            if not line.startswith(" ") and line.strip() and current_pages:
                sections.append((current_caption, list(current_pages)))
                current_pages = []
                in_toctree = False
    if current_pages:
        sections.append((current_caption, list(current_pages)))
    return sections


def api_summary_lines() -> list[str]:
    if not INVENTORY.is_file():
        return []
    data = json.loads(INVENTORY.read_text(encoding="utf-8"))
    feats = data.get("features", [])
    abi = sum(1 for f in feats if str(f.get("feature_id", "")).startswith("abi."))
    params = sum(1 for f in feats if str(f.get("feature_id", "")).startswith("params."))
    sec = sum(1 for f in feats if str(f.get("feature_id", "")).startswith("catalog.section."))
    kw = sum(
        1
        for f in feats
        if str(f.get("feature_id", "")).startswith(("catalog.cpmd.", "catalog.dft."))
    )
    return [
        f"The C ABI exposes {len(data.get('abi_symbols', []))} functions; "
        f"the feature table lists {abi} ABI rows, {params} parameter fields, "
        f"{sec} deck sections, and {kw} catalog keywords.",
        "Authoritative API tables: reference/c-abi, reference/contracts, "
        "reference/python, reference/cpmd-options.",
    ]


def build_llms_txt() -> str:
    desc = (
        "C ABI for embedding OpenCPMD in OmniPotentRPC tools. Hosts pass "
        "Cap'n Proto CPMDParams once, ForceInput per geometry step, and receive "
        "native C results or PotentialResult messages for rgpot and other drivers."
    )
    parts = [
        "# cpmdc",
        "",
        f"> {desc}",
        "",
        "Machine-readable inventory: schema/inventory/cpmd_features.json.",
        "Human docs source: docs/orgmode/ (export to Sphinx with pixi run -e docs docbld).",
        "",
    ]
    parts.extend(api_summary_lines())
    parts.append("")
    parts.append(f"AGENTS.md at the repository root describes layout, build, test, and docs for agents.")
    parts.append(f"Full plain text of every page: {SITE}/llms-full.txt (committed as llms-full.txt).")
    parts.append("")
    for caption, pages in discover_sections():
        parts.append(f"## {caption}")
        parts.append("")
        for org in pages:
            title = page_title(org)
            url = repo_rel_doc_url(org)
            parts.append(f"- [{title}]({url})")
        parts.append("")
    parts.append("## Repository")
    parts.append("")
    parts.append(f"- [Source](https://github.com/OmniPotentRPC/cpmdc): MIT licensed cpmdc")
    parts.append(f"- [Schema]({SITE}/reference/cpmdparams-schema): Potentials.capnp CPMDParams")
    parts.append("")
    return "\n".join(parts)


def build_llms_full() -> str:
    chunks: list[str] = [
        "# cpmdc documentation (full text export)",
        "",
        f"Site: {SITE}",
        "",
    ]
    for path in sorted(ORG_ROOT.rglob("*.org")):
        if path.name == "changelog.org":
            continue
        rel = path.relative_to(ORG_ROOT).as_posix()
        chunks.append(f"## File: docs/orgmode/{rel}")
        chunks.append("")
        chunks.append(org_plaintext(path.read_text(encoding="utf-8")))
        chunks.append("")
    return "\n".join(chunks)


def write_outputs() -> tuple[str, str]:
    summary = build_llms_txt()
    full = build_llms_full()
    SUMMARY_PATH.write_text(summary, encoding="utf-8", newline="\n")
    FULL_PATH.write_text(full, encoding="utf-8", newline="\n")
    return summary, full


def check_outputs() -> list[str]:
    failures: list[str] = []
    expected_summary = build_llms_txt()
    expected_full = build_llms_full()
    for path, expected in (
        (SUMMARY_PATH, expected_summary),
        (FULL_PATH, expected_full),
    ):
        if not path.is_file():
            failures.append(f"{path}: missing; run scripts/generate_llms_txt.py")
            continue
        actual = path.read_text(encoding="utf-8")
        if actual != expected:
            failures.append(
                f"{path}: stale (run python3 scripts/generate_llms_txt.py to refresh)"
            )
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when llms.txt or llms-full.txt are not up to date",
    )
    args = parser.parse_args()
    if args.check:
        failures = check_outputs()
        if failures:
            print("\n".join(failures), file=sys.stderr)
            return 1
        print("llms.txt and llms-full.txt are current")
        return 0
    write_outputs()
    print(f"wrote {SUMMARY_PATH}")
    print(f"wrote {FULL_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())