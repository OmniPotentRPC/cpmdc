#!/usr/bin/env python3
"""Generate llms.txt and llms-full.txt from docs/orgmode and the public inventory.

Follows https://llmstxt.org/ for the summary file. The full file concatenates plain
text derived from every non-orphan Org page under docs/orgmode/.
"""

from __future__ import annotations

import argparse
import difflib
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
REGEN_COMMAND = "python3 scripts/generate_llms_txt.py"
DEFAULT_DIFF_HEAD_LINES = 20

TOCTREE_CAPTION = re.compile(r":caption:\s*(.+)")
TOCTREE_ENTRY = re.compile(r"^\s+(\S+)\s*$")
TITLE_LINE = re.compile(r"^#\+TITLE:\s*(.+)", re.I)
HEADING_LINE = re.compile(r"^(\*+)\s+(.*\S.*)?$")
WRAPPED_ORG_LINK = re.compile(
    r"(?<!\[)\[\[([^\[\]\n]+)\]\[((?:(?!\n[ \t]*\n)[^\[\]])+)\]\](?!\])"
)


def normalize_wrapped_org_links(text: str) -> str:
    """Join well-formed ``[[target][text]]`` links that wrap across lines.

    Only the description may wrap, and never across a blank line.
    Anything opening with ``[[[`` or closing with ``]]]``, such as the Cog
    markers ``[[[cog``, ``]]]`` and ``[[[end]]]``, is left as written.
    """

    def join(match: re.Match[str]) -> str:
        description = re.sub(r"\s*\n\s*", " ", match.group(2))
        return f"[[{match.group(1)}][{description}]]"

    return WRAPPED_ORG_LINK.sub(join, text)


def repo_rel_doc_url(org_path: Path, org_root: Path = ORG_ROOT) -> str:
    rel = org_path.relative_to(org_root).with_suffix("").as_posix()
    return f"{SITE}/{rel}"


def page_title(org_path: Path) -> str:
    for line in org_path.read_text(encoding="utf-8").splitlines():
        match = TITLE_LINE.match(line.strip())
        if match:
            return match.group(1).strip()
    return org_path.stem.replace("-", " ").title()


def org_plaintext(text: str) -> str:
    """Rough Org to plain text for llms-full.txt (not a full exporter)."""
    text = normalize_wrapped_org_links(text)
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
        heading_match = HEADING_LINE.match(line)
        if heading_match:
            level = len(heading_match.group(1))
            heading = (heading_match.group(2) or "").strip()
            lines.append("#" * min(level, 6) + " " + heading)
            continue
        line = re.sub(r"~\([^)]+\)~", "", line)
        line = re.sub(r"~([^~]+)~", r"`\1`", line)
        line = re.sub(r"(?<!\[)\[\[([^][]+)\]\[([^]]*)\]\](?!\])", r"\2 (\1)", line)
        line = re.sub(r"(?<!\[)\[\[([^][]+)\]\](?!\])", r"\1", line)
        lines.append(line)
    return "\n".join(lines).strip() + "\n"


def discover_sections(org_root: Path = ORG_ROOT) -> list[tuple[str, list[Path]]]:
    """Sections and pages from index.org toctree blocks."""
    index = (org_root / "index.org").read_text(encoding="utf-8")
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
                org = org_root / f"{stem}.org"
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


def api_summary_lines(root: Path = ROOT) -> list[str]:
    inventory = root / "schema" / "inventory" / "cpmd_features.json"
    if not inventory.is_file():
        raise FileNotFoundError(
            "missing feature inventory: schema/inventory/cpmd_features.json"
        )
    data = json.loads(inventory.read_text(encoding="utf-8"))
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


def build_llms_txt(root: Path = ROOT) -> str:
    org_root = root / "docs" / "orgmode"
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
    parts.extend(api_summary_lines(root))
    parts.append("")
    parts.append(
        "AGENTS.md at the repository root describes layout, build, test, and docs for agents."
    )
    parts.append(
        f"Full plain text of every page: {SITE}/llms-full.txt (committed as llms-full.txt)."
    )
    parts.append("")
    for caption, pages in discover_sections(org_root):
        parts.append(f"## {caption}")
        parts.append("")
        for org in pages:
            title = page_title(org)
            url = repo_rel_doc_url(org, org_root)
            parts.append(f"- [{title}]({url})")
        parts.append("")
    parts.append("## Repository")
    parts.append("")
    parts.append("- [Source](https://github.com/OmniPotentRPC/cpmdc): MIT licensed cpmdc")
    parts.append(f"- [Schema]({SITE}/reference/cpmdparams-schema): Potentials.capnp CPMDParams")
    parts.append("")
    return "\n".join(parts)


def build_llms_full(root: Path = ROOT) -> str:
    org_root = root / "docs" / "orgmode"
    if not org_root.is_dir():
        raise FileNotFoundError("missing docs source directory: docs/orgmode")
    chunks: list[str] = [
        "# cpmdc documentation (full text export)",
        "",
        f"Site: {SITE}",
        "",
    ]
    org_files = sorted(
        org_root.rglob("*.org"),
        key=lambda path: path.relative_to(org_root).as_posix(),
    )
    for path in org_files:
        if path.name == "changelog.org":
            continue
        text = path.read_text(encoding="utf-8")
        if "\n:orphan:\n" in text:
            continue
        rel = path.relative_to(org_root).as_posix()
        chunks.append(f"## File: docs/orgmode/{rel}")
        chunks.append("")
        chunks.append(org_plaintext(text))
        chunks.append("")
    return "\n".join(chunks).rstrip("\n") + "\n"


def write_outputs(root: Path = ROOT) -> tuple[str, str]:
    summary = build_llms_txt(root)
    full = build_llms_full(root)
    with (root / "llms.txt").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(summary)
    with (root / "llms-full.txt").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(full)
    return summary, full


def format_unified_diff_head(
    rel_path: str,
    actual: str,
    expected: str,
    max_lines: int = DEFAULT_DIFF_HEAD_LINES,
) -> str:
    diff = list(
        difflib.unified_diff(
            actual.splitlines(keepends=True),
            expected.splitlines(keepends=True),
            fromfile=f"a/{rel_path}",
            tofile=f"b/{rel_path}",
        )
    )
    if not diff:
        return "(content differs in line endings or trailing whitespace)"
    head = diff[:max_lines]
    if len(diff) > max_lines:
        head.append(f"... ({len(diff) - max_lines} more diff lines)\n")
    return "".join(head).rstrip("\n")


def check_outputs(
    root: Path = ROOT, max_diff_lines: int = DEFAULT_DIFF_HEAD_LINES
) -> list[str]:
    failures: list[str] = []
    expected_summary = build_llms_txt(root)
    expected_full = build_llms_full(root)
    for rel, path, expected in (
        ("llms.txt", root / "llms.txt", expected_summary),
        ("llms-full.txt", root / "llms-full.txt", expected_full),
    ):
        if not path.is_file():
            failures.append(f"{rel}: missing; run `{REGEN_COMMAND}` to regenerate")
            continue
        actual = path.read_bytes().decode("utf-8")
        if actual != expected:
            diff_head = format_unified_diff_head(
                rel, actual, expected, max_lines=max_diff_lines
            )
            failures.append(
                f"{rel}: stale; run `{REGEN_COMMAND}` to regenerate\n{diff_head}"
            )
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when llms.txt or llms-full.txt are not up to date",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="repository root containing docs/orgmode and schema/inventory",
    )
    parser.add_argument(
        "--diff-lines",
        type=int,
        default=DEFAULT_DIFF_HEAD_LINES,
        help="maximum unified diff lines to print per stale file",
    )
    args = parser.parse_args()
    if args.diff_lines < 1:
        parser.error("--diff-lines must be positive")
    root = args.root.resolve()
    if args.check:
        failures = check_outputs(root, max_diff_lines=args.diff_lines)
        if failures:
            print("\n".join(failures), file=sys.stderr)
            return 1
        print("llms.txt and llms-full.txt are current")
        return 0
    write_outputs(root)
    print("wrote llms.txt")
    print("wrote llms-full.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
