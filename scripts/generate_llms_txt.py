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


# Org emphasis borders, from org-element's default emphasis components:
# pre is "-", "/", whitespace, or ('"{ and post adds "/" and the closing
# punctuation. A slash can separate two code spans (~POSDATA~/~GETFORCE~).
_LITERAL_BEGIN = re.compile(r"(?i)^#\+begin_(src|example|verse|export)\b")
_LITERAL_END = re.compile(r"(?i)^#\+end_(src|example|verse|export)\b")
# Dollar display math is included here. A line-wise scan leaves the closer
# on a later line, and a code marker inside the open block becomes a span.
_SHIELDED_MATH = re.compile(
    r"\$\$.*?\$\$|\\\(.*?\\\)|\\\[.*?\\\]|\\begin\{([A-Za-z*]+)\}.*?\\end\{\1\}",
    re.DOTALL,
)
_LINK_WITH_DESC = re.compile(r"(?<!\[)\[\[([^][]+)\]\[([^]]*)\]\](?!\])")
_LINK_BARE = re.compile(r"(?<!\[)\[\[([^][]+)\]\](?!\])")
_MATH_TOKEN = re.compile(r"\x00M(\d+)\x00")


def _is_emphasis_pre(char: str) -> bool:
    return char.isspace() or char in "-/('\"{"


def _is_emphasis_post(char: str) -> bool:
    return char.isspace() or char in "-/.,:;!?'\")}[]\\"


def _emphasis_end(line: str, start: int, marker: str) -> int | None:
    """Closing index of an Org code or verbatim span, or None.

    Only the next marker can close the span. A later marker would pull the
    first marker into the body, which is how ``=1`` inside math gets spliced
    into a following verbatim span.
    """
    if start >= len(line) or line[start] != marker:
        return None
    if start > 0 and not _is_emphasis_pre(line[start - 1]):
        return None
    end = line.find(marker, start + 1)
    if end < 0:
        return None
    body = line[start + 1 : end]
    if not body or body[0].isspace() or body[-1].isspace() or "\x00" in body:
        return None
    if end != len(line) - 1 and not _is_emphasis_post(line[end + 1]):
        return None
    return end


def _backslash_math_end(line: str, start: int) -> int | None:
    if line.startswith(r"\(", start):
        end = line.find(r"\)", start + 2)
        return None if end < 0 else end + 2
    if line.startswith(r"\[", start):
        end = line.find(r"\]", start + 2)
        return None if end < 0 else end + 2
    if line.startswith(r"\begin{", start):
        match = re.match(r"\\begin\{([A-Za-z*]+)\}", line[start:])
        if match is None:
            return None
        closer = "\\end{" + match.group(1) + "}"
        end = line.find(closer, start + match.end())
        return None if end < 0 else end + len(closer)
    return None


def _dollar_math_end(line: str, start: int) -> int | None:
    if start >= len(line) or line[start] != "$":
        return None
    if start > 0 and line[start - 1] == "$":
        return None
    if start + 1 < len(line) and line[start + 1] == "$":
        end = line.find("$$", start + 2)
        return None if end < 0 else end + 2
    end = line.find("$", start + 1)
    if end < 0:
        return None
    return end + 1


def _shield_math(text: str, vault: list[str]) -> str:
    def stash(match: re.Match[str]) -> str:
        vault.append(match.group(0))
        return f"\x00M{len(vault) - 1}\x00"

    return _SHIELDED_MATH.sub(stash, text)


def _restore_math(text: str, vault: list[str]) -> str:
    def restore(match: re.Match[str]) -> str:
        return vault[int(match.group(1))]

    return _MATH_TOKEN.sub(restore, text)


def _skip_token(line: str, start: int) -> int | None:
    if line[start] != "\x00":
        return None
    end = line.find("\x00", start + 1)
    if end < 0:
        return None
    return end + 1


def _link_match(line: str, start: int) -> re.Match[str] | None:
    if start > 0 and line[start - 1] == "[":
        return None
    described = _LINK_WITH_DESC.match(line, start)
    if described is not None:
        return described
    return _LINK_BARE.match(line, start)


def _convert_inline(line: str) -> str:
    """Convert inline markup, keeping link targets and code interiors literal.

    Code and verbatim are recognized before links, so a code span that shows
    an Org link keeps the brackets. Math is copied before either of those.
    """
    out: list[str] = []
    index = 0
    while index < len(line):
        token_end = _skip_token(line, index)
        if token_end is not None:
            out.append(line[index:token_end])
            index = token_end
            continue
        math_end = _backslash_math_end(line, index)
        if math_end is None:
            math_end = _dollar_math_end(line, index)
        if math_end is not None:
            out.append(line[index:math_end])
            index = math_end
            continue
        if line[index] in "~=":
            end = _emphasis_end(line, index, line[index])
            if end is not None:
                out.append("`" + line[index + 1 : end] + "`")
                index = end + 1
                continue
        link = _link_match(line, index)
        if link is not None:
            if link.lastindex == 2:
                out.append(f"{_convert_inline(link.group(2))} ({link.group(1)})")
            else:
                out.append(link.group(1))
            index = link.end()
            continue
        out.append(line[index])
        index += 1
    return "".join(out)


def org_plaintext(text: str) -> str:
    """Rough Org to plain text for llms-full.txt (not a full exporter).

    Backslash math, display math, dollar math, and LaTeX environments are
    copied unchanged, so ``=1`` and ``~x~`` inside them are not inline markup.
    Org code and verbatim become backtick spans. Source, example, verse, and
    export blocks are copied as written, including trailing spaces.
    """
    text = normalize_wrapped_org_links(text)
    math_vault: list[str] = []
    text = _shield_math(text, math_vault)
    lines: list[str] = []
    in_literal = False
    for raw in text.splitlines():
        if in_literal:
            if _LITERAL_END.match(raw.rstrip()):
                in_literal = False
                continue
            lines.append(_restore_math(raw, math_vault))
            continue
        line = raw.rstrip()
        if _LITERAL_BEGIN.match(line):
            in_literal = True
            continue
        if line.startswith("#+"):
            continue
        heading_match = HEADING_LINE.match(line)
        if heading_match:
            level = len(heading_match.group(1))
            heading = (heading_match.group(2) or "").strip()
            heading = _restore_math(_convert_inline(heading), math_vault)
            lines.append("#" * min(level, 6) + " " + heading)
            continue
        lines.append(_restore_math(_convert_inline(line), math_vault))
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
