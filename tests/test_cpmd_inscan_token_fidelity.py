#!/usr/bin/env python3
"""Check typed render literals against the reader-scoped OpenCPMD snapshot.

Each ``render_*_section`` line is accepted when that section's readers would
accept it. ``keyword_contains`` matches whole words (``and`` requires every
word, ``alias`` is an alternative). ``INDEX`` matches a prefix of a word.
A line whose every word is a literal that section tests is also accepted.
The check reads ``schema/inventory/opencpmd_reader_tokens.json`` and does
not need an OpenCPMD tree.

``OPENCPMD_READER_TREE`` (or ``--opencpmd-tree``) rebuilds the snapshot from
that tree and fails if it differs from the committed file.

Typed carriers the owning reader would not accept as a line are xfail with
reason ``P03 typed-carrier-fidelity``. No other literal may fail.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from opencpmd_reader_tokens import (  # noqa: E402
    canonical_snapshot,
    extract_tree,
    load_snapshot,
    section_tokens,
    tested_literals,
)

XFAIL_REASON = "P03 typed-carrier-fidelity"
# Phrases the typed renderer emits that the section reader does not test.
XFAIL = {
    "VDW": {"GRIMME", "S6", "RADIUS"},
    "RESP": {"WEIGHT", "HYPERBOLIC", "BETA"},
    "EXTE": {"EFIELD"},
    "VECTORS": {"NEWORTHO", "OVERLAP"},
    "PATH": {"SPRING"},
    "CPMD": {"SHAKE"},
}

# Typed section renderers. Helpers they call for fixed deck text are inlined.
SECTION_RENDERERS = {
    "render_system_section_with_cell": "SYSTEM",
    "render_cpmd_section": "CPMD",
    "render_dft_section": "DFT",
    "render_atoms_section": "ATOMS",
    "render_vdw_section": "VDW",
    "render_prop_section": "PROP",
    "render_linres_section": "LINRES",
    "render_pimd_section": "PIMD",
    "render_path_section": "PATH",
    "render_tddft_section": "TDDFT",
    "render_resp_section": "RESP",
    "render_exte_section": "EXTE",
    "render_vectors_section": "VECTORS",
}

# Copy the caller's text. Their string arguments are the emissions, not a
# second keyword table hidden in the helper.
_PRIMITIVE_APPEND = {
    "append_fmt",
    "append_text",
    "append_capn_text",
    "append_f64_list_line",
    "append_i32_list_line",
    "append_i32_list_inline",
    "append_f64_list_inline",
    "append_slice",
    "append_directives",
    "append_set_payload",
    "append_set_directives_for_section",
    "append_text_keyword_arg",
}
_EMITTERS = ("append_text", "append_fmt", "append_text_keyword_arg")
_PRINTF = re.compile(
    r"%[-+ #0]*\d*(?:\.\d+)?(?:hh|h|ll|l|L|z|j|t)?[diouxXfFeEgGaAcspn]"
)
_FUNC_HEAD = re.compile(
    r"^(?:static\s+)?(?:const\s+)?[\w][\w\s\*]*?\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(",
    re.MULTILINE,
)


def _decode_c_string(raw: str) -> str:
    out: list[str] = []
    i = 0
    escapes = {"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"', "'": "'"}
    while i < len(raw):
        if raw[i] == "\\" and i + 1 < len(raw):
            out.append(escapes.get(raw[i + 1], raw[i + 1]))
            i += 2
            continue
        out.append(raw[i])
        i += 1
    return "".join(out)


def _c_lex(text: str):
    """Yield ``('code', chunk)`` and ``('string', decoded, line)`` pieces."""
    i = 0
    line = 1
    code: list[str] = []

    def flush() -> None:
        if code:
            yield_buf.append(("code", "".join(code)))
            code.clear()

    yield_buf: list[tuple] = []
    while i < len(text):
        if text.startswith("/*", i):
            flush()
            i += 2
            while i < len(text) and not text.startswith("*/", i):
                if text[i] == "\n":
                    line += 1
                i += 1
            i = min(len(text), i + 2)
            continue
        if text.startswith("//", i):
            flush()
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        if text[i] == '"':
            flush()
            i += 1
            raw: list[str] = []
            start_line = line
            while i < len(text) and text[i] != '"':
                if text[i] == "\\" and i + 1 < len(text):
                    raw.append(text[i : i + 2])
                    if text[i + 1] == "\n":
                        line += 1
                    i += 2
                    continue
                if text[i] == "\n":
                    line += 1
                raw.append(text[i])
                i += 1
            i = min(len(text), i + 1)
            yield_buf.append(("string", _decode_c_string("".join(raw)), start_line))
            continue
        if text[i] == "\n":
            line += 1
        code.append(text[i])
        i += 1
    flush()
    return yield_buf


def _function_bodies(source: str) -> dict[str, str]:
    matches = list(_FUNC_HEAD.finditer(source))
    bodies: dict[str, str] = {}
    for index, match in enumerate(matches):
        name = match.group(1)
        header_end = source.find("{", match.end())
        semi = source.find(";", match.end())
        if header_end < 0 or (semi != -1 and semi < header_end):
            continue
        end = matches[index + 1].start() if index + 1 < len(matches) else len(source)
        bodies[name] = source[match.start() : end]
    return bodies


def _find_c_calls(code: str, name: str) -> list[list[str]]:
    """Return argument texts for each call of ``name`` in comment-free code."""
    calls: list[list[str]] = []
    target = name
    i = 0
    while True:
        j = code.find(target, i)
        if j < 0:
            break
        if j > 0 and (code[j - 1].isalnum() or code[j - 1] == "_"):
            i = j + 1
            continue
        end = j + len(target)
        if end < len(code) and (code[end].isalnum() or code[end] == "_"):
            i = j + 1
            continue
        k = end
        while k < len(code) and code[k].isspace():
            k += 1
        if k >= len(code) or code[k] != "(":
            i = j + 1
            continue
        args, close = _split_c_args(code, k)
        calls.append(args)
        i = close
    return calls


def _split_c_args(code: str, open_paren: int) -> tuple[list[str], int]:
    args: list[str] = []
    depth = 0
    start = open_paren + 1
    i = start
    quote = False
    while i < len(code):
        char = code[i]
        if quote:
            if char == "\\" and i + 1 < len(code):
                i += 2
                continue
            if char == '"':
                quote = False
            i += 1
            continue
        if char == '"':
            quote = True
            i += 1
            continue
        if char == "(":
            depth += 1
        elif char == ")":
            if depth == 0:
                args.append(code[start:i])
                return args, i + 1
            depth -= 1
        elif char == "," and depth == 0:
            args.append(code[start:i])
            start = i + 1
        i += 1
    args.append(code[start:])
    return args, i


def _first_string(arg: str) -> str | None:
    pieces = [item for item in _c_lex(arg) if item[0] == "string"]
    if not pieces:
        return None
    return pieces[0][1]


def normalize_emission(raw: str) -> str | None:
    """Turn one emitted C string into the keyword line the reader would see."""

    def repl(match: re.Match[str]) -> str:
        spec = match.group(0)
        prev = raw[match.start() - 1] if match.start() else ""
        if spec[-1] in "diuoxX" and prev.isalpha():
            return "1"
        return ""

    cleaned = _PRINTF.sub(repl, raw)
    collapsed = " ".join(cleaned.split())
    if not collapsed or collapsed[0] in {"&", "!"}:
        return None
    letters = [char for char in collapsed if char.isalpha()]
    if not letters or any(char.islower() for char in letters):
        return None
    return collapsed


def _emitter_literals(body: str) -> list[str]:
    """C string arguments of append calls, in source order.

    Adjacent appends belong to one deck line until an emitted newline.
    The reader tests that line, not each C fragment.
    """
    pieces: list[tuple[str, int]] = []
    for kind, *rest in _c_lex(body):
        if kind == "code":
            pieces.append(("code", rest[0]))
        else:
            pieces.append(("string", '"' + rest[0].replace("\\", "\\\\").replace('"', '\\"') + '"'))
    code = "".join(text for _, text in pieces)
    events: list[tuple[int, str]] = []
    for name in _EMITTERS:
        target = name
        i = 0
        while True:
            j = code.find(target, i)
            if j < 0:
                break
            if j > 0 and (code[j - 1].isalnum() or code[j - 1] == "_"):
                i = j + 1
                continue
            end = j + len(target)
            if end < len(code) and (code[end].isalnum() or code[end] == "_"):
                i = j + 1
                continue
            k = end
            while k < len(code) and code[k].isspace():
                k += 1
            if k >= len(code) or code[k] != "(":
                i = j + 1
                continue
            args, close = _split_c_args(code, k)
            literal = None
            for arg in args:
                literal = _first_string(arg)
                if literal is not None:
                    break
            if literal is not None:
                events.append((j, literal))
            i = close
    events.sort(key=lambda item: item[0])
    lines: list[str] = []
    pending: list[str] = []
    for _, literal in events:
        pending.append(literal)
        if "\n" in literal:
            phrase = normalize_emission("".join(pending))
            if phrase is not None:
                lines.append(phrase)
            pending = []
    if pending:
        phrase = normalize_emission("".join(pending))
        if phrase is not None:
            lines.append(phrase)
    return lines


def _called_helpers(body: str, bodies: dict[str, str]) -> list[str]:
    code = "".join(rest[0] for kind, *rest in _c_lex(body) if kind == "code")
    names: list[str] = []
    for name in bodies:
        if name in _PRIMITIVE_APPEND or name in SECTION_RENDERERS:
            continue
        if not name.startswith("append_"):
            continue
        if _find_c_calls(code, name):
            names.append(name)
    return names


def render_phrases(params_c: Path) -> dict[str, set[str]]:
    source = params_c.read_text(encoding="utf-8")
    bodies = _function_bodies(source)
    phrases: dict[str, set[str]] = {section: set() for section in set(SECTION_RENDERERS.values())}
    for func, section in SECTION_RENDERERS.items():
        if func not in bodies:
            raise SystemExit(f"missing renderer {func}")
        pending = [func]
        seen: set[str] = set()
        while pending:
            name = pending.pop()
            if name in seen or name not in bodies:
                continue
            seen.add(name)
            body = bodies[name]
            phrases[section].update(_emitter_literals(body))
            pending.extend(_called_helpers(body, bodies))
    return phrases


def _words(phrase: str) -> list[str]:
    words: list[str] = []
    for raw in phrase.upper().split():
        # ``MODE=1`` and ``TOL=`` are the keyword plus an inline value.
        word = raw.split("=", 1)[0].strip(".,;:")
        if word:
            words.append(word)
    return words


def _has_word(words: list[str], literal: str | None) -> bool:
    if not literal:
        return False
    return literal.upper() in words


def token_accepts(phrase: str, token: dict) -> bool:
    words = _words(phrase)
    if token.get("kind") == "index":
        literal = str(token.get("literal") or "").upper().strip()
        # INDEX abbreviations are prefixes (``REPL`` / ``REPLICA``). A proper
        # prefix shorter than 4 characters is not a keyword (``HYP`` matches
        # the EPR flag, not ``HYPERBOLIC``). Interior hits such as ``PRIN``
        # inside ``SPRING`` do not count.
        if not literal:
            return False
        for word in words:
            if word == literal:
                return True
            if len(literal) >= 4 and word.startswith(literal):
                return True
            for piece in word.split("-"):
                if piece == literal:
                    return True
                if len(literal) >= 4 and piece.startswith(literal):
                    return True
        return False
    keyword_ok = _has_word(words, token.get("keyword"))
    if keyword_ok and token.get("and"):
        keyword_ok = _has_word(words, token.get("and"))
    alias_ok = _has_word(words, token.get("alias"))
    if token.get("but_not") and _has_word(words, token.get("but_not")):
        return False
    return keyword_ok or alias_ok


def phrase_accepted(phrase: str, tokens: list[dict]) -> bool:
    return any(token_accepts(phrase, token) for token in tokens)


def words_are_tested(phrase: str, literals: set[str]) -> bool:
    """True when every word is a literal this section's readers test.

    ``and`` and ``alias`` companions count. ``S6`` does not count as
    ``S6GRIMME``, and a line the reader would still reject stays rejected
    when the caller uses :func:`phrase_accepted` instead.
    """
    words = _words(phrase)
    if not words:
        return False
    return all(word in literals for word in words)


def _check_drift(root: Path, snapshot_path: Path, tree: Path | None) -> list[str]:
    if tree is None:
        env = os.environ.get("OPENCPMD_READER_TREE")
        if not env:
            return []
        tree = Path(env)
    if not tree.is_dir():
        return [f"OPENCPMD_READER_TREE is not a directory: {tree}"]
    fresh = extract_tree(tree)
    if not fresh.get("opencpmd_sha"):
        return [f"{tree} has no git HEAD; pass a checkout, not a copy of src/"]
    committed = load_snapshot(snapshot_path)
    if canonical_snapshot(fresh) == canonical_snapshot(committed):
        print(f"OK: snapshot matches {tree} @ {fresh['opencpmd_sha']}")
        return []
    return [
        "snapshot drift: schema/inventory/opencpmd_reader_tokens.json "
        f"differs from {tree} @ {fresh['opencpmd_sha']}"
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpmdc-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--opencpmd-tree",
        type=Path,
        default=None,
        help="rebuild the snapshot from this OpenCPMD tree and fail on drift",
    )
    args = parser.parse_args(argv)
    root = args.cpmdc_root
    params_c = root / "src" / "cpmdc_params.c"
    snapshot_path = root / "schema" / "inventory" / "opencpmd_reader_tokens.json"
    if not params_c.is_file() or not snapshot_path.is_file():
        print(f"missing renderer or snapshot under {root}", file=sys.stderr)
        return 2
    snapshot = load_snapshot(snapshot_path)
    phrases = render_phrases(params_c)
    failures: list[str] = []
    missing_xfail: list[str] = []
    xpass: list[str] = []
    xfails: list[str] = []
    ok = 0
    for section in sorted(phrases):
        tokens = section_tokens(snapshot, section)
        literals = tested_literals(snapshot, section)
        if not tokens and section not in XFAIL:
            # EXTE's reader tests no keyword; other sections must have some.
            if section != "EXTE":
                failures.append(f"{section}: snapshot has no reader tokens")
        seen = phrases[section]
        expected = XFAIL.get(section, set())
        for phrase in sorted(expected):
            if phrase not in seen:
                missing_xfail.append(f"{section} {phrase}")
        for phrase in sorted(seen):
            # The P03 carriers are judged by whether the reader would accept
            # the line. Other lines also pass when every word is a literal
            # that section tests.
            accepted = phrase_accepted(phrase, tokens)
            if phrase in expected:
                if accepted:
                    xpass.append(f"{section} {phrase}")
                else:
                    xfails.append(f"XFAIL {section} {phrase} {XFAIL_REASON}")
                continue
            if accepted or words_are_tested(phrase, literals):
                ok += 1
                continue
            failures.append(f"FAIL {section} {phrase}")
    drift = _check_drift(root, snapshot_path, args.opencpmd_tree)
    for line in xfails:
        print(line)
    for line in failures + missing_xfail + xpass + drift:
        print(line, file=sys.stderr)
    print(
        f"phrases_ok={ok} xfail={len(xfails)} "
        f"fail={len(failures)} missing_xfail={len(missing_xfail)} "
        f"xpass={len(xpass)} drift={len(drift)}"
    )
    if failures or missing_xfail or xpass or drift:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
