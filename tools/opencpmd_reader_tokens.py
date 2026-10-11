#!/usr/bin/env python3
"""Reader-scoped OpenCPMD keyword snapshot.

For each input section, find the procedure that calls ``inscan('&SECTION')``
and the readers that procedure delegates to. Fortran comments are removed
first. String literals count only when ``keyword_contains(line, ...)`` or
``INDEX(line, ...)`` tests them. A literal in a comment, a ``WRITE``, or
another section stays out of this section.

The snapshot records the OpenCPMD source SHA and ``EXTRACTOR_VERSION``.
Refresh it with this program. The fidelity test reads the snapshot and does
not need an OpenCPMD tree.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

EXTRACTOR_VERSION = "1"

_FORTRAN_SUFFIXES = {".f", ".f90", ".f95", ".f03", ".f08"}
_PROC_START = re.compile(
    r"""^\s*
    (?:(?:recursive|pure|elemental|impure)\s+)*
    (?:
        (?:integer|real|complex|logical|double\s+precision|character|type)
        \s*(?:\([^)\n]*\))?\s+
    )?
    (subroutine|function)\s+([A-Za-z_][A-Za-z0-9_]*)
    """,
    re.IGNORECASE | re.VERBOSE,
)
_END_PROC = re.compile(r"^\s*end\s+(subroutine|function)\b", re.IGNORECASE)
_END_BARE = re.compile(r"^\s*end\s*$", re.IGNORECASE)
_BLOCK_OPEN = re.compile(
    r"^\s*(?:if\b.*\bthen\b|do\b|select\s+case\b|where\b|forall\b|associate\b|block\b|interface\b)",
    re.IGNORECASE,
)
_BLOCK_CLOSE = re.compile(
    r"^\s*end\s+(?:if|do|select|where|forall|associate|block|interface)\b",
    re.IGNORECASE,
)
_SECTION_LITERAL = re.compile(r"&([A-Za-z][A-Za-z0-9_]*)\s*$")
_CALL_NAME = re.compile(r"\bcall\s+([A-Za-z_][A-Za-z0-9_]*)\b", re.IGNORECASE)
_NAMED_ARG = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", re.IGNORECASE)
_IDENT = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*)")


@dataclass
class LogicalLine:
    text: str
    columns: list[int]


@dataclass
class Token:
    kind: str
    line: int
    keyword: str | None = None
    and_word: str | None = None
    alias: str | None = None
    but_not: str | None = None
    literal: str | None = None

    def as_json(self) -> dict[str, object]:
        out: dict[str, object] = {"kind": self.kind, "line": self.line}
        if self.kind == "keyword_contains":
            out["keyword"] = self.keyword or ""
            if self.and_word:
                out["and"] = self.and_word
            if self.alias:
                out["alias"] = self.alias
            if self.but_not:
                out["but_not"] = self.but_not
        else:
            out["literal"] = self.literal or ""
        return out

    def sort_key(self) -> tuple:
        return (
            self.line,
            self.kind,
            self.keyword or "",
            self.and_word or "",
            self.alias or "",
            self.but_not or "",
            self.literal or "",
        )


@dataclass
class Procedure:
    name: str
    file: str
    start_line: int
    body: list[LogicalLine] = field(default_factory=list)
    tokens: list[Token] = field(default_factory=list)
    calls: list[tuple[int, str]] = field(default_factory=list)
    inscans: list[tuple[int, str, int]] = field(default_factory=list)
    # inscans: (logical index, section, source line)


def _strip_fortran_comment(line: str) -> str:
    out: list[str] = []
    quote: str | None = None
    i = 0
    while i < len(line):
        char = line[i]
        if quote is not None:
            out.append(char)
            if char == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    out.append(line[i + 1])
                    i += 2
                    continue
                quote = None
            i += 1
            continue
        if char in {"'", '"'}:
            quote = char
            out.append(char)
            i += 1
            continue
        if char == "!":
            break
        out.append(char)
        i += 1
    return "".join(out).rstrip()


def _continuation_ampersand(code: str) -> int | None:
    """Index of a free-form continuation '&' at the end of the line, if any."""
    stripped = code.rstrip()
    if not stripped.endswith("&"):
        return None
    quote: str | None = None
    i = 0
    amp: int | None = None
    while i < len(stripped):
        char = stripped[i]
        if quote is not None:
            if char == quote:
                if i + 1 < len(stripped) and stripped[i + 1] == quote:
                    i += 2
                    continue
                quote = None
            i += 1
            continue
        if char in {"'", '"'}:
            quote = char
            i += 1
            continue
        if char == "&":
            amp = i
        i += 1
    if amp is None or stripped[amp + 1 :].strip():
        return None
    return amp


def fortran_logical_lines(text: str) -> list[LogicalLine]:
    """Join free-form continuations after comments are removed.

    ``columns[i]`` is the 1-based physical line that produced ``text[i]``.
    """
    physical = text.splitlines()
    logical: list[LogicalLine] = []
    pending_text: list[str] = []
    pending_cols: list[int] = []
    pending = False
    for number, raw in enumerate(physical, start=1):
        code = _strip_fortran_comment(raw)
        if pending:
            code = code.lstrip()
            if code.startswith("&"):
                code = code[1:]
        amp = _continuation_ampersand(code)
        if amp is not None:
            code = code[:amp]
            pending_text.append(code)
            pending_cols.extend([number] * len(code))
            pending = True
            continue
        pending_text.append(code)
        pending_cols.extend([number] * len(code))
        logical.append(LogicalLine("".join(pending_text), pending_cols))
        pending_text = []
        pending_cols = []
        pending = False
    if pending_text:
        logical.append(LogicalLine("".join(pending_text), pending_cols))
    return logical


def _parse_args(code: str, open_paren: int) -> tuple[list[tuple[str, int]], int]:
    """Split a parenthesised list. Each arg is ``(text, index of its start)``."""
    args: list[tuple[str, int]] = []
    depth = 0
    start = open_paren + 1
    i = start
    quote: str | None = None
    while i < len(code):
        char = code[i]
        if quote is not None:
            if char == quote:
                if i + 1 < len(code) and code[i + 1] == quote:
                    i += 2
                    continue
                quote = None
            i += 1
            continue
        if char in {"'", '"'}:
            quote = char
            i += 1
            continue
        if char == "(":
            depth += 1
        elif char == ")":
            if depth == 0:
                args.append((code[start:i], start))
                return args, i + 1
            depth -= 1
        elif char == "," and depth == 0:
            args.append((code[start:i], start))
            start = i + 1
        i += 1
    args.append((code[start:], start))
    return args, i


def _find_calls(code: str, name: str) -> list[tuple[int, list[tuple[str, int]]]]:
    target = name.lower()
    low = code.lower()
    found: list[tuple[int, list[tuple[str, int]]]] = []
    i = 0
    while True:
        j = low.find(target, i)
        if j < 0:
            break
        if j > 0 and (low[j - 1].isalnum() or low[j - 1] == "_"):
            i = j + 1
            continue
        end = j + len(target)
        if end < len(low) and (low[end].isalnum() or low[end] == "_"):
            i = j + 1
            continue
        k = end
        while k < len(code) and code[k].isspace():
            k += 1
        if k >= len(code) or code[k] != "(":
            i = j + 1
            continue
        # Skip the call when it sits inside a string. The logical line has
        # already lost comments, so a quote walk from the start is enough.
        if _index_in_string(code, j):
            i = j + 1
            continue
        args, close = _parse_args(code, k)
        found.append((j, args))
        i = close
    return found


def _index_in_string(code: str, index: int) -> bool:
    quote: str | None = None
    i = 0
    while i < index and i < len(code):
        char = code[i]
        if quote is not None:
            if char == quote:
                if i + 1 < len(code) and code[i + 1] == quote:
                    i += 2
                    continue
                quote = None
            i += 1
            continue
        if char in {"'", '"'}:
            quote = char
        i += 1
    return quote is not None


def _string_literals(expr: str, expr_start: int) -> list[tuple[str, int]]:
    """Return ``(decoded, absolute index)`` for each Fortran string literal."""
    found: list[tuple[str, int]] = []
    i = 0
    while i < len(expr):
        if expr[i] not in {"'", '"'}:
            i += 1
            continue
        quote = expr[i]
        start = i
        i += 1
        chars: list[str] = []
        while i < len(expr):
            if expr[i] == quote:
                if i + 1 < len(expr) and expr[i + 1] == quote:
                    chars.append(quote)
                    i += 2
                    continue
                i += 1
                break
            chars.append(expr[i])
            i += 1
        found.append(("".join(chars), expr_start + start))
    return found


def _joined_literal(expr: str, expr_start: int) -> tuple[str, int] | None:
    parts = _string_literals(expr, expr_start)
    if not parts:
        return None
    return "".join(part for part, _ in parts), parts[0][1]


# ``previous_line`` is the keyword line after the reader has already pulled
# the next record into ``line``. ``next_words`` is the current word of that
# line. Both are keyword tests, not output text.
_LINE_VARS = {"line", "previous_line", "next_words"}


def _is_line_expr(expr: str) -> bool:
    match = _IDENT.match(expr)
    return bool(match) and match.group(1).lower() in _LINE_VARS


def _line_of(columns: list[int], index: int) -> int:
    if not columns:
        return 1
    if index < 0:
        return columns[0]
    if index >= len(columns):
        return columns[-1]
    return columns[index]


def _keyword_token(
    args: list[tuple[str, int]], columns: list[int]
) -> Token | None:
    if not args or not _is_line_expr(args[0][0]):
        return None
    positional: list[tuple[str, int]] = []
    named: dict[str, tuple[str, int]] = {}
    for text, start in args:
        match = _NAMED_ARG.match(text)
        if match:
            named[match.group(1).lower()] = (text[match.end() :], start + match.end())
        else:
            positional.append((text, start))
    roles = {
        "keyword": positional[1] if len(positional) > 1 else None,
        "and": positional[2] if len(positional) > 2 else None,
        "alias": positional[3] if len(positional) > 3 else None,
        "but_not": positional[4] if len(positional) > 4 else None,
    }
    for role in ("keyword", "and", "alias", "but_not", "cut_at"):
        if role in named:
            roles[role] = named[role]
    decoded: dict[str, str | None] = {}
    line = 1
    for role in ("keyword", "and", "alias", "but_not"):
        slot = roles.get(role)
        if slot is None:
            decoded[role] = None
            continue
        joined = _joined_literal(slot[0], slot[1])
        if joined is None:
            decoded[role] = None
            continue
        decoded[role] = joined[0]
        if role == "keyword":
            line = _line_of(columns, joined[1])
    if not decoded.get("keyword") and not decoded.get("alias"):
        return None
    return Token(
        kind="keyword_contains",
        line=line,
        keyword=decoded.get("keyword"),
        and_word=decoded.get("and"),
        alias=decoded.get("alias"),
        but_not=decoded.get("but_not"),
    )


def _index_token(args: list[tuple[str, int]], columns: list[int]) -> Token | None:
    if len(args) < 2 or not _is_line_expr(args[0][0]):
        return None
    joined = _joined_literal(args[1][0], args[1][1])
    if joined is None or not any(ch.isalpha() for ch in joined[0]):
        return None
    return Token(
        kind="index",
        line=_line_of(columns, joined[1]),
        literal=joined[0],
    )


def _delimiter_token(args: list[tuple[str, int]], columns: list[int]) -> Token | None:
    """``index_of_delimiter(line, keyword, delim)`` tests ``keyword`` on the line."""
    if len(args) < 2 or not _is_line_expr(args[0][0]):
        return None
    joined = _joined_literal(args[1][0], args[1][1])
    if joined is None or not any(ch.isalpha() for ch in joined[0]):
        return None
    return Token(
        kind="keyword_contains",
        line=_line_of(columns, joined[1]),
        keyword=joined[0],
    )


def _section_literal(args: list[tuple[str, int]]) -> str | None:
    for text, start in args:
        joined = _joined_literal(text, start)
        if joined is None:
            continue
        match = _SECTION_LITERAL.fullmatch(joined[0].strip())
        if match:
            return match.group(1).upper()
    return None


def _call_names(code: str) -> list[str]:
    names: list[str] = []
    for match in _CALL_NAME.finditer(code):
        if _index_in_string(code, match.start()):
            continue
        names.append(match.group(1).lower())
    return names


def _parse_procedures(relative: str, text: str) -> list[Procedure]:
    logical = fortran_logical_lines(text)
    procedures: list[Procedure] = []
    stack: list[tuple[Procedure, int]] = []
    for line in logical:
        start = _PROC_START.match(line.text)
        if start:
            # A procedure header can appear inside a procedure (CONTAINS).
            # Block depth does not include the procedure itself.
            proc = Procedure(
                name=start.group(2).lower(),
                file=relative,
                start_line=_line_of(line.columns, start.start(2)),
            )
            procedures.append(proc)
            stack.append((proc, 0))
            continue
        if not stack:
            continue
        proc, depth = stack[-1]
        if _END_PROC.match(line.text) or (_END_BARE.match(line.text) and depth == 0):
            stack.pop()
            continue
        if _BLOCK_CLOSE.match(line.text):
            stack[-1] = (proc, max(0, depth - 1))
        elif _BLOCK_OPEN.match(line.text):
            stack[-1] = (proc, depth + 1)
        proc.body.append(line)
        for _, args in _find_calls(line.text, "keyword_contains"):
            token = _keyword_token(args, line.columns)
            if token is not None:
                proc.tokens.append(token)
        for _, args in _find_calls(line.text, "index"):
            token = _index_token(args, line.columns)
            if token is not None:
                proc.tokens.append(token)
        for _, args in _find_calls(line.text, "index_of_delimiter"):
            token = _delimiter_token(args, line.columns)
            if token is not None:
                proc.tokens.append(token)
        for _, args in _find_calls(line.text, "inscan"):
            section = _section_literal(args)
            if section is None:
                continue
            source = _line_of(line.columns, 0)
            proc.inscans.append((len(proc.body) - 1, section, source))
        for name in _call_names(line.text):
            proc.calls.append((len(proc.body) - 1, name))
    return procedures


def _fortran_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if any(part.startswith(".") for part in path.relative_to(root).parts):
            continue
        suffixes = {path.suffix.lower()}
        name = path.name
        if ".f90" in name.lower() or ".f" in suffixes or suffixes & _FORTRAN_SUFFIXES:
            # OpenCPMD names many readers ``*.mod.F90``.
            if path.suffix.lower() in _FORTRAN_SUFFIXES or name.lower().endswith(
                (".f90", ".f", ".f95")
            ):
                files.append(path)
    return sorted(files)


def _relative(root: Path, path: Path) -> str:
    return path.relative_to(root).as_posix()


def _regions(proc: Procedure) -> list[tuple[str, int, int, int]]:
    """``(section, body_start, body_end, source_line)`` for each inscan."""
    if not proc.inscans:
        return []
    regions: list[tuple[str, int, int, int]] = []
    for index, (body_index, section, source) in enumerate(proc.inscans):
        end = proc.inscans[index + 1][0] if index + 1 < len(proc.inscans) else len(proc.body)
        regions.append((section, body_index, end, source))
    return regions


def _tokens_in(proc: Procedure, start: int, end: int) -> list[Token]:
    if not proc.body:
        return []
    selected: list[Token] = []
    # Tokens were collected in source order while walking the body. Re-scan
    # the slice so a token that belongs to a later inscan stays there.
    for line in proc.body[start:end]:
        for _, args in _find_calls(line.text, "keyword_contains"):
            token = _keyword_token(args, line.columns)
            if token is not None:
                selected.append(token)
        for _, args in _find_calls(line.text, "index"):
            token = _index_token(args, line.columns)
            if token is not None:
                selected.append(token)
        for _, args in _find_calls(line.text, "index_of_delimiter"):
            token = _delimiter_token(args, line.columns)
            if token is not None:
                selected.append(token)
    return selected


def _calls_in(proc: Procedure, start: int, end: int) -> list[tuple[int, str]]:
    return [(index, name) for index, name in proc.calls if start <= index < end]


def _guard_literals(proc: Procedure, call_index: int) -> list[str]:
    for index in range(call_index, -1, -1):
        line = proc.body[index]
        literals: list[str] = []
        for _, args in _find_calls(line.text, "keyword_contains"):
            token = _keyword_token(args, line.columns)
            if token is None:
                continue
            for word in (token.keyword, token.and_word, token.alias, token.but_not):
                if word:
                    literals.append(word)
        for _, args in _find_calls(line.text, "index"):
            token = _index_token(args, line.columns)
            if token is not None and token.literal:
                literals.append(token.literal)
        if literals:
            return literals
        if index != call_index and _CALL_NAME.search(line.text):
            break
        stripped = line.text.strip().lower()
        if stripped.startswith("do ") or stripped.startswith("subroutine "):
            break
    return []


def _git_sha(root: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    sha = out.stdout.strip()
    return sha or None


def extract_tree(root: Path, sha: str | None = None) -> dict:
    """Build the reader-scoped token snapshot for an OpenCPMD source tree."""
    root = root.resolve()
    procedures: list[Procedure] = []
    for path in _fortran_files(root):
        text = path.read_text(encoding="utf-8", errors="replace")
        procedures.extend(_parse_procedures(_relative(root, path), text))
    by_name: dict[str, list[Procedure]] = {}
    for proc in procedures:
        by_name.setdefault(proc.name, []).append(proc)
    direct = {id(proc): bool(proc.tokens) for proc in procedures}

    def resolve(name: str, caller_file: str) -> list[Procedure]:
        candidates = by_name.get(name.lower(), [])
        same = [proc for proc in candidates if proc.file == caller_file]
        if len(same) == 1:
            return same
        if len(same) > 1:
            return same
        if len(candidates) == 1:
            return candidates
        tested = [proc for proc in candidates if direct[id(proc)]]
        if len(tested) == 1:
            return tested
        return []

    memo: dict[int, bool] = {}

    def reaches_reader(proc: Procedure, depth: int, seen: set[int]) -> bool:
        key = id(proc)
        if key in memo and depth >= 3:
            return memo[key]
        if direct[key]:
            memo[key] = True
            return True
        if depth <= 0 or key in seen:
            return False
        seen.add(key)
        for _, name in proc.calls:
            for target in resolve(name, proc.file):
                if reaches_reader(target, depth - 1, seen):
                    memo[key] = True
                    return True
        return False

    def owned_sections(proc: Procedure) -> set[str]:
        return {section for _, section, _ in proc.inscans}

    sections: dict[str, dict] = {}

    def ensure(section: str) -> dict:
        return sections.setdefault(
            section, {"owners": [], "readers": []}
        )

    def reader_node(
        proc: Procedure,
        tokens: list[Token],
        call_line: int | None,
        guard: list[str],
        section: str,
        depth: int,
        seen: set[tuple[str, str]],
    ) -> dict:
        node: dict = {
            "file": proc.file,
            "procedure": proc.name,
            "tokens": tokens,
        }
        if call_line is not None:
            node["call_line"] = call_line
        if guard:
            node["guard"] = guard
        children = []
        # A procedure that owns an inscan is an owner, not a nested reader.
        if depth < 6 and not proc.inscans:
            for index, name in proc.calls:
                children.extend(
                    _delegate(proc, index, name, section, depth + 1, seen)
                )
        if children:
            node["readers"] = children
        return node

    def _delegate(
        caller: Procedure,
        call_index: int,
        name: str,
        section: str,
        depth: int,
        seen: set[tuple[str, str]],
    ) -> list[dict]:
        nodes: list[dict] = []
        for target in resolve(name, caller.file):
            key = (target.file, target.name)
            if key in seen or depth > 6:
                continue
            if owned_sections(target):
                continue
            if not reaches_reader(target, 3, set()):
                continue
            seen.add(key)
            tokens = list(target.tokens)
            guard = _guard_literals(caller, call_index)
            call_line = _line_of(caller.body[call_index].columns, 0)
            nodes.append(
                reader_node(target, tokens, call_line, guard, section, depth, seen)
            )
        return nodes

    for proc in procedures:
        for section, start, end, source in _regions(proc):
            bucket = ensure(section)
            tokens = _tokens_in(proc, start, end)
            seen = {(proc.file, proc.name)}
            children: list[dict] = []
            for index, name in _calls_in(proc, start, end):
                children.extend(_delegate(proc, index, name, section, 1, seen))
            owner = {
                "file": proc.file,
                "procedure": proc.name,
                "line": source,
                "tokens": tokens,
            }
            if children:
                owner["readers"] = children
            bucket["owners"].append(owner)

    def sort_tokens(tokens: list[Token]) -> list[dict]:
        return [token.as_json() for token in sorted(tokens, key=lambda item: item.sort_key())]

    def sort_reader(node: dict) -> dict:
        out = {
            "file": node["file"],
            "procedure": node["procedure"],
            "tokens": sort_tokens(node.get("tokens", [])),
        }
        if "line" in node:
            out["line"] = node["line"]
        if "call_line" in node:
            out["call_line"] = node["call_line"]
        if node.get("guard"):
            out["guard"] = list(node["guard"])
        readers = node.get("readers") or []
        if readers:
            out["readers"] = [
                sort_reader(child)
                for child in sorted(
                    readers,
                    key=lambda child: (
                        child["file"],
                        child["procedure"],
                        child.get("call_line", 0),
                    ),
                )
            ]
        return out

    encoded = {
        "extractor_version": EXTRACTOR_VERSION,
        "opencpmd_sha": sha or _git_sha(root) or "",
        "sections": {},
    }
    for name in sorted(sections):
        bucket = sections[name]
        owners = [
            sort_reader(owner)
            for owner in sorted(
                bucket["owners"],
                key=lambda owner: (owner["file"], owner["line"], owner["procedure"]),
            )
        ]
        encoded["sections"][name] = {"owners": owners}
    return encoded


def canonical_snapshot(snapshot: dict) -> str:
    return json.dumps(snapshot, indent=2, sort_keys=True) + "\n"


def write_snapshot(snapshot: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical_snapshot(snapshot), encoding="utf-8")


def load_snapshot(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def iter_tokens(node: dict):
    yield from node.get("tokens", [])
    for child in node.get("readers", []):
        yield from iter_tokens(child)


def section_tokens(snapshot: dict, section: str) -> list[dict]:
    """Every token tested by the section's owners and the readers they call."""
    bucket = snapshot.get("sections", {}).get(section.upper())
    if bucket is None:
        return []
    tokens: list[dict] = []
    for owner in bucket.get("owners", []):
        tokens.extend(iter_tokens(owner))
    return tokens


def tested_literals(snapshot: dict, section: str) -> set[str]:
    """Literals the section's readers accept, uppercased.

    ``and`` and ``alias`` companions count. ``but_not`` is an exclusion, so
    it does not. A string that only appears in a comment or a ``WRITE``
    does not.
    """
    found: set[str] = set()
    for token in section_tokens(snapshot, section):
        if token.get("kind") == "index":
            literal = token.get("literal") or ""
            if literal:
                found.add(literal.upper())
            continue
        for key in ("keyword", "and", "alias"):
            word = token.get(key)
            if word:
                found.add(str(word).upper())
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path, help="OpenCPMD source tree")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="write the snapshot JSON here",
    )
    parser.add_argument(
        "--sha",
        help="source SHA to record (default: git rev-parse HEAD in the tree)",
    )
    parser.add_argument(
        "--section",
        help="print tested literals for one section and exit",
    )
    args = parser.parse_args(argv)
    if not args.tree.is_dir():
        print(f"not a directory: {args.tree}", file=sys.stderr)
        return 2
    snapshot = extract_tree(args.tree, sha=args.sha)
    if not snapshot["opencpmd_sha"]:
        print("pass --sha; the tree has no git HEAD", file=sys.stderr)
        return 2
    if args.section:
        literals = sorted(tested_literals(snapshot, args.section))
        print(f"{args.section.upper()}: {len(literals)} literals")
        for literal in literals:
            print(literal)
        return 0
    if args.output is None:
        print("pass -o or --section", file=sys.stderr)
        return 2
    write_snapshot(snapshot, args.output)
    nsections = len(snapshot["sections"])
    ntokens = sum(
        len(section_tokens(snapshot, name)) for name in snapshot["sections"]
    )
    print(
        f"wrote {args.output} sections={nsections} tokens={ntokens} "
        f"sha={snapshot['opencpmd_sha']} extractor={EXTRACTOR_VERSION}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
