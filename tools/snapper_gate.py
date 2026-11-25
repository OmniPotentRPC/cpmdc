#!/usr/bin/env python3
"""Run snapper --check and fail on the diagnostic classes the tree enforces.

A file that would reformat, a ``fused`` line (two sentences on one line), and
a region snapper refused to touch because its render check failed always
fail.
Other classes fail when named in ``--fail-on``; the rest are reported with a
count per class and per file, never dropped.

snapper-fmt 0.11.10 pairs bare ``/``, ``+``, ``*`` and ``_`` across sentences
as emphasis and never breaks before a link, so ``long`` and ``wrap`` stay
advisory: the formatter's own output carries them and a manual wrap would
fail ``--check``.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import subprocess
import sys

ALWAYS = ("reformat", "fused", "render")
RENDER_WARNING = "render check failed"


def _snapper(args: argparse.Namespace, *extra: str) -> list[str]:
    return [
        args.snapper,
        "--native",
        "--check",
        "--max-width",
        str(args.max_width),
        "--config",
        args.config,
        *extra,
    ]


def _render_failures(args: argparse.Namespace, files: list[str]) -> list[tuple[str, str]]:
    found = []
    for name in files:
        result = subprocess.run(
            _snapper(args, "--", name), capture_output=True, text=True, check=False
        )
        for line in result.stderr.splitlines():
            if RENDER_WARNING in line:
                found.append((name, line.strip()))
    return found


def run(args: argparse.Namespace) -> int:
    fail_on = set(ALWAYS) | {k for k in args.fail_on.split(",") if k}
    result = subprocess.run(
        _snapper(args, "--output-format", "json", "--", *args.files),
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode not in (0, 1):
        sys.stderr.write(result.stdout + result.stderr)
        return result.returncode
    try:
        report = json.loads(result.stdout or "[]")
    except json.JSONDecodeError:
        sys.stderr.write(result.stdout + result.stderr)
        print("snapper --check did not return JSON", file=sys.stderr)
        return 2

    failures: list[str] = []
    kinds: collections.Counter[str] = collections.Counter()
    per_file: dict[str, collections.Counter[str]] = {}
    for entry in report:
        name = entry["file"]
        if entry.get("would_reformat"):
            kinds["reformat"] += 1
            failures.append(f"{name}: reformat: snapper would rewrap this file")
        for diag in entry.get("diagnostics", []):
            kind = diag["kind"]
            kinds[kind] += 1
            if kind in fail_on:
                failures.append(f"{name}:{diag['line']}: {kind}: {diag['excerpt']}")
            else:
                per_file.setdefault(name, collections.Counter())[kind] += 1

    if RENDER_WARNING in result.stderr:
        for name, line in _render_failures(args, args.files):
            kinds["render"] += 1
            failures.append(f"{name}: render: {line}")

    for line in failures:
        print(line)
    advisory = sorted({k for c in per_file.values() for k in c})
    if per_file:
        print(f"advisory ({', '.join(advisory)}) per file:")
        for name in sorted(per_file):
            counts = ", ".join(f"{k} {n}" for k, n in sorted(per_file[name].items()))
            print(f"  {name}: {counts}")
    summary = ", ".join(f"{k} {kinds.get(k, 0)}" for k in sorted(set(kinds) | set(ALWAYS)))
    enforced = ", ".join(sorted(fail_on))
    verdict = "FAIL" if failures else "OK"
    print(
        f"snapper {verdict}: {len(report)} files reported; {summary}; "
        f"enforced: {enforced}"
    )

    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a", encoding="utf-8") as handle:
            handle.write(f"### Snapper {verdict}\n\n| Class | Count | Enforced |\n|---|---|---|\n")
            for kind in sorted(set(kinds) | set(ALWAYS)):
                handle.write(f"| {kind} | {kinds.get(kind, 0)} | {'yes' if kind in fail_on else 'no'} |\n")
    return 1 if failures else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--snapper", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--max-width", type=int, default=80)
    parser.add_argument(
        "--fail-on",
        default="",
        help="extra comma-separated classes that fail, for example long,wrap",
    )
    parser.add_argument("files", nargs="+")
    return run(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
