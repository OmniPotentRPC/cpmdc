#!/usr/bin/env python3
"""Run non-skipping towncrier newsfragment and changelog checks.

Plain `towncrier check` defaults to comparing against `origin/main` and exits
with "Checks SKIPPED: news file changes detected." whenever `CHANGELOG.md`
differs between `HEAD` and that compare branch. Because this repository lands
changes on `develop` and release or fold branches edit `CHANGELOG.md`, this
wrapper selects the integration target (`origin/${GITHUB_BASE_REF}` or
`origin/develop`), always runs strict fragment validation and a draft build,
runs `towncrier check --compare-with ...` when `CHANGELOG.md` is untouched,
and verifies release/fold invariants when `CHANGELOG.md` changed.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STRICT_INSPECT_CODE = r"""
import json
import sys
from towncrier._builder import find_fragments
from towncrier._settings import load_config_from_options

root = sys.argv[1]
base_dir, config = load_config_from_options(root, None)
_, fragment_files = find_fragments(base_dir, config, strict=True)
print(json.dumps({
    "filename": config.filename,
    "directory": config.directory or "",
    "fragments": [path for path, _ in fragment_files],
}))
"""


def _git(root: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        capture_output=True,
        text=True,
        check=check,
    )


def _ref_exists(root: Path, ref: str) -> bool:
    return _git(root, "rev-parse", "--verify", "--quiet", ref, check=False).returncode == 0


def resolve_compare_branch(root: Path, explicit: str | None) -> str:
    if explicit:
        if not _ref_exists(root, explicit):
            raise SystemExit(f"towncrier-check: compare ref not found: {explicit}")
        return explicit

    base_ref = os.environ.get("GITHUB_BASE_REF", "").strip()
    if base_ref:
        for candidate in (f"origin/{base_ref}", base_ref):
            if _ref_exists(root, candidate):
                return candidate
        raise SystemExit(
            f"towncrier-check: GITHUB_BASE_REF={base_ref!r} could not be resolved"
        )

    for candidate in ("origin/develop", "develop", "origin/main", "main"):
        if _ref_exists(root, candidate):
            return candidate
    raise SystemExit("towncrier-check: could not resolve a base branch to compare against")


def _towncrier_python_cmd() -> list[str]:
    if importlib.util.find_spec("towncrier") is not None:
        return [sys.executable]
    if shutil.which("uvx"):
        return ["uvx", "--from", "towncrier", "python"]
    raise SystemExit("towncrier-check: neither towncrier nor uvx is available")


def _towncrier_cli_cmd() -> list[str]:
    if shutil.which("towncrier"):
        return ["towncrier"]
    if shutil.which("uvx"):
        return ["uvx", "towncrier"]
    if importlib.util.find_spec("towncrier") is not None:
        return [sys.executable, "-m", "towncrier"]
    raise SystemExit("towncrier-check: neither towncrier nor uvx is available")


def inspect_fragments_strict(root: Path) -> tuple[str, list[str]]:
    cmd = [*_towncrier_python_cmd(), "-c", STRICT_INSPECT_CODE, str(root)]
    proc = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if proc.returncode != 0:
        if proc.stdout:
            print(proc.stdout, file=sys.stderr, end="")
        if proc.stderr:
            print(proc.stderr, file=sys.stderr, end="")
        raise SystemExit(proc.returncode or 1)
    data = json.loads(proc.stdout)
    return str(data["filename"]), list(data["fragments"])


def run_draft_build(root: Path, version: str) -> str:
    cmd = [*_towncrier_cli_cmd(), "build", "--draft", "--version", version]
    proc = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    if proc.returncode != 0:
        if proc.stdout:
            print(proc.stdout, file=sys.stderr, end="")
        if proc.stderr:
            print(proc.stderr, file=sys.stderr, end="")
        raise SystemExit(proc.returncode or 1)
    return proc.stdout


def changed_files(root: Path, compare_with: str, staged: bool) -> list[str]:
    out = _git(root, "diff", "--name-only", f"{compare_with}...").stdout.splitlines()
    if staged:
        out.extend(_git(root, "diff", "--name-only", "--cached").stdout.splitlines())
    seen: list[str] = []
    for rel in out:
        rel = rel.strip()
        if rel and rel not in seen:
            seen.append(rel)
    return seen


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate towncrier fragments and changelog state without skipping."
    )
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--compare-with", default=None)
    parser.add_argument("--staged", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    version = (root / "version.txt").read_text(encoding="utf-8").strip()
    compare_with = resolve_compare_branch(root, args.compare_with)

    news_filename, fragments = inspect_fragments_strict(root)
    draft_output = run_draft_build(root, version)
    files = changed_files(root, compare_with, args.staged)

    news_abs = os.path.normpath(os.path.join(str(root), news_filename))
    files_abs = {os.path.normpath(os.path.join(str(root), rel)) for rel in files}

    if not files:
        print(f"towncrier-check: compare_with={compare_with} (no branch diff)")
        print(
            f"towncrier-check: strict fragment validation passed "
            f"({len(fragments)} pending fragment(s))"
        )
        if "No significant changes." in draft_output:
            print(f"towncrier-check: towncrier build --draft --version {version}: No significant changes.")
        return 0

    if news_abs not in files_abs:
        cmd = [*_towncrier_cli_cmd(), "check", "--compare-with", compare_with]
        if args.staged:
            cmd.append("--staged")
        proc = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
        if proc.stdout:
            print(proc.stdout, end="")
        if proc.stderr:
            print(proc.stderr, file=sys.stderr, end="")
        return proc.returncode

    print(f"towncrier-check: compare_with={compare_with}")
    print("Looking at these files:")
    print("----")
    for idx, rel in enumerate(files, start=1):
        print(f"{idx}. {rel}")
    print("----")

    changelog_text = (root / news_filename).read_text(encoding="utf-8")
    marker = f"## [{version}]"
    if marker not in changelog_text:
        print(
            f"towncrier-check: {news_filename} modified but missing section {marker!r}",
            file=sys.stderr,
        )
        return 1

    release_assert = root / "scripts" / "release_assert.py"
    if release_assert.is_file():
        proc = subprocess.run(
            [sys.executable, str(release_assert), "--require-changelog"],
            cwd=root,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            if proc.stdout:
                print(proc.stdout, file=sys.stderr, end="")
            if proc.stderr:
                print(proc.stderr, file=sys.stderr, end="")
            return proc.returncode

    if fragments:
        print(
            f"towncrier-check: {news_filename} was updated for {version}, "
            "but pending newsfragments remain:",
            file=sys.stderr,
        )
        for idx, path in enumerate(fragments, start=1):
            print(f"  {idx}. {os.path.relpath(path, root)}", file=sys.stderr)
        return 1

    if "No significant changes." not in draft_output:
        print(
            "towncrier-check: towncrier build --draft still reports pending changes",
            file=sys.stderr,
        )
        return 1

    print(
        f"towncrier-check: verified {news_filename} contains {marker}, "
        f"0 pending fragments remain, and towncrier build --draft --version {version} "
        "reports no significant changes."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
