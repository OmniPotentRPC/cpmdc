#!/usr/bin/env python3
"""Verify llms.txt / llms-full.txt freshness, determinism, and actionable diffs."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGEN_COMMAND = "python3 scripts/generate_llms_txt.py"


def _copy_generator_inputs(source: Path, dest: Path) -> None:
    shutil.copytree(source / "docs" / "orgmode", dest / "docs" / "orgmode")
    (dest / "schema" / "inventory").mkdir(parents=True)
    shutil.copy2(
        source / "schema" / "inventory" / "cpmd_features.json",
        dest / "schema" / "inventory" / "cpmd_features.json",
    )
    (dest / "scripts").mkdir(parents=True)
    shutil.copy2(
        source / "scripts" / "generate_llms_txt.py",
        dest / "scripts" / "generate_llms_txt.py",
    )


def _run_generator(tree: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(tree / "scripts" / "generate_llms_txt.py"), *args],
        cwd=tree,
        capture_output=True,
        text=True,
        check=False,
    )


def test_deterministic_generation(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)

        first = _run_generator(tree)
        if first.returncode != 0:
            raise AssertionError(
                f"first generation failed: {first.stderr}{first.stdout}"
            )
        first_summary = (tree / "llms.txt").read_bytes()
        first_full = (tree / "llms-full.txt").read_bytes()

        second = _run_generator(tree)
        if second.returncode != 0:
            raise AssertionError(
                f"second generation failed: {second.stderr}{second.stdout}"
            )
        second_summary = (tree / "llms.txt").read_bytes()
        second_full = (tree / "llms-full.txt").read_bytes()

        if first_summary != second_summary:
            raise AssertionError("llms.txt output differed across two runs")
        if first_full != second_full:
            raise AssertionError("llms-full.txt output differed across two runs")

        checked = _run_generator(tree, "--check")
        if checked.returncode != 0:
            raise AssertionError(
                f"--check failed immediately after generation: {checked.stderr}"
            )

        abs_markers = (str(tree), str(tree.resolve()), str(root), str(root.resolve()))
        for label, payload in (
            ("llms.txt", first_summary.decode("utf-8")),
            ("llms-full.txt", first_full.decode("utf-8")),
            ("stdout", first.stdout),
        ):
            for marker in abs_markers:
                if marker in payload:
                    raise AssertionError(
                        f"{label} leaked absolute path {marker!r}"
                    )
            if (
                "## File: docs/orgmode/tutorials/quickstart.org" in payload
                or 'http-equiv="refresh"' in payload
                or "cpmdc.rgoswami.me/tutorials/quickstart" in payload
            ):
                raise AssertionError(
                    f"{label} leaked orphan redirect page tutorials/quickstart"
                )


def test_stale_failure_is_actionable(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generated = _run_generator(tree)
        if generated.returncode != 0:
            raise AssertionError(
                f"setup generation failed: {generated.stderr}{generated.stdout}"
            )

        abs_markers = (str(tree), str(tree.resolve()), str(root), str(root.resolve()))
        for rel in ("llms-full.txt", "llms.txt"):
            target = tree / rel
            original = target.read_text(encoding="utf-8")
            stale_tail = "".join(f"stale line {idx}\n" for idx in range(30))
            target.write_text(original + stale_tail, encoding="utf-8")

            checked = _run_generator(tree, "--check")
            if checked.returncode == 0:
                raise AssertionError(f"--check unexpectedly passed for stale {rel}")
            msg = checked.stderr
            for required in (
                rel,
                REGEN_COMMAND,
                f"--- a/{rel}",
                f"+++ b/{rel}",
                "@@",
                "more diff lines",
            ):
                if required not in msg:
                    raise AssertionError(
                        f"stale failure for {rel} missing {required!r}:\n{msg}"
                    )
            for marker in abs_markers:
                if marker in msg:
                    raise AssertionError(
                        f"stale failure for {rel} leaked absolute path {marker!r}:\n{msg}"
                    )
            target.write_text(original, encoding="utf-8")

        (tree / "llms-full.txt").unlink()
        missing = _run_generator(tree, "--check")
        if missing.returncode == 0:
            raise AssertionError("--check unexpectedly passed when llms-full.txt missing")
        for required in ("llms-full.txt", REGEN_COMMAND):
            if required not in missing.stderr:
                raise AssertionError(
                    f"missing-file failure missing {required!r}:\n{missing.stderr}"
                )
        for marker in abs_markers:
            if marker in missing.stderr:
                raise AssertionError(
                    f"missing-file failure leaked absolute path {marker!r}:\n{missing.stderr}"
                )


def check_repo_outputs(root: Path) -> int:
    result = _run_generator(root, "--check")
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return 1
    return 0


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_llms_txt.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    test_deterministic_generation(root)
    test_stale_failure_is_actionable(root)
    return check_repo_outputs(root)


if __name__ == "__main__":
    raise SystemExit(main())
