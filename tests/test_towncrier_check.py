#!/usr/bin/env python3
"""Verify scripts/towncrier_check.py runs non-skipping fragment and changelog checks."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "towncrier_check.py"


def _git(cwd: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout.strip()


def _commit(cwd: Path, message: str) -> None:
    _git(cwd, "add", "-A")
    _git(
        cwd,
        "-c",
        "user.name=Test User",
        "-c",
        "user.email=test@example.com",
        "-c",
        "commit.gpgsign=false",
        "-c",
        "core.hooksPath=/dev/null",
        "commit",
        "-m",
        message,
    )


def _init_repo(root: Path) -> None:
    _git(root, "init", "-b", "main")
    (root / "version.txt").write_text("0.7.0\n", encoding="utf-8")
    shutil.copy2(ROOT / "towncrier.toml", root / "towncrier.toml")
    (root / "CHANGELOG.md").write_text(
        "# Changelog\n\n<!-- towncrier release notes start -->\n",
        encoding="utf-8",
    )
    (root / "docs" / "newsfragments").mkdir(parents=True)
    (root / "docs" / "newsfragments" / "README.md").write_text(
        "# News fragments\n", encoding="utf-8"
    )
    (root / "src.txt").write_text("base\n", encoding="utf-8")
    _commit(root, "init main")
    _git(root, "update-ref", "refs/remotes/origin/main", "HEAD")

    # Advance develop with a CHANGELOG.md change relative to origin/main,
    # reproducing why default `towncrier check` (comparing origin/main) skips.
    _git(root, "checkout", "-b", "develop")
    (root / "CHANGELOG.md").write_text(
        "# Changelog\n\n<!-- towncrier release notes start -->\n\n## [0.7.0] - 2026-11-24\n",
        encoding="utf-8",
    )
    _commit(root, "release 0.7.0 on develop")
    _git(root, "update-ref", "refs/remotes/origin/develop", "HEAD")


def _run_check(root: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged_env = dict(os.environ)
    merged_env.pop("GITHUB_BASE_REF", None)
    if env:
        merged_env.update(env)
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--root", str(root)],
        cwd=root,
        env=merged_env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_feature_branch_missing_fragment_fails() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp) / "repo"
        repo.mkdir()
        _init_repo(repo)
        _git(repo, "checkout", "-b", "feat/no-frag", "develop")
        (repo / "src.txt").write_text("changed\n", encoding="utf-8")
        _commit(repo, "change without fragment")

        res = _run_check(repo)
        assert res.returncode != 0, res.stdout + res.stderr
        assert "Checks SKIPPED" not in res.stdout
        assert "No new newsfragments found" in res.stdout


def test_feature_branch_with_fragment_passes_despite_main_changelog_diff() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp) / "repo"
        repo.mkdir()
        _init_repo(repo)
        _git(repo, "checkout", "-b", "feat/with-frag", "develop")
        (repo / "src.txt").write_text("changed\n", encoding="utf-8")
        (repo / "docs" / "newsfragments" / "+example.fixed.md").write_text(
            "Fix example.\n", encoding="utf-8"
        )
        _commit(repo, "change with fragment")

        res = _run_check(repo)
        assert res.returncode == 0, res.stdout + res.stderr
        assert "Checks SKIPPED" not in res.stdout
        assert "Found:" in res.stdout


def test_invalid_fragment_name_fails_even_when_changelog_changes() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp) / "repo"
        repo.mkdir()
        _init_repo(repo)
        _git(repo, "checkout", "-b", "docs/bad-frag", "develop")
        (repo / "CHANGELOG.md").write_text(
            "# Changelog\n\n<!-- towncrier release notes start -->\n\n## [0.7.0] - 2026-11-24\n\n- entry\n",
            encoding="utf-8",
        )
        (repo / "docs" / "newsfragments" / "badname.md").write_text(
            "bad\n", encoding="utf-8"
        )
        _commit(repo, "invalid fragment name")

        res = _run_check(repo)
        assert res.returncode != 0
        assert "Invalid news fragment name" in res.stderr


def test_changelog_fold_fails_with_unconsumed_fragment_and_passes_when_folded() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp) / "repo"
        repo.mkdir()
        _init_repo(repo)
        (repo / "docs" / "newsfragments" / "+pending.fixed.md").write_text(
            "Pending fix.\n", encoding="utf-8"
        )
        _commit(repo, "pending fragment on develop")
        _git(repo, "update-ref", "refs/remotes/origin/develop", "HEAD")

        _git(repo, "checkout", "-b", "docs/fold", "develop")
        (repo / "CHANGELOG.md").write_text(
            "# Changelog\n\n<!-- towncrier release notes start -->\n\n## [0.7.0] - 2026-11-24\n\n- Pending fix.\n",
            encoding="utf-8",
        )
        _commit(repo, "partial fold leaving fragment")

        bad = _run_check(repo)
        assert bad.returncode != 0
        assert "pending newsfragments remain" in bad.stderr

        (repo / "docs" / "newsfragments" / "+pending.fixed.md").unlink()
        _commit(repo, "consume pending fragment")

        good = _run_check(repo, env={"GITHUB_BASE_REF": "develop"})
        assert good.returncode == 0, good.stdout + good.stderr
        assert "Checks SKIPPED" not in good.stdout
        assert "0 pending fragments remain" in good.stdout


def main() -> int:
    test_feature_branch_missing_fragment_fails()
    test_feature_branch_with_fragment_passes_despite_main_changelog_diff()
    test_invalid_fragment_name_fails_even_when_changelog_changes()
    test_changelog_fold_fails_with_unconsumed_fragment_and_passes_when_folded()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
