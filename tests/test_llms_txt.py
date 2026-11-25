#!/usr/bin/env python3
"""Verify llms.txt / llms-full.txt freshness, determinism, and actionable diffs."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import importlib.util
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


def _load_generator_module(tree: Path):
    spec = importlib.util.spec_from_file_location(
        "generate_llms_txt", tree / "scripts" / "generate_llms_txt.py"
    )
    if spec is None or spec.loader is None:
        raise AssertionError("could not import generator script")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_org_plaintext_wrapped_links_and_fix_labels(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = (
            "* Heading\n"
            "*Fix:* keep this as prose.\n"
            "See [[file:../reference/pwscale.org][Plane-wave\n"
            "scaffold reference]] today.\n"
        )
        rendered = generator.org_plaintext(sample)
        if "# Fix:* keep this as prose." in rendered:
            raise AssertionError(f"*Fix:* was misparsed as heading:\n{rendered}")
        if "[[" in rendered:
            raise AssertionError(f"wrapped link was not rewritten:\n{rendered}")
        if "Plane-wave scaffold reference (file:../reference/pwscale.org)" not in rendered:
            raise AssertionError(f"wrapped link text lost after rewrite:\n{rendered}")


def test_org_plaintext_keeps_cog_blocks(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        cog_open = (
            "# [[[cog\n"
            "# import re\n"
            '# cog.out(re.sub(r"^v\\S+", "v1", cog.previous))\n'
            "# ]]]\n"
        )
        sample = (
            cog_open
            + "#+begin_src text\n"
            "v0 ok\n"
            "#+end_src\n"
            "# [[[end]]]\n"
            "See [[file:a.org][wrapped\n"
            "link]] after the block.\n"
        )
        normalized = generator.normalize_wrapped_org_links(sample)
        if not normalized.startswith(cog_open):
            raise AssertionError(f"Cog opening block was rewritten:\n{normalized}")
        if "# [[[end]]]\n" not in normalized:
            raise AssertionError(f"Cog end marker was rewritten:\n{normalized}")
        rendered = generator.org_plaintext(sample)
        for marker in ("# [[[cog", "# import re", "# ]]]", "# [[[end]]]"):
            if f"\n{marker}\n" not in f"\n{rendered}\n":
                raise AssertionError(f"Cog line {marker!r} lost:\n{rendered}")
        if "wrapped link (file:a.org)" not in rendered:
            raise AssertionError(f"link after Cog block not rewritten:\n{rendered}")


def test_normalize_joins_only_well_formed_links(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        cases = {
            "[[file:a.org][two\n  lines]]": "[[file:a.org][two lines]]",
            "[[file:a.org][across\n\nparagraphs]]": None,
            "[[file:a.org\n]]": None,
            "# [[[cog\n# x = 1\n# ]]]\n": None,
            "x [[[a][b\nc]]] y": None,
        }
        for sample, want in cases.items():
            got = generator.normalize_wrapped_org_links(sample)
            if got != (sample if want is None else want):
                raise AssertionError(f"normalize changed {sample!r} to {got!r}")


TESTS = (
    test_deterministic_generation,
    test_stale_failure_is_actionable,
    test_org_plaintext_wrapped_links_and_fix_labels,
    test_org_plaintext_keeps_cog_blocks,
    test_normalize_joins_only_well_formed_links,
)


def check_all_tests_listed() -> None:
    defined = {
        name
        for name, value in globals().items()
        if name.startswith("test_") and callable(value)
    }
    listed = {test.__name__ for test in TESTS}
    if defined != listed:
        raise AssertionError(
            f"test functions missing from TESTS: {sorted(defined - listed)}"
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
    check_all_tests_listed()
    for test in TESTS:
        test(root)
    return check_repo_outputs(root)


if __name__ == "__main__":
    raise SystemExit(main())
