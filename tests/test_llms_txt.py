#!/usr/bin/env python3
"""Verify llms.txt / llms-full.txt freshness, determinism, and actionable diffs."""

from __future__ import annotations

import re
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


_MATH_SPAN = re.compile(r"\\\(.*?\\\)|\\\[.*?\\\]", re.DOTALL)


def markup_corruption(text: str) -> list[str]:
    """Backticks inside math spans, and lines with an odd backtick count."""
    hits: list[str] = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if line.count("`") % 2 == 1:
            hits.append(f"unbalanced backticks line {lineno}: {line}")
    for match in _MATH_SPAN.finditer(text):
        if "`" not in match.group(0):
            continue
        lineno = text.count("\n", 0, match.start()) + 1
        hits.append(f"backtick inside math line {lineno}: {match.group(0)}")
    return hits


def test_markup_corruption_flags_synthetic_input(root: Path = ROOT) -> None:
    del root
    math_hits = markup_corruption(r"pre \(\epsilon_{rf}`1\) post")
    if not any(hit.startswith("backtick inside math") for hit in math_hits):
        raise AssertionError(f"synthetic math corruption was not flagged: {math_hits}")
    odd_hits = markup_corruption("one ` tick")
    if not any(hit.startswith("unbalanced backticks") for hit in odd_hits):
        raise AssertionError(f"synthetic odd backtick was not flagged: {odd_hits}")
    clean = r"value \(\epsilon_{rf}=1\) and \(~x~\) plus `code`"
    clean_hits = markup_corruption(clean)
    if clean_hits:
        raise AssertionError(f"clean sample flagged: {clean_hits}")


def test_committed_llms_have_no_markup_corruption(root: Path = ROOT) -> None:
    for name in ("llms.txt", "llms-full.txt"):
        text = (root / name).read_text(encoding="utf-8")
        hits = markup_corruption(text)
        if hits:
            preview = "\n".join(hits[:20])
            raise AssertionError(f"{name} markup corruption:\n{preview}")


def test_org_plaintext_preserves_math_spans(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = "\n".join(
            [
                r"Keep \(\epsilon_{rf}=1\) next to =one= and \(~x~\).",
                r"Also $\epsilon_{rf}=1,~x~$ and \[E = 1,~x~\].",
                r"\begin{equation}\epsilon_{rf}=1,~x~\end{equation}",
                "Outside ~x~ and shape ~(N, 3)~.",
                "See ~[[file:../howto/mpi.org][text]]~.",
                "E0 = 85 and sigma = 4 stays.",
                "#+begin_src awk",
                r"NF==6 && $2 ~ /E[-+]/ && $1 ~ /^[0-9]+$/",
                "#+end_src",
            ]
        )
        rendered = generator.org_plaintext(sample)
        required = [
            r"\(\epsilon_{rf}=1\)",
            r"\(~x~\)",
            r"$\epsilon_{rf}=1,~x~$",
            r"\[E = 1,~x~\]",
            r"\begin{equation}\epsilon_{rf}=1,~x~\end{equation}",
            "`one`",
            "`x`",
            "`(N, 3)`",
            "`[[file:../howto/mpi.org][text]]`",
            "E0 = 85 and sigma = 4 stays.",
            r"NF==6 && $2 ~ /E[-+]/ && $1 ~ /^[0-9]+$/",
        ]
        for piece in required:
            if piece not in rendered:
                raise AssertionError(f"missing {piece!r} in:\n{rendered}")
        hits = markup_corruption(rendered)
        if hits:
            raise AssertionError(f"rendered math was corrupted:\n{hits}\n{rendered}")


TESTS = (
    test_deterministic_generation,
    test_stale_failure_is_actionable,
    test_org_plaintext_wrapped_links_and_fix_labels,
    test_org_plaintext_keeps_cog_blocks,
    test_normalize_joins_only_well_formed_links,
    test_markup_corruption_flags_synthetic_input,
    test_committed_llms_have_no_markup_corruption,
    test_org_plaintext_preserves_math_spans,
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
