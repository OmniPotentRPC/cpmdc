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


_MATH_SPAN = re.compile(
    r"\$\$.*?\$\$|\\\(.*?\\\)|\\\[.*?\\\]|\\begin\{([A-Za-z*]+)\}.*?\\end\{\1\}",
    re.DOTALL,
)


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
        sample = (
            r"Keep \(\epsilon_{rf}=1\) next to =one= and \(~x~\)."
            "\n"
            r"Also $\epsilon_{rf}=1,~x~$ and \[E = 1,~x~\]."
            "\n"
            r"\begin{equation}\epsilon_{rf}=1,~x~\end{equation}"
            "\n"
            "Outside ~x~ and shape ~(N, 3)~."
            "\n"
            "See ~[[file:../howto/mpi.org][text]]~."
            "\n"
            "E0 = 85 and sigma = 4 stays."
            "\n"
            "#+begin_src awk"
            "\n"
            r"NF==6 && $2 ~ /E[-+]/ && $1 ~ /^[0-9]+$/"
            "\n"
            "#+end_src"
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


def test_org_plaintext_keeps_parenthesized_code(root: Path = ROOT) -> None:
    """A code span whose body is a parenthesized group is copied through.

    The body may contain commas and spaces. Several such spans may sit on
    one line, separated by a comma. Deleting the whole span leaves the
    surrounding punctuation and a blank where the vector or index pair was.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = (
            "row: ~(lx, 0, 0)~, ~(xy, ly, 0)~, ~(xz, yz, lz)~ end\n"
            "pairs ~(1, 3)~, ~(2, 6)~,\n"
            "and ~(5, 7)~ as the symmetrized terms\n"
            "shape ~(N, 3)~ and empty ~()~ beside ~(h[0], 0, 0)~.\n"
            r"Math \(-dE/d(\mathrm{strain})\) stays with ~(xlo, ylo, zlo)~."
        )
        rendered = generator.org_plaintext(sample)
        required = (
            "`(lx, 0, 0)`, `(xy, ly, 0)`, `(xz, yz, lz)`",
            "`(1, 3)`, `(2, 6)`",
            "`(5, 7)`",
            "`(N, 3)`",
            "`()`",
            "`(h[0], 0, 0)`",
            "`(xlo, ylo, zlo)`",
            r"\(-dE/d(\mathrm{strain})\)",
        )
        for piece in required:
            if piece not in rendered:
                raise AssertionError(f"missing {piece!r} in:\n{rendered}")
        if "row: , , " in rendered or "pairs , ," in rendered:
            raise AssertionError(f"parenthesized code was removed:\n{rendered}")
        page = (
            root / "docs" / "orgmode" / "reference" / "lammps.org"
        ).read_text(encoding="utf-8")
        page_text = generator.org_plaintext(page)
        page_bits = (
            "one lattice vector per row: `(lx, 0, 0)`, `(xy, ly, 0)`, `(xz, yz, lz)`",
            "and the pairs `(1, 3)`, `(2, 6)`,",
            "and `(5, 7)` as the symmetrized",
        )
        for piece in page_bits:
            if piece not in page_text:
                raise AssertionError(f"lammps page lost {piece!r} in:\n{page_text}")
        full = (root / "llms-full.txt").read_text(encoding="utf-8")
        for piece in page_bits:
            if piece not in full:
                raise AssertionError(f"llms-full.txt lost {piece!r}")
        blanks = (
            "one lattice vector per row: , , ",
            "and the pairs , ,",
            "and  as the symmetrized",
        )
        for piece in blanks:
            if piece in full:
                raise AssertionError(f"llms-full.txt still has blank span {piece!r}")


def test_org_plaintext_keeps_multiline_math_blocks(root: Path = ROOT) -> None:
    """Math blocks that span lines are copied unchanged.

    A line-wise walk turns a code marker inside an open dollar block into
    a backtick span. Multiline dollar math, bracket math, and LaTeX
    environments keep the marker, including spaces on a math line.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = (
            "before ~y~\n"
            "$$\n"
            "~x~  \n"
            "$$\n"
            "$$~x~$$\n"
            "\\[\n"
            "~x~\n"
            "\\]\n"
            "\\begin{equation}\n"
            "~x~\n"
            "\\end{equation}\n"
            "\\begin{align*}\n"
            "~x~\n"
            "\\end{align*}\n"
            "after ~z~\n"
        )
        rendered = generator.org_plaintext(sample)
        required = (
            "$$\n~x~  \n$$",
            "$$~x~$$",
            "\\[\n~x~\n\\]",
            "\\begin{equation}\n~x~\n\\end{equation}",
            "\\begin{align*}\n~x~\n\\end{align*}",
            "`y`",
            "`z`",
        )
        for piece in required:
            if piece not in rendered:
                raise AssertionError(f"missing {piece!r} in:\n{rendered}")
        if "`x`" in rendered:
            raise AssertionError(f"math code marker became a span:\n{rendered}")


def test_org_plaintext_keeps_literal_trailing_whitespace(root: Path = ROOT) -> None:
    """Source and example blocks keep trailing spaces and tabs.

    Prose still drops trailing spaces. The block markers themselves are
    not part of the plain text.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = (
            "before\n"
            "#+begin_src text\n"
            "keep two  \n"
            "tab\t\n"
            "   \n"
            "#+end_src\n"
            "#+begin_example\n"
            "ex   \n"
            "#+end_example\n"
            "#+BEGIN_SRC python\n"
            "upper  \n"
            "#+END_SRC\n"
            "prose trailing   \n"
        )
        rendered = generator.org_plaintext(sample)
        required = (
            "keep two  \n",
            "tab\t\n",
            "   \n",
            "ex   \n",
            "upper  \n",
        )
        for piece in required:
            if piece not in rendered:
                raise AssertionError(f"missing {piece!r} in:\n{rendered}")
        if "prose trailing   " in rendered:
            raise AssertionError(f"prose kept trailing spaces:\n{rendered}")
        if "prose trailing\n" not in rendered:
            raise AssertionError(f"prose line missing:\n{rendered}")
        if "#+begin_" in rendered or "#+end_" in rendered or "#+BEGIN_" in rendered:
            raise AssertionError(f"block marker leaked:\n{rendered}")


def test_org_plaintext_keeps_literal_whitespace_when_block_ends_page(
    root: Path = ROOT,
) -> None:
    """A literal block that ends the page keeps spaces and a tab.

    Blank lines after the block are dropped. The last literal line stays
    byte for byte, including when nothing follows it.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        sample = (
            "intro\n"
            "#+begin_src text\n"
            "last  \t\n"
            "#+end_src\n"
            "\n"
            "\n"
        )
        rendered = generator.org_plaintext(sample)
        expected = "intro\nlast  \t\n"
        if rendered != expected:
            raise AssertionError(f"page end changed: {rendered!r}")


def test_markup_corruption_flags_dollar_math_and_environments(
    root: Path = ROOT,
) -> None:
    """Dollar math and LaTeX environments count as math spans.

    A backtick inside those spans is corruption, the same as inside
    backslash math. A clean span, and a code span outside it, is not.
    """
    del root
    flagged = (
        "$$`x`$$",
        "$$\n`x`\n$$",
        "\\begin{equation}`x`\\end{equation}",
        "\\begin{align*}\n`x`\n\\end{align*}",
    )
    for sample in flagged:
        hits = markup_corruption(sample)
        if not any(hit.startswith("backtick inside math") for hit in hits):
            raise AssertionError(
                f"math corruption was not flagged: {sample!r} -> {hits}"
            )
    clean = (
        "$$~x~$$",
        "$$\n~x~\n$$",
        "\\begin{equation}~x~\\end{equation}",
        "\\(~x~\\)",
        "\\[~x~\\]",
        "prose `code`",
    )
    for sample in clean:
        hits = markup_corruption(sample)
        if hits:
            raise AssertionError(f"clean sample flagged: {sample!r} -> {hits}")


def test_emphasis_borders_match_org(root: Path = ROOT) -> None:
    """Prematch and postmatch match Org's emphasis components.

    A slash is neither, so a slash beside a marker stays inside that one
    object. Two command names with a space on each side of the slash are
    two objects. A closing bracket is not a postmatch.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tree = Path(tmp) / "repo"
        _copy_generator_inputs(root, tree)
        generator = _load_generator_module(tree)
        pre_ok = ("-", " ", "\t", "(", "'", '"', "{")
        pre_no = ("/", "]", ")", "}", "[", "\\", ".", ",", "a")
        post_ok = (
            "-",
            " ",
            "\t",
            ".",
            ",",
            ":",
            "!",
            "?",
            ";",
            "'",
            '"',
            ")",
            "}",
            "\\",
            "[",
        )
        post_no = ("/", "]", "(", "{", "a", "1")
        for char in pre_ok:
            if not generator._is_emphasis_pre(char):
                raise AssertionError(f"{char!r} should be a prematch")
        for char in pre_no:
            if generator._is_emphasis_pre(char):
                raise AssertionError(f"{char!r} is not an Org prematch")
        for char in post_ok:
            if not generator._is_emphasis_post(char):
                raise AssertionError(f"{char!r} should be a postmatch")
        for char in post_no:
            if generator._is_emphasis_post(char):
                raise AssertionError(f"{char!r} is not an Org postmatch")
        one = generator.org_plaintext("next ~build/~ and ~/path/to/x~.")
        if "`build/`" not in one or "`/path/to/x`" not in one:
            raise AssertionError(f"slash inside a span was split:\n{one}")
        glued = generator.org_plaintext("pipelined ~POSDATA~/~GETFORCE~;")
        if glued != "pipelined `POSDATA~/~GETFORCE`;\n":
            raise AssertionError(f"slash between markers was not one object:\n{glued}")
        spaced = generator.org_plaintext("pipelined ~POSDATA~ / ~GETFORCE~;")
        if "`POSDATA` / `GETFORCE`" not in spaced:
            raise AssertionError(f"spaced names were not two spans:\n{spaced}")
        page = (
            root / "docs" / "orgmode" / "howto" / "socket-drivers.org"
        ).read_text(encoding="utf-8")
        page_text = generator.org_plaintext(page)
        if "`POSDATA` / `GETFORCE`" not in page_text:
            raise AssertionError(f"socket-drivers page lost the names:\n{page_text}")
        if "POSDATA=/GETFORCE" in page_text or "`POSDATA`/`GETFORCE`" in page_text:
            raise AssertionError(f"socket-drivers page still joins the names:\n{page_text}")
        full = (root / "llms-full.txt").read_text(encoding="utf-8")
        protocol = [
            line
            for line in full.splitlines()
            if "pipelined" in line and "GETFORCE" in line
        ]
        if len(protocol) != 1:
            raise AssertionError(f"expected one protocol line, got {protocol!r}")
        line = protocol[0]
        if "`POSDATA` / `GETFORCE`" not in line:
            raise AssertionError(f"llms-full.txt lost the spaced names:\n{line}")
        if "~" in line or "POSDATA=/" in line or "`POSDATA`/`GETFORCE`" in line:
            raise AssertionError(f"stray marker on the protocol line:\n{line}")


TESTS = (
    test_deterministic_generation,
    test_stale_failure_is_actionable,
    test_org_plaintext_wrapped_links_and_fix_labels,
    test_org_plaintext_keeps_cog_blocks,
    test_normalize_joins_only_well_formed_links,
    test_markup_corruption_flags_synthetic_input,
    test_committed_llms_have_no_markup_corruption,
    test_org_plaintext_preserves_math_spans,
    test_org_plaintext_keeps_parenthesized_code,
    test_org_plaintext_keeps_multiline_math_blocks,
    test_org_plaintext_keeps_literal_trailing_whitespace,
    test_org_plaintext_keeps_literal_whitespace_when_block_ends_page,
    test_markup_corruption_flags_dollar_math_and_environments,
    test_emphasis_borders_match_org,
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
