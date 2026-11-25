#!/usr/bin/env python3
"""tools/check_snapper.sh: --check enforces classes, --write stays in the tree."""

from __future__ import annotations

import importlib.util
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Stands in for snapper-fmt: --check reports files that contain REWRAP as
# would_reformat, a FUSED, LONG or WRAP line as that diagnostic, and RENDER
# as a render-check warning on stderr; stdin mode replaces REWRAP.
FAKE_SNAPPER = """#!{python}
import json, sys
args = sys.argv[1:]
if "--check" in args:
    files = args[args.index("--") + 1:]
    report = []
    for name in files:
        with open(name, encoding="utf-8") as handle:
            text = handle.read()
        diags = [
            {{"kind": kind.lower(), "line": n, "excerpt": line}}
            for n, line in enumerate(text.splitlines(), 1)
            for kind in ("FUSED", "LONG", "WRAP")
            if kind in line
        ]
        if "RENDER" in text:
            sys.stderr.write("snapper: render check failed; left lines 1-1 unchanged\\n")
        report.append({{"file": name, "would_reformat": "REWRAP" in text,
                        "diagnostics": diags}})
    if "json" in args:
        print(json.dumps(report))
    sys.exit(1 if any(r["would_reformat"] for r in report) else 0)
sys.stdout.write(sys.stdin.read().replace("REWRAP", "wrapped"))
"""

NEEDS_WRAP = "Keep this. REWRAP\n"
WRAPPED = "Keep this. wrapped\n"
OUTSIDE = "Outside the tree. REWRAP\n"


def _load_helper(tree: Path):
    spec = importlib.util.spec_from_file_location(
        "snapper_write", tree / "tools" / "snapper_write.py"
    )
    if spec is None or spec.loader is None:
        raise AssertionError("could not import tools/snapper_write.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _make_tree(source: Path, base: Path) -> tuple[Path, Path]:
    tree = base / "repo"
    outside = base / "outside"
    (tree / "tools").mkdir(parents=True)
    (tree / "docs").mkdir()
    outside.mkdir()
    for rel in (
        "tools/check_snapper.sh",
        "tools/snapper_gate.py",
        "tools/snapper_write.py",
        ".snapperrc.toml",
    ):
        shutil.copy2(source / rel, tree / rel)
    (tree / "docs" / "page.md").write_text(NEEDS_WRAP, encoding="utf-8")
    (tree / "docs" / "page.md").chmod(0o640)
    (tree / "docs" / "other.md").write_text(NEEDS_WRAP, encoding="utf-8")
    (outside / "secret.md").write_text(OUTSIDE, encoding="utf-8")
    (tree / "docs" / "link.md").symlink_to(outside / "secret.md")
    (tree / "docs" / "linkdir").symlink_to(outside, target_is_directory=True)
    subprocess.run(["git", "init", "-q", str(tree)], check=True)
    return tree, outside


def _fake_path(base: Path) -> dict[str, str]:
    bindir = base / "bin"
    bindir.mkdir()
    fake = bindir / "snapper-fmt"
    fake.write_text(FAKE_SNAPPER.format(python=sys.executable), encoding="utf-8")
    fake.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{bindir}{os.pathsep}{env.get('PATH', '')}"
    env["PYTHON"] = sys.executable
    return env


def _run_script(
    tree: Path, env: dict[str, str], mode: str, *paths: str
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(tree / "tools" / "check_snapper.sh"), mode, *paths],
        cwd=tree,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def _write(tree: Path, env: dict[str, str], *paths: str) -> subprocess.CompletedProcess[str]:
    return _run_script(tree, env, "--write", *paths)


def _check(tree: Path, env: dict[str, str], *paths: str) -> subprocess.CompletedProcess[str]:
    return _run_script(tree, env, "--check", *paths)


def _expect(cond: bool, message: str, result: subprocess.CompletedProcess[str] | None = None) -> None:
    if not cond:
        detail = f"\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}" if result else ""
        raise AssertionError(message + detail)


def test_write_refuses_links_and_outside_paths(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        tree, outside = _make_tree(root, base)
        env = _fake_path(base)
        for bad in (
            "docs/link.md",
            "docs/linkdir/secret.md",
            str(outside / "secret.md"),
            "../outside/secret.md",
        ):
            result = _write(tree, env, "docs/other.md", bad)
            _expect(result.returncode != 0, f"--write accepted {bad}", result)
            _expect("refusing to write" in result.stderr, f"no refusal for {bad}", result)
            _expect(
                (outside / "secret.md").read_text(encoding="utf-8") == OUTSIDE,
                f"--write {bad} changed a file outside the tree",
                result,
            )
            _expect(
                (tree / "docs" / "other.md").read_text(encoding="utf-8") == NEEDS_WRAP,
                f"--write {bad} still rewrote the valid path",
                result,
            )


def test_write_rewraps_in_place(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        tree, _ = _make_tree(root, base)
        env = _fake_path(base)
        page = tree / "docs" / "page.md"
        result = _write(tree, env, "docs/page.md")
        _expect(result.returncode == 0, "--write failed on a regular file", result)
        _expect(page.read_text(encoding="utf-8") == WRAPPED, "file was not rewrapped", result)
        _expect(not page.is_symlink(), "rewrapped file became a symlink")
        _expect(stat.S_IMODE(page.stat().st_mode) == 0o640, "file mode was not kept")
        leftovers = [p.name for p in page.parent.iterdir() if ".snapper-" in p.name]
        _expect(not leftovers, f"temporary files left behind: {leftovers}")


def test_replace_refuses_swapped_path(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        tree, outside = _make_tree(root, base)
        helper = _load_helper(tree)
        root_fd = os.open(tree, os.O_RDONLY | os.O_DIRECTORY)
        cwd = os.getcwd()
        try:
            os.chdir(tree)
            target = helper.open_target(root_fd, tree, "docs/page.md")
            page = tree / "docs" / "page.md"
            page.unlink()
            page.symlink_to(outside / "secret.md")
            try:
                helper.replace_target(target, WRAPPED.encode())
            except helper.RefusedPath:
                pass
            else:
                raise AssertionError("replace_target wrote over a swapped path")
            os.close(target.dir_fd)
        finally:
            os.chdir(cwd)
            os.close(root_fd)
        _expect(
            (outside / "secret.md").read_text(encoding="utf-8") == OUTSIDE,
            "swapped symlink target was written",
        )
        _expect(page.is_symlink(), "swapped symlink was replaced")
        leftovers = [p.name for p in page.parent.iterdir() if ".snapper-" in p.name]
        _expect(not leftovers, f"temporary files left behind: {leftovers}")


def test_check_fails_on_enforced_classes(root: Path = ROOT) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        tree, _ = _make_tree(root, base)
        env = _fake_path(base)
        docs = tree / "docs"
        (docs / "advisory.md").write_text("A LONG line.\nA WRAP line.\n", encoding="utf-8")
        (docs / "fused.md").write_text("Two FUSED. Sentences.\n", encoding="utf-8")
        (docs / "render.md").write_text("RENDER me.\n", encoding="utf-8")

        result = _check(tree, env, "docs/advisory.md")
        _expect(result.returncode == 0, "advisory classes failed the check", result)
        _expect("docs/advisory.md: long 1, wrap 1" in result.stdout,
                "advisory diagnostics were not reported per file", result)
        _expect("long 1" in result.stdout and "wrap 1" in result.stdout
                and "snapper OK" in result.stdout, "summary lacks advisory counts", result)

        for name, kind in (("fused.md", "fused"), ("render.md", "render"), ("page.md", "reformat")):
            result = _check(tree, env, f"docs/{name}")
            _expect(result.returncode != 0, f"{kind} did not fail the check", result)
            _expect(f"docs/{name}" in result.stdout and f"{kind}:" in result.stdout,
                    f"{kind} failure was not named", result)
            _expect("snapper FAIL" in result.stdout, f"{kind} summary not FAIL", result)

        strict = dict(env, SNAPPER_FAIL_ON="long")
        result = _check(tree, strict, "docs/advisory.md")
        _expect(result.returncode != 0, "SNAPPER_FAIL_ON=long did not fail", result)
        _expect("docs/advisory.md:1: long:" in result.stdout, "long line not named", result)

        summary = base / "summary.md"
        result = _check(tree, dict(env, GITHUB_STEP_SUMMARY=str(summary)), "docs/fused.md")
        _expect("| fused | 1 | yes |" in summary.read_text(encoding="utf-8"),
                "step summary lacks the fused count", result)


def test_real_snapper_write_and_fused(root: Path = ROOT) -> None:
    snapper = shutil.which("snapper-fmt") or shutil.which("snapper")
    if snapper is None:
        print("snapper-fmt not on PATH; real snapper comparison not run")
        return
    text = (
        "One sentence here. Another sentence that runs long enough to need a "
        "wrap past eighty columns.\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        tree, _ = _make_tree(root, base)
        (tree / "docs" / "page.md").write_text(text, encoding="utf-8")
        reference = base / "reference.md"
        reference.write_text(text, encoding="utf-8")
        subprocess.run(
            [snapper, "--native", "--in-place", "--max-width", "80",
             "--config", str(tree / ".snapperrc.toml"), str(reference)],
            check=True,
            capture_output=True,
        )
        env = dict(os.environ, PYTHON=sys.executable)
        result = _write(tree, env, "docs/page.md")
        _expect(result.returncode == 0, "--write failed with real snapper", result)
        got = (tree / "docs" / "page.md").read_text(encoding="utf-8")
        want = reference.read_text(encoding="utf-8")
        _expect(got == want, f"--write differs from --in-place:\n{got!r}\n{want!r}")
        _expect(got != text, "real snapper left the sample unchanged")

        # snapper-fmt 0.11.10 keeps this paragraph as is and flags two lines fused.
        fused = tree / "docs" / "fused.org"
        fused.write_text(
            "".join((
                'GROMACS can integrate on cpmdc energies and forces in process, for a QM group\n'
                'that is part of the system or all of it.\n'
                "~libcpmdc_gromacs~ implements the libcp2k C API that GROMACS' CP2K QM/MM interface calls, and answers each call with one ~cpmdc_session_calculate_au()~\n"
                'of the [[file:../reference/c-abi.org][C ABI]]. GROMACS links it in place of\n'
                'CP2K, so ~mdrun~ needs no extra process and no socket. Energies and forces need\n'
                "no change to GROMACS' source; the QM virial, and so the pressure and NPT, needs\n"
                "a short patch to GROMACS' QM/MM module (see [[*Pressure and NPT][Pressure and NPT]]).\n"
            )),
            encoding="utf-8",
        )
        result = _check(tree, env, "docs/fused.org")
        _expect(result.returncode != 0, "real snapper fused line passed the check", result)
        _expect("docs/fused.org:" in result.stdout and ": fused:" in result.stdout,
                "real snapper fused line was not named", result)


TESTS = (
    test_write_refuses_links_and_outside_paths,
    test_write_rewraps_in_place,
    test_replace_refuses_swapped_path,
    test_check_fails_on_enforced_classes,
    test_real_snapper_write_and_fused,
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


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_check_snapper.py SOURCE_ROOT", file=sys.stderr)
        return 2
    if shutil.which("git") is None:
        print("git not on PATH", file=sys.stderr)
        return 77
    root = Path(sys.argv[1]).resolve()
    check_all_tests_listed()
    for test in TESTS:
        test(root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
