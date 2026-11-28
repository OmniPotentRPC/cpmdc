#!/usr/bin/env python3
"""Install odd directory layouts and check pkg-config paths stay in the prefix.

`cpmdc.pc`, `cpmdc-restart.pc`, `cpmdc-damp.pc`, and `libcp2k.pc` derive `prefix` and
`includedir` from `${pcfiledir}`. `-Dlibdir=.`, a nested `libdir`, and an
absolute `includedir` under the prefix must resolve inside that prefix, and
the installed file must not record the absolute prefix path.
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

PC_NAMES = ("cpmdc.pc", "cpmdc-restart.pc", "cpmdc-damp.pc", "libcp2k.pc")
_VAR = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")


def run(cmd, **kw):
    print("+", shlex.join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def inside(path: Path, root: Path) -> bool:
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def expand_pc(path: Path) -> dict[str, str]:
    values = {"pcfiledir": str(path.parent)}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        values[key.strip()] = value.strip()

    def subst(text: str) -> str:
        for _ in range(12):
            def repl(match: re.Match[str]) -> str:
                name = match.group(1)
                if name not in values:
                    raise SystemExit(f"{path}: unknown ${{{name}}}")
                return values[name]

            nxt = _VAR.sub(repl, text)
            if nxt == text:
                return text
            text = nxt
        raise SystemExit(f"{path}: variable expansion did not finish")

    return {key: subst(value) for key, value in values.items()}


def find_pcs(prefix: Path) -> dict[str, Path]:
    found: dict[str, Path] = {}
    for path in prefix.rglob("*.pc"):
        if path.name in PC_NAMES and path.name not in found:
            found[path.name] = path
    missing = [name for name in PC_NAMES if name not in found]
    if missing:
        raise SystemExit(f"missing {missing} under {prefix}")
    return found


def check_pc(path: Path, prefix: Path, includedir: Path) -> None:
    raw = path.read_text(encoding="utf-8")
    prefix_s = str(prefix)
    if prefix_s and prefix_s in raw:
        raise SystemExit(f"{path}: records the absolute prefix {prefix_s}")
    if "${pcfiledir}" not in raw:
        raise SystemExit(f"{path}: does not derive paths from ${{pcfiledir}}")
    expanded = expand_pc(path)
    for key in ("prefix", "includedir", "libdir"):
        if key not in expanded:
            raise SystemExit(f"{path}: missing {key}")
        resolved = Path(expanded[key]).resolve()
        if not inside(resolved, prefix):
            raise SystemExit(f"{path}: {key} {resolved} is outside {prefix}")
    got_prefix = Path(expanded["prefix"]).resolve()
    if got_prefix != prefix.resolve():
        raise SystemExit(f"{path}: prefix {got_prefix} is not {prefix}")
    got_inc = Path(expanded["includedir"]).resolve()
    if got_inc != includedir.resolve():
        raise SystemExit(f"{path}: includedir {got_inc} is not {includedir}")
    if not got_inc.is_dir():
        raise SystemExit(f"{path}: includedir {got_inc} is not a directory")
    print(f"pc ok: {path} prefix={got_prefix} includedir={got_inc}", flush=True)


def layouts(work: Path) -> list[tuple[str, list[str], str]]:
    absolute = (work / "prefixes" / "absolute-includedir").resolve()
    return [
        ("libdir-dot", ["-Dlibdir=."], "include"),
        ("nested-libdir", ["-Dlibdir=lib/cpmdc-nested"], "include"),
        (
            "absolute-includedir",
            ["-Dlibdir=lib", f"-Dincludedir={absolute / 'include'}"],
            "include",
        ),
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--build-root", required=True)
    parser.add_argument("--meson", default="meson")
    parser.add_argument("--cc", default=os.environ.get("CC", "cc"))
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    args = parser.parse_args()

    source = Path(args.source_dir).resolve()
    work = Path(args.build_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    inner = work / "build"
    setup_env = dict(os.environ)
    setup_env["CC"] = args.cc
    setup_env["FC"] = args.fc

    first = True
    for name, extra, include_suffix in layouts(work):
        prefix = (work / "prefixes" / name).resolve()
        if prefix.exists():
            shutil.rmtree(prefix)
        cmd = [
            args.meson,
            "setup",
            inner,
            source,
            "--prefix",
            prefix,
            "-Dwith_tests=false",
            *extra,
        ]
        if not first:
            cmd.append("--reconfigure")
        first = False
        run(cmd, env=setup_env)
        run([args.meson, "install", "-C", inner], env=setup_env)
        includedir = (prefix / include_suffix).resolve()
        cp2k_includedir = includedir / "cpmdc" / "gromacs"
        for pc_name, pc_path in find_pcs(prefix).items():
            expect = cp2k_includedir if pc_name == "libcp2k.pc" else includedir
            check_pc(pc_path, prefix, expect)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
