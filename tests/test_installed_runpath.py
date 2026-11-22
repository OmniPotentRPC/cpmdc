#!/usr/bin/env python3
"""Install several directory layouts and check the executable run path.

`cpmdc-ipi-c` must record one `$ORIGIN` entry on Linux, or `@loader_path`
on macOS, whose normalised path resolves to the directory that holds
`libcpmdc`. The layouts are the default directories, `-Dbindir=.`,
absolute `bindir` and `libdir` under the prefix, a nested `libdir`, and
paths that contain `..`.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ORIGIN_TOKENS = ("$ORIGIN", "${ORIGIN}", "@loader_path", "@executable_path")


def run(cmd, **kw):
    print("+", shlex.join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def runtime_path_entries(path: Path) -> list[str]:
    """`DT_RPATH` and `DT_RUNPATH`, or Mach-O `LC_RPATH` entries."""
    if sys.platform == "darwin":
        return _macho_rpaths(path)
    readelf = shutil.which("readelf")
    if readelf is not None:
        entries = _elf_readelf(path, readelf)
        if entries is not None:
            return entries
    patchelf = shutil.which("patchelf")
    if patchelf is not None:
        return _elf_patchelf(path, patchelf)
    raise SystemExit("cannot read run paths, readelf and patchelf are missing")


def _elf_readelf(path: Path, readelf: str) -> list[str] | None:
    out = subprocess.run(
        [readelf, "-d", str(path)],
        capture_output=True,
        text=True,
    )
    if out.returncode != 0:
        return None
    entries: list[str] = []
    for line in out.stdout.splitlines():
        if "(RPATH)" not in line and "(RUNPATH)" not in line:
            continue
        if "[" not in line or "]" not in line:
            continue
        payload = line.split("[", 1)[1].rsplit("]", 1)[0]
        entries.extend(part for part in payload.split(":") if part)
    return entries


def _elf_patchelf(path: Path, patchelf: str) -> list[str]:
    out = subprocess.run(
        [patchelf, "--print-rpath", str(path)],
        capture_output=True,
        text=True,
    )
    if out.returncode != 0:
        raise SystemExit(
            f"patchelf failed on {path}: {out.stderr.strip() or out.stdout.strip()}"
        )
    return [part for part in out.stdout.strip().split(":") if part]


def _macho_rpaths(path: Path) -> list[str]:
    otool = shutil.which("otool")
    if otool is None:
        raise SystemExit("cannot read run paths, otool is missing")
    out = subprocess.run([otool, "-l", str(path)], capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(
            f"otool failed on {path}: {out.stderr.strip() or out.stdout.strip()}"
        )
    entries: list[str] = []
    current = ""
    for line in out.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("cmd "):
            current = stripped.split(None, 1)[1]
            continue
        if current == "LC_RPATH" and stripped.startswith("path "):
            value = stripped.split(None, 1)[1]
            entries.append(value.split(" (", 1)[0].strip())
    return entries


def origin_suffix(entry: str) -> str | None:
    for token in ORIGIN_TOKENS:
        if entry == token:
            return "."
        prefix = token + "/"
        if entry.startswith(prefix):
            return entry[len(prefix) :]
    return None


def expand_origin(entry: str, origin: Path) -> Path | None:
    suffix = origin_suffix(entry)
    if suffix is None:
        return None
    if suffix == ".":
        return origin.resolve()
    return (origin / suffix).resolve()


def find_one(prefix: Path, name: str) -> Path:
    found = [path for path in prefix.rglob(name) if path.is_file() and not path.is_symlink()]
    if len(found) != 1:
        raise SystemExit(f"expected one {name} under {prefix}, found {found}")
    return found[0]


def find_libdir(prefix: Path) -> Path:
    found = []
    for path in prefix.rglob("*"):
        if not path.is_file() or path.is_symlink():
            continue
        name = path.name
        if name.startswith("libcpmdc_"):
            continue
        if name.startswith("libcpmdc.so") or (
            name.startswith("libcpmdc.") and "dylib" in name
        ):
            found.append(path)
    parents = {path.parent.resolve() for path in found}
    if len(parents) != 1:
        raise SystemExit(f"expected one libcpmdc directory under {prefix}, found {found}")
    return parents.pop()


def assert_runpath(prefix: Path, libdir_suffix: str | None) -> None:
    libdir = find_libdir(prefix)
    if libdir_suffix is not None and libdir != (prefix / libdir_suffix).resolve():
        raise SystemExit(f"{libdir} is not {prefix / libdir_suffix}")
    exe = find_one(prefix, "cpmdc-ipi-c")
    entries = runtime_path_entries(exe)
    origin_entries = []
    for entry in entries:
        suffix = origin_suffix(entry)
        if suffix is None:
            continue
        if os.path.normpath(suffix) != suffix:
            raise SystemExit(f"{exe}: run path {entry!r} is not normalised")
        expanded = expand_origin(entry, exe.parent)
        origin_entries.append((entry, expanded))
    if not origin_entries:
        raise SystemExit(f"{exe}: no origin-relative run path in {entries}")
    if not any(expanded == libdir for _, expanded in origin_entries):
        raise SystemExit(
            f"{exe}: {origin_entries} does not resolve to {libdir} (entries {entries})"
        )
    print(f"run path ok: {exe} -> {origin_entries[0][0]} == {libdir}", flush=True)


def layouts(work: Path) -> list[tuple[str, list[str], str | None, str]]:
    absolute = (work / "prefixes" / "absolute").resolve()
    return [
        ("default", [], None, "bin"),
        ("bindir-dot", ["-Dbindir=."], None, "."),
        (
            "absolute",
            [f"-Dbindir={absolute}/bin/", f"-Dlibdir={absolute}/lib/"],
            "lib",
            "bin",
        ),
        (
            "nested-libdir",
            ["-Dlibdir=lib/x86_64-linux-gnu"],
            "lib/x86_64-linux-gnu",
            "bin",
        ),
        (
            "parent-segments",
            [
                "-Dbindir=bin/../bin",
                "-Dlibdir=lib/x86_64-linux-gnu/../x86_64-linux-gnu",
            ],
            "lib/x86_64-linux-gnu",
            "bin",
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
    for name, extra, libdir_suffix, bindir in layouts(work):
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
        exe = find_one(prefix, "cpmdc-ipi-c")
        expected_dir = prefix if bindir == "." else (prefix / bindir).resolve()
        if exe.parent.resolve() != expected_dir:
            raise SystemExit(f"{exe} is not in {expected_dir}")
        assert_runpath(prefix, libdir_suffix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
