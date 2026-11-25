#!/usr/bin/env python3
"""Rewrap tracked prose with snapper without following links out of the tree.

Every path is opened from a handle on the repository root, one component
at a time with ``O_NOFOLLOW``, so symlinks and paths outside the root are
refused before anything is written.
The file is read through that handle, the formatted text goes to a
temporary file in the same directory handle, and a rename on that handle
replaces the original.
The path is never reopened by name.
"""

from __future__ import annotations

import argparse
import json
import os
import secrets
import stat
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


class RefusedPath(Exception):
    pass


@dataclass
class Target:
    rel: str
    dir_fd: int
    name: str
    st: os.stat_result
    content: bytes


def _relative_to_root(root: Path, path: str) -> str:
    full = os.path.normpath(os.path.join(os.getcwd(), path))
    try:
        rel = os.path.relpath(full, root)
    except ValueError as exc:
        raise RefusedPath(f"{path}: outside the repository") from exc
    if rel == os.curdir or rel == os.pardir or rel.startswith(os.pardir + os.sep):
        raise RefusedPath(f"{path}: outside the repository")
    return rel


def open_target(root_fd: int, root: Path, path: str) -> Target:
    rel = _relative_to_root(root, path)
    *dirs, name = rel.split(os.sep)
    dir_fd = os.dup(root_fd)
    try:
        for part in dirs:
            try:
                next_fd = os.open(
                    part,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
                    dir_fd=dir_fd,
                )
            except OSError as exc:
                raise RefusedPath(
                    f"{path}: {part} is a symlink or not a directory"
                ) from exc
            os.close(dir_fd)
            dir_fd = next_fd
        try:
            fd = os.open(
                name,
                os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC,
                dir_fd=dir_fd,
            )
        except OSError as exc:
            raise RefusedPath(f"{path}: symlink or unreadable file") from exc
        try:
            st = os.fstat(fd)
            if not stat.S_ISREG(st.st_mode):
                raise RefusedPath(f"{path}: not a regular file")
            chunks = []
            while chunk := os.read(fd, 1 << 16):
                chunks.append(chunk)
        finally:
            os.close(fd)
    except BaseException:
        os.close(dir_fd)
        raise
    return Target(rel, dir_fd, name, st, b"".join(chunks))


def replace_target(target: Target, content: bytes) -> None:
    """Swap in ``content`` with a rename inside the validated directory."""
    tmp = f".{target.name}.snapper-{secrets.token_hex(6)}"
    fd = os.open(
        tmp,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
        dir_fd=target.dir_fd,
    )
    try:
        try:
            os.fchmod(fd, stat.S_IMODE(target.st.st_mode))
            view = memoryview(content)
            while view:
                view = view[os.write(fd, view) :]
            os.fsync(fd)
        finally:
            os.close(fd)
        current = os.stat(target.name, dir_fd=target.dir_fd, follow_symlinks=False)
        if not stat.S_ISREG(current.st_mode) or (current.st_dev, current.st_ino) != (
            target.st.st_dev,
            target.st.st_ino,
        ):
            raise RefusedPath(f"{target.rel}: replaced while formatting")
        os.replace(tmp, target.name, src_dir_fd=target.dir_fd, dst_dir_fd=target.dir_fd)
    except BaseException:
        try:
            os.unlink(tmp, dir_fd=target.dir_fd)
        except FileNotFoundError:
            pass
        raise


def _snapper_base(args: argparse.Namespace) -> list[str]:
    return [
        args.snapper,
        "--native",
        "--max-width",
        str(args.max_width),
        "--config",
        args.config,
    ]


def files_to_rewrap(args: argparse.Namespace, root: Path, rels: list[str]) -> set[str]:
    result = subprocess.run(
        [*_snapper_base(args), "--check", "--output-format", "json", "--", *rels],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode not in (0, 1):
        sys.stderr.write(result.stderr)
        raise SystemExit(result.returncode)
    try:
        report = json.loads(result.stdout or "[]")
    except json.JSONDecodeError as exc:
        sys.stderr.write(result.stdout + result.stderr)
        raise SystemExit("snapper --check did not return JSON") from exc
    return {
        os.path.normpath(entry["file"])
        for entry in report
        if entry.get("would_reformat")
    }


def rewrap(args: argparse.Namespace, root: Path, target: Target) -> bool:
    result = subprocess.run(
        [*_snapper_base(args), "--stdin-filepath", target.rel],
        cwd=root,
        input=target.content,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stderr.decode(errors="replace"))
        raise SystemExit(f"snapper failed on {target.rel}")
    if result.stdout == target.content:
        return False
    replace_target(target, result.stdout)
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--snapper", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--root", required=True)
    parser.add_argument("--max-width", type=int, default=80)
    parser.add_argument("paths", nargs="+")
    args = parser.parse_args(argv)

    root = Path(os.path.realpath(args.root))
    root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    targets: list[Target] = []
    refused: list[str] = []
    try:
        for path in args.paths:
            try:
                targets.append(open_target(root_fd, root, path))
            except RefusedPath as exc:
                refused.append(str(exc))
        if refused:
            for message in refused:
                print(f"refusing to write {message}", file=sys.stderr)
            return 1
        pending = files_to_rewrap(args, root, [t.rel for t in targets])
        for target in targets:
            if target.rel in pending and rewrap(args, root, target):
                print(f"rewrapped {target.rel}")
        return 0
    except RefusedPath as exc:
        print(f"refusing to write {exc}", file=sys.stderr)
        return 1
    finally:
        for target in targets:
            os.close(target.dir_fd)
        os.close(root_fd)


if __name__ == "__main__":
    raise SystemExit(main())
