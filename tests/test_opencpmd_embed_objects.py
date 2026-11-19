#!/usr/bin/env python3
"""Fail when a patched OpenCPMD source is missing from the rebuild lists.

Every compilable source named by tools/opencpmd_*.patch must be an object
in tools/rebuild_opencpmd_embed.sh, or in the dedicated rebuild script that
replaces that one member. src/SOURCES is a makefile fragment and has no object.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_SOURCE_SUFFIXES = (".F90", ".f90", ".c", ".C")
_PREFIXES = ("a/", "b/", "c/", "w/")


def _strip_prefix(path: str) -> str:
    for prefix in _PREFIXES:
        if path.startswith(prefix):
            return path[len(prefix) :]
    return path


def patched_sources(text: str) -> set[str]:
    found: set[str] = set()
    for line in text.splitlines():
        paths: list[str] = []
        if line.startswith("diff --git "):
            parts = line.split()
            paths.extend(parts[2:4])
        elif line.startswith("--- ") or line.startswith("+++ "):
            token = line.split(maxsplit=1)[1].split("\t", 1)[0]
            if token != "/dev/null":
                paths.append(token)
        for token in paths:
            path = _strip_prefix(token)
            if path.startswith("src/"):
                found.add(path)
    return found


def object_name(source: str) -> str | None:
    base = Path(source).name
    if base == "SOURCES":
        return None
    for suffix in _SOURCE_SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)] + ".o"
    return None


def embed_objects(text: str) -> list[str]:
    match = re.search(r"(?m)^objects=\((.*?)\n\)", text, re.S)
    if match is None:
        raise AssertionError("rebuild_opencpmd_embed.sh has no objects array")
    return re.findall(r"[A-Za-z0-9_.-]+\.o", match.group(1))


def script_objects(text: str) -> set[str]:
    return set(re.findall(r"[A-Za-z0-9_.-]+\.o", text))


def main() -> int:
    repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    tools = repo / "tools"
    embed_path = tools / "rebuild_opencpmd_embed.sh"
    embed = embed_objects(embed_path.read_text(encoding="utf-8"))
    if len(embed) != len(set(embed)):
        raise AssertionError(f"duplicate objects in {embed_path.name}: {embed}")
    dedicated: dict[str, set[str]] = {}
    for path in sorted(tools.glob("rebuild_opencpmd_*.sh")):
        if path.name == "rebuild_opencpmd_embed.sh":
            continue
        dedicated[path.name] = script_objects(path.read_text(encoding="utf-8"))
    covered_elsewhere = set().union(*dedicated.values()) if dedicated else set()
    missing: list[str] = []
    seen_objects: set[str] = set()
    patches = sorted(tools.glob("opencpmd_*.patch"))
    if not patches:
        raise AssertionError(f"no patches in {tools}")
    for patch in patches:
        sources = patched_sources(patch.read_text(encoding="utf-8"))
        if not sources:
            raise AssertionError(f"{patch.name}: no source paths")
        for source in sorted(sources):
            obj = object_name(source)
            if obj is None:
                continue
            seen_objects.add(obj)
            if obj in embed or obj in covered_elsewhere:
                continue
            owners = [
                name for name, objs in dedicated.items() if obj in objs
            ]
            where = ", ".join(owners) if owners else "no rebuild script"
            missing.append(f"{obj} from {patch.name} ({source}); {where}")
    if missing:
        print("objects missing from rebuild_opencpmd_embed.sh:", file=sys.stderr)
        for line in missing:
            print(f"  {line}", file=sys.stderr)
        return 1
    print(
        f"OK: {len(seen_objects)} patched objects are in the embed array "
        f"or a dedicated rebuild script"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
