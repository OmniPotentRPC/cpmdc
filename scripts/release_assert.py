#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text()


def write(path: str, text: str) -> None:
    (ROOT / path).write_text(text)


def replace_one(pattern: str, repl: str, text: str, label: str) -> str:
    new, n = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit(f"could not update {label}")
    return new


def meson_version() -> str:
    match = re.search(r"version:\s*'([^']+)'", read("meson.build"))
    if not match:
        raise SystemExit("meson.build version not found")
    return match.group(1)


def towncrier_version() -> str:
    match = re.search(r'^version\s*=\s*"([^"]+)"', read("towncrier.toml"), re.M)
    if not match:
        raise SystemExit("towncrier.toml version not found")
    return match.group(1)


def pixi_version() -> str:
    match = re.search(r'^version\s*=\s*"([^"]+)"', read("pixi.toml"), re.M)
    if not match:
        raise SystemExit("pixi.toml version not found")
    return match.group(1)


def cmake_version() -> str | None:
    path = ROOT / "CMakeLists.txt"
    if not path.exists():
        return None
    match = re.search(r"project\s*\([^)]*VERSION\s+([0-9][^)\s]*)", path.read_text(), re.S)
    if not match:
        raise SystemExit("CMakeLists.txt version not found")
    return match.group(1)


def sync(version: str) -> None:
    write(
        "meson.build",
        replace_one(
            r"version:\s*'[^']+'",
            f"version: '{version}'",
            read("meson.build"),
            "meson version",
        ),
    )
    write(
        "towncrier.toml",
        replace_one(
            r'^version\s*=\s*"[^"]+"',
            f'version = "{version}"',
            read("towncrier.toml"),
            "towncrier version",
        ),
    )
    write(
        "pixi.toml",
        replace_one(
            r'^version\s*=\s*"[^"]+"',
            f'version = "{version}"',
            read("pixi.toml"),
            "pixi version",
        ),
    )
    cmake_path = ROOT / "CMakeLists.txt"
    if cmake_path.exists():
        write(
            "CMakeLists.txt",
            replace_one(
                r"(project\s*\(\s*cpmdc\s+VERSION\s+)[^\s)]+",
                rf"\g<1>{version}",
                read("CMakeLists.txt"),
                "CMake version",
            ),
        )


def assert_versions(expected: str | None) -> str:
    versions = {
        "meson.build": meson_version(),
        "pixi.toml": pixi_version(),
        "towncrier.toml": towncrier_version(),
    }
    cmake = cmake_version()
    if cmake is not None:
        versions["CMakeLists.txt"] = cmake
    unique = set(versions.values())
    if len(unique) != 1:
        rendered = ", ".join(f"{k}={v}" for k, v in sorted(versions.items()))
        raise SystemExit(f"version mismatch: {rendered}")
    version = unique.pop()
    if expected is not None and version != expected:
        raise SystemExit(f"version mismatch: expected {expected}, found {version}")
    return version


def require_changelog(version: str) -> None:
    marker = f"## [{version}]"
    if marker not in read("CHANGELOG.md"):
        raise SystemExit(f"CHANGELOG.md missing section for {version}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("version", nargs="?")
    parser.add_argument("--sync", dest="sync_version")
    parser.add_argument("--require-changelog", action="store_true")
    args = parser.parse_args()

    if args.sync_version:
        sync(args.sync_version)
        expected = args.sync_version
    else:
        expected = args.version

    version = assert_versions(expected)
    if args.require_changelog:
        require_changelog(args.version or version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
