#!/usr/bin/env python3
import argparse
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "version.txt"
COG = ["uvx", "--from", "cogapp==3.6.0", "cog"]
COG_FILES = [
    "pixi.toml",
    "CITATION.cff",
    "include/cpmdc.h",
    "docs/orgmode/tutorials/first-energy.org",
    "docs/source/Doxyfile_cpmdc.cfg",
]


def read_version() -> str:
    raw = VERSION_FILE.read_text()
    version = raw.strip()
    # Meson keeps surrounding whitespace that cog strips, so allow none.
    if raw != f"{version}\n" or not re.fullmatch(
        r"[0-9]+\.[0-9]+\.[0-9]+(?:[-+.][0-9A-Za-z.+-]+)?", version
    ):
        raise SystemExit(f"version.txt is not one version line: {raw!r}")
    return version


def require_no_literals() -> None:
    if "version: files('version.txt')" not in (ROOT / "meson.build").read_text():
        raise SystemExit("meson.build does not read its version from version.txt")
    if re.search(r"^\s*version\s*=", (ROOT / "towncrier.toml").read_text(), re.M):
        raise SystemExit("towncrier.toml carries a version; pass --version instead")


def cog(*args: str) -> None:
    result = subprocess.run([*COG, *args, *COG_FILES], cwd=ROOT)
    if result.returncode != 0:
        raise SystemExit(
            "generated version copies differ from version.txt;"
            " run `pixi run version-sync`"
        )


def require_changelog(version: str) -> None:
    marker = f"## [{version}]"
    if marker not in (ROOT / "CHANGELOG.md").read_text():
        raise SystemExit(f"CHANGELOG.md missing section for {version}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="version.txt is the only hand-written version."
    )
    parser.add_argument("version", nargs="?")
    parser.add_argument(
        "--sync",
        action="store_true",
        help="write VERSION to version.txt if given, then regenerate every copy",
    )
    parser.add_argument("--require-changelog", action="store_true")
    args = parser.parse_args()

    require_no_literals()
    if args.sync:
        if args.version:
            VERSION_FILE.write_text(f"{args.version}\n")
        read_version()
        cog("-r")
        return 0

    version = read_version()
    if args.version is not None and version != args.version:
        raise SystemExit(f"version mismatch: expected {args.version}, found {version}")
    if args.require_changelog:
        require_changelog(version)
    cog("--check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
