#!/usr/bin/env python3
"""Install cpmdc into a prefix, then compile and run consumers against it.

Legs:
- pkg-config: cc + `pkg-config --cflags --libs cpmdc`
- CMake: `find_package(cpmdc CONFIG REQUIRED)` when cmake is available
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

CONSUMER_C = r"""
#include "cpmdc.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  const char *version = cpmdc_version();
  if (version == NULL || strstr(version, "cpmdc") == NULL) {
    fprintf(stderr, "unexpected version string\n");
    return 1;
  }
  if (cpmdc_abi_version() != CPMDC_ABI_VERSION) {
    fprintf(stderr, "ABI version mismatch\n");
    return 1;
  }
  if (cpmdc_session_create(NULL, 0) != NULL) {
    fprintf(stderr, "NULL params must not create a session\n");
    return 1;
  }
  printf("installed-consumer ok: %s abi=%d\n", version, cpmdc_abi_version());
  return 0;
}
"""

CONSUMER_CMAKE = """
cmake_minimum_required(VERSION 3.24)
project(consumer C)
find_package(cpmdc CONFIG REQUIRED)
add_executable(consumer main.c)
target_link_libraries(consumer PRIVATE cpmdc::cpmdc)
"""


def run(cmd, **kw):
    print("+", shlex.join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def lib_dirs(prefix: Path):
    return [d for d in (prefix / "lib", prefix / "lib64") if d.is_dir()] + list(
        (prefix / "lib").glob("x86_64*")
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--build-root", required=True)
    parser.add_argument("--install-prefix", required=True)
    parser.add_argument("--meson", default="meson")
    parser.add_argument("--cc", default=os.environ.get("CC", "cc"))
    parser.add_argument("--pkg-config", default="pkg-config")
    args = parser.parse_args()

    source = Path(args.source_dir).resolve()
    build_root = Path(args.build_root).resolve()
    prefix = Path(args.install_prefix).resolve()
    for stale in (build_root, prefix):
        if stale.exists():
            shutil.rmtree(stale)
    build_root.mkdir(parents=True)

    inner_build = build_root / "package-build"
    run(
        [
            args.meson,
            "setup",
            inner_build,
            source,
            "--prefix",
            prefix,
            "-Dwith_tests=false",
        ]
    )
    run([args.meson, "install", "-C", inner_build])

    pc_dirs = [d / "pkgconfig" for d in lib_dirs(prefix) if (d / "pkgconfig").is_dir()]
    if not pc_dirs:
        raise SystemExit(f"no pkgconfig dir under {prefix}")
    env = dict(os.environ)
    env["PKG_CONFIG_PATH"] = os.pathsep.join(str(d) for d in pc_dirs)

    consumer_dir = build_root / "consumer"
    consumer_dir.mkdir()
    (consumer_dir / "main.c").write_text(CONSUMER_C, encoding="utf-8")

    cflags = subprocess.run(
        [args.pkg_config, "--cflags", "cpmdc"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.split()
    libs = subprocess.run(
        [args.pkg_config, "--libs", "cpmdc"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.split()

    consumer_bin = consumer_dir / "consumer"
    run([args.cc, consumer_dir / "main.c", "-o", consumer_bin] + cflags + libs)

    run_env = dict(os.environ)
    run_env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(d) for d in lib_dirs(prefix)] + [run_env.get("LD_LIBRARY_PATH", "")]
    )
    run([consumer_bin], env=run_env)

    cmake = shutil.which("cmake")
    if cmake is None:
        print("cmake not found; skipping find_package leg", flush=True)
        return 0
    cmake_dir = build_root / "cmake-consumer"
    cmake_dir.mkdir()
    (cmake_dir / "main.c").write_text(CONSUMER_C, encoding="utf-8")
    (cmake_dir / "CMakeLists.txt").write_text(CONSUMER_CMAKE, encoding="utf-8")
    cmake_build = cmake_dir / "build"
    run(
        [
            cmake,
            "-S",
            cmake_dir,
            "-B",
            cmake_build,
            f"-DCMAKE_PREFIX_PATH={prefix}",
            f"-DCMAKE_C_COMPILER={args.cc}",
        ]
    )
    run([cmake, "--build", cmake_build])
    run([cmake_build / "consumer"], env=run_env)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
