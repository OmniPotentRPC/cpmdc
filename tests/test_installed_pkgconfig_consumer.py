#!/usr/bin/env python3
"""Install cpmdc into a prefix, then compile and run consumers against it.

Checks the installed prefix in place and after moving it to a relocated
directory, under plain `pkg-config` and `pkg-config --define-prefix`.
After the move, an installed `.so` or `.dylib` fails the run when its
`RUNPATH` or `RPATH` points at the build tree or the original prefix.

Legs:
- pkg-config: cc + `pkg-config --cflags --libs cpmdc` and `cpmdc-restart`
- libcp2k: cc + `pkg-config --cflags --libs libcp2k` from
  LIBDIR/cpmdc/gromacs/pkgconfig, the file GROMACS -DGMX_CP2K=ON reads, run
  without the prefix on LD_LIBRARY_PATH
- CMake: `find_package(cpmdc CONFIG REQUIRED)` when cmake is available

The nested `meson setup` runs with CC and FC set to the compilers the outer
build was configured with (`--cc`, `--fc`), so a compiler launcher in the
host environment does not reach it.
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
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

CONSUMER_RESTART_C = r"""
#include "cpmdc_restart.h"

#include <stdio.h>

int main(void) {
  if (cpmdc_restart_read_path(NULL, NULL, 0) != NULL) {
    fprintf(stderr, "NULL path must not read a restart image\n");
    return 1;
  }
  printf("installed-restart-consumer ok\n");
  return 0;
}
"""

CONSUMER_CP2K_C = r"""
#include "libcp2k.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  char version[256] = {0};
  cp2k_get_version(version, (int)sizeof(version));
  if (strstr(version, "cpmdc") == NULL) {
    fprintf(stderr, "unexpected cp2k_get_version: %s\n", version);
    return 1;
  }
  cp2k_init();
  cp2k_finalize();
  printf("installed libcp2k consumer ok: %s\n", version);
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


def installed_shared_objects(prefix: Path):
    """Real installed shared libraries and plugins, not symlinks."""
    found = []
    for path in prefix.rglob("*"):
        if not path.is_file() or path.is_symlink():
            continue
        name = path.name
        if ".so" in name or ".dylib" in name:
            found.append(path)
    return found


def _macho_or_elf(path: Path) -> bool:
    magic = path.read_bytes()[:4]
    return magic == b"\x7fELF" or magic in (
        b"\xfe\xed\xfa\xce",
        b"\xfe\xed\xfa\xcf",
        b"\xce\xfa\xed\xfe",
        b"\xcf\xfa\xed\xfe",
    )


def runtime_path_entries(path: Path) -> list[str]:
    """`DT_RPATH` and `DT_RUNPATH`, or Mach-O `LC_RPATH` entries."""
    if sys.platform == "darwin":
        argv = ["otool", "-l", str(path)]
    else:
        argv = ["readelf", "-d", str(path)]
    try:
        out = subprocess.run(argv, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise SystemExit(f"cannot read run paths, {exc.filename} is missing") from exc
    if out.returncode != 0:
        if _macho_or_elf(path):
            raise SystemExit(
                f"{argv[0]} failed on {path}: {out.stderr.strip() or out.stdout.strip()}"
            )
        return []
    entries: list[str] = []
    if sys.platform == "darwin":
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
    for line in out.stdout.splitlines():
        if "(RPATH)" not in line and "(RUNPATH)" not in line:
            continue
        if "[" not in line or "]" not in line:
            continue
        payload = line.split("[", 1)[1].rsplit("]", 1)[0]
        entries.extend(part for part in payload.split(":") if part)
    return entries


def _origin_relative(entry: str, library: Path) -> Path | None:
    """Resolve `$ORIGIN` and `@loader_path` against `library`'s directory.

    `@executable_path` is the main executable, which this check does not
    know, so it is left unresolved.
    """
    origin = library.parent
    for token in ("$ORIGIN", "${ORIGIN}", "@loader_path"):
        if entry == token:
            return origin.resolve()
        for sep in ("/", "\\"):
            prefix = token + sep
            if entry.startswith(prefix):
                return (origin / entry[len(prefix) :]).resolve()
    return None


def entry_hits_root(entry: str, root: Path, library: Path) -> bool:
    root_resolved = root.resolve()
    expanded = _origin_relative(entry, library)
    if expanded is not None:
        return expanded == root_resolved or root_resolved in expanded.parents
    root_s = str(root)
    if root_s and root_s in entry:
        return True
    candidate = Path(entry)
    if not candidate.is_absolute():
        return False
    resolved = candidate.resolve()
    return resolved == root_resolved or root_resolved in resolved.parents


def _assert_origin_resolution() -> None:
    """`$ORIGIN` and `@loader_path` are judged from the library directory."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp).resolve()
        libdir = root / "reloc" / "lib"
        libdir.mkdir(parents=True)
        library = libdir / "libcpmdc.so"
        library.write_bytes(b"")
        build = root / "build"
        build.mkdir()
        if entry_hits_root("$ORIGIN", build, library):
            raise SystemExit("library $ORIGIN was treated as the build tree")
        if entry_hits_root("@loader_path", build, library):
            raise SystemExit("library @loader_path was treated as the build tree")
        if not entry_hits_root("$ORIGIN/../../build", build, library):
            raise SystemExit("$ORIGIN into the build tree was ignored")
        if not entry_hits_root("@loader_path/../../build", build, library):
            raise SystemExit("@loader_path into the build tree was ignored")
        if not entry_hits_root(str(build / "lib"), build, library):
            raise SystemExit("absolute build tree entry was ignored")
        if entry_hits_root("@executable_path/../../build", build, library):
            raise SystemExit("@executable_path was resolved against the library")


def reject_fixed_runpaths(tree: Path, forbidden: list[Path]) -> None:
    libraries = installed_shared_objects(tree)
    if not libraries:
        raise SystemExit(f"no installed shared libraries under {tree}")
    for library in libraries:
        for entry in runtime_path_entries(library):
            for root in forbidden:
                if entry_hits_root(entry, root, library):
                    raise SystemExit(
                        f"{library}: run path {entry!r} points into {root} "
                        "(build tree or original prefix)"
                    )


def libcp2k_loader_flags() -> list[str]:
    """Make libcp2k resolve `libcpmdc` through the library, not the executable.

    `--as-needed` drops the unused direct `-lcpmdc`. `--enable-new-dtags`
    stores a `RUNPATH`, which the loader does not apply to dependencies of
    `libcpmdc_gromacs`.
    """
    if sys.platform.startswith("linux"):
        return ["-Wl,--as-needed", "-Wl,--enable-new-dtags"]
    return []


def lib_dirs(prefix: Path):
    """Directories that may hold the installed library and its pkg-config file.

    Debian multiarch uses lib/x86_64-linux-gnu and lib/aarch64-linux-gnu.
    A plain prefix uses lib or lib64.
    """
    found = [d for d in (prefix / "lib", prefix / "lib64") if d.is_dir()]
    lib = prefix / "lib"
    if lib.is_dir():
        found.extend(
            path
            for path in lib.iterdir()
            if path.is_dir()
            and (path.name.startswith("x86_64") or path.name.startswith("aarch64"))
        )
    return found


def check_prefix(
    prefix: Path,
    work_dir: Path,
    pkg_config: str,
    cc: str,
    extra_link: list[str],
    extra_run_dirs: list[str],
    extra_pc_args: list[str],
) -> None:
    pc_dirs = [d / "pkgconfig" for d in lib_dirs(prefix) if (d / "pkgconfig").is_dir()]
    if not pc_dirs:
        raise SystemExit(f"no pkgconfig dir under {prefix}")
    env = dict(os.environ)
    env["PKG_CONFIG_PATH"] = os.pathsep.join(str(d) for d in pc_dirs)

    consumer_dir = work_dir / "consumer"
    consumer_dir.mkdir(parents=True)
    (consumer_dir / "main.c").write_text(CONSUMER_C, encoding="utf-8")

    cflags = subprocess.run(
        [pkg_config, *extra_pc_args, "--cflags", "cpmdc"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.split()
    libs = subprocess.run(
        [pkg_config, *extra_pc_args, "--libs", "cpmdc"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.split()

    consumer_bin = consumer_dir / "consumer"
    run(
        [cc, consumer_dir / "main.c", "-o", consumer_bin]
        + cflags
        + libs
        + extra_link
    )

    run_env = dict(os.environ)
    run_env["LD_LIBRARY_PATH"] = os.pathsep.join(
        [str(d) for d in lib_dirs(prefix)]
        + extra_run_dirs
        + [run_env.get("LD_LIBRARY_PATH", "")]
    )
    run([consumer_bin], env=run_env)

    restart_dir = work_dir / "restart-consumer"
    restart_dir.mkdir()
    (restart_dir / "main.c").write_text(CONSUMER_RESTART_C, encoding="utf-8")
    restart_flags = subprocess.run(
        [pkg_config, *extra_pc_args, "--cflags", "--libs", "cpmdc-restart"],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    ).stdout.split()
    restart_bin = restart_dir / "consumer"
    run([cc, restart_dir / "main.c", "-o", restart_bin] + restart_flags)
    run([restart_bin])

    ipi_client = prefix / "bin" / "cpmdc-ipi-c"
    if os.name == "posix" and not ipi_client.is_file():
        raise SystemExit(f"cpmdc-ipi-c not installed under {prefix}")
    if ipi_client.is_file():
        run([ipi_client, "--help"], env=run_env)

    cp2k_pcs = sorted(prefix.rglob("cpmdc/gromacs/pkgconfig/libcp2k.pc"))
    if not cp2k_pcs:
        raise SystemExit(
            f"libcp2k.pc not installed under {prefix}/LIBDIR/cpmdc/gromacs/pkgconfig"
        )
    cp2k_env = dict(os.environ)
    cp2k_env["PKG_CONFIG_PATH"] = str(cp2k_pcs[0].parent)
    cp2k_flags = []
    for query in ("--cflags", "--libs"):
        cp2k_flags += subprocess.run(
            [pkg_config, *extra_pc_args, query, "libcp2k"],
            check=True,
            capture_output=True,
            text=True,
            env=cp2k_env,
        ).stdout.split()
    if not any(flag.startswith("-Wl,-rpath,") for flag in cp2k_flags):
        raise SystemExit(f"libcp2k.pc Libs carry no run path: {cp2k_flags}")
    cp2k_dir = work_dir / "libcp2k-consumer"
    cp2k_dir.mkdir()
    (cp2k_dir / "main.c").write_text(CONSUMER_CP2K_C, encoding="utf-8")
    cp2k_bin = cp2k_dir / "consumer"
    run(
        [cc, cp2k_dir / "main.c", "-o", cp2k_bin]
        + libcp2k_loader_flags()
        + cp2k_flags
        + extra_link
    )
    bare_env = dict(os.environ)
    bare_env["LD_LIBRARY_PATH"] = os.pathsep.join(extra_run_dirs)
    run([cp2k_bin], env=bare_env)

    cmake = shutil.which("cmake")
    if cmake is None:
        print("cmake not found; skipping find_package leg", flush=True)
        return
    cmake_dir = work_dir / "cmake-consumer"
    if cmake_dir.exists():
        return
    cmake_dir.mkdir()
    (cmake_dir / "main.c").write_text(CONSUMER_C, encoding="utf-8")
    (cmake_dir / "CMakeLists.txt").write_text(CONSUMER_CMAKE, encoding="utf-8")
    cmake_build = cmake_dir / "build"
    config_files = list(prefix.rglob("cpmdcConfig.cmake"))
    if not config_files:
        raise SystemExit(f"cpmdcConfig.cmake not installed under {prefix}")
    run(
        [
            cmake,
            "-S",
            cmake_dir,
            "-B",
            cmake_build,
            f"-DCMAKE_PREFIX_PATH={prefix}",
            f"-Dcpmdc_DIR={config_files[0].parent}",
            f"-DCMAKE_C_COMPILER={cc}",
            "-DCMAKE_EXE_LINKER_FLAGS=" + " ".join(extra_link),
        ]
    )
    run([cmake, "--build", cmake_build])
    run([cmake_build / "consumer"], env=run_env)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--build-root", required=True)
    parser.add_argument("--install-prefix", required=True)
    parser.add_argument("--meson", default="meson")
    parser.add_argument("--cc", default=os.environ.get("CC", "cc"))
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--pkg-config", default="pkg-config")
    args = parser.parse_args()
    _assert_origin_resolution()

    source = Path(args.source_dir).resolve()
    build_root = Path(args.build_root).resolve()
    prefix = Path(args.install_prefix).resolve()
    relocated = build_root / "relocated-prefix"
    for stale in (build_root, prefix):
        if stale.exists():
            shutil.rmtree(stale)
    build_root.mkdir(parents=True)

    # The nested setup compiles with the compilers of the outer build. A
    # launcher such as sccache in the host's CC does not reach it.
    setup_env = dict(os.environ)
    setup_env["CC"] = args.cc
    setup_env["FC"] = args.fc
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
            "-Dwith_lammps_plugin=disabled",
        ],
        env=setup_env,
    )
    run([args.meson, "install", "-C", inner_build], env=setup_env)

    # libcpmdc.so carries the Fortran embed shell; point the linker and the
    # runtime at the toolchain's libgfortran when cc lives in a conda-style env.
    toolchain_lib = Path(args.cc).resolve().parents[1] / "lib"
    extra_link = []
    extra_run_dirs = []
    if toolchain_lib.is_dir():
        extra_link = [f"-L{toolchain_lib}", f"-Wl,-rpath-link,{toolchain_lib}"]
        extra_run_dirs = [str(toolchain_lib)]

    check_prefix(
        prefix,
        build_root / "in-place",
        args.pkg_config,
        args.cc,
        extra_link,
        extra_run_dirs,
        [],
    )

    shutil.move(str(prefix), str(relocated))
    reject_fixed_runpaths(relocated, [prefix, inner_build])
    for label, pc_args in (("plain", []), ("define-prefix", ["--define-prefix"])):
        check_prefix(
            relocated,
            build_root / f"relocated-{label}",
            args.pkg_config,
            args.cc,
            extra_link,
            extra_run_dirs,
            pc_args,
        )
    shutil.move(str(relocated), str(prefix))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
