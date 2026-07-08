#!/usr/bin/env python3
"""Run capnpc-fortran and place potentials_capnp.f90 at the Meson @OUTPUT@ path."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 7:
        print(
            "usage: capnpc_fortran_compile.py CAPNP CAPNPC_FORTRAN SRC_PREFIX "
            "WORK_DIR INPUT OUTPUT",
            file=sys.stderr,
        )
        return 2
    capnp, plugin, src_prefix, work_dir, input_path, output_path = sys.argv[1:]
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            capnp,
            "compile",
            "--src-prefix=" + src_prefix,
            "-o" + plugin,
            input_path,
        ],
        check=True,
        cwd=str(work),
    )
    generated = work / "potentials_capnp.f90"
    if not generated.is_file():
        matches = list(work.glob("*_capnp.f90"))
        if not matches:
            print("capnpc-fortran produced no *_capnp.f90 in", work, file=sys.stderr)
            return 1
        generated = matches[0]
    dest = Path(output_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if generated.resolve() != dest.resolve():
        shutil.copy2(generated, dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
