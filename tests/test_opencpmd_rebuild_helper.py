#!/usr/bin/env python3
"""Verify that the OpenCPMD helper rebuilds real dependency targets."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    helper = repo / "tools" / "rebuild_opencpmd_embed.sh"

    with tempfile.TemporaryDirectory() as scratch:
        root = Path(scratch) / "cpmd-root"
        (root / "obj").mkdir(parents=True)
        (root / "lib").mkdir()
        (root / "Makefile").write_text("# test fixture\n")

        log = Path(scratch) / "make-args.json"
        fake_make = Path(scratch) / "make"
        fake_make.write_text(
            "#!/usr/bin/env python3\n"
            "import json, os, sys\n"
            "open(os.environ['MAKE_LOG'], 'w').write(json.dumps(sys.argv[1:]))\n"
        )
        fake_make.chmod(0o755)

        # One member per patched OpenCPMD source file, plus scex_utils for
        # module order.
        members = [
            "error_handling.mod.o",
            "scex_utils.mod.o",
            "embed_ctrl.mod.o",
            "rinitwf_driver.mod.o",
            "initrun_driver.mod.o",
            "rwfopt_utils.mod.o",
            "updwf_utils.mod.o",
            "rkpnt_utils.mod.o",
            "timer.mod.o",
            "c_mem_utils.o",
            "broyden_utils.mod.o",
            "calc_alm_utils.mod.o",
            "chksym_utils.mod.o",
            "detdof_utils.mod.o",
            "drhov_utils.mod.o",
            "ehpsi_utils.mod.o",
            "fftnew_utils.mod.o",
            "forcedr_driver.mod.o",
            "hfx_utils.mod.o",
            "initclust_utils.mod.o",
            "k_odiis_utils.mod.o",
            "k_pcgrad_utils.mod.o",
            "mixing_g_utils.mod.o",
            "mixing_r_utils.mod.o",
            "moverho_utils.mod.o",
            "nlccset_utils.mod.o",
            "numpw_utils.mod.o",
            "phfac_utils.mod.o",
            "pw_hfx.mod.o",
            "qvan2_utils.mod.o",
            "recpnew_utils.mod.o",
            "rhodiis_utils.mod.o",
            "rnlsmd_utils.mod.o",
            "rpiiint_utils.mod.o",
            "rwswap_utils.mod.o",
            "setbasis_utils.mod.o",
            "symtrz_utils.mod.o",
            "tauofr_utils.mod.o",
            "testex_utils.mod.o",
            "updrho_utils.mod.o",
            "vdw_utils.mod.o",
            "vpsi_utils.mod.o",
            "embed_teardown.mod.o",
        ]
        for name in members:
            (root / "obj" / name).write_bytes(b"\0")

        env = os.environ.copy()
        env["MAKE"] = str(fake_make)
        env["MAKE_LOG"] = str(log)
        subprocess.run(
            [str(helper), str(root), "3"],
            cwd=repo,
            env=env,
            check=True,
        )

        actual = json.loads(log.read_text())
        expected = [
            "-C",
            str(root / "obj"),
            "-f",
            str(root / "Makefile"),
            "-j1",
            *members,
        ]
        if actual != expected:
            raise AssertionError(f"make arguments {actual!r}, expected {expected!r}")
        listed = subprocess.run(
            ["ar", "t", str(root / "lib" / "libcpmd.a")],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.split()
        if listed != members:
            raise AssertionError(f"archive members {listed!r}, expected {members!r}")

    print("OK: OpenCPMD helper rebuilds patched objects and archive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
