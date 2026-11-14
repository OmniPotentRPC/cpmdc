#!/usr/bin/env python3
"""Run the Python tests and the ASE example against an installed cpmdc.

Usage: test_installed_wheel.py SOURCE_ROOT

The source tree supplies only the test inputs: test_python_api.py, the
pytest suite in tests/python (which needs pytest, as the ``test`` extra
installs), and the example. ``cpmdc`` must import from site-packages and
load the library and schema that the wheel carries, so the library
variables are cleared and the steps run in a scratch directory. ``capnp``
must be on PATH to encode the parameter messages.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_CLEARED = ("CPMDC_LIBRARY", "RGPOT_CPMDC_ENGINE", "RGPOT_CPMD_ENGINE", "CPMDC_SCHEMA")


def _encode(capnp: str, schema: Path, source: Path, target: Path) -> None:
    done = subprocess.run(
        [capnp, "encode", str(schema), "CPMDParams"],
        input=source.read_bytes(),
        stdout=subprocess.PIPE,
        check=True,
    )
    target.write_bytes(done.stdout)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_installed_wheel.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1]).resolve()
    capnp = shutil.which("capnp")
    if capnp is None:
        print("capnp is not on PATH", file=sys.stderr)
        return 2
    env = {key: value for key, value in os.environ.items() if key not in _CLEARED}
    env.pop("PYTHONPATH", None)
    schema = root / "schema" / "Potentials.capnp"
    with tempfile.TemporaryDirectory(prefix="cpmdc-wheel-") as scratch:
        work = Path(scratch)
        params = work / "water_blyp_params.bin"
        periodic = work / "water_blyp_periodic_params.bin"
        _encode(capnp, schema, root / "tests" / "water_blyp_params.capnp.txt", params)
        _encode(
            capnp, schema, root / "tests" / "water_blyp_periodic_params.capnp.txt", periodic
        )
        subprocess.run(
            [
                sys.executable,
                str(root / "tests" / "test_python_api.py"),
                "--installed",
                str(root),
                "-",
                str(params),
                str(periodic),
            ],
            cwd=work,
            env=env,
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                str(root / "tests" / "python"),
            ],
            cwd=work,
            env={
                **env,
                "CPMDC_TEST_INSTALLED": "1",
                "CPMDC_TEST_PARAMS": str(params),
                "CPMDC_TEST_PERIODIC_PARAMS": str(periodic),
            },
            check=True,
        )
        example = subprocess.run(
            [sys.executable, str(root / "examples" / "python" / "ase_energy.py")],
            cwd=work,
            env={**env, "CPMDC_PARAMS": str(params)},
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )
    lines = example.stdout.splitlines()
    print("\n".join(lines))
    if not lines or not lines[0].startswith("energy_ev=") or "natoms=3" not in lines:
        raise AssertionError(f"unexpected ase_energy.py output: {example.stdout!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
