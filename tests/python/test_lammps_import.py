"""The optional LAMMPS frontend guard skips only genuinely absent modules."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def test_broken_rgpot_lammps_fails_loudly(tmp_path):
    package = tmp_path / "rgpot"
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "lammps.py").write_text("import rgpot._core\n", encoding="utf-8")
    here = Path(__file__).resolve().parent
    env = {
        **os.environ,
        "PYTHONPATH": os.pathsep.join([str(tmp_path), str(here)]),
    }
    done = subprocess.run(
        [
            sys.executable,
            "-c",
            "from lammps_import import import_rgpot_lammps; import_rgpot_lammps()",
        ],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert done.returncode != 0
    assert "rgpot._core" in done.stderr
