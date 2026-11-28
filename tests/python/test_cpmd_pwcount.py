"""bench/cpmd_pwcount.py prints the plane-wave count and the two FFT meshes."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("cpmd_pwcount", ROOT / "bench" / "cpmd_pwcount.py")
pw = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pw)


def test_the_840_cell_prints_the_count_and_the_meshes(capsys):
    assert pw.main() == 0
    assert capsys.readouterr().out.splitlines() == [
        "1139.60 Angstrom^3, 100 Ry: 65002 plane waves",
        "shipped FFT mesh: 128",
        "FFTW mesh: 126",
    ]
