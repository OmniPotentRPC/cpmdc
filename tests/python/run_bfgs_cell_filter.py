#!/usr/bin/env python3
"""OpenCPMD branch of test_bfgs_moves_the_cell_through_a_cell_filter.

Meson runs the whole ``tests/python`` tree in one interpreter. OpenCPMD keeps
wavefunctions between calculator sessions in that process, which can stall cell
BFGS before ASE reports convergence even after ``park_orbitals``. A fresh
interpreter matches the isolated physics check.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if os.environ.get("CPMDC_TEST_INSTALLED") != "1":
    sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests" / "python"))

import numpy as np
import pytest

from conftest import park_orbitals
from cpmdc.ase import CPMDCCalculator

# Import after path setup; normal pytest discovery does not match this filename.
from test_ase_calculator import CUBE, relax_cell, water  # noqa: E402

def test_bfgs_cell_filter(periodic_params: Path, params: Path, opencpmd: bool) -> None:
    if not opencpmd:
        pytest.skip("the reference well has no OpenCPMD orbitals")
    park_orbitals(params)
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=17)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        converged, energies = relax_cell(atoms, steps=80)
    vol_ratio = atoms.get_volume() / np.linalg.det(CUBE)
    stress_max = float(np.abs(atoms.get_stress()[:3]).max())
    assert energies[-1] < energies[0]
    assert np.all(np.diff(energies[:5]) < 0.0)
    assert converged, f"vol_ratio={vol_ratio:.4f} stress_max={stress_max:.2e}"
    assert 0.9 < vol_ratio < 1.1
    assert stress_max < 1.0e-3


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
