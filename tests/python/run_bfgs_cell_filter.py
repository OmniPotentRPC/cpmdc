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
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests" / "python"))

import numpy as np

from conftest import park_orbitals
from cpmdc.ase import CPMDCCalculator

# Import after path setup; script is not collected by pytest.
from test_ase_calculator import CUBE, relax_cell, water  # noqa: E402

periodic_params = Path(os.environ["CPMDC_TEST_PERIODIC_PARAMS"])
cluster_params = Path(os.environ["CPMDC_TEST_PARAMS"])

park_orbitals(cluster_params)
atoms = water(cell=CUBE, pbc=True)
atoms.rattle(0.05, seed=17)
with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
    atoms.calc = calc
    converged, energies = relax_cell(atoms, steps=80)
vol_ratio = atoms.get_volume() / np.linalg.det(CUBE)
stress_max = float(np.abs(atoms.get_stress()[:3]).max())
if not (energies[-1] < energies[0]):
    raise SystemExit("energy did not fall")
if not np.all(np.diff(energies[:5]) < 0.0):
    raise SystemExit("first five steps did not lower energy")
if not converged:
    raise SystemExit(
        f"BFGS did not converge: vol_ratio={vol_ratio:.4f} "
        f"stress_max={stress_max:.2e}"
    )
if not (0.9 < vol_ratio < 1.1):
    raise SystemExit(f"volume ratio {vol_ratio:.4f} outside 10% band")
if stress_max >= 1.0e-3:
    raise SystemExit(f"stress_max {stress_max:.2e} too large")