"""cpmdc.lammps.FixExternal against the conventions of CPMDCCalculator.

The evaluation step of ``FixExternal`` is exercised without a LAMMPS
object, so these cases run on the default build whether or not the
``lammps`` Python module is installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from ase import Atoms
from ase.stress import voigt_6_to_full_3x3_stress
from conftest import WATER

from cpmdc import CPMDCError
from cpmdc.ase import CPMDCCalculator

TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def water(**kwargs: Any) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


@pytest.fixture
def tolerance(opencpmd) -> float:
    return 1e-5 if opencpmd else 0.0


def fix_without_lammps(params: Path, *, periodic: bool = False) -> Any:
    """A FixExternal that can run its evaluation step without a LAMMPS object."""
    from cpmdc.lammps import FixExternal, unit_factors

    fix = object.__new__(FixExternal)
    fix._blob = Path(params).read_bytes()
    fix._library = None
    fix._session = None
    fix._composition = None
    fix._forces = np.zeros(0)
    fix._virial = np.zeros(9)
    fix.periodic = periodic
    fix.length_per_bohr, fix.energy_per_hartree = unit_factors("metal")
    return fix


def test_the_lammps_step_agrees_with_the_ase_step_in_metal_units(periodic_params, tolerance):
    fix = fix_without_lammps(periodic_params, periodic=True)
    positions = np.array(WATER)
    rows = np.array(TRICLINIC)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    tags = np.array([1, 2, 3])
    energy, forces, virial = fix._evaluate(tags, positions, numbers, rows)
    atoms = water(cell=rows, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        assert energy == pytest.approx(atoms.get_potential_energy(), rel=1e-12, abs=tolerance)
        np.testing.assert_allclose(forces, atoms.get_forces(), rtol=1e-9, atol=1e-12 + tolerance)
        volume = atoms.get_volume()
        np.testing.assert_allclose(
            virial, -volume * voigt_6_to_full_3x3_stress(atoms.get_stress()), rtol=1e-6, atol=1e-9
        )
    fix._session.close()


def test_the_lammps_step_refuses_a_non_finite_result_like_ase_does(params):
    fix = fix_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    tags = np.array([1, 2, 3])
    with pytest.raises(CPMDCError, match="non-finite"):
        fix._evaluate(tags, np.array(WATER) * 1e200, numbers, np.eye(3))
    fix._evaluate(tags, np.array(WATER), numbers, np.eye(3))
    fix._session.close()


def test_the_lammps_step_reopens_after_a_failed_new_composition(params, monkeypatch):
    import cpmdc.lammps as adapter

    fix = fix_without_lammps(params)
    tags = np.array([1, 2, 3])
    water_numbers = np.array([8, 1, 1], dtype=np.int32)
    energy = fix._evaluate(tags, np.array(WATER), water_numbers, np.eye(3))[0]
    opened = adapter.Session

    def refuse(*args, **kwargs):
        raise CPMDCError("no session")

    monkeypatch.setattr(adapter, "Session", refuse)
    with pytest.raises(CPMDCError, match="no session"):
        fix._evaluate(tags, np.array(WATER), np.array([8, 8, 1], dtype=np.int32), np.eye(3))
    monkeypatch.setattr(adapter, "Session", opened)
    assert fix._evaluate(tags, np.array(WATER), water_numbers, np.eye(3))[0] == energy
    fix._session.close()
