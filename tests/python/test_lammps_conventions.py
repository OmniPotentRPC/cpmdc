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
from lammps_import import import_rgpot_lammps

if import_rgpot_lammps() is None:
    pytest.skip("rgpot.lammps is not installed", allow_module_level=True)

TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def water(**kwargs: Any) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


@pytest.fixture
def tolerance(opencpmd) -> float:
    return 1e-5 if opencpmd else 0.0


def potential_without_lammps(params: Path, *, periodic: bool = False) -> Any:
    """The cpmdc potential adapter without an rgpot LAMMPS driver."""
    from cpmdc.lammps import _SessionPotential

    potential = _SessionPotential(Path(params).read_bytes(), None, periodic)
    potential.open()
    return potential


def test_the_lammps_step_agrees_with_the_ase_step_in_metal_units(periodic_params, tolerance):
    potential = potential_without_lammps(periodic_params, periodic=True)
    positions = np.array(WATER)
    rows = np.array(TRICLINIC)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    result = potential.evaluate(positions, numbers, rows)
    atoms = water(cell=rows, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        assert result["energy"] == pytest.approx(
            atoms.get_potential_energy(), rel=1e-12, abs=tolerance
        )
        np.testing.assert_allclose(
            result["forces"],
            atoms.get_forces(),
            rtol=1e-9,
            atol=1e-12 + tolerance,
        )
        np.testing.assert_allclose(
            result["stress"],
            voigt_6_to_full_3x3_stress(atoms.get_stress()),
            rtol=1e-6,
            atol=1e-9,
        )
    potential.close()


def test_the_lammps_step_refuses_a_non_finite_result_like_ase_does(
    params, opencpmd, cold_start
):
    """The reference well returns an infinite energy for this geometry.

    OpenCPMD stops in RATOM while it reads the coordinates of the deck, so
    no result exists and the library gives its own reason, which the
    ``cold_start`` fixture makes the one to expect. The callback halts
    LAMMPS on any error.
    """
    potential = potential_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    reason = r"^CPMD stopgm during embed SCF$" if opencpmd else "non-finite"
    with pytest.raises(CPMDCError, match=reason):
        potential.evaluate(np.array(WATER) * 1e200, numbers, np.eye(3))
    potential.evaluate(np.array(WATER), numbers, np.eye(3))
    potential.close()


@pytest.mark.parametrize("name", ["energy", "forces", "virial"])
def test_the_lammps_step_refuses_a_non_finite_value_from_the_library(
    params, monkeypatch, name
):
    import cpmdc.lammps as adapter

    def poisoned(self, positions, numbers, cell=None, *, forces, virial):
        forces[:] = 0.0
        virial[:] = 0.0
        energy = 0.0
        if name == "energy":
            energy = float("nan")
        else:
            {"forces": forces, "virial": virial}[name][0] = float("inf")
        return energy, forces, virial

    monkeypatch.setattr(adapter.Session, "calculate_au", poisoned)
    potential = potential_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    with pytest.raises(CPMDCError, match=f"non-finite {name}"):
        potential.evaluate(np.array(WATER), numbers, np.eye(3))
    potential.close()


def test_the_lammps_step_reopens_after_a_failed_new_composition(params, monkeypatch):
    import cpmdc.lammps as adapter

    potential = potential_without_lammps(params)
    water_numbers = np.array([8, 1, 1], dtype=np.int32)
    energy = potential.evaluate(np.array(WATER), water_numbers, np.eye(3))["energy"]
    opened = adapter.Session

    def refuse(*args, **kwargs):
        raise CPMDCError("no session")

    monkeypatch.setattr(adapter, "Session", refuse)
    with pytest.raises(CPMDCError, match="no session"):
        potential.evaluate(
            np.array(WATER), np.array([8, 8, 1], dtype=np.int32), np.eye(3)
        )
    monkeypatch.setattr(adapter, "Session", opened)
    assert (
        potential.evaluate(np.array(WATER), water_numbers, np.eye(3))["energy"]
        == energy
    )
    potential.close()
