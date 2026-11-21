"""cpmdc.lammps.FixExternal against the conventions of CPMDCCalculator.

The evaluation step of ``FixExternal`` is exercised on an object built by
its real constructor against a stand-in for the LAMMPS object, so these
cases need ``rgpot.lammps`` but not the ``lammps`` Python module.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import pytest
from ase import Atoms
from ase.stress import voigt_6_to_full_3x3_stress
from conftest import WATER
from lammps_import import import_rgpot_lammps

from cpmdc import CPMDCError
from cpmdc.ase import CPMDCCalculator

if import_rgpot_lammps() is None:
    pytest.skip("rgpot.lammps is not installed", allow_module_level=True)

TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def water(**kwargs: Any) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


@pytest.fixture
def tolerance(opencpmd) -> float:
    return 1e-5 if opencpmd else 0.0


class StandInLammps:
    """The slice of the LAMMPS Python object that FixExternal.__init__ reads."""

    SETTINGS: ClassVar[dict[str, int]] = {
        "world_size": 1,
        "IMGMASK": 1023,
        "IMGMAX": 512,
        "IMGBITS": 10,
        "IMG2BITS": 20,
    }

    def __init__(self, periodic: bool, units: str = "metal") -> None:
        self.periodic = periodic
        self.units = units
        self.commands: list[str] = []
        self.callbacks: list[Any] = []

    def extract_box(self) -> tuple[Any, ...]:
        return ([0.0] * 3, [1.0] * 3, 0.0, 0.0, 0.0, [self.periodic] * 3, 0)

    def extract_global(self, name: str) -> Any:
        return {"units": self.units, "ntypes": 2}[name]

    def extract_setting(self, name: str) -> int:
        return self.SETTINGS[name]

    def command(self, text: str) -> None:
        self.commands.append(text)

    def set_fix_external_callback(self, fix_id: str, callback: Any, caller: Any) -> None:
        self.callbacks.append((fix_id, callback, caller))


def fix_without_lammps(params: Path, *, periodic: bool = False) -> Any:
    """A FixExternal from its own constructor, with a stand-in LAMMPS object."""
    from cpmdc.lammps import FixExternal

    return FixExternal(StandInLammps(periodic), params, "O H")


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
    fix.close()


def test_the_lammps_step_refuses_a_non_finite_result_like_ase_does(
    params, opencpmd, cold_start
):
    """The reference well returns an infinite energy for this geometry.

    OpenCPMD stops in RATOM while it reads the coordinates of the deck, so
    no result exists and the library gives its own reason, which the
    ``cold_start`` fixture makes the one to expect. The callback halts
    LAMMPS on any error.
    """
    fix = fix_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    tags = np.array([1, 2, 3])
    reason = r"^CPMD stopgm during embed SCF$" if opencpmd else "non-finite"
    with pytest.raises(CPMDCError, match=reason):
        fix._evaluate(tags, np.array(WATER) * 1e200, numbers, np.eye(3))
    fix._evaluate(tags, np.array(WATER), numbers, np.eye(3))
    fix.close()


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
    fix = fix_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    with pytest.raises(CPMDCError, match=f"non-finite {name}"):
        fix._evaluate(np.array([1, 2, 3]), np.array(WATER), numbers, np.eye(3))
    fix.close()


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
    fix.close()


def test_a_constructed_fix_evaluates_through_the_session_adapter(params):
    """The evaluation needs nothing the constructor did not set up.

    An earlier ``_evaluate`` read attributes that only a previous driver
    initialised. Tests that set those attributes by hand on an object made
    without ``__init__`` kept passing while every LAMMPS step failed with an
    AttributeError, so the object here comes from the real constructor.
    """
    fix = fix_without_lammps(params)
    try:
        assert fix._session is fix._adapter.session
        assert fix._session is not None
        tags = np.array([1, 2, 3])
        numbers = np.array([8, 1, 1], dtype=np.int32)
        energy, forces, virial = fix._evaluate(tags, np.array(WATER), numbers, np.eye(3))
        assert np.isfinite(energy)
        assert forces.shape == (3, 3)
        assert virial.shape == (3, 3)
        np.testing.assert_array_equal(virial, 0.0)
    finally:
        fix.close()


def test_the_lammps_step_refuses_atom_ids_that_are_not_one_to_n(params):
    fix = fix_without_lammps(params)
    numbers = np.array([8, 1, 1], dtype=np.int32)
    try:
        with pytest.raises(CPMDCError, match="atom IDs 1..N"):
            fix._evaluate(np.array([1, 2, 4]), np.array(WATER), numbers, np.eye(3))
    finally:
        fix.close()
