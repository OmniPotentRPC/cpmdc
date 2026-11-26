"""Electrostatic embedding in the LAMMPS adapter, without a LAMMPS installation.

A fake ``Session`` records what the adapter sends, so the tests cover the QM
and MM split, the force assembly, the clearing of stale charges, and the width
bound. The constructed ``FixExternal`` is driven through its callback with a
stand-in LAMMPS object.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any, ClassVar

import numpy as np
import pytest
from lammps_import import import_rgpot_lammps

if import_rgpot_lammps() is None:
    pytest.skip("rgpot.lammps is not installed", allow_module_level=True)

import cpmdc.lammps as adapter
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, CPMDCError
from cpmdc.lammps import FixExternal


class FakeSession:
    """Records the charge sets and answers with fixed forces."""

    made: ClassVar[list[FakeSession]] = []

    def __init__(self, blob: bytes, *, library: Any = None) -> None:
        self.blob = blob
        self.charge_calls: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        self.qm_calls: list[tuple[np.ndarray, np.ndarray]] = []
        self.n_mm = 0
        self.closed = False
        FakeSession.made.append(self)

    def set_external_charges(self, positions, charges, widths=None) -> None:
        self.charge_calls.append(
            (np.array(positions), np.array(charges), np.array(widths))
        )
        self.n_mm = len(charges)

    def calculate_au(self, positions, numbers, cell, *, forces, virial):
        assert cell is None
        self.qm_calls.append((np.array(positions), np.array(numbers)))
        forces[:] = np.arange(1, forces.size + 1, dtype=float)
        virial[:] = 0.0
        return 2.0, None, None

    def external_charge_forces_au(self, n_charges):
        assert n_charges == self.n_mm
        return np.arange(1, 3 * n_charges + 1, dtype=float) * 10.0

    def close(self) -> None:
        self.closed = True


@pytest.fixture(autouse=True)
def fake_session(monkeypatch):
    FakeSession.made = []
    monkeypatch.setattr(adapter, "Session", FakeSession)
    return FakeSession


class StandInLammps:
    """The slice of the LAMMPS Python object the fix reads, with settable atoms."""

    SETTINGS: ClassVar[dict[str, int]] = {
        "world_size": 1,
        "IMGMASK": 1023,
        "IMGMAX": 512,
        "IMGBITS": 10,
        "IMG2BITS": 20,
    }

    def __init__(self, periodic: bool = False) -> None:
        self.periodic = periodic
        self.callbacks: list[Any] = []
        self.atoms = {
            "type": np.array([1, 2, 3, 3]),
            "q": np.array([0.0, 0.0, -0.4, 0.4]),
        }
        self.energy = None
        self.numpy = SimpleNamespace(extract_atom=lambda name: self.atoms.get(name))

    def extract_box(self):
        return ([0.0] * 3, [10.0] * 3, 0.0, 0.0, 0.0, [self.periodic] * 3, 0)

    def extract_global(self, name):
        return {"units": "metal", "ntypes": 3}[name]

    def extract_setting(self, name):
        return self.SETTINGS[name]

    def command(self, text) -> None:
        pass

    def set_fix_external_callback(self, fix_id, callback, caller) -> None:
        self.callbacks.append(callback)

    def set_internal_variable(self, name, value) -> None:
        pass

    def fix_external_set_energy_global(self, fix_id, energy) -> None:
        self.energy = energy

    def fix_external_set_virial_global(self, fix_id, virial) -> None:
        self.virial = list(virial)


def _step(fix, positions, natoms=None):
    natoms = len(positions) if natoms is None else natoms
    forces = np.zeros((natoms, 3))
    tags = np.arange(1, natoms + 1)
    fix._callback(None, 0, natoms, tags, np.array(positions[:natoms], dtype=float), forces)
    assert fix.error is None, fix.error
    return forces


POSITIONS = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]
ELEMENTS = "O H He"


def test_potential_splits_qm_and_mm_and_returns_both_forces():
    potential = adapter._SessionPotential(b"p", None, False)
    potential.mm_width_angstrom = 0.25
    potential.open()
    potential.pending_charges = np.array([0.0, -0.4, 0.4])
    potential.pending_qm = np.array([True, False, False])
    positions = np.array(POSITIONS[:3])
    result = potential.evaluate(positions, np.array([8, 1, 1], dtype=np.int32), np.eye(3))
    session = potential.session
    np.testing.assert_allclose(session.qm_calls[0][0], positions[:1] / BOHR_TO_ANGSTROM)
    np.testing.assert_array_equal(session.qm_calls[0][1], [8])
    mm_positions, mm_charges, mm_widths = session.charge_calls[0]
    np.testing.assert_allclose(mm_positions, positions[1:].reshape(-1))
    np.testing.assert_allclose(mm_charges, [-0.4, 0.4])
    np.testing.assert_allclose(mm_widths, [0.25, 0.25])
    scale = HARTREE_TO_EV / BOHR_TO_ANGSTROM
    expected_au = np.vstack([[1.0, 2.0, 3.0], np.arange(1, 7).reshape(2, 3) * 10.0])
    np.testing.assert_allclose(result["forces"], expected_au * scale)
    assert result["energy"] == pytest.approx(2.0 * HARTREE_TO_EV)
    assert result["stress"] is None


def test_potential_refuses_mm_atoms_on_a_periodic_box():
    potential = adapter._SessionPotential(b"p", None, True)
    potential.open()
    potential.pending_charges = np.array([0.0, 0.5])
    potential.pending_qm = np.array([True, False])
    with pytest.raises(CPMDCError, match="periodic electrostatic embedding"):
        potential.evaluate(
            np.zeros((2, 3)), np.array([8, 1], dtype=np.int32), np.eye(3) * 8.0
        )


def test_potential_refuses_a_step_without_qm_atoms():
    potential = adapter._SessionPotential(b"p", None, False)
    potential.open()
    potential.pending_charges = np.array([0.5, 0.5])
    potential.pending_qm = np.array([False, False])
    with pytest.raises(CPMDCError, match="no QM atoms"):
        potential.evaluate(
            np.zeros((2, 3)), np.array([8, 1], dtype=np.int32), np.eye(3)
        )


def _delete_mm_atoms(lmp):
    """The two MM atoms leave the system; the two QM atoms stay."""
    lmp.atoms["type"] = lmp.atoms["type"][:2]
    lmp.atoms["q"] = lmp.atoms["q"][:2]


def test_fix_clears_stale_charges_when_the_mm_set_goes_empty():
    lmp = StandInLammps()
    fix = FixExternal(lmp, b"p", ELEMENTS, qm_types=[1, 2], mm_width_angstrom=0.2)
    try:
        forces = _step(fix, POSITIONS)
        session = fix._session
        assert len(session.charge_calls) == 1
        assert session.n_mm == 2
        assert np.any(forces[2:] != 0.0)

        _delete_mm_atoms(lmp)
        forces = _step(fix, POSITIONS, natoms=2)
        assert fix._session is session
        assert len(session.charge_calls) == 2
        cleared_positions, cleared_charges, cleared_widths = session.charge_calls[1]
        assert cleared_positions.size == 0
        assert cleared_charges.size == 0
        assert cleared_widths.size == 0
        assert session.n_mm == 0
        assert len(session.qm_calls[1][1]) == 2
        assert forces.shape == (2, 3)

        _step(fix, POSITIONS, natoms=2)
        assert len(session.charge_calls) == 2
    finally:
        fix.close()


def test_fix_sends_the_charges_again_after_they_were_cleared():
    lmp = StandInLammps()
    fix = FixExternal(lmp, b"p", ELEMENTS, qm_types=[1, 2])
    try:
        _step(fix, POSITIONS)
        session = fix._session
        full_types, full_q = lmp.atoms["type"].copy(), lmp.atoms["q"].copy()
        _delete_mm_atoms(lmp)
        _step(fix, POSITIONS, natoms=2)
        lmp.atoms["type"], lmp.atoms["q"] = full_types, full_q
        _step(fix, POSITIONS)
        assert len(session.charge_calls) == 3
        assert session.n_mm == 2
        np.testing.assert_allclose(session.charge_calls[2][1], [-0.4, 0.4])
    finally:
        fix.close()


def test_fix_starts_a_new_session_without_charges_when_the_composition_changes():
    lmp = StandInLammps()
    fix = FixExternal(lmp, b"p", ELEMENTS, qm_types=[1, 2])
    try:
        _step(fix, POSITIONS)
        first = fix._session
        lmp.atoms["type"] = np.array([1, 1, 1, 1])
        lmp.atoms["q"] = np.zeros(4)
        _step(fix, POSITIONS)
        second = fix._session
        assert second is not first
        assert first.closed
        assert second.charge_calls == []
        assert len(second.qm_calls[0][1]) == 4
    finally:
        fix.close()


def test_fix_without_qm_types_never_sends_charges():
    lmp = StandInLammps()
    fix = FixExternal(lmp, b"p", ELEMENTS)
    try:
        _step(fix, POSITIONS)
        assert fix._session.charge_calls == []
        assert len(fix._session.qm_calls[0][1]) == 4
    finally:
        fix.close()


def test_fix_reports_missing_charges():
    lmp = StandInLammps()
    lmp.atoms.pop("q")
    fix = FixExternal(lmp, b"p", ELEMENTS, qm_types=[1, 2])
    try:
        fix._callback(None, 0, 4, np.arange(1, 5), np.array(POSITIONS), np.zeros((4, 3)))
        assert isinstance(fix.error, CPMDCError)
        assert "provide charges" in str(fix.error)
    finally:
        fix.error = None
        fix.close()


@pytest.mark.parametrize("width", [0.0, -1.0, 1.0e-13, float("nan"), float("inf")])
def test_fix_rejects_an_unusable_mm_width(width):
    with pytest.raises(ValueError, match="mm_width_angstrom"):
        FixExternal(StandInLammps(), b"p", ELEMENTS, qm_types=[1, 2], mm_width_angstrom=width)


def test_fix_accepts_the_smallest_width_the_library_allows():
    fix = FixExternal(StandInLammps(), b"p", ELEMENTS, qm_types=[1, 2], mm_width_angstrom=1.0e-12)
    fix.close()


def test_fix_rejects_empty_qm_types():
    with pytest.raises(ValueError, match="qm_types"):
        FixExternal(StandInLammps(), b"p", ELEMENTS, qm_types=[])
