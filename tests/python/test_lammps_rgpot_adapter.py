"""The cpmdc session adapter used by rgpot's generic LAMMPS frontend."""

from __future__ import annotations

import numpy as np
import pytest

from lammps_import import import_rgpot_lammps

if import_rgpot_lammps() is None:
    pytest.skip("rgpot.lammps is not installed", allow_module_level=True)

import cpmdc.lammps
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, CPMDCError


def test_session_adapter_converts_rgpot_protocol(monkeypatch):
    made = []

    class FakeSession:
        def __init__(self, blob, *, library=None):
            self.blob = blob
            self.library = library
            self.closed = False
            self.inputs = []
            made.append(self)

        def calculate_au(self, positions, numbers, cell, *, forces, virial):
            self.inputs.append((positions.copy(), numbers.copy(), cell.copy()))
            forces[:] = np.arange(forces.size) + 1.0
            virial[:] = np.arange(9) + 2.0
            return 3.0, None, None

        def close(self):
            self.closed = True

    monkeypatch.setattr(cpmdc.lammps, "Session", FakeSession)
    potential = cpmdc.lammps._SessionPotential(b"params", "libcpmdc.so", True)
    potential.open()
    positions = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    numbers = np.array([8, 1], dtype=np.int32)
    cell = np.diag([8.0, 9.0, 10.0])
    result = potential.evaluate(positions, numbers, cell)

    received_positions, received_numbers, received_cell = made[0].inputs[0]
    np.testing.assert_allclose(received_positions, positions / BOHR_TO_ANGSTROM)
    np.testing.assert_array_equal(received_numbers, numbers)
    np.testing.assert_allclose(received_cell, cell / BOHR_TO_ANGSTROM)
    assert result["energy"] == 3.0 * HARTREE_TO_EV
    np.testing.assert_allclose(
        result["forces"],
        np.arange(1.0, 7.0).reshape(2, 3)
        * HARTREE_TO_EV
        / BOHR_TO_ANGSTROM,
    )
    np.testing.assert_allclose(
        result["stress"],
        -np.arange(2.0, 11.0).reshape(3, 3)
        * HARTREE_TO_EV
        / np.linalg.det(cell),
    )

    potential.evaluate(positions, numbers[::-1], cell)
    assert len(made) == 1
    potential.evaluate(positions, np.array([8, 8], dtype=np.int32), cell)
    assert len(made) == 2
    assert made[0].closed


@pytest.mark.parametrize("name", ["energy", "forces", "virial"])
def test_session_adapter_rejects_non_finite_results(monkeypatch, name):
    class FakeSession:
        def __init__(self, *_args, **_kwargs):
            pass

        def calculate_au(self, positions, numbers, cell, *, forces, virial):
            forces[:] = 0.0
            virial[:] = 0.0
            energy = 0.0
            if name == "energy":
                energy = float("nan")
            else:
                {"forces": forces, "virial": virial}[name][0] = float("inf")
            return energy, None, None

        def close(self):
            pass

    monkeypatch.setattr(cpmdc.lammps, "Session", FakeSession)
    potential = cpmdc.lammps._SessionPotential(b"params", None, False)
    with pytest.raises(CPMDCError, match=f"non-finite {name}"):
        potential.evaluate(
            np.array([[1.0, 2.0, 3.0]]),
            np.array([1], dtype=np.int32),
            np.eye(3),
        )


class _FillSession:
    """A session that fills forces and virial with chosen values."""

    def __init__(self, force=1.0, virial=2.0, energy=3.0):
        self.values = (force, virial, energy)
        self.cells = []

    def calculate_au(self, positions, numbers, cell, *, forces, virial):
        self.cells.append(cell)
        forces[:] = self.values[0]
        virial[:] = self.values[1]
        return self.values[2], None, None

    def close(self):
        pass


def _potential(session, periodic=True):
    potential = cpmdc.lammps._SessionPotential(b"params", None, periodic)
    potential.session = session
    return potential


def test_session_adapter_volume_is_the_absolute_determinant():
    rng = np.random.default_rng(7)
    for cell in (np.diag([8.0, 9.0, 10.0]), rng.normal(size=(3, 3)), -np.eye(3)):
        assert cpmdc.lammps._volume(cell) == pytest.approx(abs(np.linalg.det(cell)), rel=1e-13)


def test_session_adapter_accepts_large_finite_results():
    potential = _potential(_FillSession(force=1e200, virial=1e200), periodic=False)
    result = potential.evaluate(
        np.array([[1.0, 2.0, 3.0]]), np.array([1], dtype=np.int32), np.eye(3)
    )
    assert np.all(np.isfinite(result["forces"] / 1e170))


def test_session_adapter_sends_no_cell_without_periodicity():
    session = _FillSession()
    potential = _potential(session, periodic=False)
    result = potential.evaluate(
        np.array([[1.0, 2.0, 3.0]]), np.array([1], dtype=np.int32), np.eye(3)
    )
    assert session.cells == [None]
    assert result["stress"] is None


@pytest.mark.parametrize("style", ["metal", "real", "electron"])
def test_session_adapter_units_follow_the_lammps_style(style):
    from rgpot.lammps import unit_factors as rgpot_unit_factors

    length_per_angstrom, energy_per_ev = rgpot_unit_factors(style)
    length_per_bohr, energy_per_hartree = cpmdc.lammps.unit_factors(style)
    session = _FillSession()
    potential = _potential(session)
    potential.set_units(length_per_angstrom, energy_per_ev, length_per_bohr, energy_per_hartree)
    positions = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    cell = np.diag([8.0, 9.0, 10.0])
    result = potential.evaluate(positions, np.array([8, 1], dtype=np.int32), cell)

    angstrom_per_bohr = length_per_bohr / length_per_angstrom
    native_per_hartree = energy_per_hartree / energy_per_ev
    np.testing.assert_allclose(session.cells[0], cell / angstrom_per_bohr, rtol=1e-12)
    assert result["energy"] == pytest.approx(3.0 * native_per_hartree, rel=1e-12)
    np.testing.assert_allclose(
        result["forces"],
        np.full((2, 3), native_per_hartree / angstrom_per_bohr),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        result["stress"], np.full((3, 3), -2.0 * native_per_hartree / 720.0), rtol=1e-12
    )
