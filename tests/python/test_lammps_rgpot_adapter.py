"""The cpmdc session adapter used by rgpot's generic LAMMPS frontend."""

from __future__ import annotations

import numpy as np
import pytest

from lammps_import import import_rgpot_lammps

if import_rgpot_lammps() is None:
    pytest.skip("rgpot.lammps is not installed", allow_module_level=True)

import cpmdc.lammps
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV


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

    assert len(made) == 1
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


def test_electron_scales_recover_lammps_lengths(monkeypatch):
    """Electron units invert rgpot's length factor back to the LAMMPS values."""
    from rgpot.lammps import unit_factors

    received = {}

    class FakeSession:
        def __init__(self, blob, *, library=None):
            del blob, library

        def calculate_au(self, positions, numbers, cell, *, forces, virial):
            received["positions"] = positions.copy()
            received["cell"] = cell.copy()
            forces[:] = 2.0
            virial[:] = 3.0
            return 4.0, None, None

        def close(self):
            return None

    monkeypatch.setattr(cpmdc.lammps, "Session", FakeSession)
    length_per_angstrom, energy_per_ev = unit_factors("electron")
    potential = cpmdc.lammps._SessionPotential(b"params", None, True)
    potential.set_units(length_per_angstrom, energy_per_ev, 1.0, 1.0)
    potential.open()
    lammps_positions = np.array([[1.25, -0.5, 0.75]])
    lammps_cell = np.diag([8.0, 9.0, 10.0])
    result = potential.evaluate(
        lammps_positions / length_per_angstrom,
        np.array([1], dtype=np.int32),
        lammps_cell / length_per_angstrom,
    )
    np.testing.assert_allclose(received["positions"], lammps_positions)
    np.testing.assert_allclose(received["cell"], lammps_cell)
    assert result["energy"] == pytest.approx(4.0 / energy_per_ev)
    np.testing.assert_allclose(
        result["forces"],
        np.full((1, 3), 2.0) * length_per_angstrom / energy_per_ev,
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


def test_session_adapter_returns_forces_that_do_not_alias_its_buffer():
    potential = _potential(_FillSession(force=1.0))
    positions = np.array([[1.0, 2.0, 3.0]])
    numbers = np.array([1], dtype=np.int32)
    first = potential.evaluate(positions, numbers, np.eye(3))["forces"].copy()
    potential.session.values = (5.0, 2.0, 3.0)
    result = potential.evaluate(positions, numbers, np.eye(3))
    np.testing.assert_array_equal(first, np.full((1, 3), 1.0 * HARTREE_TO_EV / BOHR_TO_ANGSTROM))
    assert not np.shares_memory(result["forces"], potential._buffer)


def test_session_adapter_metal_energy_and_forces_match_the_division_route():
    """Metal units divide by BOHR_TO_ANGSTROM once, so no step rounds differently."""
    rng = np.random.default_rng(11)
    forces = rng.normal(size=(5, 3))

    class Session(_FillSession):
        def calculate_au(self, positions, numbers, cell, *, forces, virial):
            self.positions = positions.copy()
            forces[:] = self.forces.ravel()
            virial[:] = 0.0
            return 0.123456789, None, None

    session = Session()
    session.forces = forces
    potential = _potential(session)
    positions = rng.normal(size=(5, 3))
    result = potential.evaluate(positions, np.arange(1, 6, dtype=np.int32), np.eye(3) * 9.0)
    np.testing.assert_array_equal(session.positions, positions / BOHR_TO_ANGSTROM)
    assert result["energy"] == 0.123456789 * HARTREE_TO_EV
    np.testing.assert_array_equal(
        result["forces"], forces * (HARTREE_TO_EV / BOHR_TO_ANGSTROM)
    )


@pytest.mark.parametrize("name", ["energy", "forces", "virial"])
def test_session_adapter_names_the_non_finite_array_after_the_joint_check(name):
    class Session(_FillSession):
        def calculate_au(self, positions, numbers, cell, *, forces, virial):
            forces[:] = 1.0
            virial[:] = 1.0
            energy = 1.0
            if name == "energy":
                energy = float("nan")
            elif name == "forces":
                forces[0] = float("inf")
            else:
                virial[3] = float("nan")
            return energy, None, None

    with pytest.raises(cpmdc.CPMDCError, match=f"non-finite {name}"):
        _potential(Session()).evaluate(
            np.array([[1.0, 2.0, 3.0]]), np.array([1], dtype=np.int32), np.eye(3)
        )


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
