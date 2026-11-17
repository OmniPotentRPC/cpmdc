"""The cpmdc session adapter used by rgpot's generic LAMMPS frontend."""

from __future__ import annotations

import numpy as np

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
