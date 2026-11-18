"""LAMMPS electrostatic adapter logic without a LAMMPS installation."""

from __future__ import annotations

import numpy as np

from cpmdc import BOHR_TO_ANGSTROM
from cpmdc.lammps import FixExternal


class _Session:
    def __init__(self) -> None:
        self.external = None

    def set_external_charges(self, positions, charges, widths) -> None:
        self.external = (
            np.asarray(positions),
            np.asarray(charges),
            np.asarray(widths),
        )

    def calculate_au(self, positions, numbers, cell, *, forces, virial):
        assert cell is None
        forces[:] = [1.0, 2.0, 3.0]
        virial[:] = 0.0
        return 7.0, forces, virial

    def external_charge_forces_au(self, n_charges):
        assert n_charges == 1
        return np.array([4.0, 5.0, 6.0])


def test_fix_external_splits_qm_and_mm_and_returns_both_forces():
    driver = FixExternal.__new__(FixExternal)
    driver.periodic = False
    driver.length_per_bohr = BOHR_TO_ANGSTROM
    driver.energy_per_hartree = 1.0
    driver.mm_width_angstrom = 0.25
    driver._session = _Session()
    driver._composition = None
    driver._forces = np.zeros(0)
    driver._virial = np.zeros(9)
    driver._library = None
    driver._blob = b""
    tags = np.array([2, 1])
    positions = np.array([[2.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    energy, forces, virial = driver._evaluate(
        tags,
        positions,
        np.array([1, 8], dtype=np.int32),
        np.zeros((3, 3)),
        np.array([-0.4, 0.0]),
        np.array([False, True]),
    )
    assert energy == 7.0
    np.testing.assert_allclose(
        forces,
        np.array([[4.0, 5.0, 6.0], [1.0, 2.0, 3.0]]) / BOHR_TO_ANGSTROM,
    )
    np.testing.assert_allclose(driver._session.external[0], [2.0, 0.0, 0.0])
    np.testing.assert_allclose(driver._session.external[1], [-0.4])
    np.testing.assert_allclose(driver._session.external[2], [0.25])
    ordered_positions = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]) / BOHR_TO_ANGSTROM
    ordered_forces = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    np.testing.assert_allclose(
        virial,
        0.5
        * (
            ordered_positions.T @ ordered_forces
            + ordered_forces.T @ ordered_positions
        ),
    )
