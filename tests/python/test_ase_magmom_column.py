"""One magnetic moment per atom includes the column ASE stores.

ASE 3.29 ``set_initial_magnetic_moments`` keeps a column of shape ``(N, 1)``
as one number per atom. A non-collinear moment is shape ``(N, 3)``.
"""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

from cpmdc.ase import _initial_moment


def water() -> Atoms:
    return Atoms(
        "OHH",
        positions=[
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ],
    )


def test_a_one_component_magmom_column_is_collinear():
    atoms = water()
    atoms.set_initial_magnetic_moments(np.array([[1.0], [0.0], [1.0]]))
    assert atoms.get_initial_magnetic_moments().shape == (3, 1)
    assert _initial_moment(atoms) == 2


def test_a_non_integral_magmom_column_names_the_sum():
    atoms = water()
    atoms.set_initial_magnetic_moments(np.array([[0.2], [0.2], [0.0]]))
    with pytest.raises(ValueError, match="not an integer"):
        _initial_moment(atoms)


def test_a_zero_magmom_column_is_not_a_moment():
    atoms = water()
    atoms.set_initial_magnetic_moments(np.zeros((3, 1)))
    assert _initial_moment(atoms) == 0


def test_a_non_collinear_moment_still_raises():
    atoms = water()
    atoms.set_initial_magnetic_moments(
        [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )
    with pytest.raises(ValueError, match="non-collinear"):
        _initial_moment(atoms)


def test_zero_non_collinear_moments_are_not_a_moment():
    atoms = water()
    atoms.set_initial_magnetic_moments(np.zeros((3, 3)))
    assert _initial_moment(atoms) == 0
