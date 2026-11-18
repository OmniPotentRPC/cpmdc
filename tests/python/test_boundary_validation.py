"""Arguments are checked at the Python boundary, with a message that names the argument."""

from __future__ import annotations

from array import array

import pytest
from conftest import WATER

import cpmdc
from cpmdc import CPMDCError, Session

BOHR = cpmdc.BOHR_TO_ANGSTROM
WATER_BOHR = [[c / BOHR for c in row] for row in WATER]


@pytest.fixture
def session(params):
    with Session(params) as handle:
        yield handle


@pytest.mark.parametrize(
    ("numbers", "error", "fragment"),
    [
        ([8.5, 1, 1], ValueError, "atomic_numbers must be an integer"),
        ([True, 1, 1], TypeError, "not bool"),
        (["O", "H", "H"], TypeError, "got str"),
        ([None, 1, 1], TypeError, "got NoneType"),
        ([0, 1, 1], ValueError, "1..118"),
        ([8, 1, 119], ValueError, "1..118"),
        ([8, 1, 2**40], ValueError, "1..118"),
        ([float("nan"), 1, 1], ValueError, "must be an integer"),
    ],
)
def test_atomic_numbers_are_checked(session, numbers, error, fragment):
    with pytest.raises(error, match=fragment):
        session.energy_forces(WATER, numbers)
    with pytest.raises(error, match=fragment):
        session.calculate_au(WATER_BOHR, numbers)


def test_whole_floats_are_atomic_numbers(session):
    done = session.energy_forces(WATER, [8.0, 1.0, 1.0])
    assert len(done.forces_h_bohr) == 3


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_au_positions_and_cell_must_be_finite(session, value):
    bad = [list(row) for row in WATER_BOHR]
    bad[1][2] = value
    with pytest.raises(ValueError, match="positions_bohr must be finite"):
        session.calculate_au(bad, [8, 1, 1])
    with pytest.raises(ValueError, match="positions_bohr must be finite"):
        session.calculate_au(array("d", [c for row in bad for c in row]), [8, 1, 1])
    cell = [20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, value]
    with pytest.raises(ValueError, match="cell_bohr must be finite"):
        session.calculate_au(WATER_BOHR, [8, 1, 1], cell)


def test_au_values_must_be_numbers(session):
    with pytest.raises(TypeError, match="positions_bohr must contain real numbers, got str"):
        session.calculate_au([["0", "0", "0"]] * 3, [8, 1, 1])


def test_au_buffer_atomic_numbers_are_checked(session):
    with pytest.raises(ValueError, match=r"atomic_numbers must be in 1\.\.118"):
        session.calculate_au(WATER_BOHR, array("i", [8, 1, 200]))


@pytest.mark.parametrize(
    ("key", "error", "fragment"),
    [
        (1.5, ValueError, "key must be an integer"),
        (True, TypeError, "not bool"),
        ("3", TypeError, "got str"),
        (None, TypeError, "got NoneType"),
        (2**63, ValueError, "signed 64-bit"),
        (-(2**63) - 1, ValueError, "signed 64-bit"),
    ],
)
def test_orbital_keys_are_checked(session, key, error, fragment):
    with pytest.raises(error, match=fragment):
        session.select_orbitals(key)
    with pytest.raises(error, match=fragment):
        session.copy_orbitals(key, 1)
    with pytest.raises(error, match=fragment):
        session.copy_orbitals(1, key)


def test_orbital_key_extremes_reach_the_library(session):
    for key in (-(2**63), 2**63 - 1, 3.0):
        session.select_orbitals(key)


def test_a_rejected_argument_does_not_poison_the_session(session):
    with pytest.raises(ValueError):
        session.calculate_au(WATER_BOHR, [0, 1, 1])
    energy, forces, _ = session.calculate_au(WATER_BOHR, [8, 1, 1])
    assert len(forces) == 9
    assert energy == energy


def test_a_closed_session_still_reports_closed(params):
    handle = Session(params)
    handle.close()
    with pytest.raises(CPMDCError, match="closed"):
        handle.select_orbitals(1)
