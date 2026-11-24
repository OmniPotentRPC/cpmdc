"""Input packing of Session: every accepted form gives the same step.

The packing helpers take a fast route for C-contiguous ``float64`` buffers
and for lists and tuples of plain ``float`` and ``int``. Anything else goes
through the checked loops. These tests hold the two routes to one result and
to one set of refusals.
"""

from __future__ import annotations

import math
from array import array

import numpy as np
import pytest
from conftest import WATER

import cpmdc
from cpmdc import Session

NUMBERS = [8, 1, 1]


def _flat(rows):
    return [float(coord) for row in rows for coord in row]


def _frames():
    """The geometry of water in each container the packing accepts."""
    array_rows = np.array(WATER)
    readonly = array_rows.copy()
    readonly.setflags(write=False)
    return {
        "lists": [list(row) for row in WATER],
        "tuples": tuple(tuple(row) for row in WATER),
        "ints": [[0, 0, 0], [0, 1, 0], [0, -1, 0]],
        "ndarray": array_rows,
        "readonly ndarray": readonly,
        "fortran ndarray": np.asfortranarray(array_rows),
        "reversed ndarray": array_rows[::-1][::-1],
        "float32 ndarray": array_rows.astype(np.float32),
        "numpy scalars": [[np.float64(coord) for coord in row] for row in WATER],
        "bools": [[True, False, True]] * 3,
        "strings": [["0.0", "0.0", "0.1173"]] * 3,
    }


@pytest.mark.parametrize("name", sorted(_frames()))
def test_positions_pack_to_the_same_doubles(name):
    rows = _frames()[name]
    count, flat = cpmdc._positions(rows)
    assert count == 3
    assert flat.typecode == "d"
    assert flat.tolist() == _flat(np.asarray(rows, dtype=float) if name == "strings" else rows)


def test_positions_do_not_alias_the_caller_array():
    rows = np.array(WATER)
    _, flat = cpmdc._positions(rows)
    rows[0, 0] = 99.0
    assert flat[0] == WATER[0][0]


@pytest.mark.parametrize(
    ("rows", "error", "text"),
    [
        ([], ValueError, "at least one atom"),
        ([[0.0, 1.0]], ValueError, "position 0 must have three coordinates"),
        ([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0, 3.0]], ValueError, "position 1 must have three coordinates"),
        ([[0.0, math.nan, 2.0]], ValueError, "positions must be finite"),
        ([[0.0, math.inf, 2.0]], ValueError, "positions must be finite"),
        (np.array([[0.0, math.nan, 2.0]]), ValueError, "positions must be finite"),
        (np.array([[0.0, -math.inf, 2.0]]), ValueError, "positions must be finite"),
        (np.zeros((0, 3)), ValueError, "at least one atom"),
        (np.zeros(6), TypeError, ""),
        ([[None, 1.0, 2.0]], TypeError, ""),
    ],
)
def test_positions_refuse_what_the_checked_loop_refuses(rows, error, text):
    with pytest.raises(error, match=text or None):
        cpmdc._positions(rows)


def test_positions_accept_large_finite_values_whose_sum_overflows():
    big = [[1.0e308, 1.0e308, 1.0e308], [1.0e308, 1.0e308, 1.0e308]]
    for rows in (big, np.array(big)):
        count, flat = cpmdc._positions(rows)
        assert count == 2
        assert flat.tolist() == [1.0e308] * 6


@pytest.mark.parametrize(
    "numbers",
    [
        [8, 1, 1],
        (8, 1, 1),
        [8.0, 1.0, 1.0],
        np.array([8, 1, 1], dtype=np.int32),
        np.array([8, 1, 1], dtype=np.int64),
        array("i", [8, 1, 1]),
        [np.int64(8), 1, 1],
    ],
)
def test_numbers_pack_to_int32(numbers):
    packed = cpmdc._numbers(numbers, 3)
    assert packed.typecode == "i"
    assert packed.tolist() == NUMBERS


@pytest.mark.parametrize(
    ("numbers", "error", "text"),
    [
        ([0, 1, 1], ValueError, "atomic numbers must be in 1..118, got 0"),
        ([8, 1, 119], ValueError, "atomic numbers must be in 1..118, got 119"),
        ([8, -1, 1], ValueError, "got -1"),
        ([8, 1], ValueError, "one entry per atom"),
        ([8, 1, 1, 1], ValueError, "one entry per atom"),
        ([8, 1.5, 1], ValueError, "must be an integer"),
        ([True, 1, 1], TypeError, "not bool"),
        (["8", 1, 1], TypeError, "must be an integer"),
        (np.array([8, 0, 1], dtype=np.int32), ValueError, "got 0"),
        (np.array([8, 300, 1], dtype=np.int32), ValueError, "got 300"),
        ([1 << 40, 1, 1], ValueError, "got 1099511627776"),
        ([[8], [1], [1]], TypeError, "must be an integer"),
        ([(8,), (1,), (1,)], TypeError, "must be an integer"),
    ],
)
def test_numbers_refuse_what_the_checked_loop_refuses(numbers, error, text):
    with pytest.raises(error, match=text):
        cpmdc._numbers(numbers, 3)


def test_force_input_is_little_endian_whatever_the_host(monkeypatch):
    """A big-endian host holds the swapped bytes; the message is still little-endian."""
    flat = array("d", [0.5, -1.25, 2.0])
    numbers = array("i", [8, 1, 1])
    expected = cpmdc._force_input_message(flat, numbers, None)
    monkeypatch.setattr(cpmdc.sys, "byteorder", "big")
    swapped = array("d", flat)
    swapped.byteswap()
    swapped_numbers = array("i", numbers)
    swapped_numbers.byteswap()
    assert cpmdc._little_endian(swapped) == flat.tobytes()
    assert cpmdc._little_endian(swapped_numbers) == numbers.tobytes()
    assert expected == cpmdc._force_input_message(swapped, swapped_numbers, None)


def _frames_for_evaluation():
    rows = np.array(WATER)
    return [rows, rows.tolist(), tuple(map(tuple, rows.tolist()))]


def test_every_container_gives_the_same_evaluation(periodic_params):
    cell = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    with Session(periodic_params) as session:
        reference = session.energy_forces(WATER, NUMBERS, cell)
        for frame in _frames_for_evaluation():
            for numbers in (NUMBERS, np.array(NUMBERS, dtype=np.int32), tuple(NUMBERS)):
                for box in (cell, np.array(cell), np.array(cell).reshape(-1)):
                    done = session.energy_forces(frame, numbers, box)
                    assert done.energy_h == reference.energy_h
                    assert done.forces_h_bohr == reference.forces_h_bohr
                    assert done.stress_h_bohr3 == reference.stress_h_bohr3


def test_forces_are_rows_of_python_floats(params):
    with Session(params) as session:
        done = session.energy_forces(WATER, NUMBERS)
    assert len(done.forces_h_bohr) == 3
    for row in done.forces_h_bohr:
        assert type(row) is list
        assert len(row) == 3
        assert all(type(value) is float for value in row)


def test_calculate_au_takes_every_container(params):
    bohr = np.array(WATER) / cpmdc.BOHR_TO_ANGSTROM
    flat = bohr.reshape(-1)
    with Session(params) as session:
        energy, forces, virial = session.calculate_au(flat, np.array(NUMBERS, dtype=np.int32))
        for positions, numbers in (
            (flat.tolist(), NUMBERS),
            (tuple(flat.tolist()), tuple(NUMBERS)),
            (bohr.tolist(), NUMBERS),
            (array("d", flat.tolist()), array("i", NUMBERS)),
        ):
            again = session.calculate_au(positions, numbers)
            assert again[0] == energy
            assert list(again[1]) == list(forces)
            assert list(again[2]) == list(virial)


@pytest.mark.parametrize(
    ("positions", "numbers", "error", "text"),
    [
        ([0.0, math.nan, 0.0] * 3, NUMBERS, ValueError, "positions_bohr must be finite"),
        ([0.0] * 8, NUMBERS, ValueError, "positions_bohr must have 9 values, got 8"),
        ([0.0] * 9, [8, 1, 0], ValueError, "atomic numbers must be in 1..118, got 0"),
        ([0.0] * 9, [8, 1, 1.5], ValueError, "must be an integer"),
        ([True] * 9, NUMBERS, TypeError, "positions_bohr must contain real numbers"),
        (["1.0"] * 9, NUMBERS, TypeError, "positions_bohr must contain real numbers"),
    ],
)
def test_calculate_au_refuses_bad_lists(params, positions, numbers, error, text):
    with Session(params) as session:
        with pytest.raises(error, match=text):
            session.calculate_au(positions, numbers)
