"""The C and Python cell classifiers decide every cell the same way.

They also compute the same volume, bit for bit, by LU with partial
pivoting on the cell scaled by a power of two.

``cpmdc._rank`` and ``src/cpmdc_cell_rank.c`` run the same one-sided
Jacobi iteration with the same binary64 operations, so they must agree on
every bit pattern, including cells whose smallest singular value sits
within rounding of the ``numpy.linalg.matrix_rank`` cutoff, where two
different SVD routines can disagree. Meson builds ``cell_rank_probe`` from
the C file and names it in ``CPMDC_CELL_RANK_PROBE``; otherwise
``build/cell_rank_probe`` of this checkout is used, and the tests skip when
neither exists.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from conftest import ROOT

from cpmdc._rank import _FULL_RANK, _cell_state, _cell_volume

_EPS = sys.float_info.epsilon
_TINY = sys.float_info.min
_HUGE = sys.float_info.max


def _probe_path() -> Path | None:
    named = os.environ.get("CPMDC_CELL_RANK_PROBE")
    candidate = Path(named) if named else ROOT / "build" / "cell_rank_probe"
    return candidate if candidate.is_file() else None


PROBE = _probe_path()
pytestmark = pytest.mark.skipif(PROBE is None, reason="cell_rank_probe is not built")


def _c_results(cells: list[np.ndarray]) -> list[tuple[int, int, float]]:
    """``(state, volume return code, volume)`` from the C helpers."""
    assert PROBE is not None
    text = "".join(
        " ".join(float(value).hex() for value in np.asarray(cell, dtype=float).ravel()) + "\n"
        for cell in cells
    )
    done = subprocess.run(
        [str(PROBE)], input=text, capture_output=True, text=True, check=True, timeout=60
    )
    results = []
    for line in done.stdout.splitlines():
        state, rc, volume = line.split()
        results.append((int(state), int(rc), float.fromhex(volume)))
    return results


def _c_states(cells: list[np.ndarray]) -> list[int]:
    return [state for state, _, _ in _c_results(cells)]


def _python_volume(cell: np.ndarray) -> tuple[int, float]:
    try:
        volume = _cell_volume(cell)
    except ValueError:
        return -1, 0.0
    return (0, volume) if math.isfinite(volume) else (-1, volume)


def _assert_volumes_match(cells: list[np.ndarray]) -> None:
    for cell, (_, rc, volume) in zip(cells, _c_results(cells), strict=True):
        assert (rc, volume) == _python_volume(cell)


def _rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    return q * np.sign(np.diag(r))


def _rotated(rng: np.random.Generator, ratio: float) -> np.ndarray:
    """A rotated cell whose smallest singular value is ``ratio`` cutoffs."""
    top = 10.0 ** rng.uniform(-1.0, 2.0)
    middle = top * 10.0 ** rng.uniform(-6.0, 0.0)
    bottom = ratio * 3.0 * _EPS * top
    return _rotation(rng) @ np.diag([top, middle, bottom]) @ _rotation(rng)


def _numpy_ratio(cell: np.ndarray) -> float:
    singular = np.linalg.svd(cell, compute_uv=False)
    return float(singular[-1] / (singular[0] * 3.0 * _EPS))


# A seeded rotated cell that numpy.linalg.matrix_rank ranks 3, with the
# smallest singular value close to its cutoff, that the C estimate used to
# classify as rank deficient while cpmdc._rank accepted it.
_REPORTED_CLASS = np.array(
    [
        [float.fromhex(h) for h in ("0x1.5bedaba4ff1ccp+4", "-0x1.e7182853e964fp+2", "-0x1.7d099f8b90e2dp+3")],
        [float.fromhex(h) for h in ("-0x1.e113f63bcb85bp+4", "0x1.516296171d335p+3", "0x1.07a385a05a1c0p+4")],
        [float.fromhex(h) for h in ("0x1.e328d6a33c332p+5", "-0x1.521d04262f37ap+4", "-0x1.088994a9772c0p+5")],
    ]
)


def test_cell_that_split_the_layers_is_full_rank_in_both() -> None:
    assert np.linalg.matrix_rank(_REPORTED_CLASS) == 3
    # LAPACK implementations need not return identical singular values this
    # close to the cutoff. The exact cell bits and the shared classifier are
    # the cross-platform contract.
    assert 1.0 < _numpy_ratio(_REPORTED_CLASS) < 1.25
    assert _cell_state(_REPORTED_CLASS) == _FULL_RANK
    assert _c_states([_REPORTED_CLASS]) == [_FULL_RANK]


@pytest.mark.parametrize(
    ("seed", "low", "high", "straddles"),
    [
        (110, 0.9, 1.0, False),
        (111, 1.0, 1.1, False),
        (112, 1.1, 1.2, False),
        (113, 0.5, 1.5, True),
    ],
    ids=["just-below", "just-above", "reported-band", "wide"],
)
def test_rotated_cells_near_the_cutoff_get_one_answer(
    seed: int, low: float, high: float, straddles: bool
) -> None:
    rng = np.random.default_rng(seed)
    cells = [_rotated(rng, rng.uniform(low, high)) for _ in range(2000)]
    python = [_cell_state(cell) for cell in cells]
    assert _c_states(cells) == python
    _assert_volumes_match(cells)
    # QR and matrix multiplication are supplied by the platform BLAS, so a
    # nominal band wholly on one side of the cutoff need not produce both
    # answers on every architecture. The wide band genuinely straddles it.
    if straddles:
        assert {state == _FULL_RANK for state in python} == {True, False}


@pytest.mark.parametrize(("ratio", "full_rank"), [(0.25, False), (4.0, True)])
def test_cells_well_away_from_the_cutoff_follow_matrix_rank(ratio: float, full_rank: bool) -> None:
    rng = np.random.default_rng(114)
    cells = [_rotated(rng, ratio) for _ in range(500)]
    for cell, state in zip(cells, _c_states(cells), strict=True):
        assert bool(np.linalg.matrix_rank(cell) == 3) is full_rank
        assert (state == _FULL_RANK) is full_rank
        assert _cell_state(cell) == state


def test_power_of_two_scaling_never_changes_the_answer() -> None:
    rng = np.random.default_rng(115)
    base = [_rotated(rng, rng.uniform(0.9, 1.1)) for _ in range(300)]
    expected = [_cell_state(cell) for cell in base]
    for power in (-300, -60, 60, 300):
        scaled = [np.ldexp(cell, power) for cell in base]
        assert [_cell_state(cell) for cell in scaled] == expected
        assert _c_states(scaled) == expected


def test_named_cells_get_one_answer() -> None:
    near = 1.0 + 28.0 * _EPS
    cells = [
        np.zeros((3, 3)),
        np.diag([3.0, 4.0, 5.0]),
        np.diag([3.0, 4.0, -5.0]),
        np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
        np.diag([1.0, 1.0, 3.0 * _EPS]),
        np.diag([1.0, 1.0, 4.0 * _EPS]),
        np.diag([np.nan, 1.0, 1.0]),
        np.diag([np.inf, 1.0, 1.0]),
        np.diag([_HUGE / 2.0] * 3),
        np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]]) * 1.0e308,
        np.diag([np.nextafter(0.0, 1.0)] * 3),
        np.diag([1.0e-310] * 3),
        np.diag([3.1072325e-103] * 3),
        np.diag([2.8e-103] * 3),
        np.array([[1.0, 1.0, 1.0], [1.0, near, 1.0], [1.0, 1.0, near]]),
    ]
    assert _c_states(cells) == [_cell_state(cell) for cell in cells]
    assert [_cell_state(cell) for cell in cells] == [
        0, 1, 1, -1, -1, 1, -2, -2, 1, 1, -3, -3, 1, -3, 1,
    ]
    assert abs(np.prod(np.diag(cells[12]))) / _TINY == pytest.approx(1.348, abs=1e-3)
    _assert_volumes_match(cells)


def test_volume_of_a_cancelling_cell_is_exact() -> None:
    """The three-term determinant of this accepted cell cancels to zero."""
    near = 1.0 + 28.0 * _EPS
    cell = np.array([[1.0, 1.0, 1.0], [1.0, near, 1.0], [1.0, 1.0, near]])
    assert _cell_state(cell) == _FULL_RANK
    expected = (28.0 * _EPS) ** 2
    assert expected == pytest.approx(3.865418435582966e-29, rel=1e-15)
    assert float(np.linalg.det(cell)) == pytest.approx(expected, rel=1e-12)
    assert _cell_volume(cell) == expected
    assert _c_results([cell])[0][1:] == (0, expected)


def test_volume_above_the_largest_float_is_refused() -> None:
    cell = np.diag([_HUGE / 2.0] * 3)
    assert _cell_state(cell) == _FULL_RANK
    assert _cell_volume(cell) == math.inf
    assert _c_results([cell])[0][1:] == (-1, math.inf)
    assert _c_results([np.diag([2.0**341] * 3)])[0][1:] == (0, 2.0**1023)
    assert _cell_volume(np.diag([2.0**341] * 3)) == 2.0**1023
