"""A periodic cell is accepted when its lattice vectors have rank 3.

The cutoff is the one ``numpy.linalg.matrix_rank`` uses for a float64
matrix: a singular value counts when it is strictly greater than
``sigma_max * 3 * eps``. ``eps`` is the float64 unit roundoff. The
retired check compared the absolute volume with ``1e-8`` times the
product of the edge lengths.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms

from cpmdc import BOHR_TO_ANGSTROM
from cpmdc._rank import _cell_is_rank_deficient
from cpmdc.ase import CPMDCCalculator
from cpmdc.ipi import Driver, ProtocolError

_ASE = (
    "cpmdc takes a fully periodic cell; atoms.pbc is true but "
    "atoms.cell does not have three independent vectors"
)
_IPI = (
    "the server sent a cell without three independent vectors; send an "
    "all-zero cell or start cpmdc-ipi with --cluster for an isolated system"
)
_EPS = float(np.finfo(np.float64).eps)
# max(m, n) for the 3 by 3 cell, as numpy.linalg.matrix_rank applies it.
_RANK_ORDER = 3
# Retired volume ratio in CPMDCCalculator._periodic_cell and Driver._cell_rows.
_RETIRED_VOLUME_RATIO = 1.0e-8


def _orthorhombic() -> np.ndarray:
    return np.diag([3.0, 4.0, 5.0])


def _left_handed() -> np.ndarray:
    return np.diag([3.0, 4.0, -5.0])


def _skewed_triclinic() -> np.ndarray:
    """Fifteen-degree shear plus a tilt, well above either cutoff."""
    angle = math.radians(15.0)
    return np.array(
        [
            [4.0, 0.0, 0.0],
            [5.0 * math.cos(angle), 5.0 * math.sin(angle), 0.0],
            [0.4, -0.7, 6.0],
        ]
    )


def _skewed_below_retired_ratio() -> np.ndarray:
    """Triclinic shear between the rank cutoff and the retired volume ratio."""
    shear = math.sqrt(_RETIRED_VOLUME_RATIO * _RANK_ORDER * _EPS)
    return np.array([[1.0, 0.0, 0.0], [0.3, shear, 0.0], [0.2, 0.1, 0.8]])


def _coplanar() -> np.ndarray:
    return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])


def _zero_vector() -> np.ndarray:
    return np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]])


def _near_collinear(relative: float) -> np.ndarray:
    """Two finite-length vectors with the requested smallest singular value."""
    root_half = math.sqrt(0.5)
    rotation = np.array(
        [
            [root_half, 0.0, root_half],
            [0.0, 1.0, 0.0],
            [-root_half, 0.0, root_half],
        ]
    )
    return rotation @ np.diag([1.0, 1.0, relative])


def _near_coplanar(relative: float) -> np.ndarray:
    """Three finite-length vectors close to one plane."""
    rotation, _ = np.linalg.qr(
        np.array([[1.0, 2.0, 3.0], [2.0, -1.0, 1.0], [1.0, 1.0, -1.0]])
    )
    return rotation @ np.diag([1.0, 0.75, relative])


def _retired_volume_rejects(cell: np.ndarray) -> bool:
    lengths = np.linalg.norm(cell, axis=1)
    volume = abs(float(np.linalg.det(cell)))
    return bool(np.any(lengths <= 0.0) or volume <= _RETIRED_VOLUME_RATIO * float(np.prod(lengths)))


def _ase(cell: np.ndarray) -> np.ndarray | None:
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
    return CPMDCCalculator._periodic_cell(atoms)


def _wire(rows: np.ndarray) -> list[float]:
    flat = (np.asarray(rows, dtype=float).T / BOHR_TO_ANGSTROM).ravel(order="C")
    return [float(value) for value in flat]


def _ipi(cell: np.ndarray) -> list[list[float]] | None:
    rows = Driver._cell_rows(SimpleNamespace(cluster=False), _wire(cell))
    return rows


@pytest.mark.parametrize(
    "cell",
    [_orthorhombic(), _left_handed(), _skewed_triclinic()],
    ids=["orthorhombic", "left-handed", "skewed-triclinic"],
)
def test_full_rank_cells_are_accepted(cell: np.ndarray) -> None:
    assert np.linalg.matrix_rank(cell) == 3
    assert not _retired_volume_rejects(cell)
    assert not _cell_is_rank_deficient(cell)
    assert _ase(cell) is not None
    assert _ipi(cell) is not None


def test_skewed_cell_below_the_retired_volume_ratio_is_accepted() -> None:
    cell = _skewed_below_retired_ratio()
    assert np.linalg.matrix_rank(cell) == 3
    assert _retired_volume_rejects(cell)
    assert _ase(cell) is not None
    assert _ipi(cell) is not None


@pytest.mark.parametrize(
    "cell",
    [_coplanar(), _zero_vector(), np.zeros((3, 3))],
    ids=["coplanar", "zero-vector", "all-zero"],
)
def test_rank_deficient_cells_include_zero_vectors(cell: np.ndarray) -> None:
    assert np.linalg.matrix_rank(cell) < 3
    assert _cell_is_rank_deficient(cell)
    with pytest.raises(ValueError, match=re.escape(_ASE)):
        _ase(cell)
    if np.any(cell):
        with pytest.raises(ProtocolError, match=re.escape(_IPI)):
            _ipi(cell)


@pytest.mark.parametrize(
    "factory",
    [_near_collinear, _near_coplanar],
    ids=["near-collinear", "near-coplanar"],
)
@pytest.mark.parametrize(
    ("factor", "deficient"),
    [(0.5, True), (2.0, False)],
    ids=["below-boundary", "above-boundary"],
)
def test_near_singular_cells_cross_the_relative_tolerance(
    factory: Callable[[float], np.ndarray], factor: float, deficient: bool
) -> None:
    cell = factory(factor * _RANK_ORDER * _EPS)
    assert _cell_is_rank_deficient(cell) is deficient
    if deficient:
        with pytest.raises(ValueError, match=re.escape(_ASE)):
            _ase(cell)
        with pytest.raises(ProtocolError, match=re.escape(_IPI)):
            _ipi(cell)
    else:
        assert _ase(cell) is not None
        assert _ipi(cell) is not None


def test_singular_value_at_the_default_tolerance_is_deficient() -> None:
    cell = np.diag([1.0, 1.0, _RANK_ORDER * _EPS])
    assert _cell_is_rank_deficient(cell)


def test_default_tolerance_keeps_existing_valid_calls_compatible() -> None:
    cell = _skewed_triclinic()
    assert not _cell_is_rank_deficient(cell)
    assert not _cell_is_rank_deficient(cell, rtol=None)
    assert _cell_is_rank_deficient(cell, rtol=0.5)


@pytest.mark.parametrize("scale", [1.0e-100, 1.0, 1.0e100])
@pytest.mark.parametrize("factor", [0.5, 2.0], ids=["deficient", "full-rank"])
def test_rank_decision_is_invariant_under_scaling(scale: float, factor: float) -> None:
    cell = _near_coplanar(factor * _RANK_ORDER * _EPS)
    expected = factor < 1.0
    assert _cell_is_rank_deficient(cell * scale) is expected


def test_all_zero_wire_cell_stays_a_cluster() -> None:
    assert Driver._cell_rows(SimpleNamespace(cluster=False), [0.0] * 9) is None


@pytest.mark.parametrize("non_finite", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_cells_raise_a_clear_error(non_finite: float) -> None:
    cell = np.diag([2.0, 3.0, 4.0]).astype(float)
    cell[1, 1] = non_finite
    with pytest.raises(ValueError, match="cell must contain only finite values"):
        _cell_is_rank_deficient(cell)
    with pytest.raises(ValueError, match="cell must contain only finite values"):
        _ase(cell)
    with pytest.raises(ProtocolError, match="invalid cell: cell must contain only finite values"):
        _ipi(cell)


@pytest.mark.parametrize("rtol", [-1.0, float("nan"), float("inf")])
def test_invalid_relative_tolerance_is_rejected(rtol: float) -> None:
    with pytest.raises(ValueError, match="rtol must be finite and non-negative"):
        _cell_is_rank_deficient(np.eye(3), rtol=rtol)
