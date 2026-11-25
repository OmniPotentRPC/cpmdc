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
_ASE_VOLUME = (
    "cpmdc takes a fully periodic cell; atoms.pbc is true but "
    "atoms.cell volume is too small for binary64"
)
_IPI = (
    "the server sent a cell without three independent vectors; send an "
    "all-zero cell or start cpmdc-ipi with --cluster for an isolated system"
)
_IPI_VOLUME = "the server sent a cell whose volume is too small for binary64"
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


def _nearly_coplanar() -> np.ndarray:
    delta = 0.25 * _RANK_ORDER * _EPS
    return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, delta]])


def _short_edge() -> np.ndarray:
    short = 0.25 * _RANK_ORDER * _EPS
    return np.diag([1.0, 1.0, short])


def _near_collinear(relative: float) -> np.ndarray:
    """Two unit vectors rotated off the axes, the third ``relative`` long."""
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
    """Three vectors close to one plane, rotated off the axes."""
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
    [_coplanar(), _zero_vector(), _nearly_coplanar()],
    ids=["coplanar", "zero-vector", "nearly-coplanar"],
)
def test_rank_deficient_cells_raise_the_existing_message(cell: np.ndarray) -> None:
    assert np.linalg.matrix_rank(cell) < 3
    with pytest.raises(ValueError, match=re.escape(_ASE)):
        _ase(cell)
    with pytest.raises(ProtocolError, match=re.escape(_IPI)):
        _ipi(cell)


def test_short_edge_below_the_rank_cutoff_is_refused() -> None:
    """An orthorhombic cell can have a positive volume ratio and still be rank 2."""
    cell = _short_edge()
    assert np.linalg.matrix_rank(cell) < 3
    assert not _retired_volume_rejects(cell)
    with pytest.raises(ValueError, match=re.escape(_ASE)):
        _ase(cell)
    with pytest.raises(ProtocolError, match=re.escape(_IPI)):
        _ipi(cell)


def test_all_zero_cell_is_refused_as_a_periodic_cell() -> None:
    cell = np.zeros((3, 3))
    assert _cell_is_rank_deficient(cell)
    with pytest.raises(ValueError, match=re.escape(_ASE)):
        _ase(cell)


@pytest.mark.parametrize(
    "factory",
    [_near_collinear, _near_coplanar],
    ids=["near-collinear", "near-coplanar"],
)
@pytest.mark.parametrize(
    ("factor", "deficient"),
    [(0.5, True), (2.0, False)],
    ids=["below-cutoff", "above-cutoff"],
)
def test_rotated_near_singular_cells_on_each_side_of_the_cutoff(
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


def test_singular_value_equal_to_the_cutoff_is_refused() -> None:
    """The cutoff is strict: a value equal to it does not count."""
    assert _cell_is_rank_deficient(np.diag([1.0, 1.0, _RANK_ORDER * _EPS]))


# 1e-90 rather than 1e-100: the full-rank cell scaled by 1e-100 has a volume
# below the normal float64 range, which is refused for another reason.
@pytest.mark.parametrize("scale", [1.0e-90, 1.0, 1.0e100])
@pytest.mark.parametrize("factor", [0.5, 2.0], ids=["deficient", "full-rank"])
def test_rank_decision_does_not_depend_on_the_length_unit(scale: float, factor: float) -> None:
    cell = _near_coplanar(factor * _RANK_ORDER * _EPS)
    assert _cell_is_rank_deficient(cell * scale) is (factor < 1.0)


def test_all_zero_wire_cell_stays_a_cluster() -> None:
    assert Driver._cell_rows(SimpleNamespace(cluster=False), [0.0] * 9) is None


@pytest.mark.parametrize("non_finite", [math.nan, math.inf, -math.inf])
def test_non_finite_cell_is_not_reported_as_rank_deficient(non_finite: float) -> None:
    cell = np.diag([2.0, 3.0, 4.0]).astype(float)
    cell[1, 1] = non_finite
    assert not _cell_is_rank_deficient(cell)
    try:
        _ase(cell)
    except ValueError as exc:
        assert "independent vectors" not in str(exc)
    wire = _wire(np.diag([2.0, 3.0, 4.0]))
    wire[4] = non_finite
    try:
        Driver._cell_rows(SimpleNamespace(cluster=False), wire)
    except ProtocolError as exc:
        assert "independent vectors" not in str(exc)


def test_rank_check_requires_a_three_by_three_cell() -> None:
    assert not _cell_is_rank_deficient(((3, 0, 0), (0, 4, 0), (0, 0, 5)))
    with pytest.raises(ValueError, match=r"shape \(3, 3\).+got \(2, 2\)"):
        _cell_is_rank_deficient(np.eye(2))
    with pytest.raises(ValueError, match=r"shape \(3, 3\).+got \(9,\)"):
        _cell_is_rank_deficient([1.0] * 9)
    with pytest.raises(ValueError, match=r"shape \(3, 3\).+got \(\)"):
        _cell_is_rank_deficient(1.0)
    with pytest.raises(ValueError, match=r"shape \(3, 3\).+got \(\)"):
        _cell_is_rank_deficient(np.array(1.0))
    with pytest.raises(ValueError, match="rectangular"):
        _cell_is_rank_deficient([[1.0, 0.0, 0.0], [0.0, 1.0], [0.0, 0.0, 1.0]])
    with pytest.raises(ValueError, match=r"shape \(3, 3\)"):
        _cell_is_rank_deficient([[[1.0], [0.0], [0.0]]] * 3)
    with pytest.raises(TypeError):
        _cell_is_rank_deficient("123456789")
    with pytest.raises(TypeError):
        _cell_is_rank_deficient([[True, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(TypeError):
        _cell_is_rank_deficient([["1", 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


def test_rank_cutoff_does_not_overflow_for_a_large_cell() -> None:
    side = float(np.finfo(np.float64).max) / 2.0
    assert not _cell_is_rank_deficient(np.diag([side, side, side]))


def test_subnormal_cell_is_refused_before_its_volume_underflows() -> None:
    side = float(np.nextafter(0.0, 1.0))
    cell = np.diag([side, side, side])
    assert np.linalg.matrix_rank(cell) == 3
    assert _cell_is_rank_deficient(cell)
    with pytest.raises(ValueError, match=re.escape(_ASE_VOLUME)):
        _ase(cell)
    with pytest.raises(ProtocolError, match=re.escape(_IPI_VOLUME)):
        _ipi(cell)


def test_ipi_judges_the_volume_on_the_bohr_wire_values() -> None:
    """cpmdc-ipi-c sees Bohr; so does cpmdc-ipi, not the Angstrom it sends on."""
    tiny = float(np.finfo(np.float64).tiny)
    normal = 3.1072325e-103
    assert normal**3 / tiny == pytest.approx(1.3482698433963591, rel=1e-12)
    assert (normal * BOHR_TO_ANGSTROM) ** 3 < tiny
    wire = [normal, 0.0, 0.0, 0.0, normal, 0.0, 0.0, 0.0, normal]
    rows = Driver._cell_rows(SimpleNamespace(cluster=False), wire)
    assert rows is not None
    assert rows[0][0] == normal * BOHR_TO_ANGSTROM
    below = 2.8e-103
    assert below**3 < tiny
    with pytest.raises(ProtocolError, match=re.escape(_IPI_VOLUME)):
        Driver._cell_rows(
            SimpleNamespace(cluster=False), [below, 0.0, 0.0, 0.0, below, 0.0, 0.0, 0.0, below]
        )
