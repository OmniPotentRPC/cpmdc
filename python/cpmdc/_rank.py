"""Rank of a periodic cell.

Three lattice vectors are usable when the 3 by 3 matrix they form has
rank 3. The cutoff is the one ``numpy.linalg.matrix_rank`` uses for a
float64 matrix when ``tol`` is omitted: a singular value counts only when
it is strictly greater than ``sigma_max * max(m, n) * eps``, with ``eps``
the binary64 unit roundoff and ``max(m, n)`` equal to 3.

The singular values do not come from ``numpy.linalg.svd``. LAPACK's
absolute error, about ``eps * sigma_max``, is a third of the cutoff, so
two SVD routines can rank the same cell differently. They come instead
from a one-sided Jacobi iteration on the cell scaled by a power of two,
which ``src/cpmdc_cell_rank.c`` performs with the same binary64
operations in the same order. Only ``+``, ``-``, ``*``, ``/``, ``sqrt``,
``abs``, ``frexp`` and ``ldexp`` appear, each exact or correctly rounded,
so this module, the C library, and ``cpmdc-ipi-c`` reach the same
decision for the same bits, and scaling a cell by a power of two never
changes it.

A rank-3 cell is still unusable when its volume is smaller than the
least normal float64 value: the i-PI virial volume would otherwise
underflow. The volume is ``|det|`` by LU with partial pivoting on the
scaled cell, and the comparison with ``2**-1022`` is made on exponents,
so it neither overflows nor rounds.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from numpy.typing import ArrayLike

__all__ = []

# Values of the enum in src/cpmdc_cell_rank.h.
_VOLUME_UNDERFLOW = -3
_NONFINITE = -2
_RANK_DEFICIENT = -1
_ZERO = 0
_FULL_RANK = 1

_EPS = sys.float_info.epsilon
_MIN_EXP = sys.float_info.min_exp
_SWEEPS = 64


def _entries(cell: ArrayLike) -> list[float]:
    matrix = np.asarray(cell, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError(f"cell must have shape (3, 3), got {matrix.shape}")
    return [float(value) for value in matrix.ravel(order="C")]


def _scale_exponent(a: list[float]) -> int:
    largest = 0.0
    for value in a:
        if abs(value) > largest:
            largest = abs(value)
    return math.frexp(largest)[1]


def _rotate_columns(a: list[float], p: int, q: int) -> bool:
    alpha = 0.0
    beta = 0.0
    gamma = 0.0
    for row in range(3):
        ap = a[3 * row + p]
        aq = a[3 * row + q]
        alpha = alpha + ap * ap
        beta = beta + aq * aq
        gamma = gamma + ap * aq
    if gamma == 0.0 or abs(gamma) <= _EPS * math.sqrt(alpha) * math.sqrt(beta):
        return False
    zeta = (beta - alpha) / (2.0 * gamma)
    if abs(zeta) > 1.0:
        t = 1.0 / (abs(zeta) * (1.0 + math.sqrt(1.0 + 1.0 / (zeta * zeta))))
    else:
        t = 1.0 / (abs(zeta) + math.sqrt(1.0 + zeta * zeta))
    if zeta < 0.0:
        t = -t
    cosine = 1.0 / math.sqrt(1.0 + t * t)
    sine = cosine * t
    for row in range(3):
        ap = a[3 * row + p]
        aq = a[3 * row + q]
        a[3 * row + p] = cosine * ap - sine * aq
        a[3 * row + q] = sine * ap + cosine * aq
    return True


def _abs_det_lu(scaled: list[float]) -> float:
    m = list(scaled)
    for k in range(3):
        pivot = k
        best = abs(m[3 * k + k])
        for i in range(k + 1, 3):
            if abs(m[3 * i + k]) > best:
                pivot = i
                best = abs(m[3 * i + k])
        if best == 0.0:
            return 0.0
        if pivot != k:
            for j in range(3):
                m[3 * k + j], m[3 * pivot + j] = m[3 * pivot + j], m[3 * k + j]
        for i in range(k + 1, 3):
            factor = m[3 * i + k] / m[3 * k + k]
            for j in range(k + 1, 3):
                m[3 * i + j] = m[3 * i + j] - factor * m[3 * k + j]
    return abs(m[0] * m[4] * m[8])


def _cell_state(cell: ArrayLike) -> int:
    """The ``cpmdc_cell_rank3_state`` code of a row-major 3 by 3 cell."""
    entries = _entries(cell)
    if not all(math.isfinite(value) for value in entries):
        return _NONFINITE
    if not any(entries):
        return _ZERO
    exponent = _scale_exponent(entries)
    scaled = [math.ldexp(value, -exponent) for value in entries]
    a = list(scaled)
    for _ in range(_SWEEPS):
        rotated = _rotate_columns(a, 0, 1)
        rotated = _rotate_columns(a, 0, 2) or rotated
        rotated = _rotate_columns(a, 1, 2) or rotated
        if not rotated:
            break
    singular = []
    for col in range(3):
        total = 0.0
        for row in range(3):
            total = total + a[3 * row + col] * a[3 * row + col]
        singular.append(math.sqrt(total))
    largest = 0.0
    for value in singular:
        if value > largest:
            largest = value
    cutoff = largest * 3.0 * _EPS
    if any(not value > cutoff for value in singular):
        return _RANK_DEFICIENT
    determinant = _abs_det_lu(scaled)
    if determinant == 0.0:
        return _RANK_DEFICIENT
    if math.frexp(determinant)[1] + 3 * exponent < _MIN_EXP:
        return _VOLUME_UNDERFLOW
    return _FULL_RANK


def _cell_is_rank_deficient(cell: ArrayLike) -> bool:
    """True when a finite 3 by 3 cell is unusable as a periodic cell.

    A non-finite entry is not a rank question; the coordinate check
    reports it. A cell is unusable when a Jacobi singular value is not
    strictly greater than ``sigma_max * 3 * eps``, the
    ``numpy.linalg.matrix_rank`` cutoff, or when its volume is smaller than
    the least normal float64 value, before a virial calculation can
    underflow it.
    """
    state = _cell_state(cell)
    return state not in (_NONFINITE, _FULL_RANK)
