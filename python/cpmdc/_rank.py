"""Rank of a periodic cell.

Three lattice vectors are usable when the 3 by 3 matrix they form has
rank 3. The cutoff is the one ``numpy.linalg.matrix_rank`` uses for a
float64 matrix when ``tol`` is omitted: the singular values are those of
``numpy.linalg.svd``, and a value counts only when it is strictly
greater than ``sigma_max * max(m, n) * eps``. ``eps`` is
``numpy.finfo`` of that dtype, the binary64 unit roundoff for a cell
stored in float64. ``max(m, n)`` is 3.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

__all__ = []


def _cell_is_rank_deficient(cell: ArrayLike, *, rtol: float | None = None) -> bool:
    """Return whether a finite 3 by 3 cell has rank below 3.

    ``rtol`` is relative to the largest singular value. Its default is
    ``max(shape) * numpy.finfo(dtype).eps``, the relative part of the
    ``numpy.linalg.matrix_rank`` rule. Keeping it optional preserves
    existing one-argument callers and the default cutoff.
    """
    matrix = np.asarray(cell, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError(f"cell must have shape (3, 3), got {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("cell must contain only finite values")
    relative = (
        max(matrix.shape) * float(np.finfo(matrix.dtype).eps)
        if rtol is None
        else float(rtol)
    )
    if not np.isfinite(relative) or relative < 0.0:
        raise ValueError(f"cell rank rtol must be finite and non-negative, got {relative}")
    singular = np.linalg.svd(matrix, compute_uv=False)
    largest = float(singular[0])
    if largest == 0.0:
        return True
    return bool(singular[-1] <= largest * relative)
