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


def _cell_is_rank_deficient(cell: ArrayLike) -> bool:
    """True when a finite 3 by 3 cell has rank below 3.

    A non-finite entry is not a rank question; the coordinate check
    reports it. Singular values that are not strictly greater than
    ``sigma_max * max(shape) * numpy.finfo(dtype).eps`` do not count,
    which is the ``numpy.linalg.matrix_rank`` rule.
    """
    matrix = np.asarray(cell, dtype=np.float64)
    if not np.all(np.isfinite(matrix)):
        return False
    singular = np.linalg.svd(matrix, compute_uv=False)
    largest = float(singular[0])
    if largest == 0.0:
        return True
    cutoff = largest * max(matrix.shape) * float(np.finfo(singular.dtype).eps)
    return bool(np.any(singular <= cutoff))
