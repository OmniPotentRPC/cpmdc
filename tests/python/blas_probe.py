"""A check that NumPy's matrix product is usable as a reference.

Several references in this suite (``gromacs_case.ewald_qm_mm``, the ASE
calculator path) form small dense products with ``@``, which NumPy hands to
OpenBLAS. OpenBLAS picks a GEMM kernel from the CPU it detects. An emulator
that advertises a CPU extension it does not execute correctly, or a kernel
with a defect, makes ``@`` return wrong numbers without any error: on
aarch64 under ``qemu-user`` 8.2 with ``-cpu max`` the detected core is
``armv9sme`` and 195 of 385 product shapes (one OpenBLAS thread) are wrong
by O(1). A reference
built on such a product then disagrees with a correct library by O(1), and
the failure shows up in an unrelated assertion.

:func:`matmul_defect` compares ``@`` with ``np.einsum`` over the shapes the
suite uses. ``np.einsum`` without ``optimize`` is a plain C loop and does
not call BLAS.
"""

from __future__ import annotations

import platform
from collections.abc import Callable

import numpy as np

# (rows, inner, columns) of the products the references form, the Ewald
# lattice shifts (343, 3) @ (3, 3) among them, plus wider and taller ones.
SHAPES = ((343, 3, 3), (2, 8, 16), (3, 3, 257), (16, 16, 16), (64, 100, 8))

# Reordering a sum of k products changes it by at most about k eps times
# the sum of the absolute products, eps = 2**-53, so by 1.1e-14 of that at
# k = 100. A wrong kernel is off by about the size of the result itself, a
# relative 1e-1 to 1e1 of that scale. 1e-12 is two orders above the
# rounding bound and eleven below a wrong kernel.
THRESHOLD = 1e-12


def matmul_defect(matmul: Callable[[np.ndarray, np.ndarray], np.ndarray] = np.matmul) -> float:
    """Largest relative disagreement of ``matmul`` with a BLAS-free product.

    The error of each element is taken against the sum of the absolute
    products of its row and column, the scale rounding acts on, and the
    largest ratio over every shape in :data:`SHAPES` is returned. Both
    memory orders of the second operand are tried.
    """
    rng = np.random.default_rng(20240607)
    worst = 0.0
    for m, k, n in SHAPES:
        a = rng.standard_normal((m, k))
        b = rng.standard_normal((k, n))
        scale = np.einsum("ij,jk->ik", np.abs(a), np.abs(b))
        for operand in (b, np.asfortranarray(b)):
            exact = np.einsum("ij,jk->ik", a, operand)
            got = np.asarray(matmul(a, operand))
            worst = max(worst, float(np.max(np.abs(got - exact) / scale)))
    return worst


def environment_problem() -> str | None:
    """A message when ``@`` is wrong on this machine, otherwise ``None``."""
    defect = matmul_defect()
    if defect <= THRESHOLD:
        return None
    return (
        f"NumPy {np.__version__}'s matrix product is wrong on this machine "
        f"({platform.machine()}): it differs from a plain loop by {defect:.3g} "
        f"of the scale of the result, the limit being {THRESHOLD:g}. "
        "Every test that compares against a NumPy reference would be "
        "comparing against garbage, so the session stops here. "
        "OpenBLAS chooses its GEMM kernel from the detected CPU (OPENBLAS_VERBOSE=2 "
        "prints the core); set OPENBLAS_CORETYPE to an older core, for example "
        "ARMV8 on aarch64, "
        "or run an emulator with a CPU that has no SME, for example "
        "QEMU_CPU=max,sme=off for qemu-user."
    )
