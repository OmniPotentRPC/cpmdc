"""The BLAS preflight of conftest.py: it passes a correct product and stops a wrong one."""

from __future__ import annotations

import numpy as np
import pytest

import blas_probe


def test_numpy_matmul_agrees_with_a_plain_loop() -> None:
    assert blas_probe.matmul_defect() <= blas_probe.THRESHOLD
    assert blas_probe.environment_problem() is None


@pytest.mark.parametrize(
    "damage",
    [
        pytest.param(lambda c: np.zeros_like(c), id="zeros"),
        pytest.param(lambda c: np.roll(c, 1, axis=0), id="shifted-rows"),
        pytest.param(lambda c: c * (1.0 + 1e-9), id="relative-1e-9"),
    ],
)
def test_a_wrong_product_is_reported(damage) -> None:
    def broken(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return damage(a @ b)

    assert blas_probe.matmul_defect(broken) > blas_probe.THRESHOLD


def test_rounding_level_reordering_is_not_a_defect() -> None:
    def reordered(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.einsum("ij,jk->ik", a[:, ::-1], b[::-1, :])

    assert blas_probe.matmul_defect(reordered) <= blas_probe.THRESHOLD


def test_the_message_names_the_remedies(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(blas_probe, "matmul_defect", lambda: 2.5)
    text = blas_probe.environment_problem()
    assert text is not None
    for needle in ("OPENBLAS_CORETYPE", "OPENBLAS_VERBOSE", "QEMU_CPU=max,sme=off", "2.5"):
        assert needle in text
