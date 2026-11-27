"""bench/npt_block_average.py on series whose answers are known, and the explanation table it generated."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("npt_block_average", ROOT / "bench" / "npt_block_average.py")
nba = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(nba)

SUMMARY = ROOT / "bench" / "results" / "gromacs_npt_crescale_64w.org"
EXPLANATION = ROOT / "docs" / "orgmode" / "explanation" / "md-engines.org"


def ar1(n: int, phi: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(n)
    x = np.empty(n)
    x[0] = noise[0] / np.sqrt(1.0 - phi * phi)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + noise[i]
    return x


@pytest.mark.parametrize("phi", [0.0, 0.5, 0.9])
def test_statistical_inefficiency_of_ar1(phi):
    g = nba.statistical_inefficiency(ar1(200000, phi, 1))
    assert g == pytest.approx((1.0 + phi) / (1.0 - phi), rel=0.05)


def test_block_error_stops_when_no_level_has_min_blocks():
    with pytest.raises(SystemExit, match="no blocking level has at least 8 blocks"):
        nba.block_error(np.arange(6.0), min_blocks=8)


@pytest.mark.parametrize("phi", [0.0, 0.9])
def test_block_error_matches_the_ar1_standard_error(phi):
    n = 2**16
    exact = np.sqrt((1.0 / (1.0 - phi * phi)) * (1.0 + phi) / (1.0 - phi) / n)
    error = nba.block_error(ar1(n, phi, 2), min_blocks=64)
    assert error["blocks"] >= 64 and error["block_frames"] * error["blocks"] <= n
    assert error["sem"] == pytest.approx(exact, rel=0.2)
    assert error["sem"] == pytest.approx(error["sem_half"], rel=0.2)


def test_block_error_keeps_growing_for_a_drift():
    error = nba.block_error(np.linspace(0.0, 1.0, 4096), min_blocks=8)
    assert error["sem"] == pytest.approx(np.sqrt(2.0) * error["sem_half"], rel=0.05)


def test_plateau_of_ar1_lies_past_its_correlation_time():
    n = 2**16
    levels = nba.blocking(ar1(n, 0.9, 4))
    k = nba.plateau(levels, min_blocks=16)
    assert k is not None and 8 <= k <= 10
    exact = np.sqrt((1.0 / (1.0 - 0.81)) * 19.0 / n)
    assert levels[k][1] == pytest.approx(exact, rel=0.2)


def test_plateau_of_white_noise_takes_blocks_past_the_cube_root_of_2n():
    assert nba.plateau(nba.blocking(ar1(2**15, 0.0, 5)), min_blocks=16) == 6


def test_plateau_is_none_while_the_blocks_stay_shorter_than_the_correlation():
    assert nba.plateau(nba.blocking(np.linspace(0.0, 1.0, 4096)), min_blocks=8) is None
    assert nba.plateau(nba.blocking(ar1(4096, 0.999, 6)), min_blocks=8) is None


def test_equilibration_start_follows_the_transient():
    x = ar1(20000, 0.9, 3) + 40.0 * np.exp(-np.arange(20000) / 500.0)
    t0, _ = nba.equilibration_start(x, stride=50)
    assert 1500 <= t0 <= 6000


def test_equilibration_start_stays_in_the_first_half():
    t0, _ = nba.equilibration_start(np.linspace(1.0, 0.0, 400) ** 4, stride=1)
    assert t0 <= 200


def test_explanation_carries_the_generated_table():
    rows = [line for line in SUMMARY.read_text().splitlines() if line.startswith("|")]
    assert rows
    page = EXPLANATION.read_text()
    missing = [row for row in rows if row not in page]
    assert not missing, missing
