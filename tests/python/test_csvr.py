"""The CSVR binding matches the C library on a fixed seed."""

from __future__ import annotations

import ctypes
import math
from pathlib import Path

import pytest

from cpmdc import CPMDCError
from cpmdc.csvr import HARTREE_PER_KELVIN, Rng, scale
from cpmdc.csvr import _RngState


def test_rng_layout() -> None:
    assert ctypes.sizeof(_RngState) == 48


def test_hartree_per_kelvin_matches_the_joule_ratio() -> None:
    assert HARTREE_PER_KELVIN == 1.380649e-23 / 4.3597447222071e-18


def test_same_seed_matches_to_1e_9_and_another_seed_does_not() -> None:
    def run(seed: int) -> list[tuple[float, float]]:
        rng = Rng(seed)
        kinetic = 1.2
        out: list[tuple[float, float]] = []
        for _ in range(16):
            factor, work = scale(
                kinetic,
                5,
                rng,
                temperature_K=10.0,
                taut=2.0,
                timestep=0.5,
                kB=0.1,
            )
            out.append((factor, work))
            kinetic *= factor * factor
        return out

    first = run(5)
    second = run(5)
    for (left, left_work), (right, right_work) in zip(first, second):
        assert abs(left - right) <= 1e-9
        assert abs(left_work - right_work) <= 1e-9
    other = run(6)
    assert any(abs(left - right) > 1e-9 for (left, _), (right, _) in zip(first, other))


def test_infinite_taut_is_nve_within_1e_12() -> None:
    factor, work = scale(
        1.5,
        4,
        Rng(3),
        temperature_K=300.0,
        taut=math.inf,
        timestep=1.0,
        kB=HARTREE_PER_KELVIN,
    )
    assert abs(factor - 1.0) <= 1e-12
    assert abs(work) <= 1e-12


def test_work_restores_the_kinetic_energy() -> None:
    kinetic = 2.0
    factor, work = scale(
        kinetic,
        4,
        Rng(9),
        temperature_K=50.0,
        taut=3.0,
        timestep=0.2,
        kB=1.0,
    )
    assert abs(kinetic * factor * factor + work - kinetic) <= 1e-9 * (1.0 + kinetic)


def test_zero_temperature_decays_without_noise() -> None:
    factor, work = scale(
        2.0,
        4,
        Rng(1),
        temperature_K=0.0,
        taut=4.0,
        timestep=2.0,
        kB=1.0,
    )
    assert abs(factor - math.exp(-0.25)) <= 1e-12
    assert abs(work - 2.0 * (1.0 - math.exp(-0.5))) <= 1e-12


def test_ideal_gas_mean_stays_inside_the_sampling_bound() -> None:
    ndof = 2
    n = 4000
    kT = 1.0
    expected_mean = 0.5 * ndof * kT
    expected_var = 0.5 * ndof * kT * kT
    rng = Rng(7)
    kinetic = expected_mean
    mean = 0.0
    m2 = 0.0
    for i in range(1, n + 1):
        factor, _work = scale(
            kinetic,
            ndof,
            rng,
            temperature_K=kT,
            taut=1e-8,
            timestep=1.0,
            kB=1.0,
        )
        kinetic *= factor * factor
        delta = kinetic - mean
        mean += delta / i
        m2 += delta * (kinetic - mean)
    var = m2 / (n - 1)
    mean_limit = 8.0 * math.sqrt(expected_var / n)
    print(
        f"csvr_py ndof {ndof} mean {mean:.9f} exp_mean {expected_mean:.9f} "
        f"var {var:.9f} exp_var {expected_var:.9f}"
    )
    assert abs(mean - expected_mean) <= mean_limit
    assert abs(var - expected_var) <= 10.0 * expected_var * math.sqrt((12.0 / ndof + 2.0) / n)


def test_find_library_soname_is_passed_to_cdll(monkeypatch: pytest.MonkeyPatch) -> None:
    import cpmdc.csvr as csvr

    soname = "libcpmdc_csvr_not_a_file.so.0"
    opened: list[str] = []

    class _Symbol:
        def __init__(self) -> None:
            self.argtypes = None
            self.restype = None

    class _Loaded:
        def __init__(self) -> None:
            self.cpmdc_csvr_rng_seed = _Symbol()
            self.cpmdc_csvr_scale = _Symbol()

    def open_name(name: str) -> _Loaded:
        opened.append(name)
        return _Loaded()

    monkeypatch.setattr(csvr, "_LIB", None)
    monkeypatch.delenv("CPMDC_CSVR_LIBRARY", raising=False)
    monkeypatch.setattr(csvr, "_WHEEL_LIBS", Path("/no/such/cpmdc-csvr-wheel"))
    monkeypatch.setattr(csvr.ctypes.util, "find_library", lambda _name: soname)
    monkeypatch.setattr(csvr.ctypes, "CDLL", open_name)
    loaded = csvr._library()
    assert opened == [soname]
    assert loaded is not None
    assert not Path(soname).is_file()


def test_library_missing_a_symbol_raises_cpmdc_error(monkeypatch: pytest.MonkeyPatch) -> None:
    import cpmdc.csvr as csvr

    class _Empty:
        pass

    monkeypatch.setattr(csvr, "_LIB", None)
    monkeypatch.delenv("CPMDC_CSVR_LIBRARY", raising=False)
    monkeypatch.setattr(csvr, "_WHEEL_LIBS", Path("/no/such/cpmdc-csvr-wheel"))
    monkeypatch.setattr(csvr.ctypes.util, "find_library", lambda _name: "libcpmdc_csvr.so.9")
    monkeypatch.setattr(csvr.ctypes, "CDLL", lambda _name: _Empty())
    with pytest.raises(CPMDCError):
        csvr._library()


def test_rejects_bad_python_arguments() -> None:
    rng = Rng(2)
    with pytest.raises(TypeError):
        scale(True, 3, rng, temperature_K=1.0, taut=1.0, timestep=1.0, kB=1.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        scale(1.0, True, rng, temperature_K=1.0, taut=1.0, timestep=1.0, kB=1.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        scale(1.0, 3, "no", temperature_K=1.0, taut=1.0, timestep=1.0, kB=1.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        scale(1.0, 3, rng, 1.0, 1.0, 1.0, 1.0)  # type: ignore[misc]
    with pytest.raises(ValueError):
        scale(-1.0, 3, rng, temperature_K=1.0, taut=1.0, timestep=1.0, kB=1.0)
    with pytest.raises(ValueError):
        scale(1.0, 0, rng, temperature_K=1.0, taut=1.0, timestep=1.0, kB=1.0)
    with pytest.raises(ValueError):
        scale(1.0, 3, rng, temperature_K=1.0, taut=0.0, timestep=1.0, kB=1.0)
    with pytest.raises(ValueError):
        scale(1.0, 3, rng, temperature_K=-1.0, taut=1.0, timestep=1.0, kB=1.0)
    with pytest.raises(ValueError):
        scale(1.0, 3, rng, temperature_K=1.0, taut=1.0, timestep=-1.0, kB=1.0)
    with pytest.raises(ValueError):
        scale(1.0, 3, rng, temperature_K=1.0, taut=1.0, timestep=1.0, kB=-1.0)
    with pytest.raises(TypeError):
        Rng(True)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        Rng(-1)
    with pytest.raises(ValueError):
        Rng(1 << 64)
