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


_MASK = (1 << 64) - 1


def _u64(value: int) -> int:
    return value & _MASK


def _rotl(value: int, shift: int) -> int:
    value = _u64(value)
    return _u64((value << shift) | (value >> (64 - shift)))


def _splitmix64(mixed: list[int]) -> int:
    mixed[0] = _u64(mixed[0] + 0x9E3779B97F4A7C15)
    z = mixed[0]
    z = _u64((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9)
    z = _u64((z ^ (z >> 27)) * 0x94D049BB133111EB)
    return _u64(z ^ (z >> 31))


class _Replay:
    """Same generator the C library uses, so the formula sees its draws."""

    def __init__(self, seed: int) -> None:
        mixed = [seed]
        self.state = [_splitmix64(mixed) for _ in range(4)]
        if self.state[0] | self.state[1] | self.state[2] | self.state[3] == 0:
            self.state[0] = 0x9E3779B97F4A7C15
        self.spare = 0.0
        self.spare_ready = False

    def _xoshiro(self) -> int:
        s = self.state
        result = _u64(_rotl(_u64(s[1] * 5), 7) * 9)
        t = _u64(s[1] << 17)
        s[2] = _u64(s[2] ^ s[0])
        s[3] = _u64(s[3] ^ s[1])
        s[1] = _u64(s[1] ^ s[2])
        s[0] = _u64(s[0] ^ s[3])
        s[2] = _u64(s[2] ^ t)
        s[3] = _rotl(s[3], 45)
        return result

    def _unit(self) -> float:
        for _ in range(64):
            bits = self._xoshiro() >> 11
            unit = math.ldexp(float(bits), -53)
            if unit > 0.0:
                return unit
        raise RuntimeError("unit interval")

    def normal(self) -> float:
        if self.spare_ready:
            self.spare_ready = False
            return self.spare
        u1 = self._unit()
        u2 = self._unit()
        mag = math.sqrt(-2.0 * math.log(u1))
        ang = 2.0 * math.acos(-1.0) * u2
        self.spare = mag * math.sin(ang)
        self.spare_ready = True
        return mag * math.cos(ang)

    def _gamma(self, shape: float) -> float:
        if shape < 1.0:
            return self._gamma(shape + 1.0) * (self._unit() ** (1.0 / shape))
        d = shape - (1.0 / 3.0)
        c = 1.0 / math.sqrt(9.0 * d)
        for _ in range(1000000):
            x = 0.0
            v = 0.0
            accepted = False
            for _inner in range(1000000):
                x = self.normal()
                v = 1.0 + c * x
                if v > 0.0:
                    accepted = True
                    break
            if not accepted:
                raise RuntimeError("gamma")
            v = v * v * v
            u = self._unit()
            x2 = x * x
            if u < 1.0 - 0.0331 * x2 * x2 or math.log(u) < 0.5 * x2 + d * (
                1.0 - v + math.log(v)
            ):
                return d * v
        raise RuntimeError("gamma")

    def sum_of_squares(self, count: int) -> float:
        if count <= 0:
            return 0.0
        return 2.0 * self._gamma(0.5 * count)


def _paper_scale(
    seed: int, kinetic: float, ndof: int, temperature: float, timestep: float, taut: float
) -> tuple[float, float]:
    """Rescaling factor from the same normal draws, in the paper's form."""
    replay = _Replay(seed)
    c = math.exp(-timestep / taut)
    ratio = (0.5 * temperature) / kinetic
    normal = replay.normal()
    rest = replay.sum_of_squares(ndof - 1)
    one_c = 1.0 - c
    cross = 2.0 * math.exp(-0.5 * timestep / taut) * math.sqrt(ratio * one_c) * normal
    alpha2 = c + ratio * one_c * (normal * normal + rest) + cross
    return math.sqrt(alpha2), cross


def test_finite_tau_matches_the_rescaling_factor() -> None:
    dt = 1.0
    tau = 10.0
    kinetic = 1.25
    temperature = 0.8
    assert math.exp(-dt / tau) > 0.5
    for ndof in (1, 4):
        matched = False
        for seed in range(1, 64):
            expected, cross = _paper_scale(seed, kinetic, ndof, temperature, dt, tau)
            if abs(cross) <= 1e-3:
                continue
            factor, _work = scale(
                kinetic,
                ndof,
                Rng(seed),
                temperature_K=temperature,
                taut=tau,
                timestep=dt,
                kB=1.0,
            )
            assert abs(factor - expected) <= 1e-12
            matched = True
            break
        assert matched


def test_copied_generator_replays_the_next_scales() -> None:
    rng = Rng(17)
    kinetic = 1.5
    for _ in range(4):
        factor, _work = scale(
            kinetic, 3, rng, temperature_K=1.0, taut=10.0, timestep=1.0, kB=1.0
        )
        kinetic *= factor * factor
    other = Rng.__new__(Rng)
    other._state = _RngState()
    ctypes.memmove(
        ctypes.byref(other._state), ctypes.byref(rng._state), ctypes.sizeof(_RngState)
    )
    for _ in range(6):
        left, left_work = scale(
            kinetic, 3, rng, temperature_K=1.0, taut=10.0, timestep=1.0, kB=1.0
        )
        right, right_work = scale(
            kinetic, 3, other, temperature_K=1.0, taut=10.0, timestep=1.0, kB=1.0
        )
        assert left == right
        assert left_work == right_work
        kinetic *= left * left


def test_ndof_one_negative_root_stays_positive() -> None:
    dt = 1.0
    tau = 10.0
    kinetic = 0.05
    temperature = 1.0
    c = math.exp(-dt / tau)
    s = (1.0 - c) * (0.5 * temperature) / kinetic
    for seed in range(1, 10000):
        replay = _Replay(seed)
        signed = math.sqrt(c) + replay.normal() * math.sqrt(s)
        if signed >= 0.0:
            continue
        factor, _work = scale(
            kinetic,
            1,
            Rng(seed),
            temperature_K=temperature,
            taut=tau,
            timestep=dt,
            kB=1.0,
        )
        assert factor > 0.0
        assert abs(factor - abs(signed)) <= 1e-12
        assert factor != signed
        return
    raise AssertionError("no negative root")


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
