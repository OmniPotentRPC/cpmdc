"""Opt-in canonical sampling through velocity rescaling.

The session never calls this module. A host multiplies its own velocities
by the scale :func:`scale` returns. ``temperature_K`` and ``taut`` are the
ASE names for the target temperature and the thermostat time constant.
``taut`` and ``timestep`` share one time unit. ``kB`` converts one kelvin
into the energy unit of the kinetic energy, so an ASE kinetic energy uses
``ase.units.kB`` and a Hartree kinetic energy uses
:data:`HARTREE_PER_KELVIN`.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import math
import os
from numbers import Integral, Real
from pathlib import Path

from cpmdc import CPMDCError

# SI Boltzmann constant divided by the Hartree joule value used for
# HARTREE_TO_KCAL_PER_MOL.
HARTREE_PER_KELVIN = 1.380649e-23 / 4.3597447222071e-18

__all__ = ["HARTREE_PER_KELVIN", "Rng", "scale"]

_PACKAGE_DIR = Path(__file__).resolve().parent
_WHEEL_LIBS = _PACKAGE_DIR.parent / ".cpmdc.mesonpy.libs"
_SHARED_NAMES = ("libcpmdc_csvr.so", "libcpmdc_csvr.dylib")
_LIB: ctypes.CDLL | None = None


class _RngState(ctypes.Structure):
    _fields_ = [
        ("state", ctypes.c_uint64 * 4),
        ("spare_normal", ctypes.c_double),
        ("spare_ready", ctypes.c_int),
    ]


def _candidate_paths() -> list[Path]:
    paths: list[Path] = []
    named = os.environ.get("CPMDC_CSVR_LIBRARY")
    if named:
        paths.append(Path(named))
    for name in _SHARED_NAMES:
        paths.append(_WHEEL_LIBS / name)
        paths.append(_PACKAGE_DIR / ".cpmdc.mesonpy.libs" / name)
        paths.append(Path(name))
        paths.append(Path("build") / name)
        paths.append(_PACKAGE_DIR.parents[1] / "build" / name)
    found = ctypes.util.find_library("cpmdc_csvr")
    if found:
        paths.append(Path(found))
    return paths


def _library() -> ctypes.CDLL:
    global _LIB
    if _LIB is not None:
        return _LIB
    if ctypes.sizeof(_RngState) != 48:
        raise CPMDCError("cpmdc_csvr_rng layout is not 48 bytes")
    chosen: Path | None = None
    for candidate in _candidate_paths():
        if candidate.is_file():
            chosen = candidate
            break
    if chosen is None:
        raise CPMDCError(
            "libcpmdc_csvr was not found. Set CPMDC_CSVR_LIBRARY to the shared library."
        )
    try:
        lib = ctypes.CDLL(str(chosen))
    except OSError as exc:
        raise CPMDCError(f"libcpmdc_csvr failed to load: {chosen}") from exc
    lib.cpmdc_csvr_rng_seed.argtypes = [ctypes.POINTER(_RngState), ctypes.c_uint64]
    lib.cpmdc_csvr_rng_seed.restype = ctypes.c_int
    lib.cpmdc_csvr_scale.argtypes = [
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(_RngState),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.cpmdc_csvr_scale.restype = ctypes.c_int
    _LIB = lib
    return lib


class Rng:
    """Caller-owned generator for one CSVR trajectory.

    Two generators seeded with the same integer produce the same scales.
    Do not share one generator across threads.
    """

    def __init__(self, seed: int) -> None:
        if isinstance(seed, bool) or not isinstance(seed, Integral):
            raise TypeError("seed must be an integer")
        seed_i = int(seed)
        if seed_i < 0 or seed_i > 0xFFFFFFFFFFFFFFFF:
            raise ValueError("seed must fit in an unsigned 64-bit integer")
        seed_c = ctypes.c_uint64(seed_i)
        state = _RngState()
        if _library().cpmdc_csvr_rng_seed(ctypes.byref(state), seed_c) != 0:
            raise CPMDCError("csvr seed failed")
        self._state = state


def scale(
    kinetic: float,
    ndof: int,
    rng: Rng,
    *,
    temperature_K: float,
    taut: float,
    timestep: float,
    kB: float,
) -> tuple[float, float]:
    """Return the positive velocity scale and the thermostat work.

    ``work`` is the energy removed from the ions, in the energy unit of
    ``kinetic``. The ASE accumulator ``transferred_energy`` stores the
    opposite sign.
    """
    if isinstance(kinetic, bool) or not isinstance(kinetic, Real):
        raise TypeError("kinetic must be a real number")
    if isinstance(ndof, bool) or not isinstance(ndof, Integral):
        raise TypeError("ndof must be an integer")
    if not isinstance(rng, Rng):
        raise TypeError("rng must be an Rng")
    if isinstance(temperature_K, bool) or not isinstance(temperature_K, Real):
        raise TypeError("temperature_K must be a real number")
    if isinstance(taut, bool) or not isinstance(taut, Real):
        raise TypeError("taut must be a real number")
    if isinstance(timestep, bool) or not isinstance(timestep, Real):
        raise TypeError("timestep must be a real number")
    if isinstance(kB, bool) or not isinstance(kB, Real):
        raise TypeError("kB must be a real number")
    ndof_i = int(ndof)
    if ndof_i < 1 or ndof_i > 2147483647:
        raise ValueError("ndof must be at least 1 and fit in a C int")
    ndof_c = ctypes.c_int(ndof_i)
    kinetic_f = float(kinetic)
    temperature_f = float(temperature_K)
    taut_f = float(taut)
    timestep_f = float(timestep)
    kB_f = float(kB)
    if not math.isfinite(kinetic_f) or kinetic_f <= 0.0:
        raise ValueError("kinetic must be finite and positive")
    if not math.isfinite(temperature_f) or temperature_f < 0.0:
        raise ValueError("temperature_K must be finite and non-negative")
    if not math.isfinite(timestep_f) or timestep_f < 0.0:
        raise ValueError("timestep must be finite and non-negative")
    if not (math.isfinite(taut_f) or math.isinf(taut_f)) or taut_f <= 0.0:
        raise ValueError("taut must be positive or positive infinity")
    if not math.isfinite(kB_f) or kB_f < 0.0:
        raise ValueError("kB must be finite and non-negative")
    temperature = kB_f * temperature_f
    if not math.isfinite(temperature):
        raise ValueError("kB times temperature_K must be finite")
    out_scale = ctypes.c_double()
    out_work = ctypes.c_double()
    rc = _library().cpmdc_csvr_scale(
        kinetic_f,
        ndof_c,
        temperature,
        timestep_f,
        taut_f,
        ctypes.byref(rng._state),
        ctypes.byref(out_scale),
        ctypes.byref(out_work),
    )
    if rc != 0:
        raise CPMDCError("csvr scale rejected the arguments")
    return float(out_scale.value), float(out_work.value)
