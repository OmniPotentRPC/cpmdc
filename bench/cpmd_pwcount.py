#!/usr/bin/env python3
"""Plane-wave count and FFT mesh of the committed 840 kg/m^3 water cell.

OpenCPMD stores a gamma-point orbital as the hemisphere of reciprocal
lattice vectors whose kinetic energy is at most the cutoff, with the
cutoff in Rydberg so that energy is |G|^2. numpw then takes the real-space
mesh from alat/pi*sqrt(4*cutoff)+1/2 and the next even length fftchk
accepts: roots 2, 3 and 5 for the shipped FFT, and 7 as well for FFTW.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRESSURE = ROOT / "bench" / "results" / "pbe_water_pressure.org"
# cnst.mod.F90: fbohr = 1/0.529177210859, pi = 3.141592653589793
FBOHR = 1.0 / 0.529177210859
PI = 3.141592653589793
CUTOFF_RY = 100.0
DUAL = 4.0
CELL = re.compile(r"([\d.]+) Angstrom\^3, 840\.0 kg/m\^3")


def nint(value: float) -> int:
    """Fortran NINT for a positive value: half away from zero."""
    return int(math.floor(value + 0.5))


def factors_into(length: int, primes: tuple[int, ...]) -> bool:
    if length <= 0:
        return False
    for prime in primes:
        while length % prime == 0:
            length //= prime
    return length == 1


def next_even_mesh(minimum: int, primes: tuple[int, ...]) -> int:
    """fftchk(minimum, 2): the next even length whose prime factors are in primes."""
    length = minimum
    while not (length % 2 == 0 and factors_into(length, primes)):
        length += 1
    return length


def cell_volume_angstrom3(text: str) -> str:
    found = CELL.search(text)
    if found is None:
        raise SystemExit(f"{PRESSURE} has no 840.0 kg/m^3 cell volume")
    return found.group(1)


def plane_waves(alat_bohr: float, cutoff_ry: float) -> int:
    """Gamma hemisphere of |G|^2 <= cutoff, G in bohr^-1 and the cutoff in Rydberg."""
    scale = (2.0 * PI / alat_bohr) ** 2
    nmax = int(math.sqrt(cutoff_ry / scale)) + 2
    count = 0
    for ix in range(-nmax, nmax + 1):
        for iy in range(-nmax, nmax + 1):
            for iz in range(-nmax, nmax + 1):
                if scale * (ix * ix + iy * iy + iz * iz) > cutoff_ry + 1e-8:
                    continue
                if iz > 0 or (iz == 0 and iy > 0) or (iz == 0 and iy == 0 and ix >= 0):
                    count += 1
    return count


def report(volume_text: str) -> str:
    alat = (float(volume_text) ** (1.0 / 3.0)) * FBOHR
    minimum = nint(alat / PI * math.sqrt(DUAL * CUTOFF_RY) + 0.5)
    lines = [
        f"{volume_text} Angstrom^3, {CUTOFF_RY:.0f} Ry: {plane_waves(alat, CUTOFF_RY)} plane waves",
        f"shipped FFT mesh: {next_even_mesh(minimum, (2, 3, 5))}",
        f"FFTW mesh: {next_even_mesh(minimum, (2, 3, 5, 7))}",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    sys.stdout.write(report(cell_volume_angstrom3(PRESSURE.read_text(encoding="utf-8"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
