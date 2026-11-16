#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Stress and relaxed cell volume of one water against the plane-wave cutoff.

Each basis cutoff runs the periodic water deck
(tests/water_blyp_periodic_params.capnp.txt, BLYP, 10 Angstrom cube) twice:
as written, and with the constant-cutoff functional at E0 = FRACTION x the
basis cutoff and the library defaults A = 50 x the basis cutoff, sigma =
E0 / 10.

The stress table compares the analytic mean stress of the unrelaxed water
with dE/dV from two cells strained isotropically by STRAIN, and gives the
mean wall time of those three stress calls, which run in one session and
start from the stored orbitals. The relax table runs BFGS on a
FrechetCellFilter from the water rattled by 0.05 Angstrom (seed 17) to
fmax 0.01 eV/Angstrom or STEPS steps, and stops early when the volume falls
below a quarter of the start. One water has almost no intermolecular
stress, so a physical cell stays near 1000 Angstrom^3.

capnp must be on PATH. CPMDC_PSEUDO_DIR must name the pseudopotential
directory. OMP_NUM_THREADS applies as set.

usage: constant_cutoff_scan.py --library LIB [--cutoffs 25 40 60 80]
                               [--fraction 0.75] [--strain 0.015]
                               [--steps 80]
"""

from __future__ import annotations

import argparse
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from ase import Atoms  # noqa: E402
from ase.filters import FrechetCellFilter  # noqa: E402
from ase.optimize import BFGS  # noqa: E402

from cpmdc import CPMDCError, find_schema  # noqa: E402
from cpmdc.ase import CPMDCCalculator  # noqa: E402

WATER = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])
BOX = 10.0
SOURCE = ROOT / "tests" / "water_blyp_periodic_params.capnp.txt"


def params(cutoff: float) -> bytes:
    text = SOURCE.read_text(encoding="utf-8")
    assert text.count("cutOffRy = 40.0") == 2
    text = text.replace("cutOffRy = 40.0", f"cutOffRy = {cutoff:.1f}")
    done = subprocess.run(
        ["capnp", "encode", str(find_schema()), "CPMDParams"],
        input=text.encode("utf-8"),
        stdout=subprocess.PIPE,
        check=True,
    )
    return done.stdout


def water(rattle: bool = False) -> Atoms:
    atoms = Atoms("OHH", positions=WATER, cell=np.eye(3) * BOX, pbc=True)
    if rattle:
        atoms.rattle(0.05, seed=17)
    return atoms


def stress_row(calc: CPMDCCalculator, strain: float) -> tuple[float, float, float]:
    atoms = water()
    atoms.calc = calc
    times = []
    t0 = time.perf_counter()
    pressure = atoms.get_stress()[:3].mean()
    times.append(time.perf_counter() - t0)
    points = []
    for scale in (1.0 - strain, 1.0 + strain):
        strained = atoms.copy()
        strained.set_cell(atoms.cell.array * scale, scale_atoms=True)
        strained.calc = calc
        t0 = time.perf_counter()
        strained.get_stress()
        times.append(time.perf_counter() - t0)
        points.append((strained.get_volume(), strained.get_potential_energy()))
    (v_low, e_low), (v_high, e_high) = points
    return pressure, (e_high - e_low) / (v_high - v_low), float(np.mean(times))


class Collapsed(Exception):
    """The cell fell below a quarter of its start volume."""


def relax_row(calc: CPMDCCalculator, steps: int) -> dict[str, object]:
    atoms = water(rattle=True)
    atoms.calc = calc
    start_volume = atoms.get_volume()
    energies: list[float] = []
    status = "not converged"
    t0 = time.perf_counter()
    try:
        with BFGS(FrechetCellFilter(atoms), logfile=None) as opt:

            def watch() -> None:
                energies.append(atoms.get_potential_energy())
                if atoms.get_volume() < 0.25 * start_volume:
                    raise Collapsed

            opt.attach(watch)
            if opt.run(fmax=0.01, steps=steps):
                status = "converged"
            nsteps = opt.nsteps
    except Collapsed:
        status = "collapsed"
        nsteps = len(energies) - 1
    except CPMDCError as exc:
        status = f"failed: {str(exc).splitlines()[0][:60]}"
        nsteps = len(energies) - 1
    elapsed = time.perf_counter() - t0
    return {
        "status": status,
        "steps": nsteps,
        "volume": atoms.get_volume(),
        "pressure": atoms.get_stress()[:3].mean() if status in ("converged", "not converged") else float("nan"),
        "e_start": energies[0],
        "e_end": energies[-1],
        "rises": int(np.sum(np.diff(energies) > 0.0)),
        "seconds": elapsed,
    }


def cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--library", required=True)
    parser.add_argument("--cutoffs", type=float, nargs="+", default=[25.0, 40.0, 60.0, 80.0])
    parser.add_argument("--fraction", type=float, default=0.75)
    parser.add_argument("--strain", type=float, default=0.015)
    parser.add_argument("--steps", type=int, default=80)
    args = parser.parse_args()

    print(f"CPU: {cpu_model()}; OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '')}")
    print(f"E0 = {args.fraction:g} x cutoff; strain {args.strain:g}; at most {args.steps} BFGS steps")
    print()
    print("| Cutoff, Ry | E0, Ry | Mean stress, eV/A^3 | dE/dV, eV/A^3 "
          "| Difference, eV/A^3 | Time per stress call, s |")
    print("|-+-+-+-+-+-|")
    relax_rows = []
    for cutoff in args.cutoffs:
        blob = params(cutoff)
        for e0 in (None, args.fraction * cutoff):
            label = "-" if e0 is None else f"{e0:g}"
            with CPMDCCalculator(blob, library=args.library, constant_cutoff=e0) as calc:
                pressure, slope, seconds = stress_row(calc, args.strain)
                print(f"| {cutoff:g} | {label} | {pressure:+.3e} | {slope:+.3e} "
                      f"| {pressure - slope:+.2e} | {seconds:.1f} |", flush=True)
            with CPMDCCalculator(blob, library=args.library, constant_cutoff=e0) as calc:
                relax_rows.append((cutoff, label, relax_row(calc, args.steps)))
    print()
    print("| Cutoff, Ry | E0, Ry | Result | Steps | Final volume, A^3 | Final mean stress, eV/A^3 "
          "| Energy change, eV | Steps that raise the energy | Wall time, s |")
    print("|-+-+-+-+-+-+-+-+-|")
    for cutoff, label, row in relax_rows:
        print(f"| {cutoff:g} | {label} | {row['status']} | {row['steps']} | {row['volume']:.1f} "
              f"| {row['pressure']:+.1e} | {row['e_end'] - row['e_start']:+.3f} "
              f"| {row['rises']} | {row['seconds']:.0f} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
