#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""One water moved towards a face of the isolated cell of an OpenCPMD deck.

The cluster water deck (tests/water_blyp_params.capnp.txt) has a 10
Angstrom Hockney cell with its origin at a corner. The water starts at the
centre and is moved along AXIS until its outermost atom is GAP Angstrom
from the face at 10 Angstrom; a negative gap puts that atom past the face,
and "origin" is the water at the corner, across three faces. Each row is
one evaluation in one session, so each later row starts from stored
orbitals. Prints an org table of the energy and force changes from the
centred water, the net force, and the trace of sum_i x_i (x) f_i in Bohr
Hartree/Bohr, the CELL NONE virial of libcpmdc_gromacs.

CPMDC_RECENTER and CPMDC_RECENTER_MARGIN apply as set in the environment;
CPMDC_RECENTER=off gives the cell without recentring. CPMDC_PSEUDO_DIR must
name the pseudopotential directory.

usage: cluster_drift_scan.py --library LIB --params PARAMS [--axis 0|1|2]
                             [GAP ...]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from cpmdc import BOHR_TO_ANGSTROM, Session  # noqa: E402

WATER = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])
NUMBERS = [8, 1, 1]
BOX = 10.0


def positions(gap: str, axis: int) -> np.ndarray:
    if gap == "origin":
        return WATER.copy()
    x = WATER + BOX / 2
    if gap != "centre":
        x[:, axis] += (BOX - float(gap)) - x[:, axis].max()
    return x


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", required=True)
    parser.add_argument("--params", required=True)
    parser.add_argument("--axis", type=int, default=1, choices=(0, 1, 2))
    parser.add_argument("gaps", nargs="*", default=["3.0", "2.0", "1.0", "0.5", "-0.5", "origin"])
    args = parser.parse_args()
    print(f"CPMDC_RECENTER={os.environ.get('CPMDC_RECENTER', '')!s} axis={'xyz'[args.axis]}")
    print("| Gap, Angstrom | Energy, Hartree | dE, Hartree | max abs dF, Hartree/Bohr "
          "| abs sum F, Hartree/Bohr | tr sum x (x) f, Hartree |")
    print("|-+-+-+-+-+-|")
    with Session(args.params, library=args.library) as session:
        ref = None
        for gap in ["centre", *args.gaps]:
            x = positions(gap, args.axis)
            result = session.energy_forces(x.tolist(), NUMBERS)
            f = np.asarray(result.forces_h_bohr)
            trace = float(np.einsum("ia,ia->", x / BOHR_TO_ANGSTROM, f))
            if ref is None:
                ref = (result.energy_h, f)
            print(f"| {gap} | {result.energy_h:.10f} | {result.energy_h - ref[0]:+.2e} "
                  f"| {np.abs(f - ref[1]).max():.2e} | {np.linalg.norm(f.sum(axis=0)):.2e} "
                  f"| {trace:+.6f} |", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
