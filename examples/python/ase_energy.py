#!/usr/bin/env python3
"""Evaluate one water geometry with the ASE calculator.

Reads CPMDC_PARAMS and, when set, CPMDC_LIBRARY and CPMDC_SCHEMA.
"""

from __future__ import annotations

import os
import sys

from ase import Atoms

from cpmdc.ase import CPMDCCalculator


def main() -> int:
    params = os.environ.get("CPMDC_PARAMS")
    if not params:
        print("CPMDC_PARAMS is required", file=sys.stderr)
        return 2
    atoms = Atoms(
        "OHH",
        positions=[
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ],
    )
    atoms.calc = CPMDCCalculator(
        params,
        library=os.environ.get("CPMDC_LIBRARY"),
        schema=os.environ.get("CPMDC_SCHEMA"),
    )
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    print(f"energy_ev={energy:.12f}")
    print(f"natoms={len(atoms)}")
    print(f"fmax={max(abs(float(v)) for row in forces for v in row):.8f}")
    atoms.calc.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
