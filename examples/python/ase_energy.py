#!/usr/bin/env python3
"""Evaluate one water geometry with the ASE calculator.

Reads CPMDC_PARAMS, the path of a serialized CPMDParams file. CPMDC_LIBRARY
names libcpmdc when it is set; otherwise the library a wheel installed is
used. Prints the energy in eV, the atom count, and the largest absolute
force component in eV/Angstrom.
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
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    print(f"energy_ev={energy:.12f}")
    print(f"natoms={len(atoms)}")
    print(f"fmax={abs(forces).max():.8f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
