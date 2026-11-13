#!/usr/bin/env python3
"""Session, rgpot-shaped provider, and ASE calculator."""

from __future__ import annotations

import sys
from pathlib import Path


def _central_forces(session, positions, numbers, delta=1.0e-3):
    from cpmdc import HARTREE_TO_EV

    forces = []
    for atom in range(len(positions)):
        row = []
        for axis in range(3):
            plus = [list(item) for item in positions]
            minus = [list(item) for item in positions]
            plus[atom][axis] += delta
            minus[atom][axis] -= delta
            e_plus = session.energy_forces(plus, numbers).energy_h * HARTREE_TO_EV
            e_minus = session.energy_forces(minus, numbers).energy_h * HARTREE_TO_EV
            row.append(-(e_plus - e_minus) / (2.0 * delta))
        forces.append(row)
    return forces


def _close(left, right, tol):
    for a, b in zip(left, right):
        for x, y in zip(a, b):
            if abs(x - y) > tol:
                raise AssertionError(f"{left} != {right}")


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "usage: test_python_api.py SOURCE_ROOT LIBRARY PARAMS PERIODIC",
            file=sys.stderr,
        )
        return 2
    root, library, params, periodic = sys.argv[1:]
    sys.path.insert(0, str(Path(root) / "python"))
    schema = str(Path(root) / "schema" / "Potentials.capnp")

    from cpmdc import (
        CpmdcError,
        Session,
        forces_to_ev_per_angstrom,
    )

    positions = [
        [0.0, 0.0, 0.1173],
        [0.0, 0.7572, -0.4692],
        [0.0, -0.7572, -0.4692],
    ]
    numbers = [8, 1, 1]
    with Session(params, library=library, schema=schema) as session:
        first = session.energy_forces(positions, numbers)
        again = session.energy_forces(positions, numbers)
        if first.energy_h != again.energy_h:
            raise AssertionError("a repeated geometry changed the energy")
        if not session.version().startswith("cpmdc/"):
            raise AssertionError(session.version())
        converted = forces_to_ev_per_angstrom(first.forces_h_bohr)
        _close(_central_forces(session, positions, numbers), converted, 1.0e-6)
        try:
            session.energy_forces([], numbers)
            raise AssertionError("empty positions were accepted")
        except ValueError:
            pass
        try:
            session.energy_forces(positions, [0, 1, 1])
            raise AssertionError("atomic number 0 was accepted")
        except ValueError:
            pass
        try:
            session.energy_forces([[float("nan"), 0.0, 0.0]], [1])
            raise AssertionError("NaN positions were accepted")
        except ValueError:
            pass

    try:
        Session(b"", library=library)
        raise AssertionError("an empty params buffer was accepted")
    except CpmdcError:
        pass
    try:
        Session(b"not a message", library=library)
        raise AssertionError("a garbage params buffer was accepted")
    except CpmdcError:
        pass
    try:
        Session(params, library=str(Path(root) / "missing-libcpmdc.so"))
        raise AssertionError("a missing library was accepted")
    except CpmdcError:
        pass

    from cpmdc import Potential

    with Potential(periodic, library=library, schema=schema) as pot:
        caps = pot.caps()
        if caps != {"stress": True, "energies": False}:
            raise AssertionError(caps)
        raw = pot.evaluate(
            positions,
            numbers,
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
        )
        if "energy" not in raw or len(raw["forces"]) != 3:
            raise AssertionError(raw)
        if "stress" not in raw or len(raw["stress"]) != 9:
            raise AssertionError("periodic stress missing")
        if raw["energy"] == 0.0:
            raise AssertionError("periodic energy was zero")

    import numpy as np
    from ase import Atoms
    from ase.calculators.calculator import Calculator

    from cpmdc.ase import CPMDCCalculator

    if not issubclass(CPMDCCalculator, Calculator):
        raise AssertionError("CPMDCCalculator is not an ASE Calculator")
    atoms = Atoms("OHH", positions=positions)
    atoms.calc = CPMDCCalculator(params, library=library, schema=schema)
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    if forces.shape != (3, 3):
        raise AssertionError(forces.shape)
    if atoms.get_potential_energy() != energy:
        raise AssertionError("ASE cache missed the repeated energy")
    free = atoms.calc.get_property("free_energy", atoms)
    if free != energy:
        raise AssertionError("free_energy differs from energy")
    atoms.calc.close()

    periodic_atoms = Atoms(
        "OHH",
        positions=positions,
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )
    periodic_atoms.calc = CPMDCCalculator(periodic, library=library, schema=schema)
    stress = periodic_atoms.get_stress()
    if np.asarray(stress).shape != (6,):
        raise AssertionError(stress)
    periodic_atoms.calc.close()

    mixed = Atoms("OHH", positions=positions, cell=[10.0, 10.0, 10.0], pbc=[True, True, False])
    mixed.calc = CPMDCCalculator(params, library=library, schema=schema)
    try:
        mixed.get_potential_energy()
        raise AssertionError("a mixed periodic cell was accepted")
    except ValueError as exc:
        if "z" not in str(exc):
            raise
    mixed.calc.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
