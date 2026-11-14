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


def _installed(root: str) -> None:
    """The imported package, library, and schema come from site-packages."""
    import cpmdc

    package = Path(cpmdc.__file__).resolve().parent
    if package.is_relative_to(Path(root).resolve()):
        raise AssertionError(f"cpmdc imported from the source tree: {package}")
    library = cpmdc.find_library().resolve()
    if library.parent != package.parent / ".cpmdc.mesonpy.libs":
        raise AssertionError(f"{library} is not the library the wheel installed")
    schema = cpmdc.find_schema().resolve()
    if schema != package / "schema" / "Potentials.capnp":
        raise AssertionError(f"{schema} is not the schema the wheel installed")


def main() -> int:
    args = sys.argv[1:]
    installed = bool(args) and args[0] == "--installed"
    if installed:
        args = args[1:]
    if len(args) != 4:
        print(
            "usage: test_python_api.py [--installed] SOURCE_ROOT LIBRARY PARAMS PERIODIC\n"
            "--installed tests the importable cpmdc; LIBRARY - uses its own library",
            file=sys.stderr,
        )
        return 2
    root, library_arg, params, periodic = args
    library: str | None = None if library_arg == "-" else library_arg
    if installed:
        _installed(root)
    else:
        sys.path.insert(0, str(Path(root) / "python"))

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
    with Session(params, library=library) as session:
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

    with Potential(periodic, library=library) as pot:
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
    from ase.calculators.fd import (
        calculate_numerical_forces,
        calculate_numerical_stress,
    )

    from cpmdc.ase import CPMDCCalculator

    if not issubclass(CPMDCCalculator, Calculator):
        raise AssertionError("CPMDCCalculator is not an ASE Calculator")
    atoms = Atoms("OHH", positions=positions)
    atoms.calc = CPMDCCalculator(params, library=library)
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    if forces.shape != (3, 3):
        raise AssertionError(forces.shape)
    reference_forces = calculate_numerical_forces(atoms, eps=1.0e-4)
    if not np.allclose(forces, reference_forces, rtol=1.0e-5, atol=1.0e-8):
        raise AssertionError(f"forces {forces} ASE finite difference {reference_forces}")
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
    periodic_atoms.calc = CPMDCCalculator(periodic, library=library)
    stress = periodic_atoms.get_stress()
    if np.asarray(stress).shape != (6,):
        raise AssertionError(stress)
    reference_stress = calculate_numerical_stress(periodic_atoms, eps=1.0e-4)
    if not np.allclose(stress, reference_stress, rtol=1.0e-4, atol=1.0e-6):
        raise AssertionError(f"stress {stress} ASE finite difference {reference_stress}")
    periodic_atoms.calc.close()

    mixed = Atoms("OHH", positions=positions, cell=[10.0, 10.0, 10.0], pbc=[True, True, False])
    mixed.calc = CPMDCCalculator(params, library=library)
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
