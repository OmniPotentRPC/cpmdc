"""ASE calculator for one cpmdc session.

Units are ASE units: energy in eV, forces in eV/Angstrom, stress in
eV/Angstrom^3. Stress is Voigt order ``(xx, yy, zz, yz, xz, xy)`` and is
the Cartesian tensor the library reports, converted from Hartree/Bohr^3.
``free_energy`` equals ``energy``. The returned energy is the potential
energy of the step; the library does not separate an electronic entropy.

A structure with every periodic flag false is a cluster: the cell is not
sent. A mixed periodic flag raises ``ValueError`` and names the
non-periodic axes. The same session is reused, so the atom count and the
elemental composition stay fixed after the first successful step.

Importing this module imports ASE. ``import cpmdc`` does not.
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator, all_changes

from cpmdc import Potential

_ASE_KEYS = {
    "restart",
    "ignore_bad_restart_file",
    "label",
    "atoms",
    "directory",
}


def _pbc_flags(pbc: Any) -> list[bool]:
    raw = list(pbc)
    if len(raw) == 1:
        raw = [raw[0], raw[0], raw[0]]
    if len(raw) != 3:
        raise ValueError(f"pbc must be one flag or three, got {len(raw)}")
    return [bool(flag) for flag in raw]


def _voigt(stress: Sequence[float]) -> list[float]:
    """Row-major 3x3 to ASE Voigt (xx, yy, zz, yz, xz, xy)."""
    s = list(stress)
    return [s[0], s[4], s[8], s[5], s[2], s[1]]


class CPMDCCalculator(Calculator):
    """ASE calculator that evaluates one cpmdc session in process.

    Parameters
    ----------
    params:
        Serialized ``CPMDParams`` bytes, or a path to that file.
    library:
        Path of ``libcpmdc``. The environment search in
        :func:`cpmdc.find_library` is used when this is omitted.
    schema:
        Path of ``schema/Potentials.capnp``. Required for a periodic cell.
        ``CPMDC_SCHEMA`` is used when this is omitted.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress"]

    def __init__(
        self,
        params: Any,
        *,
        library: Any = None,
        schema: Any = None,
        **kwargs: Any,
    ) -> None:
        ase_kw = {key: value for key, value in kwargs.items() if key in _ASE_KEYS}
        unknown = [key for key in kwargs if key not in _ASE_KEYS]
        if unknown:
            raise TypeError(f"unknown calculator arguments: {sorted(unknown)}")
        self._potential = Potential(params, library=library, schema=schema)
        Calculator.__init__(self, **ase_kw)

    @property
    def potential(self) -> Potential:
        """The rgpot-shaped potential this calculator evaluates."""
        return self._potential

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy",),
        system_changes: Any = all_changes,
    ) -> None:
        """Fill ``self.results``. ``free_energy`` equals ``energy``."""
        Calculator.calculate(self, atoms, properties, system_changes)
        assert self.atoms is not None
        flags = _pbc_flags(self.atoms.pbc)
        cell: Any = None
        if any(flags):
            if not all(flags):
                names = [axis for axis, flag in zip(("x", "y", "z"), flags) if not flag]
                raise ValueError(
                    "cpmdc takes a fully periodic cell; non-periodic axes: "
                    + ", ".join(names)
                )
            cell = self.atoms.cell.array
        evaluated = self._potential.evaluate(
            self.atoms.get_positions(),
            self.atoms.get_atomic_numbers(),
            cell,
        )
        energy = float(evaluated["energy"])
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = np.asarray(evaluated["forces"], dtype=float)
        if "stress" in evaluated:
            self.results["stress"] = np.asarray(_voigt(evaluated["stress"]), dtype=float)

    def close(self) -> None:
        """Close the session."""
        self._potential.close()
