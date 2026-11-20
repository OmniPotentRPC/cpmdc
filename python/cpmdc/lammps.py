"""Couple cpmdc to LAMMPS through :mod:`rgpot.lammps`.

The compatibility :class:`FixExternal` keeps cpmdc's constructor and error
type. Its potential adapter converts an :class:`cpmdc.Session` to the standard
Python rgpot potential protocol; :class:`rgpot.lammps.FixExternal` owns the
LAMMPS callback, units, cell, atom ordering, MPI gather, and halt path.

Units: LAMMPS ``metal`` (Angstrom, eV), ``real`` (Angstrom, kcal/mol) and
``electron`` (Bohr, Hartree) are converted; any other unit style raises
``ValueError``. The cell sent is LAMMPS' restricted triclinic box with one
lattice vector per row, ``a = (lx, 0, 0)``, ``b = (xy, ly, 0)``,
``c = (xz, yz, lz)``. Positions are LAMMPS' Cartesian coordinates unwrapped
with the image flags, so an atom that crosses the box keeps a continuous
path. The box origin is not subtracted. The virial ``-dE/d(strain)`` is
LAMMPS' fix virial ``sum r f`` in energy units, symmetrized and sent as
``xx yy zz xy xz yz``.

Atom types map to elements through ``elements``: one symbol or atomic
number per type, in type order. Atom IDs must be the integers 1..N; the
step sent to cpmdc is in ID order.

Under ``mpirun`` every rank runs the same script, and LAMMPS spreads the
atoms over the ranks. World rank 0 owns the session: it gathers the
positions of every rank by atom ID, evaluates, and sends each rank the
forces of the atoms it owns. Energy and virial go to every rank, because
``fix external`` takes its global energy and virial from every rank.
Rank 0 first adopts ``MPI_COMM_SELF`` as the calculator communicator, so
an OpenCPMD build never waits for ranks that do not call it. More than one
rank needs ``mpi4py``.

A failed evaluation cannot raise through LAMMPS. The callback stores it in
:attr:`FixExternal.error`, sends zero forces, and sets a LAMMPS internal
variable that a ``fix halt`` added beside the fix turns into a LAMMPS error
at the end of that step. :meth:`FixExternal.run` raises the stored error
as :class:`cpmdc.CPMDCError`.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from rgpot.lammps import FixExternal as RgpotFixExternal

from cpmdc import (
    BOHR_TO_ANGSTROM,
    HARTREE_TO_EV,
    HARTREE_TO_KCAL_PER_MOL,
    CPMDCError,
    Session,
    _load,
    find_library,
)
from cpmdc.ipi import species_numbers

__all__ = ["UNIT_STYLES", "FixExternal", "box_rows", "type_numbers", "unit_factors"]

# LAMMPS unit style -> (its length unit per Bohr, its energy unit per Hartree).
UNIT_STYLES: dict[str, tuple[float, float]] = {
    "metal": (BOHR_TO_ANGSTROM, HARTREE_TO_EV),
    "real": (BOHR_TO_ANGSTROM, HARTREE_TO_KCAL_PER_MOL),
    "electron": (1.0, 1.0),
}


def unit_factors(style: str) -> tuple[float, float]:
    """Length per Bohr and energy per Hartree in a LAMMPS unit style."""
    try:
        return UNIT_STYLES[style]
    except KeyError:
        raise ValueError(
            f"cpmdc supports LAMMPS units metal, real, and electron; this input uses units {style!r}"
        ) from None


def box_rows(
    boxlo: Sequence[float], boxhi: Sequence[float], xy: float, yz: float, xz: float
) -> np.ndarray:
    """LAMMPS' restricted triclinic box as a 3x3 cell, one vector per row."""
    lx, ly, lz = (float(boxhi[k]) - float(boxlo[k]) for k in range(3))
    return np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]], dtype=float)


def type_numbers(elements: str | Sequence[str | int], n_types: int) -> np.ndarray:
    """Atomic number of each LAMMPS atom type, from one entry per type.

    ``elements`` is a string of symbols or numbers (``"O H"``, ``"8 1"``) or
    a sequence with one symbol or atomic number per type.
    """
    items = elements.split() if isinstance(elements, str) else list(elements)
    numbers: list[int] = []
    for item in items:
        found = species_numbers(str(item))
        if len(found) != 1:
            raise ValueError(f"{item!r} must name one element for one atom type")
        numbers.extend(found)
    if len(numbers) != n_types:
        raise ValueError(
            f"elements names {len(numbers)} atom types, but the LAMMPS box has {n_types}"
        )
    return np.asarray(numbers, dtype=np.int32)


def _adopt_self(library: str | os.PathLike[str] | None) -> None:
    """Make MPI_COMM_SELF the calculator communicator of this rank."""
    from mpi4py import MPI

    lib = _load(find_library(library))
    handle = MPI._addressof(MPI.COMM_SELF)
    if lib.cpmdc_adopt_calculator_comm(handle, MPI._sizeof(MPI.Comm), 1) != 0:
        # A library without MPI cannot adopt and needs no communicator: only
        # an OpenCPMD build would run collectives on the world.
        if lib.cpmdc_available():
            text = lib.cpmdc_last_error()
            raise CPMDCError(
                (text.decode("utf-8", "replace") if text else "")
                or "cpmdc_adopt_calculator_comm failed"
            )


class _SessionPotential:
    """Expose one cpmdc session through rgpot's Python potential protocol."""

    def __init__(
        self,
        blob: bytes,
        library: str | os.PathLike[str] | None,
        periodic: bool,
    ) -> None:
        self.blob = blob
        self.library = library
        self.periodic = periodic
        self.session: Session | None = None
        self.composition: tuple[int, ...] | None = None
        self.forces = np.zeros(0)
        self.virial = np.zeros(9)
        self._length_per_angstrom = 1.0
        self._energy_per_ev = 1.0
        self._length_per_bohr = BOHR_TO_ANGSTROM
        self._energy_per_hartree = HARTREE_TO_EV

    def open(self, adopt_self: bool = False) -> None:
        if self.session is not None:
            return
        if adopt_self:
            _adopt_self(self.library)
        self.session = Session(self.blob, library=self.library)

    def evaluate(
        self, positions: np.ndarray, numbers: np.ndarray, cell: np.ndarray
    ) -> dict[str, Any]:
        composition = tuple(sorted(int(z) for z in numbers))
        if self.session is None or (
            self.composition is not None and composition != self.composition
        ):
            if self.session is not None:
                self.session.close()
                self.session = None
            self.open()
        assert self.session is not None
        if self.forces.size != positions.size:
            self.forces = np.zeros(positions.size)
        length_scale = self._length_per_angstrom / self._length_per_bohr
        energy_scale = self._energy_per_hartree / self._energy_per_ev
        cell_au = cell * length_scale if self.periodic else None
        energy, _, _ = self.session.calculate_au(
            positions * length_scale,
            numbers,
            cell_au,
            forces=self.forces,
            virial=self.virial,
        )
        for name, value in (
            ("energy", energy),
            ("forces", self.forces),
            ("virial", self.virial),
        ):
            if not np.all(np.isfinite(value)):
                raise CPMDCError(f"the library returned a non-finite {name}")
        self.composition = composition
        forces = self.forces.reshape(positions.shape) * energy_scale * length_scale
        stress = None
        if self.periodic:
            volume = abs(float(np.linalg.det(cell)))
            stress = -self.virial.reshape(3, 3) * energy_scale / volume
        return {
            "energy": float(energy) * energy_scale,
            "forces": forces,
            "stress": stress,
        }

    def close(self) -> None:
        if self.session is not None:
            self.session.close()
            self.session = None


class FixExternal(RgpotFixExternal):
    """cpmdc compatibility wrapper around :class:`rgpot.lammps.FixExternal`."""

    def __init__(
        self,
        lmp: Any,
        params: bytes | str | os.PathLike[str],
        elements: str | Sequence[str | int],
        *,
        fix_id: str = "cpmdc",
        create: bool = True,
        library: str | os.PathLike[str] | None = None,
        comm: Any = None,
    ) -> None:
        blob = params if isinstance(params, bytes) else Path(params).read_bytes()
        periodicity = list(lmp.extract_box()[5])
        unit_factors(str(lmp.extract_global("units")))
        numbers = type_numbers(elements, int(lmp.extract_global("ntypes")))
        adapter = _SessionPotential(blob, library, all(periodicity))
        self._adapter = adapter
        try:
            super().__init__(
                lmp,
                adapter,
                [int(number) for number in numbers],
                fix_id=fix_id,
                create=create,
                comm=comm,
                virial="auto",
            )
        except RuntimeError as exc:
            raise CPMDCError(str(exc)) from exc
        except ValueError as exc:
            text = str(exc)
            marker = "non-periodic axes:"
            if marker in text:
                axes = text.split(marker, 1)[1].strip()
                raise ValueError(
                    "cpmdc takes a fully periodic box or none; non-periodic axes: "
                    + axes
                ) from exc
            raise
        length_per_bohr, energy_per_hartree = unit_factors(self.units)
        adapter._length_per_angstrom = self.length_per_angstrom
        adapter._energy_per_ev = self.energy_per_ev
        adapter._length_per_bohr = length_per_bohr
        adapter._energy_per_hartree = energy_per_hartree

        problem = ""
        if self._root:
            try:
                adapter.open(adopt_self=self._comm is not None)
            except Exception as exc:
                problem = str(exc) or type(exc).__name__
        if self._comm is not None:
            problem = self._comm.bcast(problem, root=0)
        if problem:
            super().close()
            raise CPMDCError(problem)

    @property
    def _session(self) -> Session | None:
        """Compatibility session attribute used by existing callers."""
        return self._adapter.session

    def _take_error(self) -> CPMDCError:
        error = super()._take_error()
        if isinstance(error, CPMDCError):
            return error
        return CPMDCError(f"cpmdc step failed: {error}")

    def close(self) -> None:
        """Close the session and remove the fixes this object added."""
        self._adapter.close()
        super().close()
