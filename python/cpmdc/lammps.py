"""Couple cpmdc to LAMMPS in process through ``fix external``.

:class:`FixExternal` hooks a ``fix external pf/callback`` of a ``lammps``
object from LAMMPS' own Python module (``pip install lammps``, or a LAMMPS
build with ``PYTHON`` support) to one :class:`cpmdc.Session`. Every time
LAMMPS needs forces, the callback sends the whole system to
``cpmdc_session_calculate_au()`` and hands back the energy, the forces, and
the virial, so ``thermo`` reports the cpmdc energy in ``pe`` and the
barostats of ``fix npt`` and ``fix nph`` see its pressure.

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
as :class:`cpmdc.CpmdcError`.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from cpmdc import (
    BOHR_TO_ANGSTROM,
    HARTREE_TO_EV,
    HARTREE_TO_KCAL_PER_MOL,
    CpmdcError,
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


def type_numbers(elements: str | Sequence[str | int], ntypes: int) -> np.ndarray:
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
    if len(numbers) != ntypes:
        raise ValueError(
            f"elements names {len(numbers)} atom types, but the LAMMPS box has {ntypes}"
        )
    return np.asarray(numbers, dtype=np.int32)


def _mpi_world(comm: Any, size: int) -> Any:
    if comm is None:
        try:
            from mpi4py import MPI
        except ImportError as exc:
            raise CpmdcError(
                f"LAMMPS runs on {size} MPI ranks; cpmdc.lammps needs mpi4py to "
                "gather the atoms on rank 0"
            ) from exc
        comm = MPI.COMM_WORLD
    if comm.Get_size() != size:
        raise CpmdcError(
            f"comm has {comm.Get_size()} ranks, but LAMMPS runs on {size}; pass "
            "the communicator the lammps object was created with"
        )
    return comm


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
            raise CpmdcError(
                (text.decode("utf-8", "replace") if text else "")
                or "cpmdc_adopt_calculator_comm failed"
            )


class FixExternal:
    """cpmdc as the force engine of a ``fix external`` in a LAMMPS run.

    Parameters
    ----------
    lmp:
        A ``lammps.lammps`` object with the box, atoms, and masses defined.
    params:
        Serialized ``CPMDParams`` bytes, or a path to that file.
    elements:
        One element symbol or atomic number per atom type, in type order.
    fix_id:
        ID of the ``fix external``. With ``create`` true (the default) this
        object issues ``fix <fix_id> all external pf/callback 1 1`` and
        ``fix_modify <fix_id> energy yes virial yes``; otherwise that fix
        must already exist on group ``all`` with those settings.
    library:
        Path of ``libcpmdc``; :func:`cpmdc.find_library` is used otherwise.
    comm:
        The ``mpi4py`` communicator of the lammps object, for more than one
        rank. ``MPI.COMM_WORLD`` is used when omitted.
    """

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
        self._lmp = lmp
        self.fix_id = fix_id
        self.units = str(lmp.extract_global("units"))
        self.length_per_bohr, self.energy_per_hartree = unit_factors(self.units)
        self.type_numbers = type_numbers(elements, int(lmp.extract_global("ntypes")))
        periodicity = list(lmp.extract_box()[5])
        if any(periodicity) and not all(periodicity):
            names = [axis for axis, flag in zip("xyz", periodicity) if not flag]
            raise ValueError(
                "cpmdc takes a fully periodic box or none; non-periodic axes: "
                + ", ".join(names)
            )
        self.periodic = all(periodicity)
        size = int(lmp.extract_setting("world_size"))
        self._comm = _mpi_world(comm, size) if size > 1 else None
        self._root = self._comm is None or self._comm.Get_rank() == 0
        self._imgmask = int(lmp.extract_setting("IMGMASK"))
        self._imgmax = int(lmp.extract_setting("IMGMAX"))
        self._imgbits = int(lmp.extract_setting("IMGBITS"))
        self._img2bits = int(lmp.extract_setting("IMG2BITS"))
        self._blob = params if isinstance(params, bytes) else Path(params).read_bytes()
        self._library = library
        self._session: Session | None = None
        self._composition: tuple[int, ...] | None = None
        self._forces = np.zeros(0)
        self._virial = np.zeros(9)
        self.error: BaseException | None = None
        self.energy = 0.0
        self.evaluations = 0
        problem = ""
        if self._root:
            try:
                if self._comm is not None:
                    _adopt_self(library)
                self._session = Session(self._blob, library=library)
            except Exception as exc:
                problem = str(exc) or type(exc).__name__
        if self._comm is not None:
            problem = self._comm.bcast(problem, root=0)
        if problem:
            raise CpmdcError(problem)
        self._flag = f"{fix_id}_failed"
        self._created = create
        if create:
            lmp.command(f"fix {fix_id} all external pf/callback 1 1")
            lmp.command(f"fix_modify {fix_id} energy yes virial yes")
        lmp.command(f"variable {self._flag} internal 0")
        lmp.command(
            f"fix {fix_id}_halt all halt 1 v_{self._flag} > 0.5 error hard message yes"
        )
        lmp.set_fix_external_callback(fix_id, self._callback, None)

    # --- one step ---------------------------------------------------------

    def _callback(
        self, _caller: Any, _step: int, nlocal: int, tag: Any, x: Any, f: Any
    ) -> None:
        try:
            self._step(int(nlocal), tag, x, f)
        except Exception as exc:
            if self.error is None:
                self.error = exc
            f[:] = 0.0
            self._lmp.set_internal_variable(self._flag, 1.0)

    def _unwrapped(self, nlocal: int, x: Any, rows: np.ndarray) -> np.ndarray:
        xu = np.array(x[:nlocal], dtype=float).reshape(nlocal, 3)
        if self.periodic and nlocal:
            image = np.asarray(
                self._lmp.numpy.extract_atom("image")[:nlocal], dtype=np.int64
            )
            shifts = np.empty((nlocal, 3), dtype=float)
            shifts[:, 0] = (image & self._imgmask) - self._imgmax
            shifts[:, 1] = ((image >> self._imgbits) & self._imgmask) - self._imgmax
            shifts[:, 2] = (image >> self._img2bits) - self._imgmax
            xu += shifts @ rows
        return xu

    def _evaluate(
        self, tags: np.ndarray, xu: np.ndarray, numbers: np.ndarray, rows: np.ndarray
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """Rank 0: energy, forces in tag order, and the 3x3 virial, LAMMPS units."""
        natoms = tags.size
        if natoms == 0 or tags.min() != 1 or tags.max() != natoms:
            raise CpmdcError(
                f"cpmdc needs atom IDs 1..N; LAMMPS holds {natoms} atoms with IDs "
                f"{tags.min() if natoms else '-'}..{tags.max() if natoms else '-'}"
            )
        index = tags - 1
        positions = np.empty((natoms, 3), dtype=float)
        positions[index] = xu / self.length_per_bohr
        ordered = np.empty(natoms, dtype=np.int32)
        ordered[index] = numbers
        composition = tuple(sorted(ordered.tolist()))
        if self._session is None or (
            self._composition is not None and composition != self._composition
        ):
            if self._session is not None:
                self._session.close()
            self._session = Session(self._blob, library=self._library)
        if self._forces.size != 3 * natoms:
            self._forces = np.zeros(3 * natoms)
        cell = rows / self.length_per_bohr if self.periodic else None
        energy, _, _ = self._session.calculate_au(
            positions, ordered, cell, forces=self._forces, virial=self._virial
        )
        self._composition = composition
        force_unit = self.energy_per_hartree / self.length_per_bohr
        forces = self._forces.reshape(natoms, 3)[index] * force_unit
        virial = self._virial.reshape(3, 3) * self.energy_per_hartree
        return energy * self.energy_per_hartree, forces, virial

    def _step(self, nlocal: int, tag: Any, x: Any, f: Any) -> None:
        lmp = self._lmp
        boxlo, boxhi, xy, yz, xz = lmp.extract_box()[:5]
        rows = box_rows(boxlo, boxhi, xy, yz, xz)
        tags = np.asarray(tag[:nlocal], dtype=np.int64)
        types = np.asarray(lmp.numpy.extract_atom("type")[:nlocal], dtype=np.int64)
        numbers = self.type_numbers[types - 1]
        xu = self._unwrapped(nlocal, x, rows)
        if self._comm is None:
            energy, forces, virial = self._evaluate(tags, xu, numbers, rows)
        else:
            comm = self._comm
            parts = comm.gather((tags, xu, numbers), root=0)
            outcome: tuple[Any, ...]
            shares: list[Any] | None = None
            if self._root:
                try:
                    all_tags = np.concatenate([part[0] for part in parts])
                    energy, ordered, virial = self._evaluate(
                        all_tags,
                        np.concatenate([part[1] for part in parts]),
                        np.concatenate([part[2] for part in parts]),
                        rows,
                    )
                    outcome = ("", energy, virial)
                    shares = []
                    start = 0
                    for part in parts:
                        stop = start + part[0].size
                        shares.append(ordered[start:stop])
                        start = stop
                except Exception as exc:
                    outcome = (str(exc) or type(exc).__name__, 0.0, np.zeros((3, 3)))
                    shares = [np.zeros((part[0].size, 3)) for part in parts]
            problem, energy, virial = comm.bcast(
                outcome if self._root else None, root=0
            )
            forces = comm.scatter(shares, root=0)
            if problem:
                raise CpmdcError(problem)
        f[:nlocal] = forces
        sym = 0.5 * (virial + virial.T)
        lmp.fix_external_set_energy_global(self.fix_id, float(energy))
        lmp.fix_external_set_virial_global(
            self.fix_id,
            [sym[0, 0], sym[1, 1], sym[2, 2], sym[0, 1], sym[0, 2], sym[1, 2]],
        )
        self.energy = float(energy)
        self.evaluations += 1

    # --- control ----------------------------------------------------------

    def run(self, steps: int, *keywords: str) -> None:
        """``run steps keywords...``; raise :class:`CpmdcError` for a failed step."""
        command = " ".join(["run", str(int(steps)), *keywords])
        try:
            self._lmp.command(command)
        except Exception as exc:
            if self.error is not None:
                raise self._take_error() from exc
            raise
        if self.error is not None:
            raise self._take_error()

    def _take_error(self) -> CpmdcError:
        error, self.error = self.error, None
        self._lmp.set_internal_variable(self._flag, 0.0)
        if isinstance(error, CpmdcError):
            return error
        return CpmdcError(f"cpmdc step failed: {error}")

    def close(self) -> None:
        """Close the session and remove the fixes this object added."""
        if self._session is not None:
            self._session.close()
            self._session = None
        lmp, self._lmp = self._lmp, None
        if lmp is None:
            return
        try:
            lmp.command(f"unfix {self.fix_id}_halt")
            if self._created:
                lmp.command(f"unfix {self.fix_id}")
        except Exception:
            return

    def __enter__(self) -> FixExternal:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()
