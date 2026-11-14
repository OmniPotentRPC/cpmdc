"""In-process session over libcpmdc.

Positions are Angstrom. :class:`Session` returns Hartree and Hartree/Bohr,
the native units of the C ABI. :class:`Potential` returns eV and eV/Angstrom
so an rgpot calculator can call it. One session is one handle. It is not
shared across threads: the C library keeps the active calculator in
process globals.

``CPMDC_LIBRARY``, then ``RGPOT_CPMDC_ENGINE``, then ``RGPOT_CPMD_ENGINE``,
name the shared library, which is the same order rgpot uses after an
explicit path. A wheel carries its own ``libcpmdc`` in
``.cpmdc.mesonpy.libs`` beside this package, and that copy is used when
none of the variables is set. A cell is encoded as ``ForceInput`` in
process. :func:`find_schema` returns ``schema/Potentials.capnp``, which a
wheel also carries, for encoding ``CPMDParams`` text with ``capnp encode``.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import struct
import warnings
from array import array
from pathlib import Path
from typing import Any, Sequence

# src/cpmdc.c CPMDC_HARTREE_EV and src/cpmdc_params.h CPMDC_BOHR_TO_ANGSTROM.
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903
# CODATA 2018 Hartree energy times the Avogadro constant, in thermochemical
# kilocalories (4184 J): the energy unit of LAMMPS "real" units.
HARTREE_TO_KCAL_PER_MOL = 4.3597447222071e-18 * 6.02214076e23 / 4184.0

__all__ = [
    "BOHR_TO_ANGSTROM",
    "HARTREE_TO_EV",
    "HARTREE_TO_KCAL_PER_MOL",
    "CpmdcError",
    "Evaluation",
    "Potential",
    "Session",
    "find_library",
    "find_schema",
    "forces_to_ev_per_angstrom",
    "stress_to_ev_per_angstrom3",
]


class CpmdcError(RuntimeError):
    """The C library rejected a call, or the host could not load it."""


class Evaluation:
    """One finished step.

    ``energy_h`` is Hartree. ``forces_h_bohr`` is an N by 3 list of
    Hartree/Bohr forces (the negative nuclear gradient). ``stress_h_bohr3``
    is the row-major Cartesian tensor in Hartree/Bohr^3 when the step had a
    periodic cell, and ``None`` when the library did not compute stress. It
    has the OpenCPMD sign, ``-(1/V) dE/d(strain)``, the opposite of ASE's.
    """

    def __init__(
        self,
        energy_h: float,
        forces_h_bohr: list[list[float]],
        stress_h_bohr3: list[float] | None,
    ) -> None:
        self.energy_h = energy_h
        self.forces_h_bohr = forces_h_bohr
        self.stress_h_bohr3 = stress_h_bohr3


class _CResult(ctypes.Structure):
    _fields_ = [
        ("ok", ctypes.c_int),
        ("_pad", ctypes.c_int),
        ("energy_h", ctypes.c_double),
        ("message", ctypes.c_char * 512),
    ]


class _CStress(ctypes.Structure):
    _fields_ = [
        ("valid", ctypes.c_int),
        ("_pad", ctypes.c_int),
        ("values", ctypes.c_double * 9),
    ]


_PACKAGE_DIR = Path(__file__).resolve().parent
# meson-python's directory for libraries a wheel installs from libdir.
_WHEEL_LIBS = _PACKAGE_DIR.parent / ".cpmdc.mesonpy.libs"
_SHARED_NAMES = ("libcpmdc.so", "libcpmdc.dylib")


def _library_names() -> tuple[str, ...]:
    return ("CPMDC_LIBRARY", "RGPOT_CPMDC_ENGINE", "RGPOT_CPMD_ENGINE")


def find_library(path: str | os.PathLike[str] | None = None) -> Path:
    """Return the shared library path.

    An explicit ``path`` wins. Otherwise the environment variables
    ``CPMDC_LIBRARY``, ``RGPOT_CPMDC_ENGINE``, and ``RGPOT_CPMD_ENGINE``
    are tried in that order, then the library a wheel installs beside
    this package, then the dynamic linker, then ``build/libcpmdc.so``
    in the working directory.
    """
    if path is not None:
        candidate = Path(path)
        if not candidate.is_file():
            raise CpmdcError(f"cpmdc library not found: {path}")
        return candidate
    for name in _library_names():
        value = os.environ.get(name)
        if value and Path(value).is_file():
            return Path(value)
    for name in _SHARED_NAMES:
        bundled = _WHEEL_LIBS / name
        if bundled.is_file():
            return bundled
    found = ctypes.util.find_library("cpmdc")
    if found:
        located = Path(found)
        if located.is_file():
            return located
    for relative in (Path("build/libcpmdc.so"), Path("build/libcpmdc.dylib")):
        if relative.is_file():
            return relative
    raise CpmdcError(
        "cpmdc library not found. Set CPMDC_LIBRARY to the shared library."
    )


def find_schema(path: str | os.PathLike[str] | None = None) -> Path:
    """Return ``schema/Potentials.capnp``.

    An explicit ``path`` wins, then ``CPMDC_SCHEMA``, then the copy
    installed inside this package by a wheel, then
    ``schema/Potentials.capnp`` in the working directory.
    """
    if path is not None:
        candidate = Path(path)
        if not candidate.is_file():
            raise CpmdcError(f"CPMD schema not found: {path}")
        return candidate
    env = os.environ.get("CPMDC_SCHEMA")
    if env and Path(env).is_file():
        return Path(env)
    bundled = _PACKAGE_DIR / "schema" / "Potentials.capnp"
    if bundled.is_file():
        return bundled
    relative = Path("schema/Potentials.capnp")
    if relative.is_file():
        return relative
    raise CpmdcError(
        "CPMD schema not found. Set CPMDC_SCHEMA to schema/Potentials.capnp."
    )


def forces_to_ev_per_angstrom(forces_h_bohr: Sequence[Sequence[float]]) -> list[list[float]]:
    """Convert Hartree/Bohr forces to eV/Angstrom."""
    scale = HARTREE_TO_EV / BOHR_TO_ANGSTROM
    return [[component * scale for component in row] for row in forces_h_bohr]


def stress_to_ev_per_angstrom3(stress_h_bohr3: Sequence[float]) -> list[float]:
    """Convert a row-major Hartree/Bohr^3 tensor to eV/Angstrom^3.

    The sign is kept; :meth:`Potential.evaluate` negates the result to
    ASE's ``(1/V) dE/d(strain)``.
    """
    if len(stress_h_bohr3) != 9:
        raise ValueError("stress must have 9 components")
    scale = HARTREE_TO_EV / (BOHR_TO_ANGSTROM ** 3)
    return [component * scale for component in stress_h_bohr3]


def _finite(value: float, what: str) -> float:
    number = float(value)
    if number != number or number in (float("inf"), float("-inf")):
        raise ValueError(f"{what} must be finite")
    return number


def _positions(positions: Sequence[Sequence[float]]) -> tuple[int, array]:
    rows = list(positions)
    if not rows:
        raise ValueError("positions must contain at least one atom")
    flat = array("d")
    for index, row in enumerate(rows):
        coords = list(row)
        if len(coords) != 3:
            raise ValueError(f"position {index} must have three coordinates")
        for coord in coords:
            flat.append(_finite(coord, "positions"))
    return len(rows), flat


def _numbers(numbers: Sequence[int], n_atoms: int) -> array:
    values = [int(number) for number in numbers]
    if len(values) != n_atoms:
        raise ValueError("atomic_numbers must have one entry per atom")
    for number in values:
        if number < 1 or number > 118:
            raise ValueError("atomic numbers must be in 1..118")
    return array("i", values)


def _cell(cell: Sequence[Sequence[float]] | Sequence[float] | None) -> list[float] | None:
    if cell is None:
        return None
    raw = list(cell)
    if len(raw) == 3 and all(not isinstance(item, (int, float)) for item in raw):
        flat: list[float] = []
        for row in raw:
            coords = list(row)  # type: ignore[arg-type]
            if len(coords) != 3:
                raise ValueError("cell must be 3x3 or 9 numbers")
            flat.extend(_finite(coord, "cell") for coord in coords)
        return flat
    if len(raw) != 9:
        raise ValueError("cell must be 3x3 or 9 numbers")
    return [_finite(value, "cell") for value in raw]  # type: ignore[arg-type]


_FORMATS = {ctypes.c_double: "d", ctypes.c_int: "i"}


def _contiguous(data: Any, ctype: Any, count: int) -> memoryview | None:
    try:
        view = memoryview(data)
    except TypeError:
        return None
    if (
        view.format.lstrip("@=") != _FORMATS[ctype]
        or view.itemsize != ctypes.sizeof(ctype)
        or not view.c_contiguous
        or view.nbytes != count * view.itemsize
    ):
        return None
    return view


def _flat(data: Any) -> list[Any]:
    flat: list[Any] = []
    for item in data:
        if hasattr(item, "__len__"):
            flat.extend(item)
        else:
            flat.append(item)
    return flat


def _in_array(data: Any, ctype: Any, count: int, what: str) -> Any:
    """``count`` values of ``ctype``; a matching writable buffer is not copied."""
    view = _contiguous(data, ctype, count)
    if view is not None:
        if view.readonly:
            return (ctype * count).from_buffer_copy(view)
        return (ctype * count).from_buffer(data)
    flat = _flat(data)
    if len(flat) != count:
        raise ValueError(f"{what} must have {count} values, got {len(flat)}")
    return (ctype * count)(*flat)


def _out_array(data: Any, count: int, what: str) -> tuple[Any, Any]:
    if data is None:
        data = array("d", bytes(8 * count))
    view = _contiguous(data, ctypes.c_double, count)
    if view is None or view.readonly:
        raise ValueError(f"{what} must be a writable C-contiguous buffer of {count} doubles")
    return data, (ctypes.c_double * count).from_buffer(data)


# ForceInput in schema/Potentials.capnp: two data words (hasCharge,
# charge, hasMultiplicity, multiplicity; all zero is the schema default)
# and five pointers (pos, atmnrs, box, lengthUnit, energyUnit).
_FORCE_INPUT_DATA_WORDS = 2
_FORCE_INPUT_POINTERS = 5
_ELEMENT_BYTE = 2
_ELEMENT_FOUR_BYTES = 4
_ELEMENT_EIGHT_BYTES = 5


def _force_input_message(flat: array, numbers: array, cell: list[float] | None) -> bytes:
    """Single-segment Cap'n Proto ``ForceInput``, laid out as ``capnp encode`` does."""
    lists: list[tuple[int, int, bytes]] = [
        (_ELEMENT_EIGHT_BYTES, len(flat), struct.pack(f"<{len(flat)}d", *flat)),
        (_ELEMENT_FOUR_BYTES, len(numbers), struct.pack(f"<{len(numbers)}i", *numbers)),
    ]
    if cell is not None:
        lists.append((_ELEMENT_EIGHT_BYTES, 9, struct.pack("<9d", *cell)))
    for text in ("angstrom", "hartree"):
        raw = text.encode("ascii") + b"\0"
        lists.append((_ELEMENT_BYTE, len(raw), raw))
    struct_word = 1
    first_pointer = struct_word + _FORCE_INPUT_DATA_WORDS
    body = bytearray(8 * (_FORCE_INPUT_DATA_WORDS + _FORCE_INPUT_POINTERS))
    pointers: dict[int, int] = {}
    next_word = first_pointer + _FORCE_INPUT_POINTERS
    slots = [0, 1, 2, 3, 4] if cell is not None else [0, 1, 3, 4]
    for slot, (size, count, payload) in zip(slots, lists):
        pointer_word = first_pointer + slot
        offset = next_word - (pointer_word + 1)
        pointers[slot] = ((offset << 2) | 1) | ((size | (count << 3)) << 32)
        padded = payload + b"\0" * (-len(payload) % 8)
        body += padded
        next_word += len(padded) // 8
    for slot, value in pointers.items():
        start = 8 * (_FORCE_INPUT_DATA_WORDS + slot)
        body[start : start + 8] = struct.pack("<Q", value)
    root = struct.pack("<Q", (_FORCE_INPUT_DATA_WORDS | (_FORCE_INPUT_POINTERS << 16)) << 32)
    segment = root + bytes(body)
    return struct.pack("<II", 0, len(segment) // 8) + segment


def _warn_schema(stacklevel: int) -> None:
    warnings.warn(
        "the schema argument is ignored and will be removed: cpmdc encodes "
        "ForceInput in process and needs no schema file",
        DeprecationWarning,
        stacklevel=stacklevel,
    )


def _load(path: Path) -> ctypes.CDLL:
    try:
        lib = ctypes.CDLL(str(path))
    except OSError as exc:
        raise CpmdcError(f"failed to load {path}") from exc
    lib.cpmdc_session_create.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    lib.cpmdc_session_create.restype = ctypes.c_void_p
    lib.cpmdc_session_destroy.argtypes = [ctypes.c_void_p]
    lib.cpmdc_session_destroy.restype = None
    lib.cpmdc_session_energy_forces.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.cpmdc_session_energy_forces.restype = _CResult
    lib.cpmdc_session_calculate_forces.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib.cpmdc_session_calculate_forces.restype = _CResult
    lib.cpmdc_session_calculate_au.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.cpmdc_session_calculate_au.restype = ctypes.c_int
    lib.cpmdc_adopt_calculator_comm.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
    lib.cpmdc_adopt_calculator_comm.restype = ctypes.c_int
    lib.cpmdc_session_last_stress.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(_CStress),
    ]
    lib.cpmdc_session_last_stress.restype = ctypes.c_int
    lib.cpmdc_last_error.argtypes = []
    lib.cpmdc_last_error.restype = ctypes.c_char_p
    lib.cpmdc_version.argtypes = []
    lib.cpmdc_version.restype = ctypes.c_char_p
    lib.cpmdc_abi_version.argtypes = []
    lib.cpmdc_abi_version.restype = ctypes.c_int
    lib.cpmdc_available.argtypes = []
    lib.cpmdc_available.restype = ctypes.c_int
    return lib


def _message(result: _CResult) -> str:
    raw = bytes(result.message).split(b"\0", 1)[0]
    return raw.decode("utf-8", "replace")


def _rows(forces: ctypes.Array, n_atoms: int) -> list[list[float]]:  # type: ignore[type-arg]
    return [
        [float(forces[3 * i]), float(forces[3 * i + 1]), float(forces[3 * i + 2])]
        for i in range(n_atoms)
    ]


class Session:
    """One ``CPMDCSession`` handle.

    ``params`` is a serialized ``CPMDParams`` message, or a path to that
    file. The session copies the bytes. Close it when the geometry series
    is finished. A change of atom count or elemental composition is an
    error from the library; start another session.

    ``schema`` is deprecated and ignored: a cell is encoded as
    ``ForceInput`` in process, so no schema file or ``capnp`` command is
    needed.
    """

    def __init__(
        self,
        params: bytes | str | os.PathLike[str],
        *,
        library: str | os.PathLike[str] | None = None,
        schema: str | os.PathLike[str] | None = None,
    ) -> None:
        if schema is not None:
            _warn_schema(stacklevel=3)
        blob = Path(params).read_bytes() if not isinstance(params, bytes) else params
        if not blob:
            raise CpmdcError("CPMDParams buffer is empty")
        self._lib = _load(find_library(library))
        buffer = ctypes.create_string_buffer(blob)
        handle = self._lib.cpmdc_session_create(buffer, len(blob))
        if not handle:
            raise CpmdcError(self._error() or "cpmdc_session_create failed")
        self._handle = handle

    def _error(self) -> str:
        text = self._lib.cpmdc_last_error()
        return text.decode("utf-8", "replace") if text else ""

    def version(self) -> str:
        """Compiled library version string."""
        text = self._lib.cpmdc_version()
        return text.decode("utf-8", "replace") if text else ""

    def abi_version(self) -> int:
        """Numeric ABI generation of the loaded library."""
        return int(self._lib.cpmdc_abi_version())

    def available(self) -> bool:
        """True when this build is linked to OpenCPMD."""
        return bool(self._lib.cpmdc_available())

    def energy_forces(
        self,
        positions: Sequence[Sequence[float]],
        atomic_numbers: Sequence[int],
        cell: Sequence[Sequence[float]] | Sequence[float] | None = None,
    ) -> Evaluation:
        """Energy, forces, and stress for one geometry.

        ``positions`` is Angstrom, shape (N, 3). ``atomic_numbers`` has one
        entry per atom, in 1..118. ``cell`` is an Angstrom 3x3 matrix, or
        nine numbers in row-major order. Omit it for an isolated system.
        A cell is sent as a ``ForceInput`` message built in process.
        A closed session raises :class:`CpmdcError`.
        """
        n_atoms, flat = _positions(positions)
        numbers = _numbers(atomic_numbers, n_atoms)
        box = _cell(cell)
        force_buf = (ctypes.c_double * (n_atoms * 3))()
        pos_buf = (ctypes.c_double * (n_atoms * 3))(*flat)
        num_buf = (ctypes.c_int * n_atoms)(*numbers)
        if not getattr(self, "_handle", None):
            raise CpmdcError("the session is closed; open a new Session")
        if box is None:
            result = self._lib.cpmdc_session_energy_forces(
                self._handle, n_atoms, pos_buf, num_buf, force_buf
            )
        else:
            encoded = _force_input_message(flat, numbers, box)
            payload = ctypes.create_string_buffer(encoded, len(encoded))
            result = self._lib.cpmdc_session_calculate_forces(
                self._handle, payload, len(encoded), force_buf, n_atoms * 3
            )
        if not result.ok:
            raise CpmdcError(_message(result) or self._error() or "evaluation failed")
        stress = _CStress()
        have = self._lib.cpmdc_session_last_stress(self._handle, ctypes.byref(stress))
        tensor = [float(stress.values[i]) for i in range(9)] if have == 0 and stress.valid else None
        return Evaluation(float(result.energy_h), _rows(force_buf, n_atoms), tensor)

    def calculate_au(
        self,
        positions_bohr: Any,
        atomic_numbers: Any,
        cell_bohr: Any = None,
        *,
        forces: Any = None,
        virial: Any = None,
    ) -> tuple[float, Any, Any]:
        """One step through ``cpmdc_session_calculate_au()``, atomic units.

        ``positions_bohr`` holds ``3 * N`` Bohr coordinates, atom-major, flat
        or as N rows. ``atomic_numbers`` has N entries. ``cell_bohr`` is
        ``None`` for a system without a cell, or nine Bohr values (or three
        rows) with one lattice vector per row. A C-contiguous ``float64`` or
        ``int32`` buffer of the right size, such as a NumPy array, is passed
        without a copy.

        Returns ``(energy_h, forces, virial)``: Hartree, ``3 * N``
        Hartree/Bohr forces ``-dE/dR``, and the nine-component row-major
        virial ``-dE/d(strain)`` in Hartree, zero without a cell or a
        stress tensor. ``forces`` and ``virial`` may name writable
        ``float64`` buffers of ``3 * N`` and nine values to fill in place;
        otherwise new ``array('d')`` objects are returned.
        """
        n_atoms = len(atomic_numbers)
        if n_atoms < 1:
            raise ValueError("atomic_numbers must contain at least one atom")
        numbers = _in_array(atomic_numbers, ctypes.c_int, n_atoms, "atomic_numbers")
        positions = _in_array(positions_bohr, ctypes.c_double, 3 * n_atoms, "positions_bohr")
        cell = None if cell_bohr is None else _in_array(cell_bohr, ctypes.c_double, 9, "cell_bohr")
        forces, force_buf = _out_array(forces, 3 * n_atoms, "forces")
        virial, virial_buf = _out_array(virial, 9, "virial")
        energy = ctypes.c_double()
        if not getattr(self, "_handle", None):
            raise CpmdcError("the session is closed; open a new Session")
        status = self._lib.cpmdc_session_calculate_au(
            self._handle, n_atoms, numbers, positions, cell, ctypes.byref(energy),
            force_buf, virial_buf,
        )
        if status != 0:
            raise CpmdcError(self._error() or "evaluation failed")
        return float(energy.value), forces, virial

    def close(self) -> None:
        """Destroy the session handle. A second close does nothing."""
        handle = getattr(self, "_handle", None)
        if handle:
            self._lib.cpmdc_session_destroy(handle)
            self._handle = None

    def __enter__(self) -> Session:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            return


class Potential:
    """rgpot-shaped provider.

    :meth:`evaluate` takes Angstrom positions, atomic numbers, and an
    optional row-major cell. It returns ``energy`` in eV, ``forces`` in
    eV/Angstrom with shape (N, 3), and ``stress`` as nine eV/Angstrom^3
    components in row-major order when the library computed a tensor.
    The stress sign is ASE's, ``(1/V) dE/d(strain)``.
    :meth:`caps` reports that stress may be present and per-atom energies
    are not. ``schema`` is deprecated and ignored, as for :class:`Session`.
    """

    def __init__(
        self,
        params: bytes | str | os.PathLike[str],
        *,
        library: str | os.PathLike[str] | None = None,
        schema: str | os.PathLike[str] | None = None,
    ) -> None:
        if schema is not None:
            _warn_schema(stacklevel=3)
        self._session = Session(params, library=library)

    def caps(self) -> dict[str, bool]:
        """Capabilities an rgpot calculator consults before a step."""
        return {"stress": True, "energies": False}

    def evaluate(
        self,
        positions: Sequence[Sequence[float]],
        numbers: Sequence[int],
        box: Sequence[Sequence[float]] | Sequence[float] | None = None,
    ) -> dict[str, Any]:
        """One step in ASE units."""
        done = self._session.energy_forces(positions, numbers, box)
        out: dict[str, Any] = {
            "energy": done.energy_h * HARTREE_TO_EV,
            "forces": forces_to_ev_per_angstrom(done.forces_h_bohr),
        }
        if done.stress_h_bohr3 is not None:
            # The C tensor uses the OpenCPMD sign, volume * sigma = -dE/d(eps).
            # rgpot and ASE want sigma = (1/V) dE/d(eps).
            out["stress"] = [
                -component
                for component in stress_to_ev_per_angstrom3(done.stress_h_bohr3)
            ]
        return out

    def close(self) -> None:
        """Close the session this provider owns."""
        self._session.close()

    def __enter__(self) -> Potential:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()
