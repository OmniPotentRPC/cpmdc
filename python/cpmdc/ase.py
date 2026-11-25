"""ASE calculator for one cpmdc session.

Importing this module registers :class:`CPMDCCalculator` under the name
``"cpmdc"`` with :func:`ase.calculators.calculator.register_calculator_class`,
so :func:`ase.calculators.calculator.get_calculator_class` returns it for
``"cpmdc"``. The calculator's ``name`` is that same string, so the name that
a trajectory or an ASE database row records is one that ``get_calculator_class``
loads.

Units are ASE units: energy in eV, forces in eV/Angstrom, stress in
eV/Angstrom^3, and dipole in e*Angstrom. Stress is Voigt order
``(xx, yy, zz, yz, xz, xy)`` and is
``(1/V) dE/d(strain)``, the ASE sign of the library tensor.
``free_energy`` equals ``energy``. The returned energy is the potential
energy of the step; the library does not separate an electronic entropy.

A structure with every periodic flag false is a cluster: the cell is not
sent, and the CPMDParams deck sets the boundary conditions. A change of
``atoms.cell`` alone is then no system change and keeps the stored step. Every flag
true sends ``atoms.cell`` as the periodic cell. The cell is rejected
when those three vectors are linearly dependent. Dependence is rank
below 3: a singular value counts only when it is strictly greater than
``sigma_max * 3 * eps``, with ``eps`` the float64 unit roundoff, the
cutoff ``numpy.linalg.matrix_rank`` uses for a 3 by 3 matrix. The
singular values come from the one-sided Jacobi iteration that the C
library also runs, so both make the same decision for the same cell. A mixed periodic flag raises ``ValueError`` and names
the non-periodic axes.

The CPMDParams message holds the method, including the charge and the
multiplicity, and ``initial_charges`` and ``initial_magmoms`` of the
``Atoms`` are read as ASE's other calculators read them. The total charge
is ``Q = sum(atoms.get_initial_charges())`` and the total moment is
``M = sum(atoms.get_initial_magnetic_moments())``, collinear only: a
non-collinear moment raises ``ValueError``. ``Q`` and ``M`` must be
integers. A sum is accepted when it lies within ``n * eps * sum(|x|)`` of
an integer, the bound on the rounding error of summing ``n`` float64
terms ``x``, with ``eps`` the float64 unit roundoff; otherwise it raises
``ValueError``. The multiplicity is ``|M| + 1``. ``cpmdc_params_method``
reports the effective value of the message. A multiplicity of 1 is the
schema default and counts as unset, whether the message left it out or
wrote 1. A charge of 0 counts as unset only when ``hasCharge`` is false.
Top-level ``hasCharge`` true marks that 0 as explicit. The sum of
``initial_charges`` is written only when that bit is false and the
effective charge is still 0. A typed ``system`` section with
``hasCharge`` true wins over the top-level field. The deck renderer
and the Fortran embed record that charge. The embed leaves
``initial_charges`` and the ``ForceInput`` charge fields unread.
The calculator applies this precedence, and reads the edited message
back, before the session opens. For each of the two settings, an
``Atoms`` value of zero changes nothing. A nonzero ``Atoms`` value is
used when the message leaves the setting unset, and passes silently
when it equals the message's value.
When it differs from a set value, ``ValueError`` names both values and
both sources and the set value is kept. A change of
``initial_charges`` or ``initial_magmoms`` evaluates the structure again,
and opens a new session when the effective charge or multiplicity
changes. Each step runs with ``directory`` as the working directory, so CPMD
reads and writes there when the message sets neither ``permanentDir``
nor ``scratchDir``.

One session serves a run of steps with the same atom count and elemental
composition; a reordering stays on it. A step with another composition
closes it and opens a new session.

Importing this module imports ASE. ``import cpmdc`` does not.
"""

from __future__ import annotations

import functools
import base64
import os
import struct
import warnings
from contextlib import AbstractContextManager, nullcontext
from typing import Any, Sequence

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import (
    CalculationFailed as _ASECalculationFailed,
    Calculator,
    CalculatorSetupError,
    Parameters,
    PropertyNotImplementedError,
    ReadError,
    all_changes,
    equal,
    external_calculators,
    register_calculator_class,
)
from ase.stress import full_3x3_to_voigt_6_stress
from ase.utils import workdir

from cpmdc import (
    ConstantCutoff,
    CPMDCError,
    MethodInfo,
    Potential,
    _warn_schema,
    params_method,
)
from cpmdc._rank import _cell_is_rank_deficient

__all__ = ["CPMDCCalculator", "CalculationFailed", "UnsupportedAtoms"]

_CALCULATOR_NAME = "cpmdc"
_NEED_PARAMS = "params is required: a CPMDParams path or bytes"
_METHOD_IN_PARAMS = (
    "the functional, the cutoff, and every other method setting come from "
    "the CPMDParams message in params; the charge and the multiplicity come "
    "from it or from atoms.initial_charges and atoms.initial_magmoms"
)
_RESTART_IGNORED = (
    "restart is ignored and will be removed: a cpmdc step writes no restart "
    "file; a CPMDParams deck with RESTART WAVEFUNCTION reads RESTART.1 on a "
    "cold start"
)
_NO_READ = (
    "cpmdc keeps no calculator file to read: a step writes none and restores "
    "no results; use parameters= for a file written by Parameters.write"
)
_NO_ATOMS = "the calculator has no atoms: set atoms.calc = calc, or pass atoms"
_NO_STRESS_CLUSTER = (
    "stress needs a fully periodic cell: set atoms.pbc = True and a cell "
    "with three independent vectors"
)
_NO_STRESS_DECK = (
    "the library computed no stress for this periodic cell: the CPMDParams "
    "deck describes an isolated system (SYMMETRY 0, POISSON SOLVER, "
    "ISOLATED MOLECULE, or CLUSTER), or CPMDC_STRESS=0 is set"
)
_NO_DIPOLE = (
    "the library computed no dipole for this geometry: the reference evaluator "
    "has no electronic dipole, and OpenCPMD recomputes one only for a periodic "
    "Gamma-point deck with DIPOLE DYNAMICS or PROP DIPOLE MOMENT BERRY"
)
_NO_OVERRIDE: tuple[int | None, int | None] = (None, None)
_FLOAT64_EPS = float(np.finfo(np.float64).eps)
# CPMDParams in schema/Potentials.capnp: ``charge`` and ``multiplicity`` are
# the two Int32 halves of the second data word. A Cap'n Proto Int32 holds its
# value XOR the schema default, so a multiplicity of 1 is stored as 0.
# ``hasCharge`` is the first bit of the third data word. A data section that
# ends before that bit stores the schema default, false.
_PARAMS_CHARGE_BYTE = 8
_PARAMS_MULTIPLICITY_BYTE = 12
_PARAMS_HAS_CHARGE_BYTE = 20
_PARAMS_MULTIPLICITY_DEFAULT = 1
_PARAMS_SCALAR_BYTES = 16
_CAPNP_WORD = 8
_INT32_MIN = -(2**31)
_INT32_MAX = 2**31 - 1
_PARAMS_INLINE_B64_KEY = "cpmdparams_base64"
_PARAMS_INLINE_B64_ERROR = (
    "params mapping must be "
    "{'cpmdparams_base64': '<standard-base64-bytes>'}"
)


class CalculationFailed(CPMDCError, _ASECalculationFailed):
    """A step that the library could not evaluate.

    It is both ASE's ``CalculationFailed`` and :class:`cpmdc.CPMDCError`,
    so either ``except`` clause catches it. The message names the working
    directory of the step and gives the library's reason.
    """


class UnsupportedAtoms(CalculatorSetupError, ValueError):
    """A structure that the calculator cannot evaluate.

    It is ASE's ``CalculatorSetupError``, which ASE raises for an ``Atoms``
    object that a calculator does not support, and a ``ValueError``, which
    this calculator raised before. The message names the ``pbc`` flags or
    the ``cell`` that cpmdc rejects.
    """


def _pbc_flags(pbc: Any) -> list[bool]:
    raw = list(np.atleast_1d(pbc))
    if len(raw) == 1:
        raw = [raw[0], raw[0], raw[0]]
    if len(raw) != 3:
        raise ValueError(f"pbc must be one flag or three, got {len(raw)}")
    return [bool(flag) for flag in raw]


def _voigt(stress: Sequence[float] | np.ndarray) -> np.ndarray:
    """Row-major 3x3 to ASE Voigt (xx, yy, zz, yz, xz, xy), shear parts averaged."""
    voigt = full_3x3_to_voigt_6_stress(np.asarray(stress, dtype=float).reshape(3, 3))
    return np.asarray(voigt, dtype=float)


def _params_bytes(params: Any) -> bytes:
    if isinstance(params, (bytes, bytearray, memoryview)):
        return bytes(params)
    with open(os.fspath(params), "rb") as handle:
        return handle.read()


def _inline_params_blob(value: Any) -> bytes:
    if not isinstance(value, dict) or set(value) != {_PARAMS_INLINE_B64_KEY}:
        raise TypeError(_PARAMS_INLINE_B64_ERROR)
    encoded = value[_PARAMS_INLINE_B64_KEY]
    if not isinstance(encoded, str):
        raise TypeError(_PARAMS_INLINE_B64_ERROR)
    try:
        return base64.b64decode(encoded.encode("ascii"), validate=True)
    except (UnicodeEncodeError, ValueError) as exc:
        raise ValueError(
            f"params[{_PARAMS_INLINE_B64_KEY!r}] is not valid standard base64"
        ) from exc


def _param_value(key: str, value: Any) -> Any:
    if value is None:
        return None
    if key == "constant_cutoff":
        cutoff = ConstantCutoff.coerce(value)
        assert cutoff is not None
        # A bare E0 stays a float, so todict() and Parameters.write keep it short.
        if cutoff.height_ry is None and cutoff.width_ry is None:
            return cutoff.cutoff_ry
        return cutoff.todict()
    if key == "params" and isinstance(value, (bytearray, memoryview)):
        return bytes(value)
    if key == "params" and isinstance(value, dict):
        return _inline_params_blob(value)
    if isinstance(value, bytes):
        return value
    path = os.fspath(value)
    # A relative path names the file under the directory of the process
    # now, so todict() and a copy opened elsewhere must keep it absolute.
    # An empty path is left for find_library or open to refuse by name.
    return os.path.abspath(path) if path else path


def _require_finite(evaluated: Any) -> dict[str, np.ndarray]:
    """The properties of ``evaluated`` as float arrays; refuse one that is NaN or infinite.

    Each property is converted once, so the caller stores the array that was
    checked.
    """
    arrays: dict[str, np.ndarray] = {}
    for name in ("energy", "forces", "stress", "dipole"):
        if name in evaluated:
            arrays[name] = np.asarray(evaluated[name], dtype=float)
            if not np.all(np.isfinite(arrays[name])):
                raise CPMDCError(f"the library returned a non-finite {name}")
    return arrays


def _integral_total(values: np.ndarray, source: str) -> int:
    """Sum ``values`` and return the integer it equals up to float64 rounding.

    Summing ``n`` float64 terms ``x`` is exact to within ``n * eps *
    sum(|x|)``, so a sum farther than that from an integer is not one.
    """
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{source} holds a non-finite value")
    total = float(np.sum(values))
    bound = values.size * _FLOAT64_EPS * float(np.sum(np.abs(values)))
    nearest = round(total)
    if abs(total - nearest) > bound:
        raise ValueError(
            f"{source} sums to {total!r}, which is not an integer within "
            f"{bound!r}; cpmdc needs an integral total charge and moment"
        )
    return int(nearest)


def _initial_charge(atoms: Atoms) -> int:
    charges = np.asarray(atoms.get_initial_charges(), dtype=float)
    return _integral_total(charges, "Atoms.initial_charges")


def _initial_moment(atoms: Atoms) -> int:
    magmoms = np.asarray(atoms.get_initial_magnetic_moments(), dtype=float)
    if magmoms.ndim == 2:
        if np.any(magmoms != 0.0):
            raise ValueError(
                "Atoms.initial_magmoms holds non-collinear moments, shape "
                f"{magmoms.shape}; cpmdc takes one collinear moment per atom"
            )
        return 0
    return _integral_total(magmoms, "Atoms.initial_magmoms")


def _params_root(blob: bytes) -> tuple[bytearray, int, int]:
    """Return ``(buffer, data_start, data_bytes)`` for the root ``CPMDParams`` struct."""
    buffer = bytearray(blob)
    where = "the CPMDParams message cannot take a charge or multiplicity from the Atoms"
    if len(buffer) < _CAPNP_WORD:
        raise CPMDCError(f"{where}: it is shorter than a Cap'n Proto segment table")
    segments = struct.unpack_from("<I", buffer, 0)[0] + 1
    table = 4 + 4 * segments
    table += -table % _CAPNP_WORD
    if len(buffer) < table:
        raise CPMDCError(f"{where}: its segment table is truncated")
    first_words = struct.unpack_from("<I", buffer, 4)[0]
    first = table
    end = first + _CAPNP_WORD * first_words
    if first_words < 1 or end > len(buffer):
        raise CPMDCError(f"{where}: its first segment is empty or truncated")
    pointer = struct.unpack_from("<Q", buffer, first)[0]
    if pointer & 3 != 0:
        raise CPMDCError(f"{where}: its root is not a struct pointer in the first segment")
    offset = (pointer >> 2) & 0x3FFFFFFF
    if offset >= 1 << 29:
        offset -= 1 << 30
    data_bytes = ((pointer >> 32) & 0xFFFF) * _CAPNP_WORD
    start = first + _CAPNP_WORD * (1 + offset)
    if start < first or start + data_bytes > end or data_bytes < _PARAMS_SCALAR_BYTES:
        raise CPMDCError(
            f"{where}: the root struct has {data_bytes} data bytes and "
            f"{_PARAMS_SCALAR_BYTES} are needed; encode it with the current schema"
        )
    return buffer, start, data_bytes


def _top_level_has_charge(blob: bytes) -> bool:
    """True when the top-level ``CPMDParams.hasCharge`` bit is set.

    A data section that ends before the bit stores the schema default, false.
    A typed ``system`` section is a separate word; the calculator reads that
    charge back from ``cpmdc_params_method`` after an edit.
    """
    buffer, start, data_bytes = _params_root(blob)
    if data_bytes <= _PARAMS_HAS_CHARGE_BYTE:
        return False
    return bool(buffer[start + _PARAMS_HAS_CHARGE_BYTE] & 1)


def _params_with_method(blob: bytes, charge: int | None, multiplicity: int | None) -> bytes:
    """``blob`` with the top-level ``charge`` or ``multiplicity`` of ``CPMDParams`` replaced.

    ``None`` leaves a field as it is. The message is edited in place at the
    word that the schema fixes, so every other byte, including a typed
    section and ``hasCharge``, stays as written. The caller checks the result
    with ``cpmdc_params_method``.
    """
    if charge is None and multiplicity is None:
        return blob
    for name, value in (("charge", charge), ("multiplicity", multiplicity)):
        if value is not None and not _INT32_MIN <= value <= _INT32_MAX:
            raise ValueError(f"{name} {value} does not fit the Int32 of CPMDParams")
    buffer, start, _data_bytes = _params_root(blob)
    if charge is not None:
        struct.pack_into("<i", buffer, start + _PARAMS_CHARGE_BYTE, charge)
    if multiplicity is not None:
        stored = multiplicity ^ _PARAMS_MULTIPLICITY_DEFAULT
        struct.pack_into("<i", buffer, start + _PARAMS_MULTIPLICITY_BYTE, stored)
    return bytes(buffer)


def _working_directory(directory: str) -> AbstractContextManager[Any]:
    # workdir() calls os.getcwd(), which raises when the current directory was deleted.
    return nullcontext() if directory == "." else workdir(directory)


class CPMDCCalculator(Calculator):
    """ASE calculator that evaluates one cpmdc session in process.

    Parameters
    ----------
    params:
        Serialized ``CPMDParams`` bytes, or a path to that file. The bytes
        are read once, when the parameter is set. A relative path is
        resolved against the working directory of the process at that
        moment and stored absolute. The message holds the method:
        functional, cutoff, charge, multiplicity, and the boundary
        conditions of a cluster.
    library:
        Path of ``libcpmdc``, stored absolute like ``params``. The
        environment search in :func:`cpmdc.find_library` is used when
        this is omitted.
    constant_cutoff:
        ``E0`` in Rydberg, a mapping of the :class:`cpmdc.ConstantCutoff`
        fields, or a ``ConstantCutoff``. The session then evaluates the
        constant-cutoff kinetic functional, whose stress carries no Pulay
        stress, for a cell filter or a barostat. ``None``, the default,
        evaluates the message as it is.

    Every keyword goes to ASE's ``Calculator.__init__``: ``directory``,
    ``label``, and ``atoms`` as ASE defines them, ``parameters`` for a
    file written by ``Parameters.write``, and ``params``, ``library``, and
    ``constant_cutoff``, which it passes to :meth:`set`. Another keyword
    raises ``TypeError``. ``schema`` and ``restart`` are deprecated and
    ignored.

    Each step runs with ``directory`` as the working directory. The
    calculator writes no files of its own there.

    The charge and the multiplicity also come from ``atoms.initial_charges``
    and ``atoms.initial_magmoms``, as the module documentation describes. The
    sum of ``initial_charges`` is used only when ``hasCharge`` is false. A
    value that differs from one that ``params`` sets, including an explicit
    charge of 0, raises ``ValueError`` and the set value is kept.
    A change of either array evaluates the structure again, and opens a new
    session when the effective charge or multiplicity changes.

    ``parameters`` holds ``params``, ``library``, and ``constant_cutoff``.
    :meth:`set` changes them, opens a new session, and clears the stored results. The
    calculator can be copied, deep-copied, and pickled; a copy opens its
    own session on its first step.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress", "dipole"]
    default_parameters: dict[str, Any] = {
        "params": None,
        "library": None,
        "constant_cutoff": None,
    }
    # Every change of the Atoms is a system change: initial_charges and
    # initial_magmoms can set the charge and the multiplicity of the session.
    ignored_changes: set[str] = set()
    discard_results_on_any_change = True
    _override: tuple[int | None, int | None] = _NO_OVERRIDE
    _method: MethodInfo | None = None

    def __init__(
        self,
        params: Any = None,
        *,
        library: Any = None,
        constant_cutoff: Any = None,
        **kwargs: Any,
    ) -> None:
        if kwargs.pop("schema", None) is not None:
            _warn_schema(stacklevel=3)
        if kwargs.pop("restart", None) is not None:
            warnings.warn(_RESTART_IGNORED, DeprecationWarning, stacklevel=2)
        self._potential: Potential | None = None
        self._blob = b""
        self._composition: bytes | None = None
        self._override = _NO_OVERRIDE
        self._method = None
        self._no_stress = ""
        self._no_dipole = ""
        # A value left as None would override a parameters= file.
        if params is not None:
            kwargs["params"] = params
        if library is not None:
            kwargs["library"] = library
        if constant_cutoff is not None:
            kwargs["constant_cutoff"] = constant_cutoff
        Calculator.__init__(self, **kwargs)
        if self.parameters.get("params") is None:
            raise TypeError(_NEED_PARAMS)

    def _get_name(self) -> str:
        return _CALCULATOR_NAME

    @property
    def potential(self) -> Potential:
        """The rgpot-shaped potential this calculator evaluates.

        It is opened on first use after :meth:`close`, a copy, or an
        unpickle.
        """
        if self._potential is None:
            self._potential = Potential(
                self._method_blob(self._override),
                library=self.parameters.get("library"),
                constant_cutoff=self.parameters.get("constant_cutoff"),
            )
        return self._potential

    def set(self, **kwargs: Any) -> dict[str, Any]:
        """Change ``params``, ``library``, or ``constant_cutoff``; return the changes.

        A change opens a new session and clears the stored results. The
        old session stays in use when the new one cannot be opened.
        ``parameters`` names a file written by ``Parameters.write``.
        ``schema`` is deprecated and ignored. Another key raises
        ``TypeError``, because the CPMDParams message holds the method.
        """
        if "parameters" in kwargs:
            loaded = Parameters.read(kwargs.pop("parameters"))
            loaded.update(kwargs)
            kwargs = dict(loaded)
        if "schema" in kwargs:
            if kwargs.pop("schema") is not None:
                _warn_schema(stacklevel=3)
        unknown = sorted(key for key in kwargs if key not in self.default_parameters)
        if unknown:
            raise TypeError(f"unknown calculator parameters: {unknown}; {_METHOD_IN_PARAMS}")
        kwargs = {key: _param_value(key, value) for key, value in kwargs.items()}
        if "params" in kwargs and kwargs["params"] is None:
            raise TypeError(_NEED_PARAMS)
        reopen = any(
            key in kwargs and not equal(kwargs[key], self.parameters.get(key))
            for key in ("params", "library", "constant_cutoff")
        )
        if reopen:
            params = kwargs.get("params", self.parameters.get("params"))
            if params is None:
                raise TypeError(_NEED_PARAMS)
            library = kwargs.get("library", self.parameters.get("library"))
            cutoff = kwargs.get("constant_cutoff", self.parameters.get("constant_cutoff"))
            blob = _params_bytes(params)
            potential = Potential(blob, library=library, constant_cutoff=cutoff)
        changed: dict[str, Any] = Calculator.set(self, **kwargs)
        if reopen:
            self._release()
            self._potential = potential
            self._blob = blob
            self._override = _NO_OVERRIDE
            self._method = None
        return changed

    def todict(self, skip_default: bool = True) -> dict[str, Any]:
        """Parameters for a trajectory or a database row.

        A ``params`` buffer given as bytes is stored as standard base64
        under ``{"cpmdparams_base64": ...}``, so
        ``CPMDCCalculator(**todict())`` rebuilds the same calculator from
        JSON data. A path is kept absolute, so
        ``CPMDCCalculator(**todict())`` also reads the same files from any
        working directory.
        """
        dct: dict[str, Any] = Calculator.todict(self, skip_default=skip_default)
        if isinstance(dct.get("params"), bytes):
            dct["params"] = {
                _PARAMS_INLINE_B64_KEY: base64.b64encode(dct["params"]).decode("ascii")
            }
        return dct

    def check_state(self, atoms: Atoms, *args: Any, **kwargs: Any) -> list[str]:
        """Return the system changes that make the stored step stale.

        The arguments after ``atoms`` go to ASE's ``Calculator.check_state``,
        so ``tol`` keeps the default of the ASE release in use. The list is
        ASE's, without ``"cell"`` while the stored structure and
        ``atoms`` are both clusters: a cluster step sends no cell, so a
        different ``atoms.cell`` gives the same step. A cell change that
        comes with a ``pbc`` change stays, because ``"pbc"`` is then in
        the list as well.
        """
        changes: list[str] = Calculator.check_state(self, atoms, *args, **kwargs)
        if (
            "cell" in changes
            and self.atoms is not None
            and not np.any(self.atoms.pbc)
            and not np.any(atoms.pbc)
        ):
            changes.remove("cell")
        return changes

    def reset(self) -> None:
        """Clear the stored results; the session stays open."""
        Calculator.reset(self)
        self._no_stress = ""
        self._no_dipole = ""

    def read(self, label: str) -> None:
        """Refuse to read ``label``: the calculator stores no results in files.

        ASE's default only sets the label, which would tell a caller that
        results were restored when none were. ``parameters`` reads a file
        that ``Parameters.write`` made.
        """
        raise ReadError(f"{_NO_READ}: {label!r}")

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] | None = None,
        system_changes: Any = all_changes,
    ) -> None:
        """Fill ``self.results`` with every property of one step.

        Energy, free energy, and forces are always computed; stress is
        computed for a fully periodic cell. The step runs with
        ``directory`` as the working directory. Asking for stress that the
        step cannot give raises ``PropertyNotImplementedError`` with the
        reason, without computing the step again. A step that the library
        rejects raises :class:`CalculationFailed`.
        """
        if not properties:
            properties = ["energy"]
        if "stress" in properties and not system_changes and self._no_stress:
            raise PropertyNotImplementedError(self._no_stress)
        if "dipole" in properties and not system_changes and self._no_dipole:
            raise PropertyNotImplementedError(self._no_dipole)
        Calculator.calculate(self, atoms, properties, system_changes)
        self.results = {}
        if self.atoms is None:
            raise ValueError(_NO_ATOMS)
        self._no_stress = ""
        self._no_dipole = ""
        override = self._method_override(self.atoms)
        if override != self._override:
            self._method_blob(override)
            self._release()
            self._override = override
        cell = self._periodic_cell(self.atoms)
        if "stress" in properties and cell is None:
            self._no_stress = _NO_STRESS_CLUSTER
            raise PropertyNotImplementedError(_NO_STRESS_CLUSTER)
        atomic_numbers = self.atoms.get_atomic_numbers()
        numbers = atomic_numbers.tolist()
        composition = np.sort(atomic_numbers).tobytes()
        if self._composition is not None and composition != self._composition:
            self._release()
        try:
            potential = self.potential
            # C latches topology_fixed before evaluating forces, so record the
            # composition before evaluate() can raise CPMDCError.
            self._composition = composition
            with _working_directory(self.directory):
                evaluated = potential.evaluate(
                    self.atoms.get_positions(),
                    numbers,
                    None if cell is None else cell.tolist(),
                )
            checked = _require_finite(evaluated)
        except CPMDCError as exc:
            where = os.path.abspath(self.directory)
            raise CalculationFailed(
                f"cpmdc could not evaluate the step in {where}: {exc}", exc.code
            ) from exc
        energy = float(checked["energy"])
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = checked["forces"]
        if "stress" in evaluated and cell is not None:
            self.results["stress"] = _voigt(checked["stress"])
        else:
            self._no_stress = _NO_STRESS_CLUSTER if cell is None else _NO_STRESS_DECK
            if "stress" in properties:
                raise PropertyNotImplementedError(self._no_stress)
        if "dipole" in evaluated:
            self.results["dipole"] = checked["dipole"]
        else:
            missing = potential._session._missing_dipole_symbol()
            self._no_dipole = missing if missing is not None else _NO_DIPOLE
            if "dipole" in properties:
                raise PropertyNotImplementedError(self._no_dipole)

    def _params_method(self) -> MethodInfo:
        """The effective method of ``params``, read once per message."""
        if self._method is None:
            self._method = params_method(self._blob, library=self.parameters.get("library"))
        return self._method

    def _method_override(self, atoms: Atoms) -> tuple[int | None, int | None]:
        """The charge and multiplicity the ``Atoms`` add to the message.

        An entry is ``None`` when the message already holds the value, or
        the ``Atoms`` carry none. A charge of 0 is already held when the
        top-level ``hasCharge`` bit is set. ``ValueError`` when a nonzero
        ``Atoms`` value differs from a value that the message sets. The
        set value is kept.
        """
        total_charge = _initial_charge(atoms)
        moment = _initial_moment(atoms)
        if total_charge == 0 and moment == 0:
            return _NO_OVERRIDE
        method = self._params_method()
        charge: int | None = None
        if total_charge != 0 and method.charge != total_charge:
            # charge != 0 is explicit on its own. hasCharge marks an explicit 0.
            # A typed system section that pins another charge is refused when
            # the edited message is read back, before a session opens.
            if method.charge != 0 or _top_level_has_charge(self._blob):
                raise ValueError(
                    f"charge {total_charge} from Atoms.initial_charges disagrees with "
                    f"charge {method.charge} from the CPMDParams message in params; "
                    "set it in one place, or make the two equal"
                )
            charge = total_charge
        multiplicity: int | None = None
        wanted = abs(moment) + 1
        if moment != 0 and method.multiplicity != wanted:
            if method.multiplicity != _PARAMS_MULTIPLICITY_DEFAULT:
                raise ValueError(
                    f"multiplicity {wanted} (moment {moment}) from Atoms.initial_magmoms "
                    f"disagrees with multiplicity {method.multiplicity} from the "
                    "CPMDParams message in params; set it in one place, or make the "
                    "two equal"
                )
            multiplicity = wanted
        return (charge, multiplicity)

    def _method_blob(self, override: tuple[int | None, int | None]) -> bytes:
        """The message a session opens with: ``params`` plus the ``Atoms`` override.

        The edited message is read back with ``cpmdc_params_method``. A
        value that a typed section or a directive of the message pins is not
        the one the top-level field sets, and that raises ``ValueError``.
        """
        charge, multiplicity = override
        blob = _params_with_method(self._blob, charge, multiplicity)
        if blob is self._blob:
            return blob
        applied = params_method(blob, library=self.parameters.get("library"))
        if charge is not None and applied.charge != charge:
            raise ValueError(
                f"charge {charge} from Atoms.initial_charges is not the charge "
                f"{applied.charge} that the CPMDParams message in params gives: a "
                "typed section or a directive of the message fixes it"
            )
        if multiplicity is not None and applied.multiplicity != multiplicity:
            raise ValueError(
                f"multiplicity {multiplicity} from Atoms.initial_magmoms is not the "
                f"multiplicity {applied.multiplicity} that the CPMDParams message in "
                "params gives: a typed section or a directive of the message fixes it"
            )
        return blob

    @staticmethod
    def _periodic_cell(atoms: Atoms) -> np.ndarray | None:
        """Return the periodic cell, or ``None`` when every axis is open.

        A fully periodic cell must have rank 3. :func:`cpmdc._rank._cell_is_rank_deficient`
        applies the ``numpy.linalg.matrix_rank`` cutoff: singular values
        count only above ``sigma_max * 3 * eps``, with ``eps`` the float64
        unit roundoff.
        """
        flags = _pbc_flags(atoms.pbc)
        if not any(flags):
            return None
        if not all(flags):
            names = [axis for axis, flag in zip(("x", "y", "z"), flags) if not flag]
            raise UnsupportedAtoms(
                "cpmdc takes a fully periodic cell; non-periodic axes: "
                + ", ".join(names)
            )
        cell_array: np.ndarray = atoms.cell.array
        if _cell_is_rank_deficient(cell_array):
            raise UnsupportedAtoms(
                "cpmdc takes a fully periodic cell; atoms.pbc is true but "
                "atoms.cell does not have three independent vectors"
            )
        return cell_array

    def _release(self) -> None:
        if self._potential is not None:
            self._potential.close()
        self._potential = None
        self._composition = None

    def close(self) -> None:
        """Close the session. The next step opens a new one."""
        self._release()

    def __enter__(self) -> CPMDCCalculator:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def __getstate__(self) -> dict[str, Any]:
        state = self.__dict__.copy()
        state["_potential"] = None
        state["_composition"] = None
        state["results"] = dict(self.results)
        state["parameters"] = Parameters(self.parameters)
        state.pop("get_spin_polarized", None)
        return state

    def __setstate__(self, state: dict[str, Any]) -> None:
        self.__dict__.update(state)
        self.get_spin_polarized = self._deprecated_get_spin_polarized


def _register_calculator() -> None:
    """Register :class:`CPMDCCalculator` with ASE under ``"cpmdc"`` once."""
    if external_calculators.get(_CALCULATOR_NAME) is CPMDCCalculator:
        return
    if _CALCULATOR_NAME not in external_calculators:
        register_calculator_class(_CALCULATOR_NAME, CPMDCCalculator)
    else:
        external_calculators[_CALCULATOR_NAME] = CPMDCCalculator


@functools.cache
def _make_ase_plugins() -> frozenset[Any]:
    """Return the ``ase.plugins`` set for ASE's plugin discovery.

    ASE releases without ``ase._4.plugins`` get an empty set. The calculator
    is named by dotted path, so building the set imports no cpmdc code.
    """
    try:
        from ase._4.plugins.plugin import CalculatorPlugin
    except ImportError:
        return frozenset()
    return frozenset(
        {
            CalculatorPlugin(
                name="cpmdc",
                citation="https://github.com/OmniPotentRPC/cpmdc",
                implementation="cpmdc.ase.CPMDCCalculator",
            ),
        }
    )


def __getattr__(name: str) -> Any:
    """Build ``__ase_plugins__`` on first read.

    ``ase._4.plugins`` loads every ``ase.plugins`` entry point while it is
    first imported. A module-level import of it here would run that load
    from inside this import, find ``__ase_plugins__`` unset, and leave
    ``ase._4.plugins.plugins`` without the calculator.
    """
    if name == "__ase_plugins__":
        return _make_ase_plugins()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


_register_calculator()
