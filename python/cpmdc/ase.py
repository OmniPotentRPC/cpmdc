"""ASE calculator for one cpmdc session.

Importing this module registers :class:`CPMDCCalculator` under the name
``"cpmdc"`` with :func:`ase.calculators.calculator.register_calculator_class`,
so :func:`ase.calculators.calculator.get_calculator_class` returns it for
``"cpmdc"``.

Units are ASE units: energy in eV, forces in eV/Angstrom, stress in
eV/Angstrom^3. Stress is Voigt order ``(xx, yy, zz, yz, xz, xy)`` and is
``(1/V) dE/d(strain)``, the ASE sign of the library tensor.
``free_energy`` equals ``energy``. The returned energy is the potential
energy of the step; the library does not separate an electronic entropy.

A structure with every periodic flag false is a cluster: the cell is not
sent, and the CPMDParams deck sets the boundary conditions. Every flag
true sends ``atoms.cell`` as the periodic cell; a cell without three
independent vectors raises ``ValueError``. A mixed periodic flag raises
``ValueError`` and names the non-periodic axes.

The CPMDParams message holds the method, including the charge and the
multiplicity; ``initial_charges`` and ``initial_magmoms`` are not read.
Each step runs with ``directory`` as the working directory, so CPMD
reads and writes there when the message sets neither ``permanentDir``
nor ``scratchDir``.

One session serves a run of steps with the same atom count and elemental
composition; a reordering stays on it. A step with another composition
closes it and opens a new session.

Importing this module imports ASE. ``import cpmdc`` does not.
"""

from __future__ import annotations

import os
import warnings
from contextlib import AbstractContextManager, nullcontext
from typing import Any, Sequence

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import (
    CalculationFailed as _ASECalculationFailed,
    Calculator,
    Parameters,
    PropertyNotImplementedError,
    all_changes,
    equal,
    external_calculators,
    register_calculator_class,
)
from ase.stress import full_3x3_to_voigt_6_stress
from ase.utils import workdir

from cpmdc import ConstantCutoff, CPMDCError, Potential, _warn_schema

__all__ = ["CPMDCCalculator", "CalculationFailed"]

_NEED_PARAMS = "params is required: a CPMDParams path or bytes"
_METHOD_IN_PARAMS = (
    "the charge, the multiplicity, the functional, the cutoff, and every "
    "other method setting come from the CPMDParams message in params"
)
_RESTART_IGNORED = (
    "restart is ignored and will be removed: a cpmdc step writes no restart "
    "file; a CPMDParams deck with RESTART WAVEFUNCTION reads RESTART.1 on a "
    "cold start"
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


class CalculationFailed(CPMDCError, _ASECalculationFailed):
    """A step that the library could not evaluate.

    It is both ASE's ``CalculationFailed`` and :class:`cpmdc.CPMDCError`,
    so either ``except`` clause catches it. The message names the working
    directory of the step and gives the library's reason.
    """


def _pbc_flags(pbc: Any) -> list[bool]:
    raw = list(np.atleast_1d(pbc))
    if len(raw) == 1:
        raw = [raw[0], raw[0], raw[0]]
    if len(raw) != 3:
        raise ValueError(f"pbc must be one flag or three, got {len(raw)}")
    return [bool(flag) for flag in raw]


def _voigt(stress: Sequence[float]) -> np.ndarray:
    """Row-major 3x3 to ASE Voigt (xx, yy, zz, yz, xz, xy), shear parts averaged."""
    voigt = full_3x3_to_voigt_6_stress(np.asarray(stress, dtype=float).reshape(3, 3))
    return np.asarray(voigt, dtype=float)


def _params_bytes(params: Any) -> bytes:
    if isinstance(params, (bytes, bytearray, memoryview)):
        return bytes(params)
    with open(os.fspath(params), "rb") as handle:
        return handle.read()


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
    if isinstance(value, bytes):
        return value
    return os.fspath(value)


def _require_finite(evaluated: Any) -> None:
    """Refuse a step whose energy, forces, or stress is NaN or infinite."""
    for name in ("energy", "forces", "stress"):
        if name in evaluated and not np.all(np.isfinite(np.asarray(evaluated[name], dtype=float))):
            raise CPMDCError(f"the library returned a non-finite {name}")


def _working_directory(directory: str) -> AbstractContextManager[Any]:
    # workdir() calls os.getcwd(), which raises when the current directory was deleted.
    return nullcontext() if directory == "." else workdir(directory)


class CPMDCCalculator(Calculator):
    """ASE calculator that evaluates one cpmdc session in process.

    Parameters
    ----------
    params:
        Serialized ``CPMDParams`` bytes, or a path to that file. The bytes
        are read once, when the parameter is set. The message holds the
        method: functional, cutoff, charge, multiplicity, and the
        boundary conditions of a cluster.
    library:
        Path of ``libcpmdc``. The environment search in
        :func:`cpmdc.find_library` is used when this is omitted.
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

    ``parameters`` holds ``params``, ``library``, and ``constant_cutoff``.
    :meth:`set` changes them, opens a new session, and clears the stored results. The
    calculator can be copied, deep-copied, and pickled; a copy opens its
    own session on its first step.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress"]
    default_parameters: dict[str, Any] = {
        "params": None,
        "library": None,
        "constant_cutoff": None,
    }
    # Charge and multiplicity come from CPMDParams, not from the Atoms.
    ignored_changes = {"initial_charges", "initial_magmoms"}
    discard_results_on_any_change = True

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
        self._composition: tuple[int, ...] | None = None
        self._no_stress = ""
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

    @property
    def potential(self) -> Potential:
        """The rgpot-shaped potential this calculator evaluates.

        It is opened on first use after :meth:`close`, a copy, or an
        unpickle.
        """
        if self._potential is None:
            self._potential = Potential(
                self._blob,
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
        return changed

    def todict(self, skip_default: bool = True) -> dict[str, Any]:
        """Parameters for a trajectory or a database row.

        A ``params`` buffer given as bytes is not JSON data and is left
        out; a path is kept.
        """
        dct: dict[str, Any] = Calculator.todict(self, skip_default=skip_default)
        if isinstance(dct.get("params"), bytes):
            del dct["params"]
        return dct

    def reset(self) -> None:
        """Clear the stored results; the session stays open."""
        Calculator.reset(self)
        self._no_stress = ""

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
        Calculator.calculate(self, atoms, properties, system_changes)
        self.results = {}
        if self.atoms is None:
            raise ValueError(_NO_ATOMS)
        self._no_stress = ""
        cell = self._periodic_cell(self.atoms)
        if "stress" in properties and cell is None:
            self._no_stress = _NO_STRESS_CLUSTER
            raise PropertyNotImplementedError(_NO_STRESS_CLUSTER)
        numbers = [int(number) for number in self.atoms.get_atomic_numbers()]
        composition = tuple(sorted(numbers))
        if self._composition is not None and composition != self._composition:
            self._release()
        try:
            potential = self.potential
            with _working_directory(self.directory):
                evaluated = potential.evaluate(
                    self.atoms.get_positions(),
                    numbers,
                    None if cell is None else cell.tolist(),
                )
            _require_finite(evaluated)
        except CPMDCError as exc:
            where = os.path.abspath(self.directory)
            raise CalculationFailed(f"cpmdc could not evaluate the step in {where}: {exc}") from exc
        self._composition = composition
        energy = float(evaluated["energy"])
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = np.asarray(evaluated["forces"], dtype=float)
        if "stress" in evaluated and cell is not None:
            self.results["stress"] = _voigt(evaluated["stress"])
        else:
            self._no_stress = _NO_STRESS_CLUSTER if cell is None else _NO_STRESS_DECK
            if "stress" in properties:
                raise PropertyNotImplementedError(self._no_stress)

    @staticmethod
    def _periodic_cell(atoms: Atoms) -> np.ndarray | None:
        flags = _pbc_flags(atoms.pbc)
        if not any(flags):
            return None
        if not all(flags):
            names = [axis for axis, flag in zip(("x", "y", "z"), flags) if not flag]
            raise ValueError(
                "cpmdc takes a fully periodic cell; non-periodic axes: "
                + ", ".join(names)
            )
        lengths = atoms.cell.lengths()
        if np.any(lengths <= 0.0) or atoms.cell.volume <= 1.0e-8 * float(np.prod(lengths)):
            raise ValueError(
                "cpmdc takes a fully periodic cell; atoms.pbc is true but "
                "atoms.cell does not have three independent vectors"
            )
        cell_array: np.ndarray = atoms.cell.array
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
    if external_calculators.get("cpmdc") is CPMDCCalculator:
        return
    if "cpmdc" not in external_calculators:
        register_calculator_class("cpmdc", CPMDCCalculator)
    else:
        external_calculators["cpmdc"] = CPMDCCalculator


def _make_ase_plugins() -> set[Any]:
    """Return the guarded ``ase.plugins`` set for ASE v4 plugin discovery."""
    try:
        from ase._4.plugins.plugin import CalculatorPlugin
    except ImportError:
        return set()
    return {
        CalculatorPlugin(
            name="cpmdc",
            citation="https://github.com/OmniPotentRPC/cpmdc",
            implementation="cpmdc.ase.CPMDCCalculator",
        ),
    }


_register_calculator()
__ase_plugins__: set[Any] = _make_ase_plugins()


