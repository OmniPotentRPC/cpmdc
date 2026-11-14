"""ASE calculator for one cpmdc session.

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

One session serves a run of steps with the same atom count and elemental
composition; a reordering stays on it. A step with another composition
closes it and opens a new session.

Importing this module imports ASE. ``import cpmdc`` does not.
"""

from __future__ import annotations

import os
import warnings
from typing import Any, Sequence

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import (
    Calculator,
    Parameters,
    PropertyNotImplementedError,
    all_changes,
    equal,
)

from cpmdc import Potential, _warn_schema

_ASE_KEYS = {
    "restart",
    "ignore_bad_restart_file",
    "label",
    "atoms",
    "directory",
}

_NO_STRESS_CLUSTER = (
    "stress needs a fully periodic cell: set atoms.pbc = True and a cell "
    "with three independent vectors"
)
_NO_STRESS_DECK = (
    "the library computed no stress for this periodic cell: the CPMDParams "
    "deck describes an isolated system (SYMMETRY 0, POISSON SOLVER, "
    "ISOLATED MOLECULE, or CLUSTER), or CPMDC_STRESS=0 is set"
)


def _pbc_flags(pbc: Any) -> list[bool]:
    raw = list(np.atleast_1d(pbc))
    if len(raw) == 1:
        raw = [raw[0], raw[0], raw[0]]
    if len(raw) != 3:
        raise ValueError(f"pbc must be one flag or three, got {len(raw)}")
    return [bool(flag) for flag in raw]


def _voigt(stress: Sequence[float]) -> list[float]:
    """Row-major 3x3 to ASE Voigt (xx, yy, zz, yz, xz, xy)."""
    s = list(stress)
    return [s[0], s[4], s[8], s[5], s[2], s[1]]


def _params_bytes(params: Any) -> bytes:
    if isinstance(params, (bytes, bytearray, memoryview)):
        return bytes(params)
    with open(os.fspath(params), "rb") as handle:
        return handle.read()


def _param_value(key: str, value: Any) -> Any:
    if value is None:
        return None
    if key == "params" and isinstance(value, (bytearray, memoryview)):
        return bytes(value)
    if isinstance(value, bytes):
        return value
    return os.fspath(value)


class CPMDCCalculator(Calculator):
    """ASE calculator that evaluates one cpmdc session in process.

    Parameters
    ----------
    params:
        Serialized ``CPMDParams`` bytes, or a path to that file. The bytes
        are read once, when the parameter is set.
    library:
        Path of ``libcpmdc``. The environment search in
        :func:`cpmdc.find_library` is used when this is omitted.
    schema:
        Deprecated and ignored. A cell no longer needs the schema file.

    ``label``, ``directory``, and ``atoms`` are passed to ASE's
    ``Calculator``. ``restart`` is deprecated: cpmdc reads no restart
    file. Any other keyword raises ``TypeError``.

    ``parameters`` holds ``params`` and ``library``. :meth:`set` changes
    them, opens a new session, and clears the stored results. The
    calculator can be copied, deep-copied, and pickled; a copy opens its
    own session on its first step.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress"]
    default_parameters: dict[str, Any] = {"params": None, "library": None}
    # Charge and multiplicity come from CPMDParams, not from the Atoms.
    ignored_changes = {"initial_charges", "initial_magmoms"}
    discard_results_on_any_change = True

    def __init__(
        self,
        params: Any,
        *,
        library: Any = None,
        schema: Any = None,
        **kwargs: Any,
    ) -> None:
        unknown = [key for key in kwargs if key not in _ASE_KEYS]
        if unknown:
            raise TypeError(f"unknown calculator arguments: {sorted(unknown)}")
        if schema is not None:
            _warn_schema(stacklevel=3)
        if kwargs.get("restart") is not None:
            warnings.warn(
                "restart is ignored and will be removed: cpmdc reads no restart file",
                DeprecationWarning,
                stacklevel=2,
            )
        self._potential: Potential | None = None
        self._blob = b""
        self._composition: tuple[int, ...] | None = None
        self._no_stress = ""
        Calculator.__init__(self, **kwargs)
        self.set(params=params, library=library)

    @property
    def potential(self) -> Potential:
        """The rgpot-shaped potential this calculator evaluates.

        It is opened on first use after :meth:`close`, a copy, or an
        unpickle.
        """
        if self._potential is None:
            self._potential = Potential(self._blob, library=self.parameters.get("library"))
        return self._potential

    def set(self, **kwargs: Any) -> dict[str, Any]:
        """Change ``params`` or ``library`` and return the changed keys.

        A change opens a new session and clears the stored results. The
        old session stays in use when the new one cannot be opened.
        ``schema`` is deprecated and ignored. Another key raises
        ``TypeError``.
        """
        if "parameters" in kwargs:
            loaded = Parameters.read(kwargs.pop("parameters"))
            loaded.update(kwargs)
            kwargs = dict(loaded)
        if "schema" in kwargs:
            if kwargs.pop("schema") is not None:
                _warn_schema(stacklevel=3)
        unknown = [key for key in kwargs if key not in self.default_parameters]
        if unknown:
            raise TypeError(f"unknown calculator parameters: {sorted(unknown)}")
        kwargs = {key: _param_value(key, value) for key, value in kwargs.items()}
        if "params" in kwargs and kwargs["params"] is None:
            raise TypeError("params is required: a CPMDParams path or bytes")
        reopen = any(
            key in kwargs and not equal(kwargs[key], self.parameters.get(key))
            for key in ("params", "library")
        )
        if reopen:
            params = kwargs.get("params", self.parameters.get("params"))
            library = kwargs.get("library", self.parameters.get("library"))
            blob = _params_bytes(params)
            potential = Potential(blob, library=library)
        changed = Calculator.set(self, **kwargs)
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
        dct = Calculator.todict(self, skip_default=skip_default)
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
        properties: Sequence[str] = ("energy",),
        system_changes: Any = all_changes,
    ) -> None:
        """Fill ``self.results`` with every property of one step.

        Energy, free energy, and forces are always computed; stress is
        computed for a fully periodic cell. Asking for stress that the
        step cannot give raises ``PropertyNotImplementedError`` with the
        reason, without computing the step again.
        """
        if "stress" in properties and not system_changes and self._no_stress:
            raise PropertyNotImplementedError(self._no_stress)
        Calculator.calculate(self, atoms, properties, system_changes)
        assert self.atoms is not None
        cell = self._periodic_cell(self.atoms)
        if "stress" in properties and cell is None:
            self._no_stress = _NO_STRESS_CLUSTER
            raise PropertyNotImplementedError(_NO_STRESS_CLUSTER)
        numbers = [int(number) for number in self.atoms.get_atomic_numbers()]
        composition = tuple(sorted(numbers))
        if self._composition is not None and composition != self._composition:
            self._release()
        evaluated = self.potential.evaluate(self.atoms.get_positions(), numbers, cell)
        self._composition = composition
        energy = float(evaluated["energy"])
        self.results["energy"] = energy
        self.results["free_energy"] = energy
        self.results["forces"] = np.asarray(evaluated["forces"], dtype=float)
        if "stress" in evaluated:
            self.results["stress"] = np.asarray(_voigt(evaluated["stress"]), dtype=float)
            self._no_stress = ""
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
        return atoms.cell.array

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
