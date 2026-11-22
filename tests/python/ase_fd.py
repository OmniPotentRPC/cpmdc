"""Finite-difference forces and stress on every ASE version cpmdc declares.

ASE 3.23 keeps the helpers in ``ase.calculators.test`` as ``numeric_forces``
and ``numeric_stress``, and names the step ``d``. ASE 3.24 moved them to
``ase.calculators.fd`` and renamed the step ``eps``. ASE 3.26 added
``force_consistent``.

The stress helpers agree. ``numeric_stress`` puts off-diagonal strain ``d``
on both sides of the deformation and divides by ``4 d V``.
``calculate_numerical_stress`` puts ``0.5 eps`` there and divides by
``2 eps V``. Those are the same tensor. Both index Voigt order as
``(xx, yy, zz, yz, xz, xy)``.

Before 3.26 the stress helper always reads ``free_energy`` and the force
helper always reads ``energy``. A caller that passes the other choice still
gets that hardcoded derivative. Calculators under test set the two energies
equal, so the derivative is the one the caller asked for.
"""

from __future__ import annotations

import inspect
from typing import Any, Callable

from ase import Atoms


def _module_functions() -> tuple[Callable[..., Any], Callable[..., Any]] | None:
    try:
        from ase.calculators.fd import (
            calculate_numerical_forces,
            calculate_numerical_stress,
        )
    except ImportError:
        return None
    return calculate_numerical_forces, calculate_numerical_stress


def _call(func: Callable[..., Any], atoms: Atoms, options: dict[str, Any]) -> Any:
    accepted = inspect.signature(func).parameters
    kwargs: dict[str, Any] = {}
    for key, value in options.items():
        if key in accepted:
            kwargs[key] = value
        elif key == "eps" and "d" in accepted:
            kwargs["d"] = value
    return func(atoms, **kwargs)


def calculate_numerical_forces(
    atoms: Atoms,
    eps: float,
    *,
    force_consistent: bool = False,
) -> Any:
    """Central-difference forces, eV/Angstrom.

    ``force_consistent`` is forwarded when this ASE accepts it. ASE 3.23
    through 3.25 always differentiate ``energy``.
    """

    found = _module_functions()
    if found is None:
        from ase.calculators.test import numeric_forces

        return _call(
            numeric_forces,
            atoms,
            {"eps": eps, "force_consistent": force_consistent},
        )
    return _call(
        found[0],
        atoms,
        {"eps": eps, "force_consistent": force_consistent},
    )


def calculate_numerical_stress(
    atoms: Atoms,
    eps: float,
    voigt: bool = True,
    *,
    force_consistent: bool = True,
) -> Any:
    """Central-difference stress, eV/Angstrom^3.

    ``force_consistent`` is forwarded when this ASE accepts it. ASE 3.23
    through 3.25 always differentiate ``free_energy``.
    """

    found = _module_functions()
    options = {"eps": eps, "voigt": voigt, "force_consistent": force_consistent}
    if found is None:
        from ase.calculators.test import numeric_stress

        return _call(numeric_stress, atoms, options)
    return _call(found[1], atoms, options)
