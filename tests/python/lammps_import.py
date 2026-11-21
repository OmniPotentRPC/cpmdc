"""Import the optional rgpot LAMMPS frontend without hiding broken installs."""

from __future__ import annotations

import importlib
import importlib.util
from types import ModuleType


def import_rgpot_lammps() -> ModuleType | None:
    """Return ``rgpot.lammps``, or ``None`` only when it is not installed."""
    try:
        spec = importlib.util.find_spec("rgpot.lammps")
    except ModuleNotFoundError as exc:
        if exc.name in {"rgpot", "rgpot.lammps"}:
            return None
        raise
    if spec is None:
        return None
    return importlib.import_module("rgpot.lammps")
