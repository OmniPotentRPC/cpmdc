"""Which ASE versions derive the CODATA 2018 Bohr and Hartree.

ASE builds Bohr and Hartree from ``mu0``. ASE 3.23 stores
``mu0 = 4 pi 1e-7`` in the 2018 table, the exact permeability from before
the 2019 SI revision. ASE 3.24 stores the CODATA 2018 permeability,
``1.25663706212e-6``. Only that later table reproduces the recommended
CODATA 2018 Bohr and Hartree to 1e-10. cpmdc's own constants are those
recommended values on every ASE version.
"""

from __future__ import annotations

import math
from typing import Any

from ase.units import create_units

# CODATA 2018 recommended permeability, Henry per metre.
CODATA2018_MU0 = 1.25663706212e-6


def units_2018() -> Any:
    return create_units("2018")


def derived_matches_recommended(units: Any | None = None) -> bool:
    """True when this ASE's 2018 ``mu0`` is the CODATA 2018 value."""

    table = units_2018() if units is None else units
    mu0 = float(table["_mu0"])
    pre2019 = 4.0 * math.pi * 1.0e-7
    return abs(mu0 - CODATA2018_MU0) < abs(mu0 - pre2019)
