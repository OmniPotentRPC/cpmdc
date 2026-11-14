"""Inputs and units for running GROMACS on cpmdc forces.

GROMACS (checked: 2024.6, 2025.5, 2026.4), double precision, configured
with ``-DGMX_CP2K=ON`` and ``CP2K_DIR``
pointing at the directory of cpmdc's ``pkgconfig/libcp2k.pc``, link
``libcpmdc_gromacs`` in place of CP2K. Its QM/MM module then calls
``cpmdc_session_calculate_au()`` once per MD step for the QM group, in
process, and adds the energy to the ``Quantum En.`` term. The QM group can
be the whole system.

:func:`qmi_input` writes the file grompp takes with ``-qmi`` when the mdp
sets ``qmmm-cp2k-qmmethod = INPUT``. GROMACS requires the CP2K keys it
checks (``CHARGE``, ``MULTIPLICITY``, ``RUN_TYPE ENERGY_FORCE``,
``COORD_FILE_FORMAT PDB``, ``CHARGE_EXTENDED TRUE``, ``COORD_FILE_NAME``);
the ``&CPMDC`` section carries the cpmdc settings that
``libcpmdc_gromacs`` reads.

Units: GROMACS works in nm, kJ/mol, and ps; its QM/MM module converts to
Bohr and Hartree with its own constants, which :data:`BOHR_TO_NM` and
:data:`HARTREE_TO_KJ_PER_MOL` reproduce (``src/gromacs/math/units.h``,
CODATA 2018).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

__all__ = [
    "BOHR_TO_NM",
    "HARTREE_TO_KJ_PER_MOL",
    "HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM",
    "qmi_input",
    "write_qmi",
]

# GROMACS' c_bohr2Nm, and c_hartree2Kj * c_avogadro built as units.h does:
# 2 R_inf h c, with h in kJ/mol ps and c in nm/ps.
BOHR_TO_NM = 0.0529177210903
_AVOGADRO = 6.02214076e23
_PLANCK = 6.62607015e-34 * _AVOGADRO / (1e-12 * 1e3)
_HARTREE_TO_KJ = 2.0 * 1.0973731568160e-02 * _PLANCK * 2.99792458e05 / _AVOGADRO
HARTREE_TO_KJ_PER_MOL = _HARTREE_TO_KJ * _AVOGADRO
HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM = HARTREE_TO_KJ_PER_MOL / BOHR_TO_NM

_CELLS = ("PERIODIC", "NONE")
_EMBEDDINGS = ("MECHANICAL", "NONE")


def _number(value: float) -> str:
    return repr(float(value))


def qmi_input(
    params: str | os.PathLike[str],
    *,
    cell: str | None = None,
    embedding: str | None = None,
    qm_charges: Sequence[float] | None = None,
    ewald_rcut: float | None = None,
    ewald_precision: float | None = None,
    qm_ranks: str | int | None = None,
    charge: int = 0,
    multiplicity: int = 1,
) -> str:
    """Text of the grompp ``-qmi`` file for ``libcpmdc_gromacs``.

    Parameters
    ----------
    params:
        Path of the serialized ``CPMDParams`` message, as mdrun will open it
        (relative paths are relative to mdrun's working directory).
    cell:
        ``"PERIODIC"`` sends the GROMACS box to cpmdc as the cell;
        ``"NONE"`` sends none, so the deck's cell applies. The library
        defaults to ``PERIODIC`` when every atom is QM and ``NONE`` otherwise.
    embedding:
        ``"MECHANICAL"`` adds the QM-MM Coulomb energy of ``qm_charges`` and
        the MM charges; ``"NONE"`` adds nothing. Required when MM atoms
        carry charges.
    qm_charges:
        Classical charge of each QM atom, in atom-index order, for
        mechanical embedding.
    ewald_rcut:
        Real-space cutoff of the QM-MM Ewald sum, nm. The library default is
        the smaller of 1 nm and half the shortest box height.
    ewald_precision:
        ``erfc(alpha * rcut)``, which also sets the reciprocal-space cutoff.
        The library default is 1e-8.
    qm_ranks:
        Under an MPI GROMACS, ``"ALL"`` (the library default) runs the QM
        calculation on every rank mdrun gives the interface, as one cpmdc
        calculator; ``1`` runs it on rank 0 alone while the others wait.
    charge, multiplicity:
        The values GROMACS reads. GROMACS warns when ``charge`` differs from
        the classical charge it removed from the QM atoms; the ``CPMDParams``
        deck sets the charge and multiplicity cpmdc uses.
    """
    lines = ["&CPMDC", f"  PARAMS {os.fspath(params)}"]
    if cell is not None:
        cell = cell.upper()
        if cell not in _CELLS:
            raise ValueError(f"cell is PERIODIC or NONE, not {cell!r}")
        lines.append(f"  CELL {cell}")
    if embedding is not None:
        embedding = embedding.upper()
        if embedding not in _EMBEDDINGS:
            raise ValueError(f"embedding is MECHANICAL or NONE, not {embedding!r}")
        lines.append(f"  EMBEDDING {embedding}")
    if qm_charges is not None:
        if embedding != "MECHANICAL":
            raise ValueError("qm_charges needs embedding='MECHANICAL'")
        values = [_number(q) for q in qm_charges]
        for start in range(0, len(values), 8):
            lines.append("  QM_CHARGES " + " ".join(values[start : start + 8]))
    elif embedding == "MECHANICAL":
        raise ValueError("embedding='MECHANICAL' needs qm_charges")
    if ewald_rcut is not None:
        lines.append(f"  EWALD_RCUT {_number(ewald_rcut)}")
    if ewald_precision is not None:
        lines.append(f"  EWALD_PRECISION {_number(ewald_precision)}")
    if qm_ranks is not None:
        ranks = str(qm_ranks).upper()
        if ranks not in ("ALL", "1"):
            raise ValueError(f"qm_ranks is 'ALL' or 1, not {qm_ranks!r}")
        lines.append(f"  QM_RANKS {ranks}")
    lines.append("&END CPMDC")
    head = [
        "&GLOBAL",
        "  RUN_TYPE ENERGY_FORCE",
        "&END GLOBAL",
        "&FORCE_EVAL",
        "  METHOD CPMDC",
    ]
    tail = [
        "  &DFT",
        f"    CHARGE {int(charge)}",
        f"    MULTIPLICITY {int(multiplicity)}",
        "  &END DFT",
        "  &SUBSYS",
        "    &TOPOLOGY",
        "      COORD_FILE_NAME gromacs.pdb",
        "      COORD_FILE_FORMAT PDB",
        "      CHARGE_EXTENDED TRUE",
        "    &END TOPOLOGY",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    return "\n".join(head + ["  " + line for line in lines] + tail) + "\n"


def write_qmi(path: str | os.PathLike[str], params: str | os.PathLike[str], **options) -> Path:
    """Write :func:`qmi_input` to ``path`` and return it."""
    target = Path(path)
    target.write_text(qmi_input(params, **options), encoding="utf-8")
    return target
