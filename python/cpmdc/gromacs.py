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
``libcpmdc_gromacs`` reads. ``links`` adds the ``&LINK`` sections that cap
a covalent bond between a QM atom and an MM atom with a hydrogen in the QM
calculation; GROMACS writes the same sections into the input it generates.

For a PBE or BLYP input generated from an mdp and stored in a tpr, set
``CPMDC_GROMACS_PARAMS`` when running mdrun. Its effective functional,
charge, and multiplicity must match the generated input. An explicit
``&CPMDC`` section takes precedence, and the charge and multiplicity of its
message must match the first ``CHARGE`` and ``MULTIPLICITY`` lines of the
file, which GROMACS reads whichever section holds them.

Units: GROMACS works in nm, kJ/mol, and ps; its QM/MM module converts to
Bohr and Hartree with its own constants, which :data:`BOHR_TO_NM` and
:data:`HARTREE_TO_KJ_PER_MOL` reproduce (``src/gromacs/math/units.h``,
CODATA 2018).
"""

from __future__ import annotations

import math
import os
from numbers import Real
from pathlib import Path
from typing import Any, Sequence

from cpmdc import _integer

__all__ = [
    "BOHR_TO_NM",
    "HARTREE_TO_KJ_PER_MOL",
    "HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM",
    "qmi_input",
    "write_qmi",
]

# GROMACS' c_bohr2Nm, and c_hartree2Kj * c_avogadro built as units.h does:
# 2 R_inf h c, with h in kJ/mol ps and c in nm/ps.
# src/gromacs/math/units.h (CODATA 2018) stores h in kJ mol^-1 ps, so a
# joule is 10^3 of those energy units and a second is 10^12 picoseconds.
BOHR_TO_NM = 0.0529177210903
_AVOGADRO = 6.02214076e23
_JOULES_PER_KILOJOULE = 10.0**3
_SECONDS_PER_PICOSECOND = 10.0**-12
_PLANCK = 6.62607015e-34 * _AVOGADRO / (_SECONDS_PER_PICOSECOND * _JOULES_PER_KILOJOULE)
_HARTREE_TO_KJ = 2.0 * 1.0973731568160e-02 * _PLANCK * 2.99792458e05 / _AVOGADRO
HARTREE_TO_KJ_PER_MOL = _HARTREE_TO_KJ * _AVOGADRO
HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM = HARTREE_TO_KJ_PER_MOL / BOHR_TO_NM

_CELLS = ("PERIODIC", "NONE")
_EMBEDDINGS = ("ELECTROSTATIC", "MECHANICAL", "NONE")
# The library refuses a Gaussian narrower than 1e-12 Angstrom.
_MIN_MM_WIDTH_NM = 1.0e-13


def _path_text(value: Any, what: str) -> str:
    if isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError(f"{what} must be a path string or path-like object, got {type(value).__name__}")
    text = os.fspath(value)
    if not isinstance(text, str):
        raise TypeError(f"{what} must resolve to text, got {type(text).__name__}")
    if not text:
        raise ValueError(f"{what} must not be empty")
    return text


def _number(value: Any, what: str, *, positive: bool = False) -> str:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{what} must be a real number, got {type(value).__name__}")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite, got {value!r}")
    if positive and number <= 0.0:
        raise ValueError(f"{what} must be positive, got {value!r}")
    return repr(number)


def _params_file_token(path: str) -> str:
    """PARAMS text stored for grompp.

    mdrun prints that file as a printf format after it substitutes the PDB
    path, so each ``%`` is written as ``%%`` and printf restores the path.
    """
    return path.replace("%", "%%")


def _gromacs_word_char(ch: str) -> bool:
    """Word characters of GROMACS ``replaceAllWords``: alnum, ``-`` and ``_``."""
    return ch.isalnum() or ch in "-_"


def _gromacs_word_occurs(token: str, text: str) -> bool:
    start = 0
    while True:
        found = text.find(token, start)
        if found < 0:
            return False
        end = found + len(token)
        before_ok = found == 0 or not _gromacs_word_char(text[found - 1])
        after_ok = end == len(text) or not _gromacs_word_char(text[end])
        if before_ok and after_ok:
            return True
        start = found + 1


def _coord_file_name(params_token: str) -> str:
    """A ``COORD_FILE_NAME`` that is not a whole word of ``params_token``.

    GROMACS replaces every whole word equal to that value with ``%s``
    before printf. ``gromacs.pdb`` is used when it does not occur.
    """
    number = 0
    while number < 1000:
        token = "gromacs.pdb" if number == 0 else f"gromacs-{number}.pdb"
        if not _gromacs_word_occurs(token, params_token):
            return token
        number += 1
    raise ValueError(
        "params path contains every COORD_FILE_NAME token qmi_input can choose"
    )


def qmi_input(
    params: str | os.PathLike[str],
    *,
    cell: str | None = None,
    embedding: str | None = None,
    qm_charges: Sequence[float] | None = None,
    mm_width: float | None = None,
    ewald_rcut: float | None = None,
    ewald_precision: float | None = None,
    epsilon_r: float | None = None,
    qm_ranks: str | int | None = None,
    links: Sequence[tuple[int, int]] | None = None,
    link_alpha: float | None = None,
    charge: int = 0,
    multiplicity: int = 1,
) -> str:
    """Text of the grompp ``-qmi`` file for ``libcpmdc_gromacs``.

    Parameters
    ----------
    params:
        Path of the serialized ``CPMDParams`` message, as mdrun will open it
        (relative paths are relative to mdrun's working directory). A ``%``
        in the path is written as ``%%`` so mdrun's printf of the input
        restores it. ``COORD_FILE_NAME`` is ``gromacs.pdb`` unless that
        whole word occurs in the written path, and then ``gromacs-1.pdb``,
        ``gromacs-2.pdb``, and so on, because GROMACS replaces every whole
        word equal to that name.
    cell:
        ``"PERIODIC"`` sends the GROMACS box to cpmdc as the cell;
        ``"NONE"`` sends none, so the deck's cell applies. The library
        defaults to ``PERIODIC`` when every atom is QM and ``NONE`` otherwise.
    embedding:
        ``"ELECTROSTATIC"`` lets Gaussian-smeared MM charges polarize the QM
        density, ``"MECHANICAL"`` adds the classical QM-MM Coulomb energy,
        and ``"NONE"`` adds nothing. Required when MM atoms carry charges.
    qm_charges:
        Classical charge of each QM atom, in atom-index order, for
        mechanical embedding.
    mm_width:
        Gaussian width in nm for electrostatic embedding (default 0.03 nm).
    ewald_rcut:
        Real-space cutoff of the QM-MM Ewald sum, nm. The library default is
        the smaller of 1 nm and half the shortest box height.
    ewald_precision:
        ``erfc(alpha * rcut)``, which also sets the reciprocal-space cutoff.
        The library default is ``1e-8``, assigned by
        ``tools/gromacs/cpmdc_gromacs.c`` when the key is omitted.
    epsilon_r:
        Relative dielectric constant (dimensionless, positive). The library
        default is ``1``; GROMACS does not forward ``.mdp`` ``epsilon-r``
        through ``libcp2k``, so pass the same value here when ``epsilon-r``
        differs from 1. Under ``"MECHANICAL"`` it divides the classical QM-MM
        Coulomb energy, forces, and virial by ``epsilon_r``; under
        ``"ELECTROSTATIC"`` each MM charge is passed as ``q / epsilon_r`` so
        the QM density and forces see the screened external charges. The QM
        self-energy is not scaled.
    qm_ranks:
        Under an MPI GROMACS, ``"ALL"`` (the library default) runs the QM
        calculation on every rank mdrun gives the interface, as one cpmdc
        calculator; ``1`` runs it on rank 0 alone while the others wait.
    links:
        The QM-MM covalent bonds to cap, as ``(qm_index, mm_index)`` pairs of
        1-based GROMACS atom numbers, QM atom first. Each becomes a ``&LINK``
        section of ``&QMMM``, and the library adds a hydrogen to the QM
        calculation on the bond, whose force it hands back to the two atoms.
    link_alpha:
        ``ALPHA_IMOMM`` of every link: the hydrogen sits at the QM atom plus
        the bond vector divided by ``link_alpha``. The library default is
        CP2K's 1.38. It must exceed 1, and needs ``links``.
    charge, multiplicity:
        The values GROMACS reads. GROMACS warns when ``charge`` differs from
        the classical charge it removed from the QM atoms. The ``CPMDParams``
        deck sets the charge and multiplicity cpmdc uses, and
        ``libcpmdc_gromacs`` stops ``mdrun`` when they differ from these two
        values, so pass the deck's charge and multiplicity for a deck that is
        not neutral and singlet. They are written in ``&DFT``, but GROMACS
        and the library read the first ``CHARGE`` and ``MULTIPLICITY`` line
        of the file in any section, ``&CPMDC`` included, so a line added
        above ``&DFT`` replaces them.
    """
    charge_value = _integer(charge, "charge")
    multiplicity_value = _integer(multiplicity, "multiplicity")
    if multiplicity_value < 1:
        raise ValueError(f"multiplicity must be 1 or more, got {multiplicity!r}")
    params_token = _params_file_token(_path_text(params, "params"))
    lines = ["&CPMDC", f"  PARAMS {params_token}"]
    if cell is not None:
        cell = cell.upper()
        if cell not in _CELLS:
            raise ValueError(f"cell is PERIODIC or NONE, not {cell!r}")
        lines.append(f"  CELL {cell}")
    if embedding is not None:
        embedding = embedding.upper()
        if embedding not in _EMBEDDINGS:
            raise ValueError(
                f"embedding is ELECTROSTATIC, MECHANICAL or NONE, not {embedding!r}"
            )
        lines.append(f"  EMBEDDING {embedding}")
    if embedding == "ELECTROSTATIC" and cell == "PERIODIC":
        raise ValueError("periodic electrostatic embedding is unsupported")
    if qm_charges is not None:
        if embedding != "MECHANICAL":
            raise ValueError("qm_charges needs embedding='MECHANICAL'")
        values = [_number(q, "qm_charges") for q in qm_charges]
        for start in range(0, len(values), 8):
            lines.append("  QM_CHARGES " + " ".join(values[start : start + 8]))
    elif embedding == "MECHANICAL":
        raise ValueError("embedding='MECHANICAL' needs qm_charges")
    if mm_width is not None:
        if embedding != "ELECTROSTATIC":
            raise ValueError("mm_width needs embedding='ELECTROSTATIC'")
        text = _number(mm_width, "mm_width", positive=True)
        if float(mm_width) < _MIN_MM_WIDTH_NM:
            raise ValueError(
                f"mm_width must be at least {_MIN_MM_WIDTH_NM:g} nm, got {mm_width!r}"
            )
        lines.append(f"  MM_WIDTH {text}")
    if ewald_rcut is not None:
        lines.append(f"  EWALD_RCUT {_number(ewald_rcut, 'ewald_rcut', positive=True)}")
    if ewald_precision is not None:
        lines.append(f"  EWALD_PRECISION {_number(ewald_precision, 'ewald_precision', positive=True)}")
    if epsilon_r is not None:
        text_eps = _number(epsilon_r, "epsilon_r", positive=True)
        if embedding not in ("MECHANICAL", "ELECTROSTATIC"):
            raise ValueError("epsilon_r needs embedding='MECHANICAL' or 'ELECTROSTATIC'")
        lines.append(f"  EPSILON_R {text_eps}")
    if qm_ranks is not None:
        ranks = str(qm_ranks).upper()
        if ranks not in ("ALL", "1"):
            raise ValueError(f"qm_ranks is 'ALL' or 1, not {qm_ranks!r}")
        lines.append(f"  QM_RANKS {ranks}")
    lines.append("&END CPMDC")
    qmmm: list[str] = []
    if links is not None:
        alpha = None if link_alpha is None else _number(link_alpha, "link_alpha", positive=True)
        if alpha is not None and float(alpha) <= 1.0:
            raise ValueError(f"link_alpha must exceed 1, got {link_alpha!r}")
        qmmm.append("  &QMMM")
        for pair in links:
            if len(pair) != 2:
                raise ValueError(f"each link is a (qm_index, mm_index) pair, got {pair!r}")
            qm_index = _integer(pair[0], "links qm_index")
            mm_index = _integer(pair[1], "links mm_index")
            if qm_index < 1 or mm_index < 1:
                raise ValueError(f"links are 1-based atom numbers, got {pair!r}")
            if qm_index == mm_index:
                raise ValueError(f"a link joins two atoms, got {pair!r}")
            qmmm += ["    &LINK", f"      QM_INDEX {qm_index}", f"      MM_INDEX {mm_index}"]
            if alpha is not None:
                qmmm.append(f"      ALPHA_IMOMM {alpha}")
            qmmm.append("    &END LINK")
        qmmm.append("  &END QMMM")
    elif link_alpha is not None:
        raise ValueError("link_alpha needs links")
    head = [
        "&GLOBAL",
        "  RUN_TYPE ENERGY_FORCE",
        "&END GLOBAL",
        "&FORCE_EVAL",
        "  METHOD CPMDC",
    ]
    tail = [
        "  &DFT",
        f"    CHARGE {charge_value}",
        f"    MULTIPLICITY {multiplicity_value}",
        "  &END DFT",
        *qmmm,
        "  &SUBSYS",
        "    &TOPOLOGY",
        f"      COORD_FILE_NAME {_coord_file_name(params_token)}",
        "      COORD_FILE_FORMAT PDB",
        "      CHARGE_EXTENDED TRUE",
        "    &END TOPOLOGY",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    return "\n".join(head + ["  " + line for line in lines] + tail) + "\n"


def write_qmi(
    path: str | os.PathLike[str],
    params: str | os.PathLike[str],
    *,
    cell: str | None = None,
    embedding: str | None = None,
    qm_charges: Sequence[float] | None = None,
    mm_width: float | None = None,
    ewald_rcut: float | None = None,
    ewald_precision: float | None = None,
    epsilon_r: float | None = None,
    qm_ranks: str | int | None = None,
    links: Sequence[tuple[int, int]] | None = None,
    link_alpha: float | None = None,
    charge: int = 0,
    multiplicity: int = 1,
) -> Path:
    """Write :func:`qmi_input` to ``path`` and return it; the options are those of :func:`qmi_input`."""
    target = Path(_path_text(path, "path"))
    text = qmi_input(
        params,
        cell=cell,
        embedding=embedding,
        qm_charges=qm_charges,
        mm_width=mm_width,
        ewald_rcut=ewald_rcut,
        ewald_precision=ewald_precision,
        epsilon_r=epsilon_r,
        qm_ranks=qm_ranks,
        links=links,
        link_alpha=link_alpha,
        charge=charge,
        multiplicity=multiplicity,
    )
    target.write_text(text, encoding="utf-8")
    return target
