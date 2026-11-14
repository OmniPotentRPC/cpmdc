"""cpmdc inside GROMACS through the libcp2k API of its QM/MM module.

Two layers. The library cases load ``libcpmdc_gromacs`` (``CPMDC_GROMACS_LIB``
or ``build/libcpmdc_gromacs.so``) with ctypes and call it as
QMMMForceProvider does; they skip when the library is not built. The GROMACS
cases run grompp and mdrun of a double-precision ``gmx_d`` built with
``-DGMX_CP2K=ON`` against cpmdc's ``libcp2k.pc`` (``CPMDC_GMX``, or ``gmx_d``
on ``PATH``) and skip with the reason when there is none. References are
``CPMDCCalculator``, which reaches the library through ``ForceInput``
instead of ``cpmdc_session_calculate_au()``, and an independent NumPy Ewald
sum for the QM-MM Coulomb energy.
"""

from __future__ import annotations

import ctypes
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from conftest import ROOT
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV
from cpmdc.ase import CPMDCCalculator
from cpmdc.gromacs import (
    BOHR_TO_NM,
    HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM,
    HARTREE_TO_KJ_PER_MOL,
    qmi_input,
)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gromacs_case as gc  # noqa: E402


def library_path() -> Path | None:
    named = os.environ.get("CPMDC_GROMACS_LIB")
    candidate = Path(named) if named else ROOT / "build" / "libcpmdc_gromacs.so"
    return candidate if candidate.is_file() else None


LIB = library_path()
GMX, GMX_PROBLEM = gc.gmx_problem() if LIB else (None, "libcpmdc_gromacs is not built")
needs_lib = pytest.mark.skipif(LIB is None, reason="libcpmdc_gromacs is not built")
needs_gmx = pytest.mark.skipif(GMX is None, reason=GMX_PROBLEM or "GROMACS")

CUBE_NM = np.diag([1.2, 1.2, 1.2])
TRICLINIC_NM = np.array([[1.2, 0.0, 0.0], [0.16, 1.13, 0.0], [-0.09, 0.12, 1.25]])
# One QM water about the origin, where the reference well is shallow, and
# three MM waters, nm. GROMACS wraps the QM hydrogens below y = 0.
QMMM_NM = np.array(
    [
        [0.0, 0.0, 0.01173],
        [0.0, 0.07572, -0.04692],
        [0.0, -0.07572, -0.04692],
        [0.30, 0.03, 0.02],
        [0.36, 0.10, 0.03],
        [0.23, 0.07, -0.03],
        [0.02, 0.05, 0.31],
        [0.07, 0.12, 0.27],
        [-0.04, 0.01, 0.25],
        [0.75, 0.80, 0.85],
        [0.81, 0.87, 0.86],
        [0.68, 0.84, 0.80],
    ]
)
QMMM_BOX_NM = np.diag([1.5, 1.5, 1.5])
QM_CHARGES = [-0.834, 0.417, 0.417]


# --- conversions --------------------------------------------------------


def test_units_are_gromacs_units():
    assert BOHR_TO_NM == pytest.approx(BOHR_TO_ANGSTROM / 10.0, rel=1e-15)
    # c_hartree2Kj * c_avogadro and c_hartreeBohr2Md of GROMACS 2025
    # math/units.h, printed with %.17g.
    assert HARTREE_TO_KJ_PER_MOL == 2625.4996394798673
    assert HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM == 49614.752589206466
    # One e^2 / Bohr is GROMACS' c_one4PiEps0 / c_bohr2Nm in kJ/mol.
    assert HARTREE_TO_KJ_PER_MOL * BOHR_TO_NM == pytest.approx(138.93545764438196, rel=1e-9)


def test_qmi_input_has_the_keys_grompp_checks(tmp_path):
    text = qmi_input(
        "w.bin", embedding="mechanical", qm_charges=[-0.8, 0.4, 0.4], cell="none",
        ewald_rcut=0.6, ewald_precision=1e-9,
    )
    for line in (
        "RUN_TYPE ENERGY_FORCE", "CHARGE 0", "MULTIPLICITY 1", "COORD_FILE_FORMAT PDB",
        "CHARGE_EXTENDED TRUE", "COORD_FILE_NAME gromacs.pdb", "PARAMS w.bin", "CELL NONE",
        "EMBEDDING MECHANICAL", "QM_CHARGES -0.8 0.4 0.4", "EWALD_RCUT 0.6",
        "EWALD_PRECISION 1e-09",
    ):
        assert line in text
    with pytest.raises(ValueError, match="needs qm_charges"):
        qmi_input("w.bin", embedding="MECHANICAL")
    with pytest.raises(ValueError, match="needs embedding='MECHANICAL'"):
        qmi_input("w.bin", qm_charges=[0.0])
    with pytest.raises(ValueError, match="PERIODIC or NONE"):
        qmi_input("w.bin", cell="XY")


# --- the library, as QMMMForceProvider calls it ---------------------------


def write_pdb(path: Path, symbols: str, qm: list[bool], charges, x_nm) -> None:
    """The PDB of QMMMInputGenerator::generateCP2KPdb (GROMACS 2022-2026)."""
    lines = []
    for i, (sym, is_qm, q, row) in enumerate(zip(symbols, qm, charges, x_nm)):
        label = " QM     1     " if is_qm else " MM     2     "
        xyz = "%7.3f %7.3f %7.3f  1.00  0.00         " % tuple(10.0 * np.asarray(row))
        lines.append(f"ATOM  {i % 100000:5d}  {sym:<3s} {label}{xyz} {sym:<3s} {q:f}")
    path.write_text("\n".join(lines) + "\n")


def setup_env(tmp: Path, qmi: str, symbols: str, qm: list[bool], charges, x_nm) -> Path:
    write_pdb(tmp / "sys.pdb", symbols, qm, charges, x_nm)
    inp = tmp / "sys.inp"
    inp.write_text(qmi.replace("gromacs.pdb", str(tmp / "sys.pdb")))
    return inp


class Shim:
    """libcpmdc_gromacs loaded with ctypes; one force environment."""

    def __init__(self, inp: Path, out: Path):
        self.lib = ctypes.CDLL(str(LIB))
        self.env = ctypes.c_int(-1)
        self.lib.cp2k_init()
        self.lib.cp2k_create_force_env(ctypes.byref(self.env), str(inp).encode(), str(out).encode())
        natom = ctypes.c_int()
        self.lib.cp2k_get_natom(self.env, ctypes.byref(natom))
        self.natoms = natom.value

    def step(self, x_bohr, box_rows_bohr):
        n = 3 * self.natoms
        pos = (ctypes.c_double * n)(*np.asarray(x_bohr, float).ravel())
        self.lib.cp2k_set_positions(self.env, pos, n)
        # QMMMForceProvider: box_d[3 * i + j] = box[j][i].
        cell = (ctypes.c_double * 9)(*np.asarray(box_rows_bohr, float).T.ravel())
        self.lib.cp2k_set_cell(self.env, cell)
        self.lib.cp2k_calc_energy_force(self.env)
        energy = ctypes.c_double()
        self.lib.cp2k_get_potential_energy(self.env, ctypes.byref(energy))
        forces = (ctypes.c_double * n)()
        self.lib.cp2k_get_forces(self.env, forces, n)
        back = (ctypes.c_double * n)()
        self.lib.cp2k_get_positions(self.env, back, n)
        assert np.array_equal(np.array(back), np.asarray(x_bohr, float).ravel())
        return energy.value, np.array(forces).reshape(-1, 3)

    def close(self):
        self.lib.cp2k_destroy_force_env(self.env)
        self.lib.cp2k_finalize()


def ase_au(params: Path, symbols: str, x_bohr, cell_bohr=None):
    """CPMDCCalculator energy and forces, Hartree and Hartree/Bohr."""
    atoms = Atoms(symbols, positions=np.asarray(x_bohr) * BOHR_TO_ANGSTROM)
    if cell_bohr is not None:
        atoms.cell = np.asarray(cell_bohr) * BOHR_TO_ANGSTROM
        atoms.pbc = True
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy() / HARTREE_TO_EV
        forces = atoms.get_forces() / (HARTREE_TO_EV / BOHR_TO_ANGSTROM)
    return energy, forces


@needs_lib
@pytest.mark.parametrize("box", [CUBE_NM, TRICLINIC_NM], ids=["cube", "triclinic"])
def test_whole_system_step_matches_ase(box, periodic_params, tmp_path):
    x = np.asarray(gc.WATERS_NM) / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    inp = setup_env(tmp_path, qmi_input(periodic_params), gc.SYMBOLS, [True] * 6, [0.0] * 6, gc.WATERS_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(x, rows)
    shim.close()
    ref_e, ref_f = ase_au(periodic_params, gc.SYMBOLS, x, rows)
    assert energy == pytest.approx(ref_e, rel=1e-12, abs=1e-14)
    np.testing.assert_allclose(forces, ref_f, rtol=1e-11, atol=1e-13)
    log = (tmp_path / "sys.out").read_text()
    assert "atoms 6 qm 6 mm 0 cell PERIODIC embedding NONE" in log


@needs_lib
def test_qm_atoms_are_unwrapped_against_the_previous_step(periodic_params, tmp_path):
    """GROMACS wraps atoms into the box; cpmdc sees a continuous path."""
    rows = CUBE_NM / BOHR_TO_NM
    x = np.asarray(gc.WATERS_NM) / BOHR_TO_NM
    inp = setup_env(tmp_path, qmi_input(periodic_params), gc.SYMBOLS, [True] * 6, [0.0] * 6, gc.WATERS_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    shim.step(x, rows)
    moved = x.copy()
    moved[1] -= [0.0, 0.0, 0.06 / BOHR_TO_NM]  # H 1 crosses z = 0
    wrapped = moved.copy()
    wrapped[1] += rows[2]
    energy, forces = shim.step(wrapped, rows)
    shim.close()
    ref_e, ref_f = ase_au(periodic_params, gc.SYMBOLS, moved, rows)
    jump_e, _ = ase_au(periodic_params, gc.SYMBOLS, wrapped, rows)
    assert abs(jump_e - ref_e) > 1e-3
    assert energy == pytest.approx(ref_e, rel=1e-12)
    np.testing.assert_allclose(forces, ref_f, rtol=1e-11, atol=1e-13)


@needs_lib
def test_cell_none_makes_the_qm_region_whole(params, tmp_path):
    rows = QMMM_BOX_NM / BOHR_TO_NM
    x = QMMM_NM / BOHR_TO_NM
    whole = x.copy()
    split = whole.copy()
    split[2] += rows[1]  # H 2 wrapped above y = 0, as GROMACS does
    qm = [True] * 3 + [False] * 9
    inp = setup_env(
        tmp_path, qmi_input(params, embedding="NONE"), "OHH" * 4, qm, [0.0] * 3 + [0.417] * 9, QMMM_NM
    )
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(split, rows)
    shim.close()
    ref_e, ref_f = ase_au(params, "OHH", whole[:3])
    assert energy == pytest.approx(ref_e, rel=1e-12)
    np.testing.assert_allclose(forces[:3], ref_f, rtol=1e-11, atol=1e-13)
    np.testing.assert_array_equal(forces[3:], 0.0)
    assert "cell NONE embedding NONE" in (tmp_path / "sys.out").read_text()


def mechanical(tmp_path, params, x_nm, box_nm, **options):
    qm = [True] * 3 + [False] * (len(x_nm) - 3)
    charges = [0.0] * 3 + [gc.TIP3P[s] for s in "OHH" * (len(x_nm) // 3 - 1)]
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES, **options)
    tmp_path.mkdir(parents=True, exist_ok=True)
    inp = setup_env(tmp_path, qmi, "OHH" * (len(x_nm) // 3), qm, charges, x_nm)
    shim = Shim(inp, tmp_path / "sys.out")
    out = shim.step(np.asarray(x_nm) / BOHR_TO_NM, np.asarray(box_nm) / BOHR_TO_NM)
    shim.close()
    return out, np.array(charges[3:])


@needs_lib
@pytest.mark.parametrize("box", [QMMM_BOX_NM, QMMM_BOX_NM @ np.array([[1, 0, 0], [0.1, 1, 0], [-0.08, 0.12, 1]])],
                         ids=["cube", "triclinic"])
def test_mechanical_embedding_matches_numpy_ewald(box, params, tmp_path):
    (energy, forces), mm_q = mechanical(tmp_path, params, QMMM_NM, box)
    x = QMMM_NM / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    rc = gc.default_rcut(rows, BOHR_TO_NM)
    coul, fq, fm = gc.ewald_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q, rows, rc, 1e-8)
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    assert energy == pytest.approx(qm_e + coul, rel=1e-11, abs=1e-12)
    np.testing.assert_allclose(forces[:3], qm_f + fq, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(forces[3:], fm, rtol=1e-10, atol=1e-12)
    words = (tmp_path / "sys.out").read_text().split()
    assert float(words[words.index("rcut_bohr") + 1]) == pytest.approx(rc, rel=1e-15)
    assert float(words[words.index("alpha_per_bohr") + 1]) == pytest.approx(
        gc.ewald_alpha(rc, 1e-8), rel=1e-15
    )
    assert float(words[words.index("coulomb_qm_mm_hartree") + 1]) == pytest.approx(coul, rel=1e-11)


@needs_lib
def test_mechanical_embedding_ewald_is_converged_and_conservative(params, tmp_path):
    """The splitting does not change the sum; the forces are -dE/dx."""
    (e1, f1), _ = mechanical(tmp_path / "a", params, QMMM_NM, QMMM_BOX_NM)
    (e2, f2), _ = mechanical(tmp_path / "b", params, QMMM_NM, QMMM_BOX_NM, ewald_rcut=0.45, ewald_precision=1e-11)
    assert e1 == pytest.approx(e2, abs=1e-8)
    np.testing.assert_allclose(f1, f2, atol=1e-7)
    h = 1e-4
    for atom, axis in [(0, 0), (2, 1), (4, 2), (9, 0)]:
        plus, minus = QMMM_NM.copy(), QMMM_NM.copy()
        plus[atom, axis] += h * BOHR_TO_NM
        minus[atom, axis] -= h * BOHR_TO_NM
        (ep, _), _ = mechanical(tmp_path / f"p{atom}{axis}", params, plus, QMMM_BOX_NM)
        (em, _), _ = mechanical(tmp_path / f"m{atom}{axis}", params, minus, QMMM_BOX_NM)
        assert -(ep - em) / (2 * h) == pytest.approx(f1[atom, axis], abs=1e-7)


@needs_lib
def test_mechanical_embedding_without_a_box_is_a_direct_sum(params, tmp_path):
    (energy, forces), mm_q = mechanical(tmp_path, params, QMMM_NM, np.zeros((3, 3)))
    x = QMMM_NM / BOHR_TO_NM
    coul, fq, fm = gc.direct_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q)
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    assert energy == pytest.approx(qm_e + coul, rel=1e-12)
    np.testing.assert_allclose(forces[:3], qm_f + fq, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(forces[3:], fm, rtol=1e-11, atol=1e-13)


SCRIPT = textwrap.dedent(
    """
    import ctypes, sys
    lib = ctypes.CDLL(sys.argv[1])
    env = ctypes.c_int(-1)
    lib.cp2k_init()
    lib.cp2k_create_force_env(ctypes.byref(env), sys.argv[2].encode(), sys.argv[3].encode())
    n = int(sys.argv[4])
    pos = (ctypes.c_double * n)(*[0.1 * i for i in range(n)])
    lib.cp2k_set_positions(env, pos, n)
    lib.cp2k_calc_energy_force(env)
    """
)


def shim_fails(tmp_path, qmi, symbols="OHH", qm=(True,) * 3, charges=(0.0,) * 3, n_el=None):
    x = np.arange(3 * len(symbols), dtype=float).reshape(-1, 3) * 0.05
    inp = setup_env(tmp_path, qmi, symbols, list(qm), list(charges), x)
    done = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(LIB), str(inp), str(tmp_path / "sys.out"),
         str(n_el if n_el is not None else 3 * len(symbols))],
        capture_output=True, text=True, timeout=120,
    )
    assert done.returncode == 1, done.stderr
    assert "cpmdc-gromacs: " in done.stderr
    return done.stderr


@needs_lib
@pytest.mark.parametrize(
    "case, message",
    [
        ("generated", "has no &CPMDC section. Set qmmm-cp2k-qmmethod = INPUT"),
        ("no params", "&CPMDC needs PARAMS"),
        ("missing params", "cannot open the CPMDParams message"),
        ("link", "has a &LINK section"),
        ("unknown", "unknown &CPMDC keyword BASIS"),
        ("electrostatic", "EMBEDDING ELECTROSTATIC is not implemented"),
        ("unset embedding", "the MM atoms carry charges: set &CPMDC EMBEDDING"),
        ("charge count", "QM_CHARGES lists 2 charges for 3 QM atoms"),
        ("vsite", "QM atom 3 is X, not an element"),
        ("atom count", "cp2k_set_positions got 6 values; the PDB file has 3 atoms"),
        ("no qm", "has no QM atoms"),
    ],
)
def test_refused_inputs_stop_with_a_message(case, message, params, tmp_path):
    good = qmi_input(params)
    kwargs: dict = {}
    if case == "generated":
        qmi = "&GLOBAL\n RUN_TYPE ENERGY_FORCE\n&END GLOBAL\n&FORCE_EVAL\n METHOD QMMM\n" \
              " &SUBSYS\n  &TOPOLOGY\n   COORD_FILE_NAME gromacs.pdb\n" \
              "  &END TOPOLOGY\n &END SUBSYS\n&END FORCE_EVAL\n"
    elif case == "no params":
        qmi = good.replace(f"PARAMS {params}", "CELL NONE")
    elif case == "missing params":
        qmi = qmi_input(tmp_path / "absent.bin")
    elif case == "link":
        qmi = good.replace("&END CPMDC", "&END CPMDC\n  &LINK\n  &END LINK")
    elif case == "unknown":
        qmi = good.replace("&END CPMDC", "  BASIS DZVP\n&END CPMDC")
    elif case == "electrostatic":
        qmi = good.replace("&END CPMDC", "  EMBEDDING ELECTROSTATIC\n&END CPMDC")
    elif case == "unset embedding":
        qmi = good
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "charge count":
        qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=[-0.834, 0.417])
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "vsite":
        qmi = good
        kwargs = {"symbols": "OHX"}
    elif case == "atom count":
        qmi = good
        kwargs = {"n_el": 6}
    else:
        qmi = good
        kwargs = {"qm": (False,) * 3}
    stderr = shim_fails(tmp_path, qmi, **kwargs)
    assert message in stderr
    assert message in (tmp_path / "sys.out").read_text()


# --- GROMACS ---------------------------------------------------------------


def whole_system(tmp_path, params, box, **mdp):
    gc.write_system(
        tmp_path, positions=gc.WATERS_NM, box=box, qm_molecules=[0, 1],
        qmi=qmi_input(params), mdp={"coulombtype": "Cut-off", "rcoulomb": 0.5, "rvdw": 0.5, **mdp},
    )
    gc.run(GMX, tmp_path)


@needs_gmx
@pytest.mark.parametrize("box", [CUBE_NM, TRICLINIC_NM], ids=["cube", "triclinic"])
def test_gromacs_whole_system_step_matches_ase(box, periodic_params, tmp_path):
    whole_system(tmp_path, periodic_params, box)
    frame = gc.read_trr(tmp_path / "md.trr")[0]
    assert frame["double"]
    np.testing.assert_allclose(frame["box"], box, atol=1e-12)
    table = gc.energies(GMX, tmp_path, "Quantum-En.", "Potential")
    x = frame["x"] / BOHR_TO_NM
    ref_e, ref_f = ase_au(periodic_params, gc.SYMBOLS, x, frame["box"] / BOHR_TO_NM)
    assert table[0, 1] == pytest.approx(ref_e * HARTREE_TO_KJ_PER_MOL, rel=1e-10)
    assert table[0, 2] == pytest.approx(table[0, 1], rel=1e-14)
    np.testing.assert_allclose(
        frame["f"], ref_f * HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM, rtol=1e-10, atol=1e-9
    )


@needs_gmx
def test_gromacs_nve_conserves_energy(periodic_params, tmp_path):
    """The reference well pulls both waters through the box corner."""
    whole_system(
        tmp_path, periodic_params, CUBE_NM, nsteps=2000, nstxout=0, nstvout=0, nstfout=0,
        nstenergy=20, **{"gen-vel": "yes", "gen-temp": 300, "gen-seed": 4928},
    )
    table = gc.energies(GMX, tmp_path, "Potential", "Total-Energy")
    swing = table[:, 1].max() - table[:, 1].min()
    assert swing > 1.0
    assert np.abs(table[:, 2] - table[0, 2]).max() / swing < 2e-3


def qmmm_system(tmp_path, params, qmi, with_qmi=True, qm_lj=False, **mdp):
    gc.write_system(
        tmp_path, positions=QMMM_NM, box=QMMM_BOX_NM, qm_molecules=[0], qmi=qmi,
        mdp={"rcoulomb": 0.7, "rvdw": 0.7, **mdp}, qm_lj=qm_lj,
    )
    return gc.run(GMX, tmp_path, check=False, with_qmi=with_qmi)


@needs_gmx
def test_gromacs_mechanical_embedding_matches_references(params, tmp_path):
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES)
    done = qmmm_system(tmp_path, params, qmi)
    assert done.returncode == 0, done.stderr
    frame = gc.read_trr(tmp_path / "md.trr")[0]
    rows = frame["box"] / BOHR_TO_NM
    x = frame["x"] / BOHR_TO_NM
    x[1:3] -= np.rint((x[1:3] - x[0]) @ np.linalg.inv(rows)) @ rows
    mm_q = np.array([gc.TIP3P[s] for s in "OHH" * 3])
    coul, fq, _ = gc.ewald_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q, rows,
                                 gc.default_rcut(rows, BOHR_TO_NM), 1e-8)
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    table = gc.energies(GMX, tmp_path, "Quantum-En.")
    assert table[0, 1] == pytest.approx((qm_e + coul) * HARTREE_TO_KJ_PER_MOL, rel=1e-10)
    # The QM atom types carry no Lennard-Jones and grompp zeroes their
    # charges, so every force on them comes from libcpmdc_gromacs.
    np.testing.assert_allclose(
        frame["f"][:3], (qm_f + fq) * HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM, rtol=1e-10, atol=1e-9
    )


@needs_gmx
def test_gromacs_qmmm_nve_conserves_energy(params, tmp_path):
    """PME for MM-MM, the cpmdc Ewald sum for QM-MM, TIP3P LJ for both."""
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES)
    done = qmmm_system(
        tmp_path, params, qmi, qm_lj=True, nsteps=2000, nstxout=0, nstvout=0, nstfout=0, nstenergy=20,
        **{"gen-vel": "yes", "gen-temp": 300, "gen-seed": 4928},
    )
    assert done.returncode == 0, done.stderr
    table = gc.energies(GMX, tmp_path, "Potential", "Total-Energy")
    swing = table[:, 1].max() - table[:, 1].min()
    assert swing > 1.0
    assert np.abs(table[:, 2] - table[0, 2]).max() / swing < 5e-3


@needs_gmx
@pytest.mark.parametrize(
    "case, message",
    [
        ("generated", "has no &CPMDC section. Set qmmm-cp2k-qmmethod = INPUT"),
        ("electrostatic", "EMBEDDING ELECTROSTATIC is not implemented"),
        ("unset embedding", "the MM atoms carry charges: set &CPMDC EMBEDDING"),
    ],
)
def test_gromacs_unsupported_configs_stop_mdrun(case, message, params, tmp_path):
    if case == "generated":
        done = qmmm_system(tmp_path, params, "", with_qmi=False,
                           **{"qmmm-cp2k-qmmethod": "PBE"})
    elif case == "electrostatic":
        done = qmmm_system(tmp_path, params, qmi_input(params).replace(
            "&END CPMDC", "  EMBEDDING ELECTROSTATIC\n  &END CPMDC"))
    else:
        done = qmmm_system(tmp_path, params, qmi_input(params))
    assert done.returncode != 0
    assert "cpmdc-gromacs: " in done.stderr
    assert message in done.stderr
