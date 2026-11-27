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
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase_fd import calculate_numerical_stress
from conftest import ROOT
import cpmdc.gromacs as cpmdc_gromacs
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Session
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
OPENCPMD = LIB is not None and bool(ctypes.CDLL(str(LIB)).cpmdc_available())
# Strain derivatives to 1e-8 of the stress need the reference evaluator;
# test_gromacs_opencpmd.py runs these cases on OpenCPMD with a tight SCF.
reference_precision = pytest.mark.skipif(
    OPENCPMD, reason="OpenCPMD is linked; test_gromacs_opencpmd.py covers it"
)


def test_private_ewald_default_is_not_exposed() -> None:
    assert not hasattr(cpmdc_gromacs, "_DEFAULT_EWALD_PRECISION")


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
        ewald_rcut=0.6, ewald_precision=1e-9, epsilon_r=2.5,
    )
    for line in (
        "RUN_TYPE ENERGY_FORCE", "CHARGE 0", "MULTIPLICITY 1", "COORD_FILE_FORMAT PDB",
        "CHARGE_EXTENDED TRUE", "COORD_FILE_NAME gromacs.pdb", "PARAMS w.bin", "CELL NONE",
        "EMBEDDING MECHANICAL", "QM_CHARGES -0.8 0.4 0.4", "EWALD_RCUT 0.6",
        "EWALD_PRECISION 1e-09", "EPSILON_R 2.5",
    ):
        assert line in text
    with pytest.raises(ValueError, match="needs qm_charges"):
        qmi_input("w.bin", embedding="MECHANICAL")
    with pytest.raises(ValueError, match="needs embedding='MECHANICAL'"):
        qmi_input("w.bin", qm_charges=[0.0])
    with pytest.raises(ValueError, match="needs embedding='MECHANICAL'"):
        qmi_input("w.bin", epsilon_r=2.5)
    with pytest.raises(ValueError, match="PERIODIC or NONE"):
        qmi_input("w.bin", cell="XY")
    electrostatic = qmi_input(
        "w.bin", embedding="electrostatic", cell="none", mm_width=0.025, epsilon_r=2.5
    )
    assert "EMBEDDING ELECTROSTATIC" in electrostatic
    assert "MM_WIDTH 0.025" in electrostatic
    assert "EPSILON_R 2.5" in electrostatic
    with pytest.raises(ValueError, match="mm_width needs embedding='ELECTROSTATIC'"):
        qmi_input("w.bin", embedding="NONE", mm_width=0.03)
    with pytest.raises(ValueError, match="mm_width needs embedding='ELECTROSTATIC'"):
        qmi_input("w.bin", mm_width=0.03)
    for bad in (0.0, -0.01, float("nan"), float("inf"), 1.0e-14):
        with pytest.raises(ValueError, match="mm_width"):
            qmi_input("w.bin", embedding="ELECTROSTATIC", cell="none", mm_width=bad)


def _gromacs_word_char(ch: str) -> bool:
    """Word characters of GROMACS replaceAllWords: alnum, '-' and '_'."""
    return ch.isalnum() or ch in "-_"


def _replace_gromacs_words(text: str, old: str, new: str) -> str:
    """Replace whole words the way GROMACS replaceAllWords does."""
    parts: list[str] = []
    start = 0
    cursor = 0
    while True:
        found = text.find(old, cursor)
        if found < 0:
            parts.append(text[start:])
            return "".join(parts)
        end = found + len(old)
        before_ok = found == 0 or not _gromacs_word_char(text[found - 1])
        after_ok = end == len(text) or not _gromacs_word_char(text[end])
        if before_ok and after_ok:
            parts.append(text[start:found])
            parts.append(new)
            start = end
            cursor = end
        else:
            cursor = found + 1


def _coord_name(text: str) -> str:
    for line in text.splitlines():
        words = line.split()
        if words and words[0] == "COORD_FILE_NAME":
            return words[1]
    raise AssertionError("COORD_FILE_NAME missing")


def _after_mdrun(text: str, pdb: str) -> str:
    """Text mdrun writes: replace the coord-file word, then printf `%s`."""
    name = _coord_name(text)
    return _replace_gromacs_words(text, name, "%s") % (pdb,)


@pytest.mark.parametrize("path", ["calc%d.bin", "100%s.bin", "half%.bin"])
def test_a_percent_in_the_params_path_survives_mdrun(path: str):
    rewritten = _after_mdrun(qmi_input(path), "/tmp/sys.pdb")
    assert any(line.split() == ["PARAMS", path] for line in rewritten.splitlines())


def test_write_qmi_passes_the_mm_width_through(tmp_path):
    target = cpmdc_gromacs.write_qmi(
        tmp_path / "sys.qmi",
        "w.bin",
        embedding="ELECTROSTATIC",
        cell="NONE",
        mm_width=0.025,
        epsilon_r=2.5,
    )
    text = target.read_text(encoding="utf-8")
    assert "EMBEDDING ELECTROSTATIC" in text
    assert "MM_WIDTH 0.025" in text
    assert "EPSILON_R 2.5" in text


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


def generated_input(functional="BLYP", charge=0, multiplicity=1) -> str:
    """Relevant parts of GROMACS' generated CP2K QM/MM input."""
    return textwrap.dedent(
        f"""\
        &GLOBAL
          RUN_TYPE ENERGY_FORCE
        &END GLOBAL
        &FORCE_EVAL
          METHOD QMMM
          &DFT
            CHARGE {charge}
            MULTIPLICITY {multiplicity}
            &XC
              &XC_FUNCTIONAL {functional}
              &END XC_FUNCTIONAL
            &END XC
          &END DFT
          &QMMM
          &END QMMM
          &SUBSYS
            &TOPOLOGY
              COORD_FILE_NAME gromacs.pdb
              COORD_FILE_FORMAT PDB
              CHARGE_EXTENDED TRUE
            &END TOPOLOGY
          &END SUBSYS
        &END FORCE_EVAL
        """
    )


@pytest.fixture(scope="session")
def pbe_params(tmp_path_factory) -> Path:
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is needed to encode the PBE parser fixture")
    source = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text()
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=source.replace("BLYP", "PBE").encode(),
        stdout=subprocess.PIPE,
        check=True,
    )
    path = tmp_path_factory.mktemp("gromacs-pbe") / "water_pbe_params.bin"
    path.write_bytes(done.stdout)
    return path


@pytest.fixture(scope="session")
def charged_params(tmp_path_factory) -> Path:
    """The water message with the top-level charge +1 and multiplicity 2."""
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is needed to encode the charged parser fixture")
    source = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text()
    charged = source.replace("charge = 0,", "charge = 1,").replace(
        "multiplicity = 1,", "multiplicity = 2,", 1
    )
    assert charged != source
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=charged.encode(),
        stdout=subprocess.PIPE,
        check=True,
    )
    path = tmp_path_factory.mktemp("gromacs-charged") / "water_charged_params.bin"
    path.write_bytes(done.stdout)
    return path


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

    def step(self, x_bohr, box_rows_bohr, rows_layout=False):
        n = 3 * self.natoms
        pos = (ctypes.c_double * n)(*np.asarray(x_bohr, float).ravel())
        self.lib.cp2k_set_positions(self.env, pos, n)
        # QMMMForceProvider: box_d[3 * i + j] = box[j][i] in GROMACS 2024 and
        # 2025, box[i][j] in GROMACS 2026.
        rows = np.asarray(box_rows_bohr, float)
        cell = (ctypes.c_double * 9)(*(rows if rows_layout else rows.T).ravel())
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

    def stress(self):
        """cp2k_get_stress_tensor: the 3x3 tensor, Hartree/Bohr^3, and availability."""
        out = (ctypes.c_double * 9)(*[np.nan] * 9)
        available = ctypes.c_int(-1)
        self.lib.cp2k_get_stress_tensor(self.env, out, ctypes.byref(available))
        # column-major: element [a][b] is out[a + 3 * b]
        return np.array(out).reshape(3, 3).T, available.value

    def close(self):
        self.lib.cp2k_destroy_force_env(self.env)
        self.lib.cp2k_finalize()


EV_PER_ANGSTROM3 = HARTREE_TO_EV / BOHR_TO_ANGSTROM**3


def ase_stress(shim_stress):
    """ASE's tensile-positive eV/Angstrom^3 from CP2K's pressure-positive Hartree/Bohr^3."""
    return -np.asarray(shim_stress) * EV_PER_ANGSTROM3


class ShimCalculator(Calculator):
    """libcpmdc_gromacs as an ASE calculator, for ASE's numerical stress."""

    implemented_properties = ["energy", "free_energy"]

    def __init__(self, shim: "Shim"):
        super().__init__()
        self.shim = shim

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        # A sheared cell is not lower triangular, which only CP2K's layout can carry.
        energy, _ = self.shim.step(self.atoms.positions / BOHR_TO_ANGSTROM,
                                   np.asarray(self.atoms.cell) / BOHR_TO_ANGSTROM, rows_layout=True)
        # No electronic entropy: free_energy is energy. ASE 3.23 through 3.25
        # always differentiate free_energy and reject force_consistent.
        electronvolts = energy * HARTREE_TO_EV
        self.results = {"energy": electronvolts, "free_energy": electronvolts}


def numerical_stress(atoms, eps=1e-6):
    """ASE's central-difference strain derivative, (1/V) dE/d(strain), eV/Angstrom^3.

    ``force_consistent=False`` asks for the energy. ASE 3.26 and newer honor
    it. ASE 3.23 through 3.25 always read ``free_energy``; the calculators
    here set that equal to the energy, and the shear formula matches.
    """
    return calculate_numerical_stress(atoms, eps=eps, voigt=False, force_consistent=False)


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
    if not OPENCPMD:
        # The reference well is anchored at the origin, so a wrapped atom
        # changes its energy; a periodic OpenCPMD deck does not see the wrap.
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
    # Both sides sum the same terms, with the same alpha, rc and kmax, so they
    # differ by rounding only. With about 1e3 terms and partial sums of at most
    # 5e-3 Ha that is at most 1e3 * 2**-53 * 5e-3 = 6e-16 Ha on the energy and
    # as much on a force element. The runs measured at most 4e-17 on aarch64
    # and on x86-64, the latter also with contracted multiply-adds (-mfma -O3).
    # The tolerances below are three orders above that bound. One image pair
    # counted on one side only would add 1e-10 to 5e-10 Ha (erfc(alpha rc) =
    # 1e-8 times |q_i q_j| / rc), two orders above them; a k vector on the
    # cutoff adds about 1e-12 Ha, so the reciprocal cutoff is checked only to
    # that level.
    assert energy == pytest.approx(qm_e + coul, rel=1e-11, abs=1e-12)
    np.testing.assert_allclose(forces[:3], qm_f + fq, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(forces[3:], fm, rtol=1e-10, atol=1e-12)
    words = (tmp_path / "sys.out").read_text().split()
    # pytest.approx(rel=r) alone still accepts an absolute 1e-12, which is a
    # relative 7e-14 at rc and 3.5e-12 at alpha, so abs=0 says what is checked.
    # The shim forms the volume by cofactors and NumPy by LU; rc = V / 2|b x c|
    # then differs by a few units of rounding, 1.5e-15 relative in the cubic and
    # triclinic boxes, and alpha = f(rc) follows with the same relative error
    # (d ln alpha / d ln rc = -1, bisection resolution 2**-60). 1e-14 is about
    # six times the measured difference.
    assert float(words[words.index("rcut_bohr") + 1]) == pytest.approx(rc, rel=1e-14, abs=0.0)
    assert float(words[words.index("alpha_per_bohr") + 1]) == pytest.approx(
        gc.ewald_alpha(rc, 1e-8), rel=1e-14, abs=0.0
    )
    assert float(words[words.index("coulomb_qm_mm_hartree") + 1]) == pytest.approx(
        coul, rel=1e-11, abs=1e-12
    )


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


@needs_lib
@pytest.mark.parametrize("epsilon_r", [2.5, 78.4])
@pytest.mark.parametrize(
    "box",
    [QMMM_BOX_NM, QMMM_BOX_NM @ np.array([[1, 0, 0], [0.1, 1, 0], [-0.08, 0.12, 1]]), np.zeros((3, 3))],
    ids=["cube", "triclinic", "direct"],
)
def test_mechanical_embedding_epsilon_r_scales_coulomb_only(box, epsilon_r, params, tmp_path):
    """EPSILON_R scales only the mechanical-embedding QM-MM Coulomb energy and forces."""
    (e1, f1), mm_q = mechanical(tmp_path / "eps1", params, QMMM_NM, box)
    (e_eps, f_eps), _ = mechanical(tmp_path / "eps", params, QMMM_NM, box, epsilon_r=epsilon_r)
    x = QMMM_NM / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    if np.any(box):
        rc = gc.default_rcut(rows, BOHR_TO_NM)
        coul_1, fq_1, fm_1 = gc.ewald_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q, rows, rc, 1e-8)
        coul_e, fq_e, fm_e = gc.ewald_qm_mm(
            x[:3], QM_CHARGES, x[3:], mm_q, rows, rc, 1e-8, epsilon_r=epsilon_r
        )
    else:
        coul_1, fq_1, fm_1 = gc.direct_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q)
        coul_e, fq_e, fm_e = gc.direct_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q, epsilon_r=epsilon_r)
    assert coul_e == pytest.approx(coul_1 / epsilon_r, rel=1e-14)
    assert e_eps == pytest.approx(qm_e + coul_e, rel=1e-11, abs=1e-12)
    assert (e_eps - qm_e) == pytest.approx((e1 - qm_e) / epsilon_r, rel=1e-12, abs=1e-13)
    np.testing.assert_allclose(f_eps[:3], qm_f + fq_e, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(f_eps[:3] - qm_f, (f1[:3] - qm_f) / epsilon_r, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(f_eps[3:], fm_e, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(f_eps[3:], f1[3:] / epsilon_r, rtol=1e-12, atol=1e-13)
    assert abs(coul_1) > 1e-4
    words = (tmp_path / "eps" / "sys.out").read_text().split()
    assert float(words[words.index("epsilon_r") + 1]) == pytest.approx(epsilon_r, rel=1e-15)
    assert float(words[words.index("coulomb_qm_mm_hartree") + 1]) == pytest.approx(coul_e, rel=1e-11)


@needs_lib
def test_mechanical_embedding_epsilon_r_is_conservative(params, tmp_path):
    """Under EPSILON_R, forces remain -dE/dx for both QM and MM atoms."""
    epsilon_r = 3.25
    (e0, f0), _ = mechanical(tmp_path / "c", params, QMMM_NM, QMMM_BOX_NM, epsilon_r=epsilon_r)
    assert np.isfinite(e0)
    h = 1e-4
    for atom, axis in [(0, 0), (2, 1), (4, 2), (9, 0)]:
        plus, minus = QMMM_NM.copy(), QMMM_NM.copy()
        plus[atom, axis] += h * BOHR_TO_NM
        minus[atom, axis] -= h * BOHR_TO_NM
        (ep, _), _ = mechanical(
            tmp_path / f"p{atom}{axis}", params, plus, QMMM_BOX_NM, epsilon_r=epsilon_r
        )
        (em, _), _ = mechanical(
            tmp_path / f"m{atom}{axis}", params, minus, QMMM_BOX_NM, epsilon_r=epsilon_r
        )
        assert -(ep - em) / (2 * h) == pytest.approx(f0[atom, axis], abs=1e-7)


@needs_lib
def test_pdb_and_qm_charges_keep_double_precision(params, tmp_path):
    """libcpmdc_gromacs parses PDB MM charges and QM_CHARGES in double precision."""
    # Include a small charge below 5e-7 e that 6-decimal %lf rounding would zero,
    # plus 15-digit perturbations beyond both 6-decimal PDB and float32 precision.
    precise_qm = [-0.834000456789012, 0.417000228394506, 0.417000228394506]
    precise_mm = [
        -0.834123456789012,
        0.417061728394506,
        0.417061728394506,
        3.456789123456e-07,
        -1.728394561728e-07,
        -1.728394561728e-07,
        -0.834000000234567,
        0.417000000117283,
        0.417000000117284,
    ]
    lines = []
    for i, (sym, is_qm, q, row) in enumerate(
        zip("OHH" * 4, [True] * 3 + [False] * 9, [0.0] * 3 + precise_mm, QMMM_NM)
    ):
        label = " QM     1     " if is_qm else " MM     2     "
        xyz = "%7.3f %7.3f %7.3f  1.00  0.00         " % tuple(10.0 * np.asarray(row))
        lines.append(f"ATOM  {i % 100000:5d}  {sym:<3s} {label}{xyz} {sym:<3s} {q:.17g}")
    pdb_hi = tmp_path / "hi" / "sys.pdb"
    pdb_hi.parent.mkdir()
    pdb_hi.write_text("\n".join(lines) + "\n")
    inp_hi = pdb_hi.parent / "sys.inp"
    inp_hi.write_text(
        qmi_input(params, embedding="MECHANICAL", qm_charges=precise_qm).replace(
            "gromacs.pdb", str(pdb_hi)
        )
    )
    shim = Shim(inp_hi, pdb_hi.parent / "sys.out")
    x = QMMM_NM / BOHR_TO_NM
    e_hi, f_hi = shim.step(x, np.zeros((3, 3)))
    shim.close()

    # Compare against the same system rounded to 6 decimals (PDB %lf) and float32.
    rounded_mm = [round(q, 6) for q in precise_mm]
    rounded_qm = [round(q, 6) for q in precise_qm]
    (tmp_path / "lo").mkdir()
    inp_lo = setup_env(
        tmp_path / "lo",
        qmi_input(params, embedding="MECHANICAL", qm_charges=rounded_qm),
        "OHH" * 4,
        [True] * 3 + [False] * 9,
        [0.0] * 3 + rounded_mm,
        QMMM_NM,
    )
    shim_lo = Shim(inp_lo, tmp_path / "lo" / "sys.out")
    e_lo, _ = shim_lo.step(x, np.zeros((3, 3)))
    shim_lo.close()

    coul_hi, fq_hi, fm_hi = gc.direct_qm_mm(x[:3], precise_qm, x[3:], precise_mm)
    coul_lo, _, _ = gc.direct_qm_mm(x[:3], rounded_qm, x[3:], rounded_mm)
    coul_f32, _, _ = gc.direct_qm_mm(
        x[:3], np.asarray(precise_qm, np.float32), x[3:], np.asarray(precise_mm, np.float32)
    )
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    assert abs(coul_hi - coul_lo) > 5e-9
    assert abs(e_hi - e_lo) > 5e-9
    assert abs(coul_hi - coul_f32) > 1e-12
    assert abs(e_hi - (qm_e + coul_f32)) > 1e-12
    assert e_hi == pytest.approx(qm_e + coul_hi, rel=1e-14, abs=1e-15)
    assert e_lo == pytest.approx(qm_e + coul_lo, rel=1e-14, abs=1e-15)
    np.testing.assert_allclose(f_hi[:3], qm_f + fq_hi, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(f_hi[3:], fm_hi, rtol=1e-13, atol=1e-15)


@needs_lib
def test_electrostatic_embedding_accepts_charged_mm_atoms(params, tmp_path):
    qm = [True] * 3 + [False] * 9
    mm_q = np.array([gc.TIP3P[s] for s in "OHH" * 3])
    charges = np.concatenate([np.zeros(3), mm_q])
    qmi = qmi_input(
        params, embedding="ELECTROSTATIC", cell="NONE", mm_width=0.03
    )
    inp = setup_env(tmp_path, qmi, "OHH" * 4, qm, charges, QMMM_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(QMMM_NM / BOHR_TO_NM, QMMM_BOX_NM / BOHR_TO_NM)
    shim.close()
    rows = QMMM_BOX_NM / BOHR_TO_NM
    x = QMMM_NM / BOHR_TO_NM
    qm_x = x[:3].copy()
    qm_x[1:] -= np.rint((qm_x[1:] - qm_x[0]) @ np.linalg.inv(rows)) @ rows
    mm_x = x[3:].copy()
    mm_x -= np.rint((mm_x - qm_x[0]) @ np.linalg.inv(rows)) @ rows
    with Session(params) as session:
        session.set_external_charges(
            mm_x.reshape(-1) * BOHR_TO_ANGSTROM, mm_q, np.full(9, 0.30)
        )
        ref = session.energy_forces(
            qm_x * BOHR_TO_ANGSTROM, [8, 1, 1]
        )
        mm_forces = np.asarray(session.external_charge_forces_au(9)).reshape(-1, 3)
    assert energy == pytest.approx(ref.energy_h, rel=1e-12)
    np.testing.assert_allclose(forces[:3], ref.forces_h_bohr, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(forces[3:], mm_forces, rtol=1e-11, atol=1e-13)


def _electrostatic_step(tmp_path, params, x_nm, box_nm, *, mm_width=0.03, epsilon_r=None):
    qm = [True] * 3 + [False] * 9
    mm_q = np.array([gc.TIP3P[s] for s in "OHH" * 3])
    charges = np.concatenate([np.zeros(3), mm_q])
    tmp_path.mkdir(exist_ok=True)
    qmi = qmi_input(
        params, embedding="ELECTROSTATIC", cell="NONE", mm_width=mm_width, epsilon_r=epsilon_r
    )
    inp = setup_env(tmp_path, qmi, "OHH" * 4, qm, charges, x_nm)
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(x_nm / BOHR_TO_NM, box_nm / BOHR_TO_NM)
    stress, available = shim.stress()
    shim.close()
    return energy, forces, stress, available, mm_q


@needs_lib
@pytest.mark.parametrize("mm_width", [0.02, 0.03])
def test_electrostatic_embedding_epsilon_r_one_matches_default(mm_width, params, tmp_path):
    """EPSILON_R=1 matches omitting EPSILON_R for energy, forces, and virial."""
    e_def, f_def, s_def, avail_def, _ = _electrostatic_step(
        tmp_path / "default", params, QMMM_NM, QMMM_BOX_NM, mm_width=mm_width
    )
    e_one, f_one, s_one, avail_one, _ = _electrostatic_step(
        tmp_path / "one", params, QMMM_NM, QMMM_BOX_NM, mm_width=mm_width, epsilon_r=1.0
    )
    assert avail_def == 1 and avail_one == 1
    assert e_one == pytest.approx(e_def, rel=1e-14, abs=1e-15)
    np.testing.assert_allclose(f_one, f_def, rtol=1e-14, atol=1e-15)
    np.testing.assert_allclose(s_one, s_def, rtol=1e-14, atol=1e-15)
    words = (tmp_path / "one" / "sys.out").read_text().split()
    assert float(words[words.index("epsilon_r") + 1]) == pytest.approx(1.0, rel=1e-15)


@needs_lib
@pytest.mark.parametrize("epsilon_r", [2.5, 78.4])
@pytest.mark.parametrize("mm_width", [0.02, 0.06])
def test_electrostatic_embedding_epsilon_r_and_mm_width(mm_width, epsilon_r, params, tmp_path):
    """Under ELECTROSTATIC, EPSILON_R passes q / epsilon_r with MM_WIDTH to cpmdc.

    The QM density sees each MM charge as q / epsilon_r while the QM
    self-energy stays unscaled. In the reference evaluator, the embedding
    term is pure quadratic polarization (-0.5 * alpha * |E|^2), so the
    embedding energy shift, QM embedding force shift, and MM-site forces scale
    by 1 / epsilon_r^2, matching Session(params) with mm_q / epsilon_r.
    """
    e_bare, f_bare = ase_au(params, "OHH", QMMM_NM[:3] / BOHR_TO_NM)
    e1, f1, _, _, mm_q = _electrostatic_step(
        tmp_path / "eps1", params, QMMM_NM, QMMM_BOX_NM, mm_width=mm_width, epsilon_r=1.0
    )
    e_eps, f_eps, _, _, _ = _electrostatic_step(
        tmp_path / "eps", params, QMMM_NM, QMMM_BOX_NM, mm_width=mm_width, epsilon_r=epsilon_r
    )
    rows = QMMM_BOX_NM / BOHR_TO_NM
    x = QMMM_NM / BOHR_TO_NM
    qm_x = x[:3].copy()
    qm_x[1:] -= np.rint((qm_x[1:] - qm_x[0]) @ np.linalg.inv(rows)) @ rows
    mm_x = x[3:].copy()
    mm_x -= np.rint((mm_x - qm_x[0]) @ np.linalg.inv(rows)) @ rows
    with Session(params) as session:
        session.set_external_charges(
            mm_x.reshape(-1) * BOHR_TO_ANGSTROM,
            mm_q / epsilon_r,
            np.full(len(mm_q), mm_width * 10.0),
        )
        ref = session.energy_forces(qm_x * BOHR_TO_ANGSTROM, [8, 1, 1])
        ref_mm_forces = np.asarray(session.external_charge_forces_au(len(mm_q))).reshape(-1, 3)

    assert e_eps == pytest.approx(ref.energy_h, rel=1e-12, abs=1e-13)
    np.testing.assert_allclose(f_eps[:3], ref.forces_h_bohr, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(f_eps[3:], ref_mm_forces, rtol=1e-11, atol=1e-13)

    assert abs(e1 - e_bare) > 1e-8
    assert np.max(np.abs(f1[3:])) > 1e-8
    assert (e_eps - e_bare) == pytest.approx((e1 - e_bare) / (epsilon_r**2), rel=1e-11, abs=1e-14)
    np.testing.assert_allclose(
        f_eps[:3] - f_bare, (f1[:3] - f_bare) / (epsilon_r**2), rtol=1e-11, atol=1e-14
    )
    np.testing.assert_allclose(f_eps[3:], f1[3:] / (epsilon_r**2), rtol=1e-11, atol=1e-14)
    words = (tmp_path / "eps" / "sys.out").read_text().split()
    assert float(words[words.index("epsilon_r") + 1]) == pytest.approx(epsilon_r, rel=1e-15)


@needs_lib
def test_electrostatic_embedding_epsilon_r_is_conservative(params, tmp_path):
    """Under ELECTROSTATIC with EPSILON_R and MM_WIDTH, forces equal -dE/dx."""
    epsilon_r = 3.25
    mm_width = 0.025
    box_nm = np.diag([1.7, 1.7, 1.7])
    e0, f0, _, _, _ = _electrostatic_step(
        tmp_path / "c", params, QMMM_NM, box_nm, mm_width=mm_width, epsilon_r=epsilon_r
    )
    e_other, _, _, _, _ = _electrostatic_step(
        tmp_path / "w", params, QMMM_NM, box_nm, mm_width=0.06, epsilon_r=epsilon_r
    )
    assert np.isfinite(e0)
    assert abs(e0 - e_other) > 1e-8
    h = 1e-4
    for atom, axis in [(0, 0), (2, 1), (4, 2), (9, 0)]:
        plus, minus = QMMM_NM.copy(), QMMM_NM.copy()
        plus[atom, axis] += h * BOHR_TO_NM
        minus[atom, axis] -= h * BOHR_TO_NM
        ep, _, _, _, _ = _electrostatic_step(
            tmp_path / f"p{atom}{axis}", params, plus, box_nm,
            mm_width=mm_width, epsilon_r=epsilon_r,
        )
        em, _, _, _, _ = _electrostatic_step(
            tmp_path / f"m{atom}{axis}", params, minus, box_nm,
            mm_width=mm_width, epsilon_r=epsilon_r,
        )
        assert -(ep - em) / (2 * h) == pytest.approx(f0[atom, axis], abs=1e-7)


@needs_lib
def test_electrostatic_embedding_with_link_applies_z1_and_distributes_link_force(
    params, tmp_path
):
    """Under ELECTROSTATIC + &LINK, the link MM host charge is zeroed (Z1) and link forces map to QM/MM hosts."""
    qm = [True] * 3 + [False] * 9
    mm_q = np.array([gc.TIP3P[s] for s in "OHH" * 3])
    charges = np.concatenate([np.zeros(3), mm_q])
    link_alpha = 1.4
    qmi = qmi_input(
        params,
        embedding="ELECTROSTATIC",
        cell="NONE",
        mm_width=0.03,
        links=[(1, 4)],
        link_alpha=link_alpha,
    )
    inp = setup_env(tmp_path, qmi, "OHH" * 4, qm, charges, QMMM_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(QMMM_NM / BOHR_TO_NM, QMMM_BOX_NM / BOHR_TO_NM)
    shim.close()

    rows = QMMM_BOX_NM / BOHR_TO_NM
    inv_rows = np.linalg.inv(rows)
    x = QMMM_NM / BOHR_TO_NM
    qm_x = x[:3].copy()
    qm_x[1:] -= np.rint((qm_x[1:] - qm_x[0]) @ inv_rows) @ rows
    link_mm = x[3].copy()
    link_mm -= np.rint((link_mm - qm_x[0]) @ inv_rows) @ rows
    h_link = qm_x[0] + (link_mm - qm_x[0]) / link_alpha
    calc_x = np.vstack([qm_x, h_link])

    mm_x = x[3:].copy()
    mm_x -= np.rint((mm_x - qm_x[0]) @ inv_rows) @ rows
    expected_q = mm_q.copy()
    expected_q[0] = 0.0

    with Session(params) as session:
        session.set_external_charges(
            mm_x.reshape(-1) * BOHR_TO_ANGSTROM, expected_q, np.full(9, 0.30)
        )
        ref = session.energy_forces(calc_x * BOHR_TO_ANGSTROM, [8, 1, 1, 1])
        mm_forces = np.asarray(session.external_charge_forces_au(9)).reshape(-1, 3)

    a = 1.0 / link_alpha
    qm_ref_forces = np.asarray(ref.forces_h_bohr, dtype=float)
    expected_forces = np.zeros_like(forces)
    expected_forces[:3] = qm_ref_forces[:3]
    expected_forces[0] += (1.0 - a) * qm_ref_forces[3]
    expected_forces[3:] = mm_forces
    expected_forces[3] += a * qm_ref_forces[3]

    assert energy == pytest.approx(ref.energy_h, rel=1e-12)
    np.testing.assert_allclose(forces, expected_forces, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(mm_forces[0], 0.0, atol=1e-15)
    assert np.linalg.norm(forces[3]) > 1e-6


# --- stress -------------------------------------------------------------------


@needs_lib
@reference_precision
@pytest.mark.parametrize("box", [CUBE_NM, TRICLINIC_NM], ids=["cube", "triclinic"])
def test_whole_system_stress_matches_ase(box, periodic_params, tmp_path):
    """CELL PERIODIC: cpmdc's stress, against CPMDCCalculator and a strain difference."""
    x = np.asarray(gc.WATERS_NM) / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    inp = setup_env(tmp_path, qmi_input(periodic_params), gc.SYMBOLS, [True] * 6, [0.0] * 6, gc.WATERS_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    shim.step(x, rows)
    stress, available = shim.stress()
    shim.close()
    assert available == 1
    atoms = Atoms(gc.SYMBOLS, positions=x * BOHR_TO_ANGSTROM, cell=rows * BOHR_TO_ANGSTROM, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        analytic = atoms.get_stress(voigt=False)
        numeric = numerical_stress(atoms)
    ours = ase_stress(stress)
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ours, analytic, rtol=1e-13, atol=1e-15 * scale)
    np.testing.assert_allclose(ours, numeric, rtol=0, atol=1e-8 * scale)
    assert "first stress request: available 1" in (tmp_path / "sys.out").read_text()


@needs_lib
def test_gromacs_2026_cell_layout_reads_the_same_box(periodic_params, tmp_path):
    """GROMACS 2026 sends the lattice vectors one after the other.

    The reference energy depends on the cell only through sum h_ij^2, which
    a transposed cell keeps; its stress, through h h^T, does not.
    """
    x = np.asarray(gc.WATERS_NM) / BOHR_TO_NM
    rows = TRICLINIC_NM / BOHR_TO_NM
    out = {}
    for layout in (False, True):
        work = tmp_path / str(layout)
        work.mkdir()
        inp = setup_env(work, qmi_input(periodic_params), gc.SYMBOLS, [True] * 6, [0.0] * 6, gc.WATERS_NM)
        shim = Shim(inp, work / "sys.out")
        out[layout] = (*shim.step(x, rows, rows_layout=layout), shim.stress()[0])
        shim.close()
    assert out[True][0] == out[False][0]
    np.testing.assert_array_equal(out[True][1], out[False][1])
    np.testing.assert_array_equal(out[True][2], out[False][2])
    atoms = Atoms(gc.SYMBOLS, positions=x * BOHR_TO_ANGSTROM, cell=rows * BOHR_TO_ANGSTROM, pbc=True)
    transposed = atoms.copy()
    transposed.set_cell(rows.T * BOHR_TO_ANGSTROM)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        ref = atoms.get_stress(voigt=False)
    with CPMDCCalculator(periodic_params) as calc:
        transposed.calc = calc
        wrong = transposed.get_stress(voigt=False)
    scale = np.abs(ref).max()
    np.testing.assert_allclose(ase_stress(out[True][2]), ref, rtol=1e-13, atol=1e-15 * scale)
    assert not np.allclose(wrong, ref, rtol=1e-13, atol=1e-15 * scale)


MECHANICAL_TRICLINIC_NM = QMMM_BOX_NM @ np.array([[1, 0, 0], [0.1, 1, 0], [-0.08, 0.12, 1]])


@needs_lib
@reference_precision
@pytest.mark.parametrize("embedding", ["MECHANICAL", "NONE"])
@pytest.mark.parametrize("box", [QMMM_BOX_NM, MECHANICAL_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_qmmm_stress_matches_ase_numerical_stress(box, embedding, params, tmp_path):
    """CELL NONE: sum x (x) f of the QM atoms plus the Ewald strain derivative.

    The reference is ASE's central difference of libcpmdc_gromacs' own
    energy under strains of the cell that scale every atom with it, as
    GROMACS' pressure coupling does. A fixed EWALD_RCUT keeps the splitting
    fixed under the strain.
    """
    mm_q = [gc.TIP3P[s] for s in "OHH" * 3] if embedding == "MECHANICAL" else [0.0] * 9
    options = {"embedding": embedding, "ewald_rcut": 0.7, "ewald_precision": 1e-12}
    if embedding == "MECHANICAL":
        options["qm_charges"] = QM_CHARGES
    inp = setup_env(tmp_path, qmi_input(params, **options), "OHH" * 4, [True] * 3 + [False] * 9,
                    [0.0] * 3 + mm_q, QMMM_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    x = QMMM_NM / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    _, forces = shim.step(x, rows)
    stress, available = shim.stress()
    assert available == 1
    atoms = Atoms("OHH" * 4, positions=x * BOHR_TO_ANGSTROM, cell=rows * BOHR_TO_ANGSTROM, pbc=True)
    atoms.calc = ShimCalculator(shim)
    numeric = numerical_stress(atoms)
    shim.close()
    ours = ase_stress(stress)
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ours, numeric, rtol=0, atol=1e-8 * scale)
    if embedding == "NONE":
        volume = abs(np.linalg.det(rows))
        moment = x[:3].T @ forces[:3]
        expected = 0.5 * (moment + moment.T) / volume
        np.testing.assert_allclose(stress, expected, rtol=1e-14, atol=1e-14 * np.abs(expected).max())


@needs_lib
@reference_precision
@pytest.mark.parametrize("box", [QMMM_BOX_NM, MECHANICAL_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_qmmm_stress_with_epsilon_r_matches_ase_numerical_stress(box, params, tmp_path):
    """EPSILON_R scales the Ewald virial by 1/epsilon_r while keeping the QM virial unscaled."""
    epsilon_r = 2.5
    mm_q = [gc.TIP3P[s] for s in "OHH" * 3]
    x = QMMM_NM / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    stresses = {}
    for label, emb, extra in (
        ("qm_only", "NONE", {}),
        ("eps1", "MECHANICAL", {"qm_charges": QM_CHARGES}),
        ("eps", "MECHANICAL", {"qm_charges": QM_CHARGES, "epsilon_r": epsilon_r}),
    ):
        work = tmp_path / label
        work.mkdir()
        qmi = qmi_input(
            params,
            embedding=emb,
            ewald_rcut=0.7,
            ewald_precision=1e-12,
            **extra,
        )
        inp = setup_env(
            work,
            qmi,
            "OHH" * 4,
            [True] * 3 + [False] * 9,
            [0.0] * 3 + (mm_q if emb == "MECHANICAL" else [0.0] * 9),
            QMMM_NM,
        )
        shim = Shim(inp, work / "sys.out")
        shim.step(x, rows)
        stress, available = shim.stress()
        assert available == 1
        stresses[label] = stress
        if label == "eps":
            atoms = Atoms(
                "OHH" * 4,
                positions=x * BOHR_TO_ANGSTROM,
                cell=rows * BOHR_TO_ANGSTROM,
                pbc=True,
            )
            atoms.calc = ShimCalculator(shim)
            numeric = numerical_stress(atoms)
        shim.close()

    ours = ase_stress(stresses["eps"])
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ours, numeric, rtol=0, atol=1e-8 * scale)
    expected = stresses["qm_only"] + (stresses["eps1"] - stresses["qm_only"]) / epsilon_r
    np.testing.assert_allclose(
        stresses["eps"],
        expected,
        rtol=1e-13,
        atol=1e-14 * np.abs(expected).max(),
    )


@needs_lib
def test_stress_is_unavailable_without_a_box_or_a_periodic_deck(params, tmp_path):
    (tmp_path / "a").mkdir()
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES)
    inp = setup_env(tmp_path / "a", qmi, "OHH" * 4, [True] * 3 + [False] * 9,
                    [0.0] * 3 + [gc.TIP3P[s] for s in "OHH" * 3], QMMM_NM)
    shim = Shim(inp, tmp_path / "a" / "sys.out")
    shim.step(QMMM_NM / BOHR_TO_NM, np.zeros((3, 3)))
    stress, available = shim.stress()
    shim.close()
    assert available == 0
    np.testing.assert_array_equal(stress, 0.0)
    assert "first stress request: available 0 (no GROMACS box)" in (tmp_path / "a" / "sys.out").read_text()

    # The cluster deck with the GROMACS box as its cell: cpmdc computes no stress.
    (tmp_path / "b").mkdir()
    inp = setup_env(tmp_path / "b", qmi_input(params, cell="PERIODIC"), gc.SYMBOLS, [True] * 6,
                    [0.0] * 6, gc.WATERS_NM)
    shim = Shim(inp, tmp_path / "b" / "sys.out")
    shim.step(np.asarray(gc.WATERS_NM) / BOHR_TO_NM, CUBE_NM / BOHR_TO_NM)
    stress, available = shim.stress()
    shim.close()
    assert available == 0
    np.testing.assert_array_equal(stress, 0.0)
    assert "cpmdc computed no stress for this cell" in (tmp_path / "b" / "sys.out").read_text()


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


def shim_fails(
    tmp_path,
    qmi,
    symbols="OHH",
    qm=(True,) * 3,
    charges=(0.0,) * 3,
    n_el=None,
    gromacs_params=...,
):
    x = np.arange(3 * len(symbols), dtype=float).reshape(-1, 3) * 0.05
    inp = setup_env(tmp_path, qmi, symbols, list(qm), list(charges), x)
    env = dict(os.environ)
    if gromacs_params is None:
        env.pop("CPMDC_GROMACS_PARAMS", None)
    elif gromacs_params is not ...:
        env["CPMDC_GROMACS_PARAMS"] = str(gromacs_params)
    done = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(LIB), str(inp), str(tmp_path / "sys.out"),
         str(n_el if n_el is not None else 3 * len(symbols))],
        capture_output=True, text=True, timeout=120, env=env,
    )
    assert done.returncode == 1, done.stderr
    assert "cpmdc-gromacs: " in done.stderr
    return done.stderr


@needs_lib
@pytest.mark.parametrize("functional", ["PBE", "BLYP"])
def test_generated_input_matches_params(functional, pbe_params, params, tmp_path, monkeypatch):
    deck = pbe_params if functional == "PBE" else params
    monkeypatch.setenv("CPMDC_GROMACS_PARAMS", str(deck))
    inp = setup_env(
        tmp_path,
        generated_input(functional),
        "OHH",
        [True] * 3,
        [0.0] * 3,
        np.zeros((3, 3)),
    )
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
def test_explicit_cpmdc_section_wins_over_params_variable(params, tmp_path, monkeypatch):
    monkeypatch.setenv("CPMDC_GROMACS_PARAMS", str(tmp_path / "absent.bin"))
    inp = setup_env(
        tmp_path,
        qmi_input(params),
        "OHH",
        [True] * 3,
        [0.0] * 3,
        np.zeros((3, 3)),
    )
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"charge": 1}, "&DFT CHARGE 1 disagrees with &CPMDC PARAMS charge 0"),
        ({"multiplicity": 2}, "&DFT MULTIPLICITY 2 disagrees with &CPMDC PARAMS multiplicity 1"),
    ],
    ids=["charge", "multiplicity"],
)
def test_explicit_input_disagreeing_with_params_is_refused(kwargs, message, params, tmp_path):
    stderr = shim_fails(tmp_path, qmi_input(params, **kwargs))
    assert message in stderr
    assert message in (tmp_path / "sys.out").read_text()


@needs_lib
def test_explicit_input_defaults_do_not_hide_a_charged_params(charged_params, tmp_path):
    stderr = shim_fails(tmp_path, qmi_input(charged_params))
    assert "&DFT CHARGE 0 disagrees with &CPMDC PARAMS charge 1" in stderr


@needs_lib
def test_explicit_input_agreeing_with_a_charged_params_is_accepted(charged_params, tmp_path):
    inp = setup_env(
        tmp_path,
        qmi_input(charged_params, charge=1, multiplicity=2),
        "OHH",
        [True] * 3,
        [0.0] * 3,
        np.zeros((3, 3)),
    )
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
@pytest.mark.parametrize(
    "deck,charge,multiplicity",
    [("neutral", 0, 1), ("charged", 1, 2)],
)
def test_capped_region_compares_the_deck_charge(
    deck, charge, multiplicity, params, charged_params, tmp_path
):
    """A link hydrogen is not counted into the charge the deck is compared with."""
    blob = params if deck == "neutral" else charged_params
    symbols = "OHHOHH"
    qm = [True, True, True, False, False, False]
    charges = [0.0] * 6
    x = np.arange(18, dtype=float).reshape(6, 3) * 0.05
    matched = tmp_path / deck
    matched.mkdir()
    qmi = qmi_input(
        blob, embedding="NONE", links=[(1, 4)], charge=charge, multiplicity=multiplicity
    )
    shim = Shim(setup_env(matched, qmi, symbols, qm, charges, x), matched / "sys.out")
    shim.close()
    recounted = charge - 1
    refused = tmp_path / f"{deck}-recount"
    refused.mkdir()
    stderr = shim_fails(
        refused,
        qmi_input(
            blob,
            embedding="NONE",
            links=[(1, 4)],
            charge=recounted,
            multiplicity=multiplicity,
        ),
        symbols=symbols,
        qm=qm,
        charges=charges,
    )
    assert f"&DFT CHARGE {recounted} disagrees with &CPMDC PARAMS charge {charge}" in stderr


@needs_lib
@pytest.mark.parametrize(
    "line,message",
    [
        ("CHARGE one", "&DFT CHARGE value one is not an integer"),
        ("CHARGE 0.5", "&DFT CHARGE value 0.5 is not an integer"),
        ("MULTIPLICITY 1.0", "&DFT MULTIPLICITY value 1.0 is not an integer"),
    ],
)
def test_explicit_input_charge_must_be_an_integer(line, message, params, tmp_path):
    key = line.split()[0]
    qmi = qmi_input(params).replace(
        "    CHARGE 0\n" if key == "CHARGE" else "    MULTIPLICITY 1\n", f"    {line}\n"
    )
    assert line in qmi
    assert message in shim_fails(tmp_path, qmi)


@needs_lib
def test_explicit_input_uses_the_first_charge_and_multiplicity_tokens(params, tmp_path):
    qmi = qmi_input(params).replace(
        "    CHARGE 0\n", "    CHARGE 0\n    CHARGE 1\n"
    ).replace(
        "    MULTIPLICITY 1\n", "    MULTIPLICITY 1\n    MULTIPLICITY 2\n"
    )
    inp = setup_env(tmp_path, qmi, "OHH", [True] * 3, [0.0] * 3, np.zeros((3, 3)))
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
def test_explicit_input_first_tokens_outside_dft_define_qm_state(charged_params, tmp_path):
    qmi = qmi_input(charged_params, charge=0, multiplicity=1).replace(
        "&FORCE_EVAL\n", "&FORCE_EVAL\n    CHARGE 1\n    MULTIPLICITY 2\n", 1
    )
    inp = setup_env(tmp_path, qmi, "OHH", [True] * 3, [0.0] * 3, np.zeros((3, 3)))
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
def test_charge_inside_cpmdc_is_compared_before_dft(params, tmp_path):
    """GROMACS reads the first CHARGE line in any section, so &CPMDC counts."""
    qmi = qmi_input(params).replace("  &CPMDC\n", "  &CPMDC\n    CHARGE 1\n", 1)
    assert qmi.index("    CHARGE 1\n") < qmi.index("  &DFT\n")
    stderr = shim_fails(tmp_path, qmi)
    assert "&DFT CHARGE 1 disagrees with &CPMDC PARAMS charge 0" in stderr


@needs_lib
def test_charge_inside_cpmdc_agreeing_with_params_is_accepted(charged_params, tmp_path):
    qmi = qmi_input(charged_params).replace(
        "  &CPMDC\n", "  &CPMDC\n    CHARGE 1\n    MULTIPLICITY 2\n", 1
    )
    assert "    CHARGE 0\n" in qmi and "    MULTIPLICITY 1\n" in qmi
    inp = setup_env(tmp_path, qmi, "OHH", [True] * 3, [0.0] * 3, np.zeros((3, 3)))
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
def test_explicit_input_charge_mismatch_uses_first_token_anywhere(params, tmp_path):
    qmi = qmi_input(params).replace("&FORCE_EVAL\n", "&FORCE_EVAL\n    CHARGE 1\n", 1)
    stderr = shim_fails(tmp_path, qmi)
    assert "CHARGE 1 disagrees with &CPMDC PARAMS charge 0" in stderr


@needs_lib
@pytest.mark.parametrize(
    "qmi,message",
    [
        (
            generated_input("PBE"),
            "generated input functional PBE disagrees with "
            "CPMDC_GROMACS_PARAMS functional BLYP",
        ),
        (
            generated_input(charge=1),
            "generated input CHARGE 1 disagrees with CPMDC_GROMACS_PARAMS charge 0",
        ),
        (
            generated_input(multiplicity=2),
            "generated input MULTIPLICITY 2 disagrees with "
            "CPMDC_GROMACS_PARAMS multiplicity 1",
        ),
    ],
    ids=["functional", "charge", "multiplicity"],
)
def test_generated_input_mismatch_is_refused(qmi, message, params, tmp_path):
    stderr = shim_fails(tmp_path, qmi, gromacs_params=params)
    assert message in stderr


@needs_lib
def test_generated_input_without_params_variable_is_refused(tmp_path):
    stderr = shim_fails(tmp_path, generated_input(), gromacs_params=None)
    assert "CPMDC_GROMACS_PARAMS is not set" in stderr


@needs_lib
def test_generated_input_cannot_invent_qm_charges(params, tmp_path):
    stderr = shim_fails(
        tmp_path,
        generated_input(),
        symbols="OHHOHH",
        qm=(True,) * 3 + (False,) * 3,
        charges=(0.0,) * 3 + (-0.834, 0.417, 0.417),
        gromacs_params=params,
    )
    assert "charged MM atoms but no QM_CHARGES" in stderr


@needs_lib
def test_generated_input_reads_charge_tokens_outside_dft(params, tmp_path):
    qmi = generated_input().replace("&FORCE_EVAL\n", "&FORCE_EVAL\n    CHARGE 1\n", 1)
    stderr = shim_fails(tmp_path, qmi, gromacs_params=params)
    assert "generated input CHARGE 1 disagrees with CPMDC_GROMACS_PARAMS charge 0" in stderr


@needs_lib
def test_generated_input_rejects_extra_tokens(params, tmp_path):
    qmi = generated_input().replace("CHARGE 0", "CHARGE 0 1")
    stderr = shim_fails(tmp_path, qmi, gromacs_params=params)
    assert "CHARGE takes one value" in stderr


@needs_lib
def test_generated_input_rejects_truncated_functional_token(params, tmp_path):
    qmi = generated_input(functional="X" * 80)
    stderr = shim_fails(tmp_path, qmi, gromacs_params=params)
    assert "&XC_FUNCTIONAL token is too long" in stderr


@needs_lib
def test_generated_input_rejects_extra_xc_functional_tokens(params, tmp_path):
    qmi = generated_input("PBE").replace("&XC_FUNCTIONAL PBE", "&XC_FUNCTIONAL PBE EXTRA")
    stderr = shim_fails(tmp_path, qmi, gromacs_params=params)
    assert "&XC_FUNCTIONAL in &DFT takes one value" in stderr


@needs_lib
def test_explicit_cpmdc_ignores_multivalue_section_headers(params, tmp_path):
    qmi = qmi_input(params).replace(
        "&SUBSYS\n",
        "&SUBSYS\n    &KIND O H\n    &END KIND\n",
    ).replace(
        "&FORCE_EVAL\n",
        "&FORCE_EVAL\n  &EXTRA a b\n  &END EXTRA\n",
    ).replace(
        "  &DFT\n",
        "  &DFT\n    &XC\n      &XC_FUNCTIONAL PBE EXTRA\n      &END XC_FUNCTIONAL\n    &END XC\n",
    )
    inp = setup_env(
        tmp_path,
        qmi,
        "OHH",
        [True] * 3,
        [0.0] * 3,
        np.zeros((3, 3)),
    )
    shim = Shim(inp, tmp_path / "sys.out")
    shim.close()


@needs_lib
@pytest.mark.parametrize(
    "case, message",
    [
        ("generated", "CPMDC_GROMACS_PARAMS is not set"),
        ("no params", "&CPMDC needs PARAMS"),
        ("missing params", "cannot open the CPMDParams message"),
        ("link", "has a &LINK section without QM_INDEX"),
        ("unknown", "unknown &CPMDC keyword BASIS"),
        ("unset embedding", "the MM atoms carry charges: set &CPMDC EMBEDDING"),
        ("epsilon zero", "&CPMDC EPSILON_R must be positive"),
        ("epsilon negative", "&CPMDC EPSILON_R must be positive"),
        ("epsilon nan", "&CPMDC EPSILON_R: nan is not a number"),
        ("epsilon none", "&CPMDC EPSILON_R needs EMBEDDING MECHANICAL"),
        ("epsilon unset", "&CPMDC EPSILON_R needs EMBEDDING MECHANICAL"),
        ("pdb nan", "atom 4 has no finite charge"),
        ("width without electrostatic", "&CPMDC MM_WIDTH needs EMBEDDING ELECTROSTATIC"),
        ("tiny width", "&CPMDC MM_WIDTH must be a finite width of at least 1e-13 nm"),
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
        qmi = generated_input()
    elif case == "no params":
        qmi = good.replace(f"PARAMS {params}", "CELL NONE")
    elif case == "missing params":
        qmi = qmi_input(tmp_path / "absent.bin")
    elif case == "link":
        qmi = good.replace("&END CPMDC", "&END CPMDC\n  &LINK\n  &END LINK")
    elif case == "unknown":
        qmi = good.replace("&END CPMDC", "  BASIS DZVP\n&END CPMDC")
    elif case == "unset embedding":
        qmi = good
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "epsilon zero":
        qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES).replace(
            "&END CPMDC", "  EPSILON_R 0\n&END CPMDC"
        )
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "epsilon negative":
        qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES).replace(
            "&END CPMDC", "  EPSILON_R -2.5\n&END CPMDC"
        )
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "epsilon nan":
        qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES).replace(
            "&END CPMDC", "  EPSILON_R nan\n&END CPMDC"
        )
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (-0.834, 0.417, 0.417)}
    elif case == "epsilon none":
        qmi = good.replace("&END CPMDC", "  EMBEDDING NONE\n  EPSILON_R 4.0\n&END CPMDC")
    elif case == "epsilon unset":
        qmi = good.replace("&END CPMDC", "  EPSILON_R 4.0\n&END CPMDC")
    elif case == "pdb nan":
        qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES)
        kwargs = {"symbols": "OHHOHH", "qm": (True,) * 3 + (False,) * 3,
                  "charges": (0.0,) * 3 + (float("nan"), 0.417, 0.417)}
    elif case == "width without electrostatic":
        qmi = good.replace("&END CPMDC", "  EMBEDDING NONE\n  MM_WIDTH 0.03\n&END CPMDC")
    elif case == "tiny width":
        qmi = good.replace(
            "&END CPMDC", "  EMBEDDING ELECTROSTATIC\n  MM_WIDTH 1e-14\n&END CPMDC"
        )
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
    stderr = shim_fails(
        tmp_path, qmi, gromacs_params=None if case == "generated" else ..., **kwargs
    )
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
def test_gromacs_generated_input_uses_params_variable(periodic_params, tmp_path):
    """GROMACS' generated route translates the QM atoms into its QM box."""
    gc.write_system(
        tmp_path,
        positions=gc.WATERS_NM,
        box=CUBE_NM,
        qm_molecules=[0, 1],
        qmi="",
        mdp={
            "coulombtype": "Cut-off",
            "rcoulomb": 0.5,
            "rvdw": 0.5,
            "qmmm-cp2k-qmmethod": "BLYP",
            "qmmm-cp2k-qmcharge": 0,
            "qmmm-cp2k-qmmultiplicity": 1,
        },
    )
    done = gc.run(
        GMX,
        tmp_path,
        with_qmi=False,
        env_extra={"CPMDC_GROMACS_PARAMS": str(periodic_params)},
    )
    assert done.returncode == 0, done.stderr
    frame = gc.read_trr(tmp_path / "md.trr")[0]
    table = gc.energies(GMX, tmp_path, "Quantum-En.")
    assert np.isfinite(table[0, 1])
    assert np.isfinite(frame["f"]).all()
    log = (tmp_path / "md_cp2k.out").read_text()
    assert "atoms 6 qm 6 mm 0 cell PERIODIC embedding NONE" in log
    backend = "opencpmd" if OPENCPMD else "reference"
    assert f"first step: backend {backend}" in log


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


def qmmm_system(
    tmp_path, params, qmi, with_qmi=True, qm_lj=False, env_extra=None, **mdp
):
    gc.write_system(
        tmp_path, positions=QMMM_NM, box=QMMM_BOX_NM, qm_molecules=[0], qmi=qmi,
        mdp={"rcoulomb": 0.7, "rvdw": 0.7, **mdp}, qm_lj=qm_lj,
    )
    return gc.run(
        GMX,
        tmp_path,
        check=False,
        with_qmi=with_qmi,
        env_extra=env_extra,
    )


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
@pytest.mark.parametrize("box", [QMMM_BOX_NM, MECHANICAL_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_gromacs_epsilon_r_matches_references_and_strain_derivative(box, params, tmp_path):
    """GROMACS epsilon-r + &CPMDC EPSILON_R: Quantum En., forces, and virial match references."""
    epsilon_r = 2.5
    qmi = qmi_input(
        params,
        embedding="MECHANICAL",
        qm_charges=QM_CHARGES,
        ewald_rcut=0.7,
        ewald_precision=1e-12,
        epsilon_r=epsilon_r,
    )
    mdp = {
        "integrator": "md-vv",
        "rcoulomb": 0.7,
        "rvdw": 0.7,
        "fourier-nx": 24,
        "fourier-ny": 24,
        "fourier-nz": 24,
        "ewald-rtol": 1e-10,
        "epsilon-r": epsilon_r,
        "nstxout": 1,
        "nstvout": 0,
        "nstfout": 1,
    }

    def run_box(work, x_nm, box_nm, *, with_lj=True):
        gc.write_system(
            work, positions=x_nm, box=box_nm, qm_molecules=[0], qmi=qmi, mdp=mdp, qm_lj=with_lj
        )
        gc.run(GMX, work)
        return gc.energies(GMX, work, "Quantum-En.", "Potential")[0]

    # Step 0 energy & forces without QM LJ so frame["f"][:3] is pure libcpmdc_gromacs.
    e_row = run_box(tmp_path / "ef", QMMM_NM, box, with_lj=False)
    frame = gc.read_trr(tmp_path / "ef" / "md.trr")[0]
    rows = frame["box"] / BOHR_TO_NM
    x = frame["x"] / BOHR_TO_NM
    x[1:3] -= np.rint((x[1:3] - x[0]) @ np.linalg.inv(rows)) @ rows
    mm_q = np.array([gc.TIP3P[s] for s in "OHH" * 3])
    coul, fq, _ = gc.ewald_qm_mm(
        x[:3], QM_CHARGES, x[3:], mm_q, rows, 0.7 / BOHR_TO_NM, 1e-12, epsilon_r=epsilon_r
    )
    qm_e, qm_f = ase_au(params, "OHH", x[:3])
    assert e_row[1] == pytest.approx((qm_e + coul) * HARTREE_TO_KJ_PER_MOL, rel=1e-10)
    np.testing.assert_allclose(
        frame["f"][:3], (qm_f + fq) * HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM, rtol=1e-10, atol=1e-9
    )

    # Full strain derivative of GROMACS' Potential energy (MM-MM PME with epsilon-r + QM + QM-MM).
    x_nm = QMMM_NM + 0.4
    run_box(tmp_path / "0", x_nm, box)
    virial, _ = virial_and_pressure(tmp_path / "0")
    h = 1e-5
    pairs = [(a, b) for a in range(3) for b in range(a, 3)]
    derivative = np.zeros((3, 3))
    for a, b in pairs:
        plus = run_box(tmp_path / f"{a}{b}+", *strained(x_nm, box, a, b, h))[2]
        minus = run_box(tmp_path / f"{a}{b}-", *strained(x_nm, box, a, b, -h))[2]
        derivative[a, b] = (plus - minus) / (2 * h)
    scale = np.abs(virial).max()
    for a, b in pairs:
        assert virial[b, a] == pytest.approx(0.5 * derivative[a, b], abs=1e-8 * scale), (a, b)
    np.testing.assert_allclose(virial, virial.T, atol=1e-9 * scale)


@needs_gmx
def test_gromacs_pdb_rounds_mm_charges_to_six_decimals(params, tmp_path):
    """Upstream GROMACS writes md_cp2k.pdb MM charges with %lf (6 decimals)."""
    mm_top_charges = (-0.8341234, 0.4170617, 0.4170617)
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES)
    gc.write_system(
        tmp_path,
        positions=QMMM_NM,
        box=QMMM_BOX_NM,
        qm_molecules=[0],
        qmi=qmi,
        mdp={"rcoulomb": 0.7, "rvdw": 0.7},
        mm_charges=mm_top_charges,
    )
    gc.run(GMX, tmp_path)
    pdb_lines = [
        line.split()[-1]
        for line in (tmp_path / "md_cp2k.pdb").read_text().splitlines()
        if line.startswith("ATOM")
    ]
    assert pdb_lines[:3] == ["0.000000", "0.000000", "0.000000"]
    assert pdb_lines[3:6] == ["-0.834123", "0.417062", "0.417062"]

    frame = gc.read_trr(tmp_path / "md.trr")[0]
    rows = frame["box"] / BOHR_TO_NM
    x = frame["x"] / BOHR_TO_NM
    x[1:3] -= np.rint((x[1:3] - x[0]) @ np.linalg.inv(rows)) @ rows
    rc = gc.default_rcut(rows, BOHR_TO_NM)
    coul_pdb, _, _ = gc.ewald_qm_mm(
        x[:3], QM_CHARGES, x[3:], np.array([float(s) for s in pdb_lines[3:]]), rows, rc, 1e-8
    )
    coul_top, _, _ = gc.ewald_qm_mm(
        x[:3], QM_CHARGES, x[3:], np.array(list(mm_top_charges) * 3), rows, rc, 1e-8
    )
    qm_e, _ = ase_au(params, "OHH", x[:3])
    q_en = gc.energies(GMX, tmp_path, "Quantum-En.")[0, 1]
    assert abs(coul_pdb - coul_top) > 1e-8
    assert q_en == pytest.approx((qm_e + coul_pdb) * HARTREE_TO_KJ_PER_MOL, rel=1e-11)
    assert abs(q_en - (qm_e + coul_top) * HARTREE_TO_KJ_PER_MOL) > 1e-5


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


TENSOR = ("XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ")
EV_TO_KJ_PER_MOL = HARTREE_TO_KJ_PER_MOL / HARTREE_TO_EV
# GROMACS' c_presfac, bar per kJ/(mol nm^3), against one eV/Angstrom^3.
EV_PER_ANGSTROM3_TO_BAR = 1.602176634e6


def virial_and_pressure(workdir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Vir-XX .. Vir-ZZ (kJ/mol) and Pres-XX .. Pres-ZZ (bar) of the first frame."""
    columns = gc.energy_columns(GMX, workdir, *(f"{q}-{t}" for q in ("Vir", "Pres") for t in TENSOR))
    return tuple(np.array([columns[f"{q}-{t}"][0] for t in TENSOR]).reshape(3, 3) for q in ("Vir", "Pres"))


@needs_gmx
@pytest.mark.parametrize("box", [CUBE_NM, TRICLINIC_NM], ids=["cube", "triclinic"])
def test_gromacs_virial_is_the_qm_stress_times_volume(box, periodic_params, tmp_path):
    """Whole-system QM: GROMACS' virial is -V sigma / 2, and only that.

    md-vv from rest has no kinetic energy at step 0, so the pressure is
    -2 Xi / V, which is ASE's stress with the sign flipped.
    """
    whole_system(tmp_path, periodic_params, box, integrator="md-vv")
    frame = gc.read_trr(tmp_path / "md.trr")[0]
    atoms = Atoms(gc.SYMBOLS, positions=frame["x"] * 10.0, cell=frame["box"] * 10.0, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        sigma = atoms.get_stress(voigt=False)
    virial, pressure = virial_and_pressure(tmp_path)
    expected = 0.5 * atoms.get_volume() * sigma * EV_TO_KJ_PER_MOL
    scale = np.abs(expected).max()
    np.testing.assert_allclose(virial, expected, rtol=1e-9, atol=1e-12 * scale)
    np.testing.assert_allclose(pressure, -sigma * EV_PER_ANGSTROM3_TO_BAR, rtol=1e-8,
                               atol=1e-11 * scale)
    log = (tmp_path / "md_cp2k.out").read_text()
    assert "first stress request: available 1" in log
    assert "returned no stress tensor" not in (tmp_path / "md.log").read_text()


def strained(x_nm, box_nm, a, b, h):
    """Every position and the box under x_a -> x_a + h x_b; a <= b keeps the box lower triangular."""
    eps = np.zeros((3, 3))
    eps[a, b] = h
    deform = np.eye(3) + eps
    return np.asarray(x_nm) @ deform.T, np.asarray(box_nm) @ deform.T


@needs_gmx
@pytest.mark.parametrize("box", [QMMM_BOX_NM, MECHANICAL_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_gromacs_qmmm_virial_is_the_strain_derivative_of_the_energy(box, params, tmp_path):
    """QM/MM with PME: Xi_ba = dE/d(eps_ab) / 2 over GROMACS' own potential energy.

    The classical terms, the QM energy and the cpmdc QM-MM Ewald sum all
    enter E; a QM virial left out of GROMACS shows as a mismatch. Fixed
    PME grids and a fixed EWALD_RCUT keep both splittings fixed under the
    strain.
    """
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES, ewald_rcut=0.7,
                    ewald_precision=1e-12)
    mdp = {"integrator": "md-vv", "rcoulomb": 0.7, "rvdw": 0.7, "fourier-nx": 24, "fourier-ny": 24,
           "fourier-nz": 24, "ewald-rtol": 1e-10, "nstxout": 0, "nstvout": 0, "nstfout": 0}

    def potential(work, x_nm, box_nm):
        gc.write_system(work, positions=x_nm, box=box_nm, qm_molecules=[0], qmi=qmi, mdp=mdp,
                        qm_lj=True)
        gc.run(GMX, work)
        return gc.energies(GMX, work, "Potential")[0, 1]

    # Inside the box: mdrun breaks a water across the boundary of a box
    # strained to a negative tilt when it makes the molecules whole.
    x_nm = QMMM_NM + 0.4
    potential(tmp_path / "0", x_nm, box)
    virial, _ = virial_and_pressure(tmp_path / "0")
    h = 1e-5
    pairs = [(a, b) for a in range(3) for b in range(a, 3)]
    derivative = np.zeros((3, 3))
    for a, b in pairs:
        plus = potential(tmp_path / f"{a}{b}+", *strained(x_nm, box, a, b, h))
        minus = potential(tmp_path / f"{a}{b}-", *strained(x_nm, box, a, b, -h))
        derivative[a, b] = (plus - minus) / (2 * h)
    scale = np.abs(virial).max()
    for a, b in pairs:
        assert virial[b, a] == pytest.approx(0.5 * derivative[a, b], abs=1e-8 * scale), (a, b)
    np.testing.assert_allclose(virial, virial.T, atol=1e-9 * scale)


@needs_gmx
@pytest.mark.parametrize(
    "case, message",
    [
        ("generated", "CPMDC_GROMACS_PARAMS is not set"),
        ("unset embedding", "the MM atoms carry charges: set &CPMDC EMBEDDING"),
    ],
)
def test_gromacs_unsupported_configs_stop_mdrun(case, message, params, tmp_path):
    if case == "generated":
        done = qmmm_system(tmp_path, params, "", with_qmi=False,
                           env_extra={"CPMDC_GROMACS_PARAMS": ""},
                           **{"qmmm-cp2k-qmmethod": "PBE"})
    else:
        done = qmmm_system(tmp_path, params, qmi_input(params))
    assert done.returncode != 0
    assert "cpmdc-gromacs: " in done.stderr
    assert message in done.stderr


# --- MPI GROMACS: domain decomposition --------------------------------------

LAUNCHER, GMX_MPI, GMX_MPI_PROBLEM = (
    gc.gmx_mpi_problem() if LIB else (None, None, "libcpmdc_gromacs is not built")
)
needs_gmx_mpi = pytest.mark.skipif(GMX_MPI is None, reason=GMX_MPI_PROBLEM or "MPI GROMACS")

# 64 waters on a 0.6 nm lattice in a 2.4 nm box. Molecule 0, QM, sits at
# the origin, where the reference well is shallow, with its hydrogens across
# x = 0 and y = 0: the periodic boundary, which is a domain boundary of every
# decomposition.
DD_BOX_NM = np.eye(3) * 2.4
DD_QM = 0


def dd_positions() -> np.ndarray:
    x = gc.water_lattice(4, 0.6)
    x[3 * DD_QM : 3 * DD_QM + 3] = [
        [0.01, 0.01, 0.01173],
        [-0.0857, 0.01, -0.0117],
        [0.0339, -0.0827, -0.0117],
    ]
    return x


def dd_run(
    tmp_path, params, ranks, qm_ranks="ALL", mdrun_args=("-npme", "0"), epsilon_r=None, **kw
):
    qmi = qmi_input(
        params,
        embedding="MECHANICAL",
        qm_charges=QM_CHARGES,
        epsilon_r=epsilon_r,
        qm_ranks=qm_ranks,
    )
    mdp = {
        "rcoulomb": 0.55,
        "rvdw": 0.55,
        "rlist": 0.55,
        "nstlist": 1,
        "verlet-buffer-tolerance": -1,
        "nsteps": 20,
        "nstxout": 10,
        "nstvout": 0,
        "nstfout": 10,
        "nstenergy": 1,
        "gen-vel": "yes",
        "gen-temp": 300,
        "gen-seed": 4928,
    }
    if epsilon_r is not None:
        mdp["epsilon-r"] = epsilon_r
    gc.write_system(
        tmp_path, positions=dd_positions(), box=DD_BOX_NM, qm_molecules=[DD_QM], qmi=qmi,
        qm_lj=True,
        # A pair list rebuilt every step with no buffer: the Verlet buffer
        # lets each decomposition miss different pairs between rebuilds. The
        # cutoff is off the lattice spacing, where rounding would decide pairs.
        mdp=mdp,
    )
    return gc.run(GMX_MPI, tmp_path, check=False, launcher=LAUNCHER, ranks=ranks,
                  mdrun_args=["-dlb", "no", *mdrun_args], timeout=600, **kw)


def dd_result(path):
    done = (path / "md.log").read_text()
    frames = gc.read_trr(path / "md.trr")
    table = gc.energies(GMX_MPI, path, "Quantum-En.", "Potential", "Total-Energy")
    tensors = gc.energy_columns(GMX_MPI, path, *(f"{q}-{t}" for q in ("Vir", "Pres") for t in TENSOR))
    return done, frames, table, tensors


@pytest.fixture(scope="module")
def single_rank(params, tmp_path_factory):
    if GMX_MPI is None:
        pytest.skip(GMX_MPI_PROBLEM)
    path = tmp_path_factory.mktemp("dd-1")
    done = dd_run(path, params, 1)
    assert done.returncode == 0, done.stderr[-3000:]
    return dd_result(path)


@needs_gmx_mpi
@pytest.mark.parametrize(
    "ranks, qm_ranks", [(2, "ALL"), (4, "ALL"), (4, "1")], ids=["2-all", "4-all", "4-root"]
)
def test_mpi_domain_decomposition_matches_one_rank(ranks, qm_ranks, single_rank, params, tmp_path):
    done = dd_run(tmp_path, params, ranks, qm_ranks)
    assert done.returncode == 0, done.stderr[-3000:]
    log, frames, table, tensors = dd_result(tmp_path)
    grid = re.search(r"Domain decomposition grid (\d+) x (\d+) x (\d+)", log)
    assert grid and int(grid.group(1)) >= 2, "no decomposition along x"
    out = (tmp_path / "md_cp2k.out").read_text()
    assert f"ranks {ranks} qm_ranks {ranks if qm_ranks == 'ALL' else 1}" in out
    _, ref_frames, ref_table, ref_tensors = single_rank
    assert table.shape == ref_table.shape
    np.testing.assert_allclose(table[:, 1], ref_table[:, 1], rtol=1e-10, atol=0)
    np.testing.assert_allclose(table[:, 2:], ref_table[:, 2:], rtol=1e-10, atol=0)
    # The QM virial enters once, from the main rank, whatever the rank count.
    for name, ref in ref_tensors.items():
        np.testing.assert_allclose(tensors[name], ref, rtol=0, atol=1e-9 * np.abs(ref).max(),
                                   err_msg=name)
    assert [f["step"] for f in frames] == [f["step"] for f in ref_frames] == [0, 10, 20]
    qm = slice(3 * DD_QM, 3 * DD_QM + 3)
    for frame, ref in zip(frames, ref_frames):
        np.testing.assert_allclose(frame["x"], ref["x"], rtol=0, atol=1e-10)
        scale = np.abs(ref["f"]).max()
        np.testing.assert_allclose(frame["f"], ref["f"], rtol=0, atol=1e-10 * scale)
        np.testing.assert_allclose(frame["f"][qm], ref["f"][qm], rtol=1e-10, atol=1e-10)


@needs_gmx_mpi
def test_mpi_domain_decomposition_with_epsilon_r(params, tmp_path):
    done_1 = dd_run(tmp_path / "r1", params, 1, epsilon_r=2.5)
    assert done_1.returncode == 0, done_1.stderr[-3000:]
    _, ref_frames, ref_table, ref_tensors = dd_result(tmp_path / "r1")

    done_4 = dd_run(tmp_path / "r4", params, 4, "ALL", epsilon_r=2.5)
    assert done_4.returncode == 0, done_4.stderr[-3000:]
    log, frames, table, tensors = dd_result(tmp_path / "r4")
    grid = re.search(r"Domain decomposition grid (\d+) x (\d+) x (\d+)", log)
    assert grid and int(grid.group(1)) >= 2, "no decomposition along x"
    assert "epsilon_r 2.5" in (tmp_path / "r4" / "md_cp2k.out").read_text()
    np.testing.assert_allclose(table[:, 1], ref_table[:, 1], rtol=1e-10, atol=0)
    np.testing.assert_allclose(table[:, 2:], ref_table[:, 2:], rtol=1e-10, atol=0)
    for name, ref in ref_tensors.items():
        np.testing.assert_allclose(
            tensors[name], ref, rtol=0, atol=1e-9 * np.abs(ref).max(), err_msg=name
        )
    for frame, ref in zip(frames, ref_frames):
        np.testing.assert_allclose(frame["x"], ref["x"], rtol=0, atol=1e-10)
        np.testing.assert_allclose(
            frame["f"], ref["f"], rtol=0, atol=1e-10 * np.abs(ref["f"]).max()
        )


@needs_gmx_mpi
def test_mpi_separate_pme_ranks_are_refused(params, tmp_path):
    # libcpmdc_gromacs runs collectives on the simulation communicator, which
    # relies on every rank of it calling the QM/MM interface.
    done = dd_run(tmp_path, params, 4, mdrun_args=("-npme", "1"))
    assert done.returncode != 0
    assert "not compatible with QMMM" in " ".join((done.stderr + done.stdout).split())
