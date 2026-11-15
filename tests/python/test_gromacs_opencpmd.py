"""The stress of libcpmdc_gromacs with OpenCPMD as the QM engine.

These cases need ``libcpmdc_gromacs`` linked to an OpenCPMD build of cpmdc
(``CPMDC_GROMACS_LIB``), ``capnp`` to encode the decks with a tight SCF, and
``CPMDC_PSEUDO_DIR``; they skip otherwise. The reference is ASE's central
difference of the library's own energy under strains of the cell, as in the
reference-evaluator cases of ``test_gromacs.py``.
"""

from __future__ import annotations

import ctypes
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from conftest import ROOT
from cpmdc import BOHR_TO_ANGSTROM
from cpmdc.ase import CPMDCCalculator
from cpmdc.gromacs import BOHR_TO_NM, qmi_input
from test_gromacs import (
    LIB,
    MECHANICAL_TRICLINIC_NM,
    QM_CHARGES,
    QMMM_BOX_NM,
    QMMM_NM,
    Shim,
    ShimCalculator,
    ase_stress,
    gc,
    numerical_stress,
    setup_env,
)


# The libcpmdc that libcpmdc_gromacs links, which Meson builds and installs
# beside it; CPMDCCalculator would otherwise search CPMDC_LIBRARY and build/.
CPMDC_LIB = Path(LIB).parent / "libcpmdc.so" if LIB else None


def opencpmd_problem() -> str | None:
    if LIB is None:
        return "libcpmdc_gromacs is not built"
    if not ctypes.CDLL(str(LIB)).cpmdc_available():
        return f"{LIB} is not linked to OpenCPMD"
    if not CPMDC_LIB.is_file():
        return f"no libcpmdc beside {LIB}"
    if shutil.which("capnp") is None:
        return "capnp is not on PATH"
    if not os.environ.get("CPMDC_PSEUDO_DIR"):
        return "CPMDC_PSEUDO_DIR is not set"
    return None


PROBLEM = opencpmd_problem()
pytestmark = pytest.mark.skipif(PROBLEM is not None, reason=PROBLEM or "")

CONVERGENCE = 1e-8
# Each cell is a new setup at the deck's cutoff, so a plane wave that the
# strain takes across the cutoff changes the basis, and the energy, in a
# step; OpenCPMD's stress is the derivative at a fixed basis. A strain of
# 1e-6 crosses none in these cells (1e-4 does in the triclinic one). The
# cluster cell does not change; there the error of the difference goes as
# the strain squared and is up to 1.2e-5 of the stress at 1e-5.
PERIODIC_STRAIN = 1e-6
CLUSTER_STRAIN = 1e-5
WATER_NM = np.array([[0.50, 0.50, 0.51173], [0.50, 0.57572, 0.45308], [0.50, 0.42428, 0.45308]])
CUBE_NM = np.eye(3)
TRICLINIC_NM = np.array([[1.0, 0.0, 0.0], [0.13, 0.94, 0.0], [-0.075, 0.1, 1.04]])


def tight(source: str, tmp: Path) -> Path:
    text = (ROOT / "tests" / source).read_text()
    text = re.sub(r"convergenceOrbitals = [0-9.eE+-]+", f"convergenceOrbitals = {CONVERGENCE:g}", text)
    text = re.sub(r"maxStep = \d+", "maxStep = 300", text)
    text = re.sub(r"maxIter = \d+", "maxIter = 300", text)
    target = tmp / (Path(source).stem + ".bin")
    target.write_bytes(subprocess.run(
        ["capnp", "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode(), stdout=subprocess.PIPE, check=True,
    ).stdout)
    return target


@pytest.fixture(scope="module")
def cluster_params(tmp_path_factory):
    return tight("water_blyp_params.capnp.txt", tmp_path_factory.mktemp("cluster"))


@pytest.fixture(scope="module")
def periodic_params(tmp_path_factory):
    return tight("water_blyp_periodic_params.capnp.txt", tmp_path_factory.mktemp("periodic"))


def shim_and_strain_stress(inp: Path, out: Path, symbols: str, x_nm, box_nm, strain):
    """cp2k_get_stress_tensor and ASE's numerical stress of the same library, eV/Angstrom^3."""
    shim = Shim(inp, out)
    x = np.asarray(x_nm) / BOHR_TO_NM
    rows = np.asarray(box_nm) / BOHR_TO_NM
    shim.step(x, rows)
    stress, available = shim.stress()
    assert available == 1
    atoms = Atoms(symbols, positions=x * BOHR_TO_ANGSTROM, cell=rows * BOHR_TO_ANGSTROM, pbc=True)
    atoms.calc = ShimCalculator(shim)
    numeric = numerical_stress(atoms, eps=strain)
    shim.close()
    return ase_stress(stress), numeric


@pytest.mark.parametrize("box", [CUBE_NM, TRICLINIC_NM], ids=["cube", "triclinic"])
def test_opencpmd_whole_system_stress_matches_ase(box, periodic_params, tmp_path):
    """CELL PERIODIC: OpenCPMD's stress tensor, through the drop-in."""
    inp = setup_env(tmp_path, qmi_input(periodic_params), "OHH", [True] * 3, [0.0] * 3, WATER_NM)
    ours, numeric = shim_and_strain_stress(inp, tmp_path / "sys.out", "OHH", WATER_NM, box,
                                           PERIODIC_STRAIN)
    atoms = Atoms("OHH", positions=WATER_NM * 10.0, cell=box * 10.0, pbc=True)
    with CPMDCCalculator(periodic_params, library=CPMDC_LIB) as calc:
        atoms.calc = calc
        analytic = atoms.get_stress(voigt=False)
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ours, analytic, rtol=0, atol=1e-7 * scale)
    # At a strain of 1e-6 the difference resolves 1e-12 Hartree of energy,
    # which is the SCF's reproducibility: the cube is 1.5e-6 of the stress off.
    np.testing.assert_allclose(ours, numeric, rtol=0, atol=3e-6 * scale)


@pytest.mark.parametrize("box", [QMMM_BOX_NM, MECHANICAL_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_opencpmd_qmmm_stress_matches_ase(box, cluster_params, tmp_path):
    """CELL NONE, mechanical embedding: sum x (x) f of the Hockney cluster and the Ewald term.

    The QM water sits at the centre of the deck's 10 Angstrom cluster box.
    """
    x_nm = QMMM_NM + 0.5
    qmi = qmi_input(cluster_params, embedding="MECHANICAL", qm_charges=QM_CHARGES, ewald_rcut=0.7,
                    ewald_precision=1e-12)
    inp = setup_env(tmp_path, qmi, "OHH" * 4, [True] * 3 + [False] * 9,
                    [0.0] * 3 + [gc.TIP3P[s] for s in "OHH" * 3], x_nm)
    ours, numeric = shim_and_strain_stress(inp, tmp_path / "sys.out", "OHH" * 4, x_nm, box,
                                           CLUSTER_STRAIN)
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ours, numeric, rtol=0, atol=3e-5 * scale)
    assert "backend opencpmd" in (tmp_path / "sys.out").read_text()
