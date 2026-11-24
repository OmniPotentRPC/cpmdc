"""Link atoms: QM-MM covalent bonds capped with a hydrogen in the QM calculation.

The library cases load ``libcpmdc_gromacs`` as ``test_gromacs.py`` does and
check the capped energy, the chain-rule distribution of the link force, the
nearest image of the MM atom, the QM-MM Coulomb sum without the charge of the
MM atom of each cut bond, and the strain derivative. The GROMACS cases need a
double-precision ``gmx_d`` with the CP2K interface: grompp writes the
``&LINK`` sections of a bonded QM/MM group itself, and mdrun's forces must be
the library's.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from cpmdc import BOHR_TO_ANGSTROM
from cpmdc.gromacs import (
    BOHR_TO_NM,
    HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM,
    HARTREE_TO_KJ_PER_MOL,
    qmi_input,
)
from test_gromacs import (
    GMX,
    GMX_MPI,
    GMX_MPI_PROBLEM,
    LAUNCHER,
    LIB,
    QM_CHARGES,
    SCRIPT,
    Shim,
    ShimCalculator,
    ase_au,
    ase_stress,
    generated_input,
    gc,
    needs_gmx,
    needs_gmx_mpi,
    needs_lib,
    numerical_stress,
    reference_precision,
    setup_env,
    strained,
    virial_and_pressure,
)

# A QM water about the origin whose oxygen (atom 1) is bonded to the oxygen
# of an MM water (atom 4), nm. The bond is 0.14 nm and crosses nothing else.
LINK_NM = np.array(
    [
        [0.0, 0.0, 0.01173],
        [0.0, 0.07572, -0.04692],
        [0.0, -0.07572, -0.04692],
        [0.12, 0.04, 0.05],
        [0.17, 0.11, 0.04],
        [0.15, -0.03, 0.10],
        [0.62, 0.55, 0.60],
        [0.68, 0.62, 0.61],
        [0.55, 0.59, 0.55],
    ]
)
LINK_SYMBOLS = "OHHOHHOHH"
LINK_BOX_NM = np.diag([1.4, 1.4, 1.4])
LINK_TRICLINIC_NM = LINK_BOX_NM @ np.array([[1, 0, 0], [0.1, 1, 0], [-0.08, 0.12, 1]])
QM_FLAGS = [True] * 3 + [False] * 6
MM_CHARGES = [gc.TIP3P[s] for s in "OHH" * 2]
ALPHA_DEFAULT = 1.38


def hydrogen(x_qm, x_mm, alpha=ALPHA_DEFAULT):
    return x_qm + (x_mm - x_qm) / alpha


def run_link(tmp_path, params, qmi, x_nm, box_nm, charges=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    charges = [0.0] * 9 if charges is None else charges
    inp = setup_env(tmp_path, qmi, LINK_SYMBOLS, QM_FLAGS, charges, x_nm)
    shim = Shim(inp, tmp_path / "sys.out")
    out = shim.step(np.asarray(x_nm) / BOHR_TO_NM, np.asarray(box_nm) / BOHR_TO_NM)
    return shim, out


def capped_reference(params, x_bohr, alpha=ALPHA_DEFAULT, box_rows=None):
    """The QM water plus its hydrogen, and the forces handed to the two bonded atoms."""
    x_qm, x_mm = x_bohr[0], x_bohr[3]
    if box_rows is not None:
        shift = np.round((x_mm - x_qm) @ np.linalg.inv(box_rows))
        x_mm = x_mm - shift @ box_rows
    cluster = np.vstack([x_bohr[:3], hydrogen(x_qm, x_mm, alpha)])
    energy, f = ase_au(params, "OHHH", cluster)
    a = 1.0 / alpha
    forces = np.zeros((len(x_bohr), 3))
    forces[:3] = f[:3]
    forces[0] += (1.0 - a) * f[3]
    forces[3] += a * f[3]
    return energy, forces


@needs_lib
@pytest.mark.parametrize("box", [LINK_BOX_NM, LINK_TRICLINIC_NM, np.zeros((3, 3))],
                         ids=["cube", "triclinic", "no-box"])
def test_link_hydrogen_is_part_of_the_qm_calculation(box, params, tmp_path):
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)])
    shim, (energy, forces) = run_link(tmp_path, params, qmi, LINK_NM, box)
    shim.close()
    x = LINK_NM / BOHR_TO_NM
    ref_e, ref_f = capped_reference(params, x, box_rows=box / BOHR_TO_NM if box.any() else None)
    assert energy == pytest.approx(ref_e, rel=1e-12, abs=1e-14)
    np.testing.assert_allclose(forces, ref_f, rtol=1e-11, atol=1e-13)
    np.testing.assert_array_equal(forces[6:], 0.0)
    assert "links 1" in (tmp_path / "sys.out").read_text()
    # A bare water would not feel the MM oxygen: the cap changes the energy.
    bare_e, _ = ase_au(params, "OHH", x[:3])
    assert abs(energy - bare_e) > 1e-6


@needs_lib
def test_link_force_distribution_sums_to_the_hydrogen_force(params, tmp_path):
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)])
    shim, (_, forces) = run_link(tmp_path, params, qmi, LINK_NM, LINK_BOX_NM)
    shim.close()
    x = LINK_NM / BOHR_TO_NM
    cluster = np.vstack([x[:3], hydrogen(x[0], x[3])])
    _, f = ase_au(params, "OHHH", cluster)
    np.testing.assert_allclose(forces[:4].sum(axis=0), f.sum(axis=0), rtol=0, atol=1e-13)
    a = 1.0 / ALPHA_DEFAULT
    np.testing.assert_allclose(forces[3], a * f[3], rtol=1e-11, atol=1e-13)


@needs_lib
def test_link_forces_are_the_gradient_of_the_energy(params, tmp_path):
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)])
    shim, (_, forces) = run_link(tmp_path / "0", params, qmi, LINK_NM, LINK_BOX_NM)
    shim.close()
    h = 1e-4
    for atom, axis in [(0, 0), (0, 2), (1, 1), (3, 0), (3, 1), (4, 2)]:
        energies = []
        for sign in (1, -1):
            x = LINK_NM.copy()
            x[atom, axis] += sign * h * BOHR_TO_NM
            shim, (e, _) = run_link(tmp_path / f"{atom}{axis}{sign}", params, qmi, x, LINK_BOX_NM)
            shim.close()
            energies.append(e)
        assert -(energies[0] - energies[1]) / (2 * h) == pytest.approx(forces[atom, axis], abs=1e-8)


@needs_lib
@pytest.mark.parametrize("box", [LINK_BOX_NM, LINK_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_link_uses_the_nearest_image_of_the_mm_atom(box, params, tmp_path):
    """GROMACS wraps the MM atom into the box; the hydrogen still sits on the bond."""
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)])
    shim, (e0, f0) = run_link(tmp_path / "a", params, qmi, LINK_NM, box)
    shim.close()
    for lattice in (box[0], -box[1], box[2] + box[0]):
        wrapped = LINK_NM.copy()
        wrapped[3] += lattice
        shim, (e1, f1) = run_link(tmp_path / f"b{lattice[0]:.3f}{lattice[1]:.3f}", params, qmi, wrapped, box)
        shim.close()
        assert e1 == pytest.approx(e0, rel=1e-13)
        np.testing.assert_allclose(f1, f0, rtol=1e-11, atol=1e-13)


@needs_lib
def test_link_alpha_moves_the_hydrogen(params, tmp_path):
    for alpha in (1.2, 1.7):
        qmi = qmi_input(params, embedding="NONE", links=[(1, 4)], link_alpha=alpha)
        shim, (energy, forces) = run_link(tmp_path / str(alpha), params, qmi, LINK_NM, LINK_BOX_NM)
        shim.close()
        ref_e, ref_f = capped_reference(params, LINK_NM / BOHR_TO_NM, alpha)
        assert energy == pytest.approx(ref_e, rel=1e-12)
        np.testing.assert_allclose(forces, ref_f, rtol=1e-11, atol=1e-13)


@needs_lib
def test_two_links_on_one_qm_atom_and_on_one_mm_atom(params, tmp_path):
    """Hydrogens are added bond by bond; forces accumulate on shared atoms."""
    links = [(1, 4), (1, 7), (2, 4)]
    qmi = qmi_input(params, embedding="NONE", links=links)
    x_nm = LINK_NM.copy()
    x_nm[6] = [-0.10, 0.05, 0.11]
    x_nm[7] = [-0.15, 0.10, 0.12]
    x_nm[8] = [-0.12, 0.00, 0.17]
    shim, (energy, forces) = run_link(tmp_path, params, qmi, x_nm, LINK_BOX_NM)
    shim.close()
    x = x_nm / BOHR_TO_NM
    caps = [hydrogen(x[q - 1], x[m - 1]) for q, m in links]
    ref_e, f = ase_au(params, "OHH" + "H" * 3, np.vstack([x[:3], *caps]))
    a = 1.0 / ALPHA_DEFAULT
    expected = np.zeros((9, 3))
    expected[:3] = f[:3]
    for (q, m), fh in zip(links, f[3:]):
        expected[q - 1] += (1.0 - a) * fh
        expected[m - 1] += a * fh
    assert energy == pytest.approx(ref_e, rel=1e-12)
    np.testing.assert_allclose(forces, expected, rtol=1e-11, atol=1e-13)


@needs_lib
@pytest.mark.parametrize("box", [LINK_BOX_NM, LINK_TRICLINIC_NM, np.zeros((3, 3))],
                         ids=["cube", "triclinic", "no-box"])
def test_mechanical_embedding_drops_the_charge_of_the_cut_bond(box, params, tmp_path):
    qmi = qmi_input(params, embedding="MECHANICAL", qm_charges=QM_CHARGES, links=[(1, 4)])
    # MM charges of the PDB: the three atoms of the bonded water, then a free one.
    charges = [0.0] * 3 + MM_CHARGES
    shim, (energy, forces) = run_link(tmp_path, params, qmi, LINK_NM, box, charges)
    shim.close()
    x = LINK_NM / BOHR_TO_NM
    mm_q = np.array(charges[3:])
    mm_q[0] = 0.0
    cell = box / BOHR_TO_NM
    if box.any():
        coul, fq, fm = gc.ewald_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q, cell,
                                      gc.default_rcut(cell, BOHR_TO_NM), 1e-8)
    else:
        coul, fq, fm = gc.direct_qm_mm(x[:3], QM_CHARGES, x[3:], mm_q)
    ref_e, ref_f = capped_reference(params, x, box_rows=cell if box.any() else None)
    assert energy == pytest.approx(ref_e + coul, rel=1e-11, abs=1e-12)
    ref_f[:3] += fq
    ref_f[3:] += fm
    np.testing.assert_allclose(forces, ref_f, rtol=1e-10, atol=1e-12)
    # The deleted charge is the only difference from keeping every MM charge.
    assert abs(coul) > 1e-6


@needs_lib
@reference_precision
@pytest.mark.parametrize("embedding", ["MECHANICAL", "NONE"])
@pytest.mark.parametrize("box", [LINK_BOX_NM, LINK_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_link_stress_matches_ase_numerical_stress(box, embedding, params, tmp_path):
    """CELL NONE: the link hydrogen strains with its two bonded atoms."""
    options = {"embedding": embedding, "ewald_rcut": 0.7, "ewald_precision": 1e-12, "links": [(1, 4)]}
    charges = [0.0] * 9
    if embedding == "MECHANICAL":
        options["qm_charges"] = QM_CHARGES
        charges = [0.0] * 3 + MM_CHARGES
    shim, (_, forces) = run_link(tmp_path, params, qmi_input(params, **options), LINK_NM, box, charges)
    x = LINK_NM / BOHR_TO_NM
    rows = box / BOHR_TO_NM
    stress, available = shim.stress()
    assert available == 1
    atoms = Atoms(LINK_SYMBOLS, positions=x * BOHR_TO_ANGSTROM, cell=rows * BOHR_TO_ANGSTROM, pbc=True)
    atoms.calc = ShimCalculator(shim)
    numeric = numerical_stress(atoms)
    shim.close()
    scale = np.abs(numeric).max()
    np.testing.assert_allclose(ase_stress(stress), numeric, rtol=0, atol=1e-8 * scale)


@needs_lib
def test_generated_input_with_link_sections(params, tmp_path, monkeypatch):
    """GROMACS' generated input names each bond in &QMMM, with no alpha."""
    monkeypatch.setenv("CPMDC_GROMACS_PARAMS", str(params))
    text = generated_input().replace(
        "  &QMMM\n  &END QMMM\n",
        "  &QMMM\n    &QM_KIND O\n      MM_INDEX 1\n    &END QM_KIND\n"
        "    &LINK\n      QM_INDEX 1\n      MM_INDEX 4\n    &END LINK\n  &END QMMM\n",
    )
    inp = setup_env(tmp_path, text, LINK_SYMBOLS, QM_FLAGS, [0.0] * 9, LINK_NM)
    shim = Shim(inp, tmp_path / "sys.out")
    energy, forces = shim.step(LINK_NM / BOHR_TO_NM, np.zeros((3, 3)))
    shim.close()
    ref_e, ref_f = capped_reference(params, LINK_NM / BOHR_TO_NM)
    assert energy == pytest.approx(ref_e, rel=1e-12)
    np.testing.assert_allclose(forces, ref_f, rtol=1e-11, atol=1e-13)
    assert "links 1" in (tmp_path / "sys.out").read_text()


@needs_lib
@pytest.mark.parametrize(
    "link, message",
    [
        ("      MM_INDEX 4\n", "&LINK section without QM_INDEX"),
        ("      QM_INDEX 1\n", "&LINK section without MM_INDEX"),
        ("      QM_INDEX 4\n      MM_INDEX 1\n", "&LINK QM_INDEX 4 is an MM atom"),
        ("      QM_INDEX 1\n      MM_INDEX 2\n", "&LINK MM_INDEX 2 is a QM atom"),
        ("      QM_INDEX 1\n      MM_INDEX 10\n", "&LINK names atom 10, but the PDB file has 9 atoms"),
        ("      QM_INDEX 0\n      MM_INDEX 4\n", "is a 1-based atom index"),
        ("      QM_INDEX 1\n      MM_INDEX 4\n      ALPHA_IMOMM 0.9\n", "ALPHA_IMOMM must exceed 1"),
        ("      QM_INDEX 1\n      MM_INDEX 4\n      LINK_TYPE PSEUDO\n", "LINK_TYPE is IMOMM"),
        ("      QM_INDEX 1\n      MM_INDEX 4\n      RADIUS 1.0\n", "unsupported &LINK keyword RADIUS"),
        ("      QM_INDEX 1\n      MM_INDEX 4\n    &END LINK\n    &LINK\n"
         "      QM_INDEX 1\n      MM_INDEX 4\n", "appears twice"),
    ],
)
def test_bad_link_sections_are_refused(link, message, params, tmp_path):
    qmi = qmi_input(params, embedding="NONE").replace(
        "  &DFT", f"  &QMMM\n    &LINK\n{link}    &END LINK\n  &END QMMM\n  &DFT", 1
    )
    stderr = shim_fails_for_links(tmp_path, qmi)
    assert message in stderr
    assert message in (tmp_path / "sys.out").read_text()


def shim_fails_for_links(tmp_path, qmi):
    inp = setup_env(tmp_path, qmi, LINK_SYMBOLS, QM_FLAGS, [0.0] * 9, LINK_NM)
    done = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(LIB), str(inp), str(tmp_path / "sys.out"), "27"],
        capture_output=True, text=True, timeout=120, env=dict(os.environ),
    )
    assert done.returncode == 1, done.stderr
    assert "cpmdc-gromacs: " in done.stderr
    return done.stderr


# --- GROMACS ----------------------------------------------------------------

BONDED_TOPOLOGY = """
[ moleculetype ]
CAP 2

[ atoms ]
1 OQ 1 CAP O1 1 0.0 15.9994
2 HQ 1 CAP H2 1 0.0 1.008
3 HQ 1 CAP H3 1 0.0 1.008
4 OQ 1 CAP O4 1 0.0 15.9994
5 HQ 1 CAP H5 1 0.0 1.008
6 HQ 1 CAP H6 1 0.0 1.008

[ bonds ]
1 2 1 0.09572 502416.0
1 3 1 0.09572 502416.0
1 4 1 0.14 0.0
4 5 1 0.09572 0.0
4 6 1 0.09572 0.0

[ angles ]
2 1 3 1 104.52 628.02
"""


def write_bonded_system(
    workdir: Path, qmi: str, x_nm=LINK_NM[:6], box_nm=LINK_BOX_NM, waters=None, **mdp
) -> None:
    """One 6-atom molecule whose first water is QM, bonded to the second, MM.

    ``waters`` adds that many TIP3P waters (three rows of positions each).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    rows = list(np.asarray(x_nm)) + ([] if waters is None else list(np.asarray(waters)))
    lines = ["bonded", str(len(rows))]
    names = ["O1", "H2", "H3", "O4", "H5", "H6"]
    for i, row in enumerate(rows):
        if i < 6:
            residue, name = "CAP", names[i]
        else:
            residue, name = "SOL", ("OW", "HW1", "HW2")[(i - 6) % 3]
        molecule = 1 if i < 6 else 2 + (i - 6) // 3
        lines.append(
            f"{molecule:5d}{residue:<5s}{name:>5s}{i + 1:5d}" + "".join(f"{v:15.10f}" for v in row)
        )
    lines.append(gc.gro_box(np.asarray(box_nm)))
    (workdir / "conf.gro").write_text("\n".join(lines) + "\n")
    molecules = "CAP 1\n" + "SOL 1\n" * ((len(rows) - 6) // 3)
    (workdir / "topol.top").write_text(
        gc.TOPOLOGY_HEAD
        + BONDED_TOPOLOGY
        + gc.MOLECULE.format(
            name="SOL", o="OW", h="HW", qo=-0.834, qh1=0.417, qh2=0.417
        )
        + f"\n[ system ]\nbonded\n\n[ molecules ]\n{molecules}"
    )
    (workdir / "index.ndx").write_text(
        "[ System ]\n" + " ".join(map(str, range(1, len(rows) + 1))) + "\n[ QM ]\n1 2 3\n"
    )
    settings = {
        "integrator": "md", "nsteps": 0, "dt": 0.00025, "cutoff-scheme": "Verlet",
        "coulombtype": "PME", "rcoulomb": 0.5, "rvdw": 0.5, "nstxout": 1, "nstvout": 1,
        "nstfout": 1, "nstenergy": 1, "nstcalcenergy": 1, "tcoupl": "no", "pcoupl": "no",
        "constraints": "none", "comm-mode": "None", "qmmm-cp2k-active": "true",
        "qmmm-cp2k-qmgroup": "QM", "qmmm-cp2k-qmmethod": "INPUT",
    }
    settings.update(mdp)
    (workdir / "md.mdp").write_text("".join(f"{k} = {v}\n" for k, v in settings.items()))
    (workdir / "cpmdc.inp").write_text(qmi)


def mdrun_forces(workdir: Path) -> np.ndarray:
    frames = gc.read_trr(workdir / "md.trr")
    return frames[0]["f"] / HARTREE_PER_BOHR_TO_KJ_PER_MOL_NM


@needs_gmx
@pytest.mark.parametrize("route", ["explicit", "generated"])
def test_gromacs_caps_a_bond_across_the_qm_boundary(route, params, tmp_path):
    """mdrun's forces are the library's; the bond between the groups has no MM force."""
    x_nm = LINK_NM[:6].copy()
    env = {}
    if route == "explicit":
        write_bonded_system(tmp_path, qmi_input(params, embedding="NONE", links=[(1, 4)]))
        with_qmi = True
    else:
        write_bonded_system(tmp_path, "", **{"qmmm-cp2k-qmmethod": "BLYP"})
        env["CPMDC_GROMACS_PARAMS"] = str(params)
        with_qmi = False
        # Its own input places the QM atoms' centre of geometry in the middle
        # of the 1 nm QM box, and GROMACS sends every atom shifted with it.
        x_nm += 0.5 - x_nm[:3].mean(axis=0)
    gc.run(GMX, tmp_path, with_qmi=with_qmi, env_extra=env)
    ref_e, ref_f = capped_reference(params, x_nm / BOHR_TO_NM)
    np.testing.assert_allclose(mdrun_forces(tmp_path), ref_f, rtol=1e-9, atol=1e-10)
    quantum = gc.energies(GMX, tmp_path, "Quantum-En.")[0, 1]
    assert quantum == pytest.approx(ref_e * HARTREE_TO_KJ_PER_MOL, rel=1e-10)
    text = (tmp_path / "md_cp2k.out").read_text()
    assert "links 1" in text
    if route == "generated":
        assert "&LINK" in (tmp_path / "md_cp2k.inp").read_text()


@needs_gmx
@pytest.mark.parametrize("box", [LINK_BOX_NM, LINK_TRICLINIC_NM], ids=["cube", "triclinic"])
def test_gromacs_virial_includes_the_link_hydrogen(box, params, tmp_path):
    """Xi_ba = dE/d(eps_ab) / 2 for a QM/MM system with a capped bond.

    Every classical term of the bonded molecule is zero, so GROMACS' potential
    is the QM energy of the capped water, and the strain moves the hydrogen
    with the bonded atoms. This needs the virial patch.
    """
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)])
    mdp = {"integrator": "md-vv", "nstxout": 0, "nstvout": 0, "nstfout": 0}

    def potential(work, x_nm, box_nm):
        write_bonded_system(work, qmi, x_nm, box_nm, **mdp)
        gc.run(GMX, work)
        return gc.energies(GMX, work, "Potential")[0, 1]

    x_nm = LINK_NM[:6] + 0.4
    potential(tmp_path / "0", x_nm, box)
    virial, _ = virial_and_pressure(tmp_path / "0")
    h = 1e-5
    pairs = [(a, b) for a in range(3) for b in range(a, 3)]
    for a, b in pairs:
        plus = potential(tmp_path / f"{a}{b}+", *strained(x_nm, box, a, b, h))
        minus = potential(tmp_path / f"{a}{b}-", *strained(x_nm, box, a, b, -h))
        scale = np.abs(virial).max()
        assert virial[b, a] == pytest.approx(0.5 * (plus - minus) / (2 * h), abs=1e-8 * scale), (a, b)
    assert np.abs(virial).max() > 0.0


# --- MPI GROMACS --------------------------------------------------------------

DD_BOX_NM = np.eye(3) * 2.4


def dd_waters() -> np.ndarray:
    """A 4 x 4 x 4 lattice of waters, less those within 0.35 nm of the bonded molecule."""
    x = gc.water_lattice(4, 0.6)
    keep = [
        m for m in range(len(x) // 3)
        if np.min(np.linalg.norm(x[3 * m : 3 * m + 3, None] - LINK_NM[None, :6], axis=2)) > 0.35
    ]
    return np.vstack([x[3 * m : 3 * m + 3] for m in keep])


def dd_link_run(tmp_path, params, ranks, qm_ranks="ALL"):
    qmi = qmi_input(params, embedding="NONE", links=[(1, 4)], qm_ranks=qm_ranks)
    # QM water and its bonded MM water sit across the periodic boundary at the
    # origin, so every decomposition splits the capped bond between domains.
    write_bonded_system(
        tmp_path, qmi, LINK_NM[:6], DD_BOX_NM, dd_waters(), **{
            "rcoulomb": 0.55, "rvdw": 0.55, "rlist": 0.55, "nstlist": 1,
            "verlet-buffer-tolerance": -1, "nsteps": 20, "nstxout": 10, "nstvout": 0,
            "nstfout": 10, "nstenergy": 1, "gen-vel": "yes", "gen-temp": 300, "gen-seed": 4928,
        },
    )
    return gc.run(GMX_MPI, tmp_path, check=False, launcher=LAUNCHER, ranks=ranks,
                  mdrun_args=["-dlb", "no", "-npme", "0"], timeout=600)


@pytest.fixture(scope="module")
def single_rank_link(params, tmp_path_factory):
    if GMX_MPI is None:
        pytest.skip(GMX_MPI_PROBLEM)
    path = tmp_path_factory.mktemp("dd-link-1")
    done = dd_link_run(path, params, 1)
    assert done.returncode == 0, done.stderr[-3000:]
    return gc.read_trr(path / "md.trr"), gc.energies(GMX_MPI, path, "Quantum-En.", "Potential")


@needs_gmx_mpi
@pytest.mark.parametrize(
    "ranks, qm_ranks", [(2, "ALL"), (4, "ALL"), (4, "1")], ids=["2-all", "4-all", "4-root"]
)
def test_mpi_link_atoms_match_one_rank(ranks, qm_ranks, single_rank_link, params, tmp_path):
    """The link hydrogen is built from rank 0's positions on every QM rank."""
    done = dd_link_run(tmp_path, params, ranks, qm_ranks)
    assert done.returncode == 0, done.stderr[-3000:]
    grid = re.search(r"Domain decomposition grid (\d+) x (\d+) x (\d+)",
                     (tmp_path / "md.log").read_text())
    assert grid and int(grid.group(1)) >= 2, "no decomposition along x"
    assert (
        f"ranks {ranks} qm_ranks {ranks if qm_ranks == 'ALL' else 1} "
        "epsilon_r 1 links 1"
    ) in (tmp_path / "md_cp2k.out").read_text()
    ref_frames, ref_table = single_rank_link
    table = gc.energies(GMX_MPI, tmp_path, "Quantum-En.", "Potential")
    np.testing.assert_allclose(table[:, 1:], ref_table[:, 1:], rtol=1e-10, atol=0)
    frames = gc.read_trr(tmp_path / "md.trr")
    assert [f["step"] for f in frames] == [f["step"] for f in ref_frames] == [0, 10, 20]
    for frame, ref in zip(frames, ref_frames):
        np.testing.assert_allclose(frame["x"], ref["x"], rtol=0, atol=1e-10)
        np.testing.assert_allclose(frame["f"], ref["f"], rtol=0, atol=1e-10 * np.abs(ref["f"]).max())
        np.testing.assert_allclose(frame["f"][:6], ref["f"][:6], rtol=1e-10, atol=1e-10)
