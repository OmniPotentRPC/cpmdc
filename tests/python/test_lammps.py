"""cpmdc inside LAMMPS: cpmdc.lammps.FixExternal and pair_style cpmdc.

LAMMPS comes from its Python module (``pip install lammps``); every test
skips with the reason when that module cannot start LAMMPS. The plugin is
``CPMDC_LAMMPS_PLUGIN`` or ``build/cpmdcplugin.so`` of this checkout, and
the plugin cases skip when neither exists. Energies, forces, and virials
are compared with ``CPMDCCalculator``, which reaches the library through
``ForceInput`` instead of ``cpmdc_session_calculate_au()``. The MPI cases
run ``lammps_case.py`` under ``mpiexec``.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from conftest import ROOT
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, CpmdcError
from cpmdc.ase import CPMDCCalculator
from cpmdc.lammps import UNIT_STYLES, FixExternal, box_rows, type_numbers, unit_factors

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lammps_case import CMDARGS, THERMO, attach, build, measure


def lammps_problem() -> str:
    try:
        import lammps
    except Exception as exc:  # ImportError, or a missing MPI library
        return f"the lammps Python module is not importable: {exc}"
    try:
        lmp = lammps.lammps(cmdargs=CMDARGS)
    except Exception as exc:
        return f"the lammps Python module cannot start LAMMPS: {exc}"
    lmp.close()
    return ""


LAMMPS_PROBLEM = lammps_problem()
pytestmark = pytest.mark.skipif(bool(LAMMPS_PROBLEM), reason=LAMMPS_PROBLEM or "LAMMPS")


def plugin_path() -> Path | None:
    named = os.environ.get("CPMDC_LAMMPS_PLUGIN")
    candidate = Path(named) if named else ROOT / "build" / "cpmdcplugin.so"
    return candidate if candidate.is_file() else None


PLUGIN = plugin_path()
ROUTES = [
    "fix",
    pytest.param(
        "plugin",
        marks=pytest.mark.skipif(PLUGIN is None, reason="cpmdcplugin.so is not built"),
    ),
]

# Two waters, Angstrom. Every atom is inside both cells below.
WATERS = [
    [1.0, 1.0, 1.1173],
    [1.0, 1.7572, 0.5308],
    [1.0, 0.2428, 0.5308],
    [5.0, 4.0, 4.1],
    [5.6, 4.7, 4.2],
    [4.3, 4.4, 3.6],
]
SYMBOLS = "OHHOHH"
TYPES = [1, 2, 2, 1, 2, 2]
MASSES = [15.999, 1.008]
CUBE = [[9.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 9.0]]
TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def per_angstrom(units: str) -> float:
    """LAMMPS length units per Angstrom."""
    return unit_factors(units)[0] / BOHR_TO_ANGSTROM


def per_ev(units: str) -> float:
    """LAMMPS energy units per eV."""
    return unit_factors(units)[1] / HARTREE_TO_EV


def spec(
    route: str, units: str, params: Path, cell=TRICLINIC, positions=WATERS, **extra
):
    scale = per_angstrom(units)
    out = {
        "route": route,
        "units": units,
        "cell": None if cell is None else (np.asarray(cell) * scale).tolist(),
        "positions": (np.asarray(positions) * scale).tolist(),
        "types": TYPES,
        "masses": MASSES,
        "elements": ["O", "H"],
        "params": str(params),
        "plugin": str(PLUGIN) if PLUGIN else "",
    }
    out.update(extra)
    return out


def run_inline(case: dict, *commands: str) -> tuple[dict, object]:
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string("\n".join(THERMO))
    for command in commands or ("run 0",):
        if driver is not None and command.startswith("run "):
            driver.run(*command.split()[1:])
        else:
            lmp.command(command)
    return measure(lmp), (lmp, driver)


def finish(handles) -> None:
    lmp, driver = handles
    if driver is not None:
        driver.close()
    lmp.close()


def ase_reference(params: Path, positions, cell) -> dict:
    """eV, eV/Angstrom, and the virial -dE/d(strain) in eV, xx yy zz xy xz yz."""
    atoms = Atoms(SYMBOLS, positions=positions)
    if cell is not None:
        atoms.cell = cell
        atoms.pbc = True
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        out = {"pe": atoms.get_potential_energy(), "forces": atoms.get_forces()}
        if cell is not None:
            w = -atoms.get_stress(voigt=False) * atoms.get_volume()
            out["virial"] = np.array(
                [w[0, 0], w[1, 1], w[2, 2], w[0, 1], w[0, 2], w[1, 2]]
            )
        else:
            out["virial"] = np.zeros(6)
    return out


def assert_matches(got: dict, ref: dict, units: str) -> None:
    energy = per_ev(units)
    force = energy / per_angstrom(units)
    assert got["pe"] == pytest.approx(ref["pe"] * energy, rel=1e-11, abs=1e-12)
    np.testing.assert_allclose(
        got["forces"], ref["forces"] * force, rtol=1e-10, atol=1e-12
    )
    np.testing.assert_allclose(
        got["virial"], ref["virial"] * energy, rtol=1e-10, atol=1e-10
    )


# --- conversions --------------------------------------------------------


def test_unit_table():
    assert set(UNIT_STYLES) == {"metal", "real", "electron"}
    assert unit_factors("electron") == (1.0, 1.0)
    assert unit_factors("metal") == (BOHR_TO_ANGSTROM, HARTREE_TO_EV)
    assert unit_factors("real")[1] == pytest.approx(627.509474063, rel=1e-11)
    for style in ("lj", "si", "cgs", "micro", "nano"):
        with pytest.raises(ValueError, match="units metal, real, and electron"):
            unit_factors(style)


def test_box_rows_is_the_lammps_h_matrix_transposed():
    rows = box_rows([1.0, 2.0, 3.0], [10.0, 10.5, 12.4], 1.2, 0.9, -0.7)
    np.testing.assert_array_equal(rows, TRICLINIC)


def test_type_numbers():
    np.testing.assert_array_equal(type_numbers("O H", 2), [8, 1])
    np.testing.assert_array_equal(type_numbers(["H", 8, "8"], 3), [1, 8, 8])
    with pytest.raises(
        ValueError, match="names 1 atom types, but the LAMMPS box has 2"
    ):
        type_numbers(["O"], 2)
    with pytest.raises(ValueError, match="one element"):
        type_numbers(["OH"], 1)
    with pytest.raises(ValueError, match="unknown element"):
        type_numbers(["Qq"], 1)


# --- one step against ASE -----------------------------------------------


@pytest.mark.parametrize("route", ROUTES)
@pytest.mark.parametrize("units", ["metal", "real", "electron"])
@pytest.mark.parametrize("cell", [CUBE, TRICLINIC], ids=["cube", "triclinic"])
def test_periodic_step_matches_ase(route, units, cell, periodic_params):
    got, handles = run_inline(spec(route, units, periodic_params, cell=cell))
    finish(handles)
    ref = ase_reference(periodic_params, WATERS, cell)
    assert np.abs(ref["virial"]).max() > 0.0
    assert_matches(got, ref, units)


@pytest.mark.parametrize("route", ROUTES)
def test_cluster_step_matches_ase(route, params):
    got, handles = run_inline(spec(route, "metal", params, cell=None))
    finish(handles)
    ref = ase_reference(params, WATERS, None)
    assert_matches(got, ref, "metal")


@pytest.mark.parametrize("route", ROUTES)
def test_types_map_to_elements_in_type_order(route, periodic_params):
    case = spec(route, "metal", periodic_params, elements="H 8")
    case["types"] = [2 if t == 1 else 1 for t in TYPES]
    case["masses"] = MASSES[::-1]
    got, handles = run_inline(case)
    finish(handles)
    assert_matches(got, ase_reference(periodic_params, WATERS, TRICLINIC), "metal")


@pytest.mark.parametrize("route", ROUTES)
def test_positions_are_unwrapped(route, periodic_params):
    """An atom moved out of the box keeps its path; LAMMPS keeps it inside."""
    shift = np.array([-2.5, 9.0, 4.0])
    got, handles = run_inline(
        spec(route, "metal", periodic_params),
        f"displace_atoms all move {shift[0]} {shift[1]} {shift[2]} units box",
        "run 0",
    )
    lmp, _ = handles
    wrapped = np.array(lmp.gather_atoms("x", 1, 3)).reshape(-1, 3)
    finish(handles)
    moved = np.asarray(WATERS) + shift
    assert np.abs(wrapped - moved).max() > 1.0
    assert_matches(got, ase_reference(periodic_params, moved, TRICLINIC), "metal")


@pytest.mark.parametrize("route", ROUTES)
def test_thermo_pressure_is_the_virial(route, periodic_params):
    got, handles = run_inline(spec(route, "metal", periodic_params))
    lmp, _ = handles
    volume = lmp.get_thermo("vol")
    pressure = lmp.get_thermo("press")
    nktv2p = lmp.extract_global("nktv2p")
    finish(handles)
    # No velocities: the pressure is the trace of the virial over 3V.
    trace = sum(got["virial"][:3])
    assert pressure == pytest.approx(trace / (3.0 * volume) * nktv2p, rel=1e-12)


# --- dynamics -----------------------------------------------------------


def drift(log: np.ndarray, column: int) -> float:
    """Largest departure of a conserved column over the swing of pe."""
    swing = log[:, 1].max() - log[:, 1].min()
    assert swing > 0.0
    return float(np.abs(log[:, column] - log[0, column]).max() / swing)


def thermo_rows(lmp, steps: int, every: int, *columns: str) -> np.ndarray:
    rows = []
    for _ in range(steps // every):
        lmp.command(f"run {every}")
        rows.append([lmp.get_thermo(name) for name in columns])
    return np.array(rows)


@pytest.mark.parametrize("route", ROUTES)
def test_nve_conserves_energy_across_the_boundary(route, periodic_params):
    """The reference well pulls both waters through the box corner."""
    case = spec(route, "metal", periodic_params, cell=CUBE)
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string(
        "\n".join(
            [
                "velocity all create 300 4928 loop geom",
                "fix nve all nve",
                "timestep 0.00025",
                "thermo_style custom step pe ke etotal",
            ]
        )
    )
    log = thermo_rows(lmp, 800, 20, "step", "pe", "ke", "etotal")
    image = np.array(lmp.gather_atoms("image", 0, 1))
    if driver is not None:
        driver.close()
    lmp.close()
    unmoved = 512 | 512 << 10 | 512 << 20
    assert np.any(image != unmoved), "no atom crossed the box"
    assert drift(log, 3) < 2e-3


@pytest.mark.parametrize("route", ROUTES)
def test_npt_conserves_the_extended_energy(route, periodic_params):
    """The barostat integrates the cpmdc virial; econserve stays put.

    The reference energy is not translation invariant, so the box dilates
    about the origin, the point the strain derivative refers to. The target
    stress sits near the reference's own so the volume stays moderate.
    """
    case = spec(route, "metal", periodic_params, cell=CUBE)
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string(
        "\n".join(
            [
                "velocity all create 300 4928 loop geom",
                "fix npt all npt temp 300 300 0.025 aniso -20000 -20000 0.25 fixedpoint 0 0 0",
                "timestep 0.00025",
                "thermo_style custom step pe vol econserve",
            ]
        )
    )
    log = thermo_rows(lmp, 400, 20, "step", "pe", "vol", "econserve")
    if driver is not None:
        driver.close()
    lmp.close()
    assert abs(log[-1, 2] - log[0, 2]) > 1e-3 * log[0, 2]
    assert drift(log, 3) < 5e-3


@pytest.mark.parametrize("route", ROUTES)
def test_virial_is_the_strain_derivative(route, periodic_params):
    """Every component LAMMPS reports equals -dE/d(strain) on a tilted cell.

    ``aniso`` and ``tri`` NPT runs on a tilted box do not isolate the shear
    terms (``aniso`` strains a tilted cell without accounting for it), so
    the six components are checked against central differences instead.
    Lower-triangular strains keep the cell in the restricted LAMMPS form.
    """
    got, handles = run_inline(spec(route, "metal", periodic_params))
    finish(handles)
    virial = np.array(got["virial"])
    cell = np.asarray(TRICLINIC)
    step = 1e-5
    voigt = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 0): 3, (2, 0): 4, (2, 1): 5}
    for (row, col), slot in voigt.items():
        energies = []
        for sign in (1.0, -1.0):
            strain = np.eye(3)
            strain[row, col] += sign * step
            case = spec(
                route,
                "metal",
                periodic_params,
                cell=cell @ strain,
                positions=np.asarray(WATERS) @ strain,
            )
            strained, handles = run_inline(case)
            finish(handles)
            energies.append(strained["pe"])
        derivative = (energies[0] - energies[1]) / (2.0 * step)
        assert -derivative == pytest.approx(
            virial[slot], abs=1e-6 * np.abs(virial).max()
        )


# --- errors -------------------------------------------------------------


@pytest.mark.parametrize("route", ROUTES)
def test_other_units_are_refused(route, periodic_params):
    case = spec(route, "metal", periodic_params)
    case["units"] = "lj"
    lmp = build(case)
    try:
        with pytest.raises(Exception, match="units metal, real, and electron"):
            attach(lmp, case)
    finally:
        lmp.close()


@pytest.mark.parametrize("route", ROUTES)
def test_mixed_periodicity_is_refused(route, periodic_params):
    case = spec(route, "metal", periodic_params, cell=CUBE)
    lmp = build(case)
    lmp.command("change_box all boundary p p f")
    try:
        with pytest.raises(Exception, match="fully periodic box"):
            attach(lmp, case)
            lmp.command("run 0")
    finally:
        lmp.close()


@pytest.mark.parametrize("route", ROUTES)
def test_wrong_element_count_is_refused(route, periodic_params):
    case = spec(route, "metal", periodic_params, elements=["O"])
    lmp = build(case)
    try:
        with pytest.raises(Exception, match="atom type"):
            attach(lmp, case)
    finally:
        lmp.close()


def test_plugin_reports_a_bad_params_file(tmp_path, periodic_params):
    if PLUGIN is None:
        pytest.skip("cpmdcplugin.so is not built")
    lmp = build(spec("plugin", "metal", periodic_params))
    try:
        lmp.command(f"plugin load {PLUGIN}")
        with pytest.raises(Exception, match="cannot read a CPMDParams message"):
            lmp.command(f"pair_style cpmdc {tmp_path / 'missing.bin'}")
        garbage = tmp_path / "garbage.bin"
        garbage.write_bytes(b"\x00" * 3)
        lmp.command(f"pair_style cpmdc {garbage}")
        lmp.command("pair_coeff * * O H")
        with pytest.raises(Exception, match="pair_style cpmdc: "):
            lmp.command("run 0")
    finally:
        lmp.close()


def test_failed_step_stops_the_run(periodic_params):
    case = spec("fix", "metal", periodic_params)
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string("fix nve all nve\nthermo_style custom step pe")
    try:
        driver.run(2)
        assert driver.evaluations == 3

        def broken(*_args, **_kwargs):
            raise CpmdcError("SCF did not converge")

        driver._session.calculate_au = broken
        with pytest.raises(CpmdcError, match="SCF did not converge"):
            driver.run(5)
        assert lmp.extract_global("ntimestep") == 3
        assert driver.error is None
    finally:
        driver.close()
        lmp.close()


def test_bad_params_raise_before_the_run(tmp_path, periodic_params):
    lmp = build(spec("fix", "metal", periodic_params))
    try:
        garbage = tmp_path / "garbage.bin"
        garbage.write_bytes(b"\x00" * 3)
        with pytest.raises(CpmdcError):
            FixExternal(lmp, garbage, "O H")
    finally:
        lmp.close()


# --- MPI ----------------------------------------------------------------


def run_case(case: dict, tmp_path: Path, ranks: int | None) -> dict:
    path = tmp_path / f"case-{case['route']}-{ranks}.json"
    path.write_text(json.dumps(case), encoding="utf-8")
    command = [sys.executable, str(HERE / "lammps_case.py"), str(path)]
    if ranks is not None:
        launcher = shutil.which("mpiexec") or shutil.which("mpirun")
        if launcher is None:
            pytest.skip("no mpiexec or mpirun on PATH")
        command = [launcher, "-n", str(ranks), *command]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(ROOT / "python"), env.get("PYTHONPATH", "")])
    )
    done = subprocess.run(
        command, capture_output=True, text=True, timeout=300, env=env, check=False
    )
    assert done.returncode == 0, done.stderr + done.stdout
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    assert len(lines) == 1, done.stdout
    return json.loads(lines[0])


@pytest.mark.parametrize("route", ROUTES)
def test_two_ranks_match_one(route, periodic_params, tmp_path):
    if route == "fix":
        pytest.importorskip("mpi4py")
    positions = [
        [0.6 + 1.4 * k, 1.0 + 0.9 * (k % 3), 1.2 + 1.1 * (k % 2)] for k in range(6)
    ]
    case = spec(
        route,
        "real",
        periodic_params,
        positions=positions,
        setup=["processors * 1 1"],
        run=[
            "velocity all create 300 77 loop geom",
            "fix nve all nve",
            "timestep 0.25",
            "run 40",
        ],
    )
    one = run_case(case, tmp_path, None)
    two = run_case(case, tmp_path, 2)
    if two["world_size"] != 2:
        pytest.skip("mpiexec did not start one two-rank MPI world that LAMMPS joins")
    assert min(two["nlocal"]) > 0, two["nlocal"]
    assert two["pe"] == pytest.approx(one["pe"], rel=1e-12)
    np.testing.assert_allclose(two["forces"], one["forces"], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(two["virial"], one["virial"], rtol=1e-10, atol=1e-10)
