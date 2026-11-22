"""cpmdc inside LAMMPS: cpmdc.lammps.FixExternal and pair_style cpmdc.

LAMMPS comes from its Python module; every test skips with the reason when
that module cannot start LAMMPS.
``rgpot.lammps`` is the frontend ``FixExternal`` subclasses. A released
rgpot wheel does not ship it, so fix-external cases skip when it is absent;
a broken installation still fails collection. Pair-style plugin cases do
not need that frontend. The plugin is ``CPMDC_LAMMPS_PLUGIN`` or
``build/cpmdcplugin.so`` of this checkout, and the plugin cases skip when
neither exists. Energies, forces, and virials are compared with
``CPMDCCalculator``, which reaches the library through ``ForceInput``
instead of ``cpmdc_session_calculate_au()``. The MPI cases run
``lammps_case.py`` under ``mpiexec``.

A test that holds only for one library is marked ``cpmd_backend`` with the
reason. Collection omits it unless ``conftest.linked_backend()`` matches,
the same bit the ``opencpmd`` fixture reads, and the fixture then supplies
that name. The other library does not collect the test and does not skip it.
"""

from __future__ import annotations

import ctypes
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from conftest import ROOT, park_orbitals
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, ConstantCutoff, CPMDCError
from cpmdc.ase import CPMDCCalculator
from cpmdc.lammps import UNIT_STYLES, FixExternal, box_rows, type_numbers, unit_factors

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lammps_case import CMDARGS, THERMO, attach, build, measure
from lammps_import import import_rgpot_lammps

RGPOT_LAMMPS = import_rgpot_lammps()
needs_rgpot_lammps = pytest.mark.skipif(
    RGPOT_LAMMPS is None, reason="rgpot.lammps is not installed"
)


def lammps_problem() -> str:
    try:
        import lammps
    except Exception as exc:  # ImportError, or a missing MPI library
        return f"the lammps Python module is not importable: {exc}"
    try:
        lmp = lammps.lammps(cmdargs=list(CMDARGS))
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
    pytest.param("fix", marks=needs_rgpot_lammps),
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


def ase_reference(params: Path, positions, cell, *, symbols=SYMBOLS, **calculator) -> dict:
    """eV, eV/Angstrom, and the virial -dE/d(strain) in eV, xx yy zz xy xz yz."""
    atoms = Atoms(symbols, positions=positions)
    if cell is not None:
        atoms.cell = cell
        atoms.pbc = True
    with CPMDCCalculator(params, **calculator) as calc:
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
def test_accepts_sparse_atom_ids(route, periodic_params):
    removed = {1, 4}
    kept = [position for index, position in enumerate(WATERS) if index not in removed]
    symbols = "".join(symbol for index, symbol in enumerate(SYMBOLS) if index not in removed)
    case = spec(route, "metal", periodic_params)
    case["after_create"] = [
        "group cpmdc_removed id 2 5",
        "delete_atoms group cpmdc_removed compress no",
    ]
    lmp = build(case)
    driver = attach(lmp, case)
    try:
        lmp.commands_string("\n".join(THERMO))
        if driver is None:
            lmp.command("run 0")
        else:
            driver.run(0)
        assert lmp.get_thermo("pe") == pytest.approx(
            ase_reference(periodic_params, kept, TRICLINIC, symbols=symbols)["pe"],
            rel=1e-11,
            abs=1e-12,
        )
        nlocal = int(lmp.extract_setting("nlocal"))
        tags = np.asarray(lmp.numpy.extract_atom("id")[:nlocal], dtype=np.int64)
        assert set(tags) == {1, 3, 4, 6}
        forces = np.asarray(lmp.numpy.extract_atom("f")[:nlocal], dtype=float)
        expected = ase_reference(periodic_params, kept, TRICLINIC, symbols=symbols)["forces"]
        np.testing.assert_allclose(forces[np.argsort(tags)], expected, rtol=1e-10, atol=1e-12)
    finally:
        if driver is not None:
            driver.close()
        lmp.close()


@pytest.mark.parametrize("route", ROUTES)
def test_mid_run_deletion_reopens_with_sparse_atom_ids(route, periodic_params):
    removed = {1, 4}
    kept = [position for index, position in enumerate(WATERS) if index not in removed]
    symbols = "".join(symbol for index, symbol in enumerate(SYMBOLS) if index not in removed)
    case = spec(route, "metal", periodic_params)
    lmp = build(case)
    driver = attach(lmp, case)
    try:
        lmp.command("thermo_style custom step pe")
        if driver is None:
            lmp.command("run 0")
        else:
            driver.run(0)
        old_session = driver._session if driver is not None else None
        lmp.commands_string(
            "group cpmdc_removed id 2 5\n"
            "delete_atoms group cpmdc_removed compress no"
        )
        if driver is None:
            lmp.command("run 0")
        else:
            driver.run(0)
            assert driver._session is not old_session
        assert lmp.get_thermo("pe") == pytest.approx(
            ase_reference(periodic_params, kept, TRICLINIC, symbols=symbols)["pe"],
            rel=1e-11,
            abs=1e-12,
        )
        nlocal = int(lmp.extract_setting("nlocal"))
        tags = np.asarray(lmp.numpy.extract_atom("id")[:nlocal], dtype=np.int64)
        assert set(tags) == {1, 3, 4, 6}
        forces = np.asarray(lmp.numpy.extract_atom("f")[:nlocal], dtype=float)
        expected = ase_reference(periodic_params, kept, TRICLINIC, symbols=symbols)["forces"]
        np.testing.assert_allclose(forces[np.argsort(tags)], expected, rtol=1e-10, atol=1e-12)
    finally:
        if driver is not None:
            driver.close()
        lmp.close()


@pytest.mark.parametrize("route", ROUTES)
def test_duplicate_atom_ids_are_refused(route, periodic_params):
    case = spec(route, "metal", periodic_params)
    lmp = build(case)
    tags = lmp.numpy.extract_atom("id")
    tags[1] = tags[0]
    driver = attach(lmp, case)
    try:
        lmp.command("thermo_style custom step pe")
        with pytest.raises(Exception, match="unique atom IDs"):
            if driver is None:
                lmp.command("run 0")
            else:
                driver.run(0)
    finally:
        if driver is not None:
            driver.close()
        lmp.close()


@pytest.mark.parametrize("route", ROUTES)
def test_types_map_to_elements_in_type_order(route, periodic_params):
    case = spec(route, "metal", periodic_params, elements="H 8")
    case["types"] = [2 if t == 1 else 1 for t in TYPES]
    case["masses"] = MASSES[::-1]
    got, handles = run_inline(case)
    finish(handles)
    assert_matches(got, ase_reference(periodic_params, WATERS, TRICLINIC), "metal")


def unwrapped_step(route: str, params: Path) -> tuple[dict, dict]:
    """LAMMPS and ASE at the two waters moved out of the triclinic box."""
    shift = np.array([-2.5, 9.0, 4.0])
    got, handles = run_inline(
        spec(route, "metal", params),
        f"displace_atoms all move {shift[0]} {shift[1]} {shift[2]} units box",
        "run 0",
    )
    lmp, _ = handles
    wrapped = np.array(lmp.gather_atoms("x", 1, 3)).reshape(-1, 3)
    finish(handles)
    moved = np.asarray(WATERS) + shift
    assert np.abs(wrapped - moved).max() > 1.0
    return got, ase_reference(params, moved, TRICLINIC)


@pytest.mark.cpmd_backend(
    "reference",
    reason=(
        "the reference forces are linear in the coordinates, so the 1e-15 "
        "Angstrom a wrap and unwrap leaves cannot move them past rtol 1e-10; "
        "OpenCPMD's SCF answers it with 1.3e-11 eV/Angstrom"
    ),
)
@pytest.mark.parametrize("route", ROUTES)
def test_positions_are_unwrapped(route, periodic_params, cpmd_backend):
    """An atom moved out of the box keeps its path; LAMMPS keeps it inside.

    The reference well is linear, so this comparison holds at the
    tolerances of the single steps.
    """
    assert cpmd_backend == "reference"
    got, ref = unwrapped_step(route, periodic_params)
    assert_matches(got, ref, "metal")


@pytest.mark.cpmd_backend(
    "opencpmd",
    reason=(
        "the SCF answers the 1e-15 Angstrom of a wrap and unwrap with changes "
        "near 1e-11 eV/Angstrom; the reference well runs the exact comparison"
    ),
)
@pytest.mark.parametrize("route", ROUTES)
def test_positions_are_unwrapped_on_plane_waves(
    route, periodic_params, cpmd_backend, record_property
):
    """The unwrapped path on OpenCPMD, at the roundoff of its SCF.

    LAMMPS wraps the atoms into the box and ``cpmdc`` unwraps them, which
    restores the path to about 1e-15 Angstrom. OpenCPMD's converged forces
    answer that with differences near 1e-11 eV/Angstrom, so the forces are
    held to 1e-9 eV/Angstrom absolute.
    """
    assert cpmd_backend == "opencpmd"
    got, ref = unwrapped_step(route, periodic_params)
    for name in ("forces", "virial"):
        miss = np.abs(np.subtract(got[name], ref[name])).max()
        record_property(f"{name}_miss", float(miss))
    record_property("pe_miss_ev", abs(got["pe"] - ref["pe"]))
    assert got["pe"] == pytest.approx(ref["pe"], rel=1e-11, abs=1e-9)
    np.testing.assert_allclose(got["forces"], ref["forces"], rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(got["virial"], ref["virial"], rtol=1e-8, atol=1e-8)


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


# A, sigma, and E0 in Ry, the order of system.constantCutoff.
CONSTANT_CUTOFF = (2000.0, 3.0, 30.0)
# The deck stops at a largest orbital gradient of 1e-5, which leaves an
# energy error that a strain of 1e-5 or a 1e-3 Angstrom step divides into
# the difference. At 1e-7 that error is below the tolerances of the
# finite differences here.
TIGHT_ORBITALS = 1.0e-7


def periodic_deck(tmp_path_factory, name: str, *, orbitals=None, system="") -> Path:
    """The periodic water message, with a tighter SCF or more of ``system``."""
    if shutil.which("capnp") is None:
        pytest.skip("capnp is not on PATH")
    text = (ROOT / "tests" / "water_blyp_periodic_params.capnp.txt").read_text()
    if orbitals is not None:
        text, count = re.subn(
            r"convergenceOrbitals = [^,]+,", f"convergenceOrbitals = {orbitals!r},", text
        )
        assert count == 1
        # A cold setup needs more than the deck's 40 iterations to reach 1e-7.
        text, count = re.subn(
            r"maxStep = 40, maxIter = 40", "maxStep = 100, maxIter = 100", text
        )
        assert count == 1
    if system:
        text, count = re.subn(
            r"(cell = \[[^\]]*\])", lambda m: f"{m.group(1)}, {system}", text
        )
        assert count == 1
    target = tmp_path_factory.mktemp(name) / f"{name}.bin"
    target.write_bytes(subprocess.run(
        ["capnp", "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode(), stdout=subprocess.PIPE, check=True,
    ).stdout)
    return target


def constant_cutoff_system() -> str:
    return "constantCutoff = [%r, %r, %r]" % CONSTANT_CUTOFF


@pytest.fixture(scope="module")
def constant_cutoff_params(tmp_path_factory) -> Path:
    """The periodic water message with system.constantCutoff."""
    return periodic_deck(
        tmp_path_factory, "water_constant_cutoff", system=constant_cutoff_system()
    )


@pytest.fixture(scope="module")
def tight_params(tmp_path_factory) -> Path:
    """The periodic water message at an orbital threshold of 1e-7."""
    return periodic_deck(tmp_path_factory, "water_tight", orbitals=TIGHT_ORBITALS)


@pytest.fixture(scope="module")
def tight_constant_cutoff_params(tmp_path_factory) -> Path:
    """The constant-cutoff message at an orbital threshold of 1e-7."""
    return periodic_deck(
        tmp_path_factory,
        "water_tight_constant_cutoff",
        orbitals=TIGHT_ORBITALS,
        system=constant_cutoff_system(),
    )


@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="the reference well has no plane waves",
)
@pytest.mark.parametrize("route", ROUTES)
def test_constant_cutoff_virial_matches_ase(
    route, periodic_params, constant_cutoff_params, cpmd_backend
):
    """system.constantCutoff reaches the virial that LAMMPS integrates.

    It is the virial of CPMDCCalculator with the same functional from its
    constant_cutoff keyword, and lies far from the fixed-basis virial of
    the plain deck, which carries the Pulay stress.
    """
    assert cpmd_backend == "opencpmd"
    got, handles = run_inline(spec(route, "metal", constant_cutoff_params, cell=CUBE))
    finish(handles)
    height, width, cutoff = CONSTANT_CUTOFF
    corrected = ase_reference(
        periodic_params, WATERS, CUBE,
        constant_cutoff=ConstantCutoff(cutoff, height_ry=height, width_ry=width),
    )
    assert_matches(got, corrected, "metal")
    fixed = ase_reference(periodic_params, WATERS, CUBE)
    assert np.mean(np.subtract(fixed["virial"], got["virial"])[:3]) < -5.0


# --- finite differences ---------------------------------------------------

VOIGT = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 0): 3, (2, 0): 4, (2, 1): 5}


def step_virial(route: str, params: Path, cell) -> np.ndarray:
    got, handles = run_inline(spec(route, "metal", params, cell=cell))
    finish(handles)
    return np.array(got["virial"])


def strain_virial(route: str, params: Path, cell, step: float) -> np.ndarray:
    """-dE/d(strain) by central differences in eV, xx yy zz xy xz yz.

    Lower-triangular strains keep the cell in the restricted LAMMPS form.
    Each strained cell is a new LAMMPS system and a new session.
    """
    cell = np.asarray(cell)
    out = np.zeros(6)
    for (row, col), slot in VOIGT.items():
        energies = []
        for sign in (1.0, -1.0):
            strain = np.eye(3)
            strain[row, col] += sign * step
            case = spec(
                route,
                "metal",
                params,
                cell=cell @ strain,
                positions=np.asarray(WATERS) @ strain,
            )
            strained, handles = run_inline(case)
            finish(handles)
            energies.append(strained["pe"])
        out[slot] = -(energies[0] - energies[1]) / (2.0 * step)
    return out


def largest_miss(derivative: np.ndarray, virial: np.ndarray) -> float:
    """Largest component of the difference, over the largest of the virial."""
    return float(np.abs(derivative - virial).max() / np.abs(virial).max())


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="the reference forces are linear, so a difference is exact there",
)
def test_forces_are_the_energy_gradient(tight_params, cpmd_backend, record_property):
    """Each force LAMMPS reports is minus the slope of the energy it reports.

    A central difference with a 1e-3 Angstrom step on each coordinate of
    the two waters, at an orbital threshold of 1e-7, met the forces to
    2.6e-5 eV/Angstrom, and a 2e-3 Angstrom step to 1.1e-4, the h^2 of the
    difference. At the deck's 1e-5 the 1e-3 step missed by 5.4e-4. The
    route does not change the numbers; the single steps match ASE on both.
    """
    assert cpmd_backend == "opencpmd"
    case = spec("fix", "metal", tight_params)
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string("\n".join(THERMO))
    step = 1e-3
    try:
        driver.run(0)
        forces = np.array(measure(lmp)["forces"])
        start = np.array(lmp.gather_atoms("x", 1, 3)).reshape(-1, 3)

        def energy(positions: np.ndarray) -> float:
            flat = positions.ravel()
            lmp.scatter_atoms("x", 1, 3, (ctypes.c_double * flat.size)(*flat))
            driver.run(0)
            return lmp.get_thermo("pe")

        difference = np.zeros_like(start)
        for atom, axis in np.ndindex(*start.shape):
            ahead, behind = start.copy(), start.copy()
            ahead[atom, axis] += step
            behind[atom, axis] -= step
            difference[atom, axis] = -(energy(ahead) - energy(behind)) / (2.0 * step)
    finally:
        driver.close()
        lmp.close()
    record_property("largest_force_ev_per_a", float(np.abs(forces).max()))
    record_property("miss_ev_per_a", float(np.abs(difference - forces).max()))
    assert np.abs(forces).max() > 1.0
    np.testing.assert_allclose(difference, forces, rtol=0.0, atol=1e-4)


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="test_virial_is_the_strain_derivative runs on the reference well",
)
def test_virial_is_the_strain_derivative_of_one_basis(
    tight_params, cpmd_backend, record_property
):
    """The plain-deck virial is -dE/d(strain) at a fixed set of plane waves.

    Every strained cell is a new setup at the deck's 40 Ry, so the
    difference holds the basis only while no plane wave crosses that
    sphere. On this triclinic cell a strain of 1e-5 keeps all of them; at
    3e-5 the difference already misses by 1.8 times the largest component.
    At an orbital threshold of 1e-7 the 1e-5 difference met the virial to
    2.5e-7 of the largest component. At the deck's 1e-5 it missed by
    1.1e-3, the SCF error over the step.
    """
    assert cpmd_backend == "opencpmd"
    virial = step_virial("fix", tight_params, TRICLINIC)
    derivative = strain_virial("fix", tight_params, TRICLINIC, 1e-5)
    record_property("miss", largest_miss(derivative, virial))
    assert largest_miss(derivative, virial) < 1e-6


def richardson(route: str, params: Path, cell, step: float) -> tuple[np.ndarray, ...]:
    """Differences at step and 2 step, and the extrapolation of the two."""
    fine = strain_virial(route, params, cell, step)
    coarse = strain_virial(route, params, cell, 2.0 * step)
    return fine, coarse, (4.0 * fine - coarse) / 3.0


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="the reference well has no plane waves",
)
def test_virial_is_the_strain_derivative_in_place(
    params, tight_constant_cutoff_params, monkeypatch, cpmd_backend, record_property
):
    """Strains of 2.5e-3 and 5e-3 on one held basis meet the virial.

    ``CPMDC_CELL_IN_PLACE`` lets a constant-cutoff deck keep the plane
    waves of its setup cell through a strain, so every difference here
    sees one basis. The misses were 3.2e-4 and 1.3e-3 of the largest
    component, a ratio of 4.0, the h^2 of a central difference, and the
    Richardson extrapolation of the two met the virial to 4.8e-6.
    """
    assert cpmd_backend == "opencpmd"
    monkeypatch.setenv("CPMDC_CELL_IN_PLACE", "0.01")
    park_orbitals(params)
    virial = step_virial("fix", tight_constant_cutoff_params, TRICLINIC)
    fine, coarse, extrapolated = richardson(
        "fix", tight_constant_cutoff_params, TRICLINIC, 2.5e-3
    )
    ratio = largest_miss(coarse, virial) / largest_miss(fine, virial)
    record_property("miss_5e-3", largest_miss(coarse, virial))
    record_property("miss_2.5e-3", largest_miss(fine, virial))
    record_property("miss_richardson", largest_miss(extrapolated, virial))
    assert 3.0 < ratio < 5.0
    assert largest_miss(extrapolated, virial) < 5e-5


@pytest.fixture(scope="module")
def reference_cell_params(tmp_path_factory) -> Path:
    """The tight message with REFERENCE CELL and MESH at TRICLINIC.

    72 72 80 is the mesh OpenCPMD picks for TRICLINIC at 40 Ry. Without
    MESH, a strain of 5e-3 already moves it to 72 72 72.
    """
    rows = ", ".join(repr(float(value)) for value in np.ravel(TRICLINIC))
    return periodic_deck(
        tmp_path_factory,
        "water_reference_cell",
        orbitals=TIGHT_ORBITALS,
        system=(
            f"referenceCell = [{rows}], referenceCellVectors = true, "
            "mesh = [72, 72, 80]"
        ),
    )


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="the reference well has no plane waves",
)
@pytest.mark.xfail(
    strict=True,
    reason=(
        "with REFERENCE CELL a cold setup builds the G vectors from the "
        "reference lattice, not the strained cell: the strain difference is "
        "+3.28 eV where the virial is -32.42 eV (xx), at a step of 1e-5 and "
        "of 5e-3 alike"
    ),
)
def test_reference_cell_holds_the_basis(
    reference_cell_params, cpmd_backend, record_property
):
    """REFERENCE CELL holds the plain deck's plane waves through a strain.

    Every strained cell keeps the 10347 plane waves of TRICLINIC and its
    mesh, so the Richardson extrapolation of strains of 2.5e-3 and 5e-3
    should meet the fixed-basis virial as the in-place differences do.
    """
    assert cpmd_backend == "opencpmd"
    virial = step_virial("fix", reference_cell_params, TRICLINIC)
    fine, coarse, extrapolated = richardson(
        "fix", reference_cell_params, TRICLINIC, 2.5e-3
    )
    record_property("virial_ev", virial.tolist())
    record_property("difference_2.5e-3_ev", fine.tolist())
    record_property("miss_5e-3", largest_miss(coarse, virial))
    record_property("miss_2.5e-3", largest_miss(fine, virial))
    record_property("miss_richardson", largest_miss(extrapolated, virial))
    assert largest_miss(extrapolated, virial) < 5e-5


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="the reference well has no plane waves",
)
def test_constant_cutoff_virial_is_the_strain_derivative(
    tight_constant_cutoff_params, monkeypatch, cpmd_backend, record_property
):
    """The constant-cutoff virial is -dE/d(strain) of new setups.

    Each strained cube is a new setup with constant_cutoff on. At a strain
    of 1e-4 no plane wave crosses the 40 Ry sphere, and the difference met
    the virial to 6.9e-7 of the largest component, the SCF error of the
    orbital threshold of 1e-7 over that step. At 5e-3 the sets
    differ, and the components between E0 + 2 sigma (36 Ry) and the cutoff
    that the functional leaves nearly empty moved the difference by 2.0e-2;
    the plain deck at 1e-3 misses by 1.3 times the largest component.
    """
    assert cpmd_backend == "opencpmd"
    monkeypatch.delenv("CPMDC_CELL_IN_PLACE", raising=False)
    virial = step_virial("fix", tight_constant_cutoff_params, CUBE)
    near = strain_virial("fix", tight_constant_cutoff_params, CUBE, 1e-4)
    record_property("miss_1e-4", largest_miss(near, virial))
    assert largest_miss(near, virial) < 1e-5
    far = strain_virial("fix", tight_constant_cutoff_params, CUBE, 5e-3)
    record_property("miss_5e-3", largest_miss(far, virial))
    assert largest_miss(far, virial) < 4e-2


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


@pytest.mark.cpmd_backend(
    "reference",
    reason=(
        "the bound is 2e-3 of the swing of pe, the fall into the well; on "
        "OpenCPMD that swing is 0.4 eV of thermal motion, and the Verlet "
        "fluctuation of etotal at this step is 8.6e-4 eV RMS, above the "
        "8e-4 eV bound"
    ),
)
@pytest.mark.parametrize("route", ROUTES)
def test_nve_conserves_energy_across_the_boundary(route, periodic_params, cpmd_backend):
    """The reference well pulls both waters through the box corner.

    etotal stays within 2e-3 of the swing of pe over the 800 steps.
    """
    assert cpmd_backend == "reference"
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


def nve_etotal(params: Path, timestep_fs: float, steps: int, every: int) -> np.ndarray:
    """Time in ps and etotal in eV of the two waters at 300 K in the cube."""
    case = spec("fix", "metal", params, cell=CUBE)
    lmp = build(case)
    driver = attach(lmp, case)
    try:
        lmp.commands_string(
            "\n".join(
                [
                    "velocity all create 300 4928 loop geom",
                    "fix nve all nve",
                    f"timestep {timestep_fs * 1e-3!r}",
                    "thermo_style custom step pe ke etotal",
                    "run 0",
                ]
            )
        )
        first = [lmp.get_thermo("step"), lmp.get_thermo("etotal")]
        log = np.vstack([first, thermo_rows(lmp, steps, every, "step", "etotal")])
    finally:
        driver.close()
        lmp.close()
    return np.column_stack([log[:, 0] * timestep_fs * 1e-3, log[:, 1]])


def drift_and_wobble(log: np.ndarray) -> tuple[float, float]:
    """Change of the straight line through etotal, and the RMS about it."""
    time, etotal = log[:, 0], log[:, 1]
    line = np.polyfit(time, etotal, 1)
    return float(line[0] * (time[-1] - time[0])), float(
        np.std(etotal - np.polyval(line, time))
    )


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="test_nve_conserves_energy_across_the_boundary runs on the reference well",
)
def test_nve_energy_error_falls_with_the_timestep_squared(
    tight_params, cpmd_backend, record_property
):
    """etotal holds, and its Verlet fluctuation is quartered by half the step.

    Two runs cover the same 50 fs from the same velocities: 100 steps of
    0.5 fs and 200 of 0.25 fs, read every 0.5 fs, at an orbital threshold
    of 1e-7. The straight line through etotal moved 8.5e-5 eV over the
    0.25 fs run. The RMS about it was 3.45e-3 eV at 0.5 fs and 8.6e-4 eV at
    0.25 fs, a ratio of 4.0: velocity Verlet's error, not SCF noise, which
    would not scale with the step.

    At the deck's 1e-5 the line falls at 2.9e-2 eV/ps, 17 times the 1.7e-3
    eV/ps at 1e-7, while the RMS stays at 8.6e-4 eV. The 800 steps of the
    reference test, 200 fs, add that fall to the fluctuation.
    """
    assert cpmd_backend == "opencpmd"
    coarse = nve_etotal(tight_params, 0.5, 100, 1)
    fine = nve_etotal(tight_params, 0.25, 200, 2)
    np.testing.assert_allclose(coarse[:, 0], fine[:, 0])
    drift, fine_wobble = drift_and_wobble(fine)
    _, coarse_wobble = drift_and_wobble(coarse)
    record_property("drift_ev", drift)
    record_property("rms_0.5fs_ev", coarse_wobble)
    record_property("rms_0.25fs_ev", fine_wobble)
    assert abs(drift) < 5e-4
    assert 3.0 < coarse_wobble / fine_wobble < 5.0


@pytest.mark.cpmd_backend(
    "reference",
    reason=(
        "the target, -20000 bar, sits near the well's own stress; the plain "
        "OpenCPMD deck sets up a new basis whenever the cell moves, and at a "
        "strain of 1e-3 its energy departs from its virial by 1.3 times the "
        "largest component"
    ),
)
@pytest.mark.parametrize("route", ROUTES)
def test_npt_conserves_the_extended_energy(route, periodic_params, cpmd_backend):
    """The barostat integrates the cpmdc virial; econserve stays put.

    The reference energy is not translation invariant, so the box dilates
    about the origin, the point the strain derivative refers to. The target
    stress sits near the reference's own so the volume stays moderate.
    """
    assert cpmd_backend == "reference"
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


@needs_rgpot_lammps
@pytest.mark.cpmd_backend(
    "opencpmd",
    reason="test_npt_conserves_the_extended_energy runs on the reference well",
)
def test_npt_conserves_the_extended_energy_at_a_constant_cutoff(
    params, tight_constant_cutoff_params, monkeypatch, cpmd_backend, record_property
):
    """An iso NPT run of the two waters on the constant-cutoff virial.

    A plain deck sets up a new basis whenever the cell moves, which moves
    its energy by the Pulay stress, so the run uses the constant-cutoff
    deck at an orbital threshold of 1e-7 and ``CPMDC_CELL_IN_PLACE=0.05``,
    which keeps the plane waves of the setup through these volumes.

    The target is 0 bar: the pressure of the two waters averaged -135 bar
    over the run. The pressure at the first step, +18 kbar, is the virial
    of the stretched O-H bonds, which swings by about 15 kbar each
    vibration; as the target, it compressed the cube from 729 to 202
    Angstrom^3 in 50 fs.

    Over 200 steps of 0.25 fs with a 0.25 ps piston the volume moved by
    3.7%, and econserve left its start by at most 6.2e-3 of the 1.85 eV
    swing of pe. NVE on the same deck and step leaves by 5.9e-3 of its
    1.89 eV swing. The constant-cutoff forces at the start reach 9.1
    eV/Angstrom, twice the 4.5 of the plain deck in this cube, so the
    waters ring harder than in the NVE test.
    """
    assert cpmd_backend == "opencpmd"
    monkeypatch.setenv("CPMDC_CELL_IN_PLACE", "0.05")
    park_orbitals(params)
    case = spec("fix", "metal", tight_constant_cutoff_params, cell=CUBE)
    lmp = build(case)
    driver = attach(lmp, case)
    try:
        lmp.commands_string(
            "\n".join(
                [
                    "velocity all create 300 4928 loop geom",
                    "fix npt all npt temp 300 300 0.025 iso 0 0 0.25",
                    "timestep 0.00025",
                    "thermo_style custom step pe vol econserve press",
                    "run 0",
                ]
            )
        )
        columns = ("step", "pe", "vol", "econserve", "press")
        first = [lmp.get_thermo(name) for name in columns]
        log = np.vstack([first, thermo_rows(lmp, 200, 2, *columns)])
    finally:
        driver.close()
        lmp.close()
    moved = (log[:, 2].max() - log[:, 2].min()) / log[0, 2]
    record_property("volume_range", float(moved))
    record_property("econserve_departure", drift(log, 3))
    record_property("pe_swing_ev", float(log[:, 1].max() - log[:, 1].min()))
    record_property("mean_pressure_bar", float(log[:, 4].mean()))
    assert moved > 1e-2
    assert drift(log, 3) < 1e-2


@pytest.mark.cpmd_backend(
    "reference",
    reason=(
        "the reference energy is quadratic, so the difference is exact; on "
        "OpenCPMD the deck's 1e-5 orbital threshold over a 1e-5 strain "
        "misses by 1.1e-3"
    ),
)
@pytest.mark.parametrize("route", ROUTES)
def test_virial_is_the_strain_derivative(route, periodic_params, cpmd_backend):
    """Every component LAMMPS reports equals -dE/d(strain) on a tilted cell.

    ``aniso`` and ``tri`` NPT runs on a tilted box do not isolate the shear
    terms (``aniso`` strains a tilted cell without accounting for it), so
    the six components are checked against central differences instead.
    The reference energy is quadratic, so the difference is exact.
    """
    assert cpmd_backend == "reference"
    virial = step_virial(route, periodic_params, TRICLINIC)
    derivative = strain_virial(route, periodic_params, TRICLINIC, 1e-5)
    assert largest_miss(derivative, virial) < 1e-6


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


PER_ATOM = [
    "compute peratom all pe/atom",
    "compute peratom all stress/atom NULL",
    "compute peratom all pe/atom pair",
    "compute peratom all centroid/stress/atom NULL",
]


@pytest.mark.parametrize("command", PER_ATOM)
def test_plugin_refuses_per_atom_energy_and_virial_before_evaluating(
    command, tmp_path, periodic_params
):
    """The energy and the virial are totals; a compute that wants a split is refused.

    The params file is garbage, so a step that reached cpmdc would fail with a
    session error instead of the refusal.
    """
    if PLUGIN is None:
        pytest.skip("cpmdcplugin.so is not built")
    garbage = tmp_path / "garbage.bin"
    garbage.write_bytes(b"\x00" * 3)
    lmp = build(spec("plugin", "metal", periodic_params))
    try:
        attach(lmp, spec("plugin", "metal", garbage))
        lmp.command(command)
        with pytest.raises(
            Exception,
            match=r"computes no per-atom energy or virial \(compute peratom\)",
        ):
            lmp.command("run 0")
        lmp.command("uncompute peratom")
        with pytest.raises(Exception, match="pair_style cpmdc: "):
            lmp.command("run 0")
    finally:
        lmp.close()


def test_plugin_runs_after_a_per_atom_compute_is_removed(periodic_params):
    if PLUGIN is None:
        pytest.skip("cpmdcplugin.so is not built")
    case = spec("plugin", "metal", periodic_params)
    lmp = build(case)
    try:
        attach(lmp, case)
        lmp.command("compute peratom all pe/atom")
        with pytest.raises(Exception, match="computes no per-atom energy or virial"):
            lmp.command("run 0")
        lmp.command("uncompute peratom")
        lmp.commands_string("\n".join(THERMO))
        lmp.command("run 0")
        assert_matches(measure(lmp), ase_reference(periodic_params, WATERS, TRICLINIC), "metal")
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


@needs_rgpot_lammps
def test_failed_step_stops_the_run(periodic_params):
    case = spec("fix", "metal", periodic_params)
    lmp = build(case)
    driver = attach(lmp, case)
    lmp.commands_string("fix nve all nve\nthermo_style custom step pe")
    try:
        driver.run(2)
        assert driver.evaluations == 3

        def broken(*_args, **_kwargs):
            raise CPMDCError("SCF did not converge")

        driver._session.calculate_au = broken
        with pytest.raises(CPMDCError, match="SCF did not converge"):
            driver.run(5)
        assert lmp.extract_global("ntimestep") == 3
        assert driver.error is None
    finally:
        driver.close()
        lmp.close()


@needs_rgpot_lammps
def test_bad_params_raise_before_the_run(tmp_path, periodic_params):
    lmp = build(spec("fix", "metal", periodic_params))
    try:
        garbage = tmp_path / "garbage.bin"
        garbage.write_bytes(b"\x00" * 3)
        with pytest.raises(CPMDCError):
            FixExternal(lmp, garbage, "O H")
    finally:
        lmp.close()


# --- MPI ----------------------------------------------------------------


def case_env() -> dict[str, str]:
    """Environment for ``lammps_case.py``, in process or under a launcher.

    Every case runs on one host. On GitHub runners with the Azure MANA NIC,
    UCX finds ``mana_0`` but cannot open it, and Open MPI ranks of the lammps
    environment then hang until the test times out. ``ob1`` over ``self``
    and ``sm`` needs no network device. MPICH ignores these variables, and
    values already in the environment win.
    """
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(ROOT / "python"), env.get("PYTHONPATH", "")])
    )
    env.setdefault("OMPI_MCA_pml", "ob1")
    env.setdefault("OMPI_MCA_btl", "self,sm")
    return env


def run_case(case: dict, tmp_path: Path, ranks: int | None) -> dict:
    path = tmp_path / f"case-{case['route']}-{ranks}.json"
    path.write_text(json.dumps(case), encoding="utf-8")
    command = [sys.executable, str(HERE / "lammps_case.py"), str(path)]
    if ranks is not None:
        launcher = shutil.which("mpiexec") or shutil.which("mpirun")
        if launcher is None:
            pytest.skip("no mpiexec or mpirun on PATH")
        command = [launcher, "-n", str(ranks), *command]
    done = subprocess.run(
        command, capture_output=True, text=True, timeout=300, env=case_env(), check=False
    )
    assert done.returncode == 0, done.stderr + done.stdout
    lines = [line for line in done.stdout.splitlines() if line.startswith("{")]
    assert len(lines) == 1, done.stdout
    return json.loads(lines[0])


def run_case_failure(case: dict, tmp_path: Path, ranks: int) -> subprocess.CompletedProcess:
    path = tmp_path / f"case-failure-{case['route']}-{ranks}.json"
    path.write_text(json.dumps(case), encoding="utf-8")
    launcher = shutil.which("mpiexec") or shutil.which("mpirun")
    if launcher is None:
        pytest.skip("no mpiexec or mpirun on PATH")
    return subprocess.run(
        [launcher, "-n", str(ranks), sys.executable, str(HERE / "lammps_case.py"), str(path)],
        capture_output=True,
        text=True,
        timeout=30,
        env=case_env(),
        check=False,
    )


@pytest.mark.parametrize("route", ROUTES)
def test_two_ranks_match_one(route, periodic_params, tmp_path):
    """One launched rank and two ranks report the same step as one process.

    The geometry is the two waters, which this deck converges. A scattered
    arrangement of the same types misses the orbital threshold in the
    deck's 40 iterations, so only the reference well can finish it, and
    that run never reaches the rank comparison.
    """
    if route == "fix":
        pytest.importorskip("mpi4py")
    case = spec(
        route,
        "real",
        periodic_params,
        setup=["processors * 1 1"],
        run=[
            "velocity all create 300 77 loop geom",
            "fix nve all nve",
            "timestep 0.25",
            "run 40",
        ],
    )
    one = run_case(case, tmp_path, None)
    launched = run_case(case, tmp_path, 1)
    two = run_case(case, tmp_path, 2)
    assert launched["world_size"] == 1, launched["world_size"]
    assert two["world_size"] == 2, (
        f"mpiexec -n 2 gave LAMMPS a world of {two['world_size']} rank(s); "
        "the LAMMPS library and the launcher use different MPI installations"
    )
    assert min(two["nlocal"]) > 0, two["nlocal"]
    assert launched["pe"] == pytest.approx(one["pe"], rel=1e-12)
    np.testing.assert_allclose(launched["forces"], one["forces"], rtol=1e-10, atol=1e-12)
    assert two["pe"] == pytest.approx(one["pe"], rel=1e-12)
    np.testing.assert_allclose(two["forces"], one["forces"], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(two["virial"], one["virial"], rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("route", ROUTES)
def test_two_rank_duplicate_id_error_does_not_deadlock(route, periodic_params, tmp_path):
    if route == "fix":
        pytest.importorskip("mpi4py")
    case = spec(
        route,
        "metal",
        periodic_params,
        setup=["processors * 1 1"],
        duplicate_nonroot_id=True,
        thermo=False,
    )
    done = run_case_failure(case, tmp_path, 2)
    assert done.returncode != 0
    assert "unique atom IDs" in done.stderr + done.stdout
