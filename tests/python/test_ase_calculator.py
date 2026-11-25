"""CPMDCCalculator against ASE's Calculator contract.

The default build evaluates the reference well: a harmonic well about the
origin, ``E = sum_i 0.5 k Z_i |r_i|^2`` plus small deck and cell terms.
It is exact to differentiate but not translation invariant.
With OpenCPMD linked the same cases run on an SCF, which converges the
energy to about 1e-7 eV, so comparisons that are exact on the reference
well take a tolerance there.
"""

from __future__ import annotations

import copy
import ctypes
import json
import os
import pickle
import re
import shutil
import struct
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import (
    CalculationFailed,
    Calculator,
    PropertyNotImplementedError,
    get_calculator_class,
    names as calculator_names,
)
from ase_fd import calculate_numerical_forces, calculate_numerical_stress
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from conftest import INSTALLED, ROOT, WATER, park_orbitals

import cpmdc
from cpmdc import (
    BOHR_TO_ANGSTROM,
    HARTREE_TO_EV,
    ConstantCutoff,
    CPMDCError,
    Potential,
    Session,
)
from cpmdc.ase import CPMDCCalculator

CUBE = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
ORTHORHOMBIC = [[10.0, 0.0, 0.0], [0.0, 11.0, 0.0], [0.0, 0.0, 12.0]]
TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def water(**kwargs) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """One entry per library evaluation by any calculator."""
    seen: list[int] = []
    evaluate = Potential.evaluate

    def counted(self, *args, **kwargs):
        seen.append(1)
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", counted)
    return seen


@pytest.fixture
def calc(params):
    with CPMDCCalculator(params) as calculator:
        yield calculator


@pytest.fixture
def periodic_calc(periodic_params):
    with CPMDCCalculator(periodic_params) as calculator:
        yield calculator


@pytest.fixture
def cwds(monkeypatch: pytest.MonkeyPatch) -> list[Path]:
    """The working directory of each library evaluation."""
    seen: list[Path] = []
    evaluate = Potential.evaluate

    def recorded(self, *args, **kwargs):
        seen.append(Path.cwd())
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", recorded)
    return seen


@pytest.fixture(scope="module")
def cluster_variant(tmp_path_factory: pytest.TempPathFactory):
    """Encode the cluster water message with one field changed."""
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("put capnp on PATH to encode a CPMDParams variant")
    source = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text(encoding="utf-8")

    def encode(old: str = "", new: str = "") -> Path:
        assert old in source
        done = subprocess.run(
            [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
            input=source.replace(old, new).encode("utf-8"),
            stdout=subprocess.PIPE,
            check=True,
        )
        target = tmp_path_factory.mktemp("variant") / "water.params.bin"
        target.write_bytes(done.stdout)
        return target

    return encode


def test_is_an_ase_calculator(calc, params):
    assert isinstance(calc, Calculator)
    assert calc.implemented_properties == ["energy", "free_energy", "forces", "stress", "dipole"]
    assert calc.default_parameters == {"params": None, "library": None, "constant_cutoff": None}
    assert calc.parameters["params"] == str(params)
    assert calc.todict() == {"params": str(params)}
    json.dumps(calc.todict())


def test_units_are_ev_and_angstrom(calc, params):
    atoms = water()
    atoms.calc = calc
    with Session(params) as session:
        native = session.energy_forces(WATER, [8, 1, 1])
    assert atoms.get_potential_energy() == pytest.approx(native.energy_h * HARTREE_TO_EV, rel=1e-14)
    scale = HARTREE_TO_EV / BOHR_TO_ANGSTROM
    expected = np.array(native.forces_h_bohr) * scale
    np.testing.assert_allclose(atoms.get_forces(), expected, rtol=1e-14)
    assert atoms.get_potential_energy(force_consistent=True) == atoms.get_potential_energy()
    assert atoms.get_forces().shape == (3, 3)
    assert atoms.get_forces().dtype == np.float64


def test_a_periodic_step_matches_session_and_potential(periodic_calc, periodic_params, opencpmd):
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.rattle(0.05, seed=19)
    atoms.calc = periodic_calc
    energy, forces, stress = atoms.get_potential_energy(), atoms.get_forces(), atoms.get_stress()
    numbers = [int(number) for number in atoms.numbers]
    with Session(periodic_params) as session:
        native = session.energy_forces(atoms.positions, numbers, atoms.cell.array)
    with Potential(periodic_params) as pot:
        provider = pot.evaluate(atoms.positions, numbers, atoms.cell.array)
    assert native.stress_h_bohr3 is not None
    tensor = -np.reshape(native.stress_h_bohr3, (3, 3)) * HARTREE_TO_EV / BOHR_TO_ANGSTROM**3
    voigt = [tensor[0, 0], tensor[1, 1], tensor[2, 2], tensor[1, 2], tensor[0, 2], tensor[0, 1]]
    scf = dict(rel=0, abs=1e-5) if opencpmd else dict(rel=1e-14)
    assert energy == pytest.approx(native.energy_h * HARTREE_TO_EV, **scf)
    assert energy == pytest.approx(provider["energy"], **scf)
    force_atol = 1e-3 if opencpmd else 0.0
    expected = np.array(native.forces_h_bohr) * HARTREE_TO_EV / BOHR_TO_ANGSTROM
    np.testing.assert_allclose(forces, expected, rtol=1e-14, atol=force_atol)
    np.testing.assert_allclose(forces, provider["forces"], rtol=1e-14, atol=force_atol)
    stress_atol = 1e-4 if opencpmd else 0.0
    np.testing.assert_allclose(stress, voigt, rtol=1e-12, atol=stress_atol)
    np.testing.assert_allclose(
        atoms.get_stress(voigt=False).ravel(), provider["stress"], rtol=1e-12, atol=stress_atol
    )


@pytest.fixture(scope="session")
def dipole_params(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Periodic water whose deck asks for DIPOLE DYNAMICS."""
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("put capnp on PATH to encode the dipole deck")
    source = (ROOT / "tests" / "water_blyp_periodic_params.capnp.txt").read_text(
        encoding="utf-8"
    )
    source = source.replace("maxIter = 40,", "maxIter = 40, dipoleDynamics = true,")
    encoded = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=source.encode("utf-8"),
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    params = tmp_path_factory.mktemp("dipole") / "water-dipole.params.bin"
    params.write_bytes(encoded)
    return params


def test_opencpmd_dipole_is_recomputed_and_reaches_ase(dipole_params, opencpmd):
    if not opencpmd:
        pytest.skip("the reference evaluator has no electronic dipole")
    params = dipole_params
    moved = np.asarray(WATER, dtype=float)
    with Session(params) as session:
        first = session.energy_forces(moved, [8, 1, 1], CUBE)
        moved[1, 0] += 0.05
        second = session.energy_forces(moved, [8, 1, 1], CUBE)
    assert first.dipole_e_bohr is not None
    assert second.dipole_e_bohr is not None
    assert not np.allclose(first.dipole_e_bohr, second.dipole_e_bohr)
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        # The independent session can converge to a slightly different point
        # within the deck's 1e-5 orbital-gradient threshold.
        np.testing.assert_allclose(
            atoms.get_dipole_moment(),
            np.asarray(first.dipole_e_bohr) * BOHR_TO_ANGSTROM,
            rtol=0,
            atol=1.0e-4,
        )


def test_opencpmd_dipole_value_is_pinned(dipole_params, opencpmd):
    """Real OpenCPMD water in a 10 Angstrom cube, BLYP, 70 Ry.

    The oxygen sits at +z and the hydrogens at -z, so the dipole, which points
    from the negative to the positive charge, points along -z. OpenCPMD gives
    (6.1e-6, 2.9e-6, -0.64866) e*Bohr. The range is wide enough for another
    BLAS and tight enough that a dipole scaled by one half or by a unit
    factor fails.
    """
    if not opencpmd:
        pytest.skip("the reference evaluator has no electronic dipole")
    flipped = np.asarray(WATER, dtype=float) * [1.0, 1.0, -1.0]
    with Session(dipole_params) as session:
        up = session.energy_forces(WATER, [8, 1, 1], CUBE)
        assert up.dipole_e_bohr is not None
        assert session.dipole_possible()
        x, y, z = up.dipole_e_bohr
        assert -0.68 < z < -0.62
        assert abs(x) < 1.0e-3
        assert abs(y) < 1.0e-3
        assert 0.62 < float(np.linalg.norm(up.dipole_e_bohr)) < 0.68

        down = session.energy_forces(flipped, [8, 1, 1], CUBE)
        assert down.dipole_e_bohr is not None
        assert 0.62 < down.dipole_e_bohr[2] < 0.68

        # calculate_au has no dipole output; the step's dipole is read after it.
        bohr = np.asarray(WATER, dtype=float).ravel() / BOHR_TO_ANGSTROM
        cell = np.asarray(CUBE, dtype=float).ravel() / BOHR_TO_ANGSTROM
        session.calculate_au(bohr, [8, 1, 1], cell)
        last = cpmdc_last_dipole_values(session)
    # A warm restart converges to within the deck's 1e-5 orbital gradient,
    # which moves the dipole by about 1e-4.
    np.testing.assert_allclose(last, up.dipole_e_bohr, rtol=0, atol=1.0e-3)


def cpmdc_last_dipole_values(session: Session) -> list[float]:
    """The dipole of the last step through the C ABI, electron Bohr."""
    from cpmdc import _CDipole

    out = _CDipole()
    assert session._lib.cpmdc_session_last_dipole(session._handle, ctypes.byref(out)) == 0
    assert out.valid
    return [float(out.values[i]) for i in range(3)]


def test_caps_report_a_dipole_only_when_the_evaluator_can_make_one(
    params, periodic_params, dipole_params, opencpmd
):
    with Potential(params) as pot:
        assert pot.caps()["dipole"] is opencpmd
    # A deck without a dipole request produces none, and caps then says so.
    with Potential(periodic_params) as pot:
        pot.evaluate(WATER, [8, 1, 1], CUBE)
        assert pot.caps() == {"stress": True, "dipole": False, "energies": False}
    if opencpmd:
        with Potential(dipole_params) as pot:
            assert pot.caps()["dipole"] is True
            assert "dipole" in pot.evaluate(WATER, [8, 1, 1], CUBE)
            assert pot.caps() == {"stress": True, "dipole": True, "energies": False}


@pytest.mark.parametrize("periodic", [False, True])
def test_forces_match_finite_differences(params, periodic_params, opencpmd, periodic):
    if periodic:
        atoms = water(cell=TRICLINIC, pbc=True)
        atoms.calc = CPMDCCalculator(periodic_params)
    else:
        atoms = water()
        atoms.calc = CPMDCCalculator(params)
    atoms.rattle(0.05, seed=3)
    forces = atoms.get_forces()
    if opencpmd:
        numeric = calculate_numerical_forces(atoms, eps=1.0e-3)
        np.testing.assert_allclose(forces, numeric, rtol=0, atol=2.0e-3)
    else:
        numeric = calculate_numerical_forces(atoms, eps=1.0e-4)
        np.testing.assert_allclose(forces, numeric, rtol=1.0e-6, atol=1.0e-9)
    atoms.calc.close()


def test_ase_gradient_test(calc, opencpmd):
    """ASE numerical forces, down to the smallest step the session can tell apart.

    ``gradient_test`` always returns its last sample, a 1e-8 Angstrom step.
    The session reports one result for geometries that agree to 1e-8
    Angstrom, so that sample is not a new geometry. ASE 3.23 updates the
    displaced coordinate in place and lands on that identity. Later ASE
    releases rebuild the coordinate from the original and can slip past it.
    The check uses ASE's own finite-difference helper from 0.1 Angstrom
    down to 1e-7 Angstrom.
    """
    if opencpmd:
        pytest.skip("steps below 1e-4 Angstrom are inside the SCF tolerance")
    atoms = water()
    atoms.rattle(0.05, seed=7)
    atoms.calc = calc
    analytic = atoms.get_forces()
    for eps in np.logspace(-1, -7, 7):
        numeric = calculate_numerical_forces(atoms, eps=float(eps))
        np.testing.assert_allclose(analytic, numeric, rtol=0, atol=1.0e-6)


@pytest.mark.parametrize("cell", [CUBE, TRICLINIC], ids=["cubic", "triclinic"])
def test_stress_matches_numerical_stress(periodic_calc, opencpmd, cell):
    atoms = water(cell=cell, pbc=True)
    atoms.rattle(0.05, seed=5)
    atoms.calc = periodic_calc
    stress = atoms.get_stress()
    assert stress.shape == (6,)
    numeric = calculate_numerical_stress(atoms, eps=1.0e-5)
    np.testing.assert_allclose(stress, numeric, rtol=1.0e-5, atol=1.0e-4 if opencpmd else 1.0e-10)
    full = atoms.get_stress(voigt=False)
    np.testing.assert_allclose(full, full.T)


def test_stress_has_the_ase_sign(periodic_calc, opencpmd):
    """sigma = (1/V) dE/d(strain): the reference energy rises on expansion."""
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    sxx = atoms.get_stress()[0]
    eps = 1.0e-5
    energies = []
    for sign in (1.0, -1.0):
        strained = atoms.copy()
        strain = np.eye(3)
        strain[0, 0] += sign * eps
        strained.set_cell(atoms.cell.array @ strain, scale_atoms=True)
        strained.calc = periodic_calc
        energies.append(strained.get_potential_energy())
    expected = (energies[0] - energies[1]) / (2.0 * eps * atoms.get_volume())
    assert sxx > 0.0
    assert sxx == pytest.approx(expected, rel=1e-2 if opencpmd else 1e-6)


def test_translation_moves_along_the_forces(calc, opencpmd):
    """A rigid shift changes the energy by the work of the mean force.

    The relation is exact for a quadratic energy and holds with both sides
    zero for a translation-invariant one.
    """
    atoms = water()
    atoms.calc = calc
    shift = np.array([0.03, -0.02, 0.05])
    e0, f0 = atoms.get_potential_energy(), atoms.get_forces()
    atoms.translate(shift)
    e1, f1 = atoms.get_potential_energy(), atoms.get_forces()
    assert e1 - e0 == pytest.approx(
        -0.5 * np.sum((f0 + f1) * shift), rel=1e-9, abs=1e-9 if opencpmd else 1e-14
    )


def test_reference_well_is_anchored_at_the_origin(calc, opencpmd):
    if opencpmd:
        pytest.skip("OpenCPMD is linked; the reference well is not")
    atoms = water()
    atoms.calc = calc
    shift = np.array([0.03, -0.02, 0.05])
    f0 = atoms.get_forces()
    atoms.translate(shift)
    per_charge = (atoms.get_forces() - f0) / atoms.numbers[:, None]
    np.testing.assert_allclose(per_charge, np.tile(per_charge[0], (3, 1)), rtol=1e-9)
    direction = per_charge[0] / np.linalg.norm(per_charge[0])
    np.testing.assert_allclose(direction, -shift / np.linalg.norm(shift), rtol=1e-9)


@pytest.mark.parametrize("periodic", [False, True])
def test_translation_invariance(params, periodic_params, opencpmd, periodic):
    if not opencpmd:
        pytest.skip("the reference well is anchored at the origin")
    if periodic:
        atoms = water(cell=CUBE, pbc=True)
        atoms.calc = CPMDCCalculator(periodic_params)
        shift = atoms.cell.array[0]
    else:
        atoms = water()
        atoms.calc = CPMDCCalculator(params)
        shift = np.array([0.3, -0.2, 0.5])
    e0, f0 = atoms.get_potential_energy(), atoms.get_forces()
    if not periodic:
        # cpmdc removes the net force of an isolated cell. A periodic deck
        # keeps the grid's, 0.09 eV/Angstrom along z for this water at 50 Ry.
        np.testing.assert_allclose(f0.sum(axis=0), 0.0, atol=1e-10)
    atoms.translate(shift)
    assert atoms.get_potential_energy() == pytest.approx(e0, abs=1e-5)
    np.testing.assert_allclose(atoms.get_forces(), f0, atol=1e-3)
    atoms.calc.close()


def test_unchanged_atoms_reuse_the_step(calc, calls):
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    atoms.get_forces()
    atoms.get_potential_energy(force_consistent=True)
    calc.get_property("energy", atoms.copy())
    assert len(calls) == 1


@pytest.mark.parametrize(
    "change",
    ["positions", "numbers", "cell", "pbc"],
)
def test_a_system_change_recomputes(periodic_calc, calls, opencpmd, change):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    first = atoms.get_potential_energy()
    if change == "positions":
        atoms.positions[1, 0] += 0.01
    elif change == "numbers":
        # The OpenCPMD deck has pseudopotentials for O and H only.
        atoms.numbers = [1, 8, 1] if opencpmd else [8, 1, 3]
    elif change == "cell":
        atoms.set_cell(np.array(CUBE) * 1.01)
    else:
        atoms.pbc = False
    assert periodic_calc.check_state(atoms) == [change]
    if opencpmd and change == "pbc":
        # Without a cell the periodic deck keeps its own 10 Angstrom cube,
        # which is CUBE, so OpenCPMD gives the same energy again.
        atoms.get_potential_energy()
    else:
        assert atoms.get_potential_energy() != first
    assert len(calls) == 2


def test_allow_calculation_false(calc, calls):
    atoms = water()
    atoms.calc = calc
    assert calc.get_property("energy", atoms, allow_calculation=False) is None
    assert calls == []
    energy = atoms.get_potential_energy()
    assert calc.get_property("energy", atoms, allow_calculation=False) == energy
    assert calc.get_property("stress", atoms, allow_calculation=False) is None
    moved = atoms.copy()
    moved.positions[0, 2] += 0.1
    assert calc.get_property("forces", moved, allow_calculation=False) is None
    assert len(calls) == 1


def test_calculation_required_follows_the_stored_step(periodic_calc, calls):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    every = ["energy", "free_energy", "forces", "stress"]
    assert periodic_calc.calculation_required(atoms, ["energy"])
    atoms.get_forces()
    assert not periodic_calc.calculation_required(atoms, every)
    for change in ("positions", "cell", "numbers"):
        moved = atoms.copy()
        if change == "positions":
            moved.positions[2, 1] += 0.01
        elif change == "cell":
            moved.set_cell(np.array(CUBE) * 0.99)
        else:
            moved.numbers = [1, 8, 1]
        assert periodic_calc.calculation_required(moved, ["energy"])
    assert len(calls) == 1


def test_reset_clears_the_step_and_keeps_the_session(calc, calls, opencpmd):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    session = calc.potential
    calc.reset()
    assert calc.results == {}
    assert calc.atoms is None
    assert calc.calculation_required(atoms, ["energy"])
    assert atoms.get_potential_energy() == pytest.approx(energy, rel=0, abs=1e-5 if opencpmd else 0)
    assert calc.potential is session
    assert len(calls) == 2


@pytest.mark.parametrize("name", ["energies", "magmom", "magmoms", "charges"])
def test_unsupported_properties_raise(calc, name):
    atoms = water()
    atoms.calc = calc
    with pytest.raises(PropertyNotImplementedError):
        calc.get_property(name, atoms)


def test_an_unavailable_dipole_raises_once(calc, calls):
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    for _ in range(2):
        with pytest.raises(PropertyNotImplementedError, match="computed no dipole"):
            atoms.get_dipole_moment()
    assert len(calls) == 1
    assert "dipole" not in calc.results


def test_a_valid_dipole_is_an_ase_vector(calc, monkeypatch):
    evaluate = Potential.evaluate

    def with_dipole(self, *args, **kwargs):
        result = evaluate(self, *args, **kwargs)
        result["dipole"] = [0.25, -0.5, 0.75]
        return result

    monkeypatch.setattr(Potential, "evaluate", with_dipole)
    atoms = water()
    atoms.calc = calc
    dipole = atoms.get_dipole_moment()
    np.testing.assert_array_equal(dipole, [0.25, -0.5, 0.75])
    assert dipole.shape == (3,)
    assert dipole.dtype == np.float64


def test_stress_of_a_cluster_raises_once(calc, calls):
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    for _ in range(2):
        with pytest.raises(PropertyNotImplementedError, match="fully periodic"):
            atoms.get_stress()
    assert len(calls) == 1


def test_stress_of_an_isolated_deck_raises_with_the_reason(calc, calls):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = calc
    for _ in range(2):
        with pytest.raises(PropertyNotImplementedError, match="isolated system"):
            atoms.get_stress()
    assert len(calls) == 1
    assert "energy" in calc.results


def test_mixed_boundaries_are_refused(calc):
    atoms = water(cell=CUBE, pbc=[True, True, False])
    atoms.calc = calc
    with pytest.raises(ValueError, match="non-periodic axes: z"):
        atoms.get_potential_energy()


@pytest.mark.parametrize(
    "cell",
    [None, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]],
    ids=["no cell", "coplanar"],
)
def test_periodic_flags_need_a_cell(periodic_calc, cell):
    atoms = water(cell=cell, pbc=True)
    atoms.calc = periodic_calc
    with pytest.raises(ValueError, match="three independent vectors"):
        atoms.get_potential_energy()


def test_a_cluster_ignores_its_cell(calc):
    boxed = water(cell=CUBE)
    boxed.calc = calc
    bare = water()
    bare.calc = copy.deepcopy(calc)
    assert boxed.get_potential_energy() == bare.get_potential_energy()


def test_set_returns_changes_and_resets(calc, params, periodic_params, calls):
    atoms = water()
    atoms.calc = calc
    first = atoms.get_potential_energy()
    assert calc.set(params=params) == {}
    assert "energy" in calc.results
    assert calc.set(params=str(periodic_params)) == {"params": str(periodic_params)}
    assert calc.results == {}
    assert atoms.get_potential_energy() != first
    blob = params.read_bytes()
    assert calc.set(params=blob) == {"params": blob}
    saved = calc.todict()
    assert saved["params"].keys() == {"cpmdparams_base64"}
    assert atoms.get_potential_energy() == first
    assert len(calls) == 3


def test_set_rejects_bad_input_and_keeps_the_session(calc, params):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    with pytest.raises(TypeError, match="unknown calculator parameters"):
        calc.set(cutoff=70.0)
    with pytest.raises(FileNotFoundError):
        calc.set(params=params.parent / "missing.bin")
    with pytest.raises(TypeError, match="params is required"):
        calc.set(params=None)
    assert calc.parameters["params"] == str(params)
    moved = atoms.copy()
    moved.calc = calc
    moved.positions[0, 0] += 0.01
    assert moved.get_potential_energy() != energy


@pytest.mark.parametrize("key, value", [("charge", -1), ("multiplicity", 3), ("xc", "PBE")])
def test_method_keywords_point_at_cpmdparams(params, key, value):
    message = rf"unknown calculator parameters: \['{key}'\].*CPMDParams"
    with pytest.raises(TypeError, match=message):
        CPMDCCalculator(params, **{key: value})
    with CPMDCCalculator(params) as calc:
        with pytest.raises(TypeError, match=message):
            calc.set(**{key: value})
        assert key not in calc.parameters


def test_params_are_required():
    with pytest.raises(TypeError, match="params is required"):
        CPMDCCalculator()


def test_keywords_go_through_calculator_init(params, tmp_path):
    with CPMDCCalculator(params=params) as calc:
        assert calc.parameters["params"] == str(params)
        saved = tmp_path / "calc.ase"
        calc.parameters.write(saved)
    with CPMDCCalculator(parameters=str(saved)) as loaded:
        assert loaded.parameters == calc.parameters
    atoms = water()
    with CPMDCCalculator(params, atoms=atoms) as attached:
        assert atoms.calc is attached
        assert np.isfinite(atoms.get_potential_energy())


def test_todict_rebuilds_an_equal_calculator(calc, opencpmd):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    with CPMDCCalculator(**calc.todict()) as rebuilt:
        assert rebuilt.todict() == calc.todict()
        other = water()
        other.calc = rebuilt
        assert other.get_potential_energy() == pytest.approx(
            energy, rel=0, abs=1e-5 if opencpmd else 0
        )


@pytest.mark.parametrize(
    "old, new",
    [("charge = 0,", "charge = 2,"), ("multiplicity = 1,", "multiplicity = 3,")],
    ids=["charge", "multiplicity"],
)
def test_charge_and_multiplicity_come_from_cpmdparams(cluster_variant, opencpmd, old, new):
    if opencpmd:
        pytest.skip("a charged or open-shell water is a separate OpenCPMD SCF")
    plain = cluster_variant()
    atoms = water()
    with CPMDCCalculator(plain) as calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        assert calc.set(params=cluster_variant(old, new)) != {}
        assert calc.results == {}
        assert atoms.get_potential_energy() != energy


def test_steps_run_in_the_calculator_directory(params, tmp_path, monkeypatch, cwds, opencpmd):
    monkeypatch.chdir(tmp_path)
    here = Path.cwd()
    atoms = water()
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy()
    moved = water()
    with CPMDCCalculator(params, directory="run") as calc:
        moved.calc = calc
        assert moved.get_potential_energy() == pytest.approx(
            energy, rel=0, abs=1e-5 if opencpmd else 0
        )
    assert cwds == [here, here / "run"]
    assert Path.cwd() == here


def test_label_names_the_directory(params, tmp_path, monkeypatch, cwds):
    monkeypatch.chdir(tmp_path)
    atoms = water()
    with CPMDCCalculator(params, label="run/water") as calc:
        assert (calc.directory, calc.prefix) == ("run", "water")
        atoms.calc = calc
        atoms.get_potential_energy()
    assert cwds == [Path.cwd() / "run"]


def test_restart_is_ignored(params):
    with pytest.warns(DeprecationWarning, match="restart is ignored"):
        calc = CPMDCCalculator(params, restart="old/run")
    assert calc.directory == "."
    assert calc.label is None
    calc.close()


def test_a_failed_step_raises_calculation_failed(params, tmp_path, monkeypatch):
    evaluate = Potential.evaluate
    broken = [True]

    def flaky(self, *args, **kwargs):
        if broken[0]:
            raise CPMDCError("SCF did not converge")
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", flaky)
    monkeypatch.chdir(tmp_path)
    here = Path.cwd()
    atoms = water()
    with CPMDCCalculator(params, directory="run") as calc:
        atoms.calc = calc
        with pytest.raises(CalculationFailed, match="SCF did not converge") as failed:
            atoms.get_potential_energy()
        assert isinstance(failed.value, CPMDCError)
        assert str(here / "run") in str(failed.value)
        assert Path.cwd() == here
        assert calc.results == {}
        broken[0] = False
        assert np.isfinite(atoms.get_potential_energy())


def test_a_step_without_atoms_raises(calc):
    with pytest.raises(ValueError, match="no atoms"):
        calc.get_potential_energy()


@pytest.mark.parametrize("how", ["copy", "deepcopy", "pickle"])
def test_copies_keep_results_and_own_a_session(periodic_calc, opencpmd, how):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    energy, stress = atoms.get_potential_energy(), atoms.get_stress()
    if how == "copy":
        clone = copy.copy(periodic_calc)
    elif how == "deepcopy":
        clone = copy.deepcopy(periodic_calc)
    else:
        clone = pickle.loads(pickle.dumps(periodic_calc))
    assert clone.parameters == periodic_calc.parameters
    assert clone.results.keys() == periodic_calc.results.keys()
    assert clone.get_property("energy", atoms, allow_calculation=False) == energy
    periodic_calc.close()
    moved = atoms.copy()
    moved.positions[0, 0] += 0.01
    moved.calc = clone
    assert moved.get_potential_energy() != energy
    np.testing.assert_allclose(clone.get_stress(atoms), stress, atol=1e-5 if opencpmd else 0)
    assert periodic_calc.results.get("energy") == energy
    clone.close()


def test_close_then_reuse(calc, opencpmd):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    calc.close()
    calc.close()
    atoms.positions[0, 0] += 0.01
    atoms.get_potential_energy()
    atoms.positions[0, 0] -= 0.01
    assert atoms.get_potential_energy() == pytest.approx(energy, rel=0, abs=1e-5 if opencpmd else 0)


def test_a_new_composition_opens_a_new_session(calc, opencpmd):
    atoms = water()
    atoms.calc = calc
    first = atoms.get_potential_energy()
    session = calc.potential
    reordered = Atoms("HOH", positions=[WATER[1], WATER[0], WATER[2]])
    reordered.calc = calc
    assert reordered.get_potential_energy() == pytest.approx(first, rel=1e-14)
    assert calc.potential is session
    smaller = Atoms("OH", positions=WATER[:2])
    smaller.calc = calc
    smaller.get_potential_energy()
    assert calc.potential is not session
    assert atoms.get_potential_energy() == pytest.approx(first, rel=0, abs=1e-5 if opencpmd else 0)


@pytest.mark.parametrize("name", ["energy", "forces"])
def test_a_refused_step_records_its_composition(calc, opencpmd, monkeypatch, name):
    """A non-finite step still leaves its composition on the session it ran in."""
    evaluate = Potential.evaluate

    def refuse_two_atoms(self, positions, numbers, box=None):
        out = evaluate(self, positions, numbers, box)
        if len(numbers) == 2:
            if name == "energy":
                out["energy"] = float("nan")
            else:
                out["forces"][0][0] = float("inf")
        return out

    monkeypatch.setattr(Potential, "evaluate", refuse_two_atoms)
    atoms = water()
    atoms.calc = calc
    first = atoms.get_potential_energy()
    session = calc.potential
    smaller = Atoms("OH", positions=WATER[:2])
    smaller.calc = calc
    with pytest.raises(CalculationFailed, match=f"non-finite {name}"):
        smaller.get_potential_energy()
    refused = calc.potential
    assert refused is not session
    assert calc._composition == np.sort(smaller.numbers).tobytes()
    atoms.positions[0, 0] += 0.01
    atoms.get_potential_energy()
    atoms.positions[0, 0] -= 0.01
    assert calc.potential is not refused
    assert calc._composition == np.sort(atoms.numbers).tobytes()
    assert atoms.get_potential_energy() == pytest.approx(first, rel=0, abs=1e-5 if opencpmd else 0)


def test_deprecated_arguments_warn(params):
    with pytest.warns(DeprecationWarning, match="schema"):
        Session(params, schema="Potentials.capnp").close()
    with pytest.warns(DeprecationWarning, match="schema"):
        Potential(params, schema="Potentials.capnp").close()
    with pytest.warns(DeprecationWarning, match="schema"):
        calc = CPMDCCalculator(params, schema="Potentials.capnp")
    with pytest.warns(DeprecationWarning, match="schema"):
        assert calc.set(schema="Potentials.capnp") == {}
    with pytest.warns(DeprecationWarning, match="restart"):
        CPMDCCalculator(params, restart="old").close()
    calc.close()


def test_bfgs_writes_a_trajectory(calc, params, tmp_path):
    atoms = water()
    atoms.rattle(0.1, seed=11)
    atoms.calc = calc
    path = tmp_path / "opt.traj"
    with BFGS(atoms, trajectory=str(path), logfile=None) as opt:
        converged = opt.run(fmax=0.01, steps=200)
    assert converged
    frames = read(path, index=":")
    assert len(frames) == opt.nsteps + 1
    assert frames[-1].get_potential_energy() < frames[0].get_potential_energy()
    assert np.linalg.norm(frames[-1].get_forces(), axis=1).max() < 0.01
    np.testing.assert_allclose(frames[-1].positions, atoms.positions)
    with Trajectory(path) as traj:
        first = traj[0]
    assert first.calc.name == "cpmdccalculator"
    assert first.calc.parameters["params"] == str(params)


def test_bfgs_trajectory_stores_periodic_stress(periodic_calc, tmp_path):
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.1, seed=13)
    atoms.calc = periodic_calc
    atoms.get_stress()
    path = tmp_path / "periodic.traj"
    with BFGS(atoms, trajectory=str(path), logfile=None) as opt:
        assert opt.run(fmax=0.01, steps=200)
    last = read(path, index=-1)
    np.testing.assert_allclose(last.get_stress(), atoms.get_stress())


def relax_cell(atoms: Atoms, steps: int) -> tuple[bool, list[float]]:
    energies: list[float] = []
    with BFGS(FrechetCellFilter(atoms), logfile=None) as opt:
        opt.attach(lambda: energies.append(atoms.get_potential_energy()))
        converged = opt.run(fmax=0.01, steps=steps)
    return converged, energies


def test_bfgs_moves_the_cell_through_a_cell_filter(params, periodic_params, opencpmd):
    """At a constant cutoff the cell relaxes to a sensible volume.

    The 40 Ry deck keeps its basis and evaluates the 30 Ry constant-cutoff
    functional. One water molecule in a 10 Angstrom box has almost no
    intermolecular stress, so on OpenCPMD the cell converges within 10% of
    its start while the energy falls. The first steps lower it every time.
    """
    if opencpmd:
        runner = ROOT / "tests" / "python" / "run_bfgs_cell_filter.py"
        env = os.environ.copy()
        env["CPMDC_TEST_PARAMS"] = str(params)
        env["CPMDC_TEST_PERIODIC_PARAMS"] = str(periodic_params)
        subprocess.run([sys.executable, str(runner)], check=True, env=env)
        return
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=17)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        converged, energies = relax_cell(atoms, steps=4)
    assert energies[-1] < energies[0]
    assert np.all(np.diff(energies[:5]) < 0.0)
    assert not np.allclose(atoms.cell.array, CUBE)


def test_opencpmd_in_process_bfgs_converges_after_parking(
    params, periodic_params, opencpmd
):
    """Cell BFGS after parking reaches the force threshold.

    Parking replaces the orbitals OpenCPMD keeps for the process. The
    rattled water at a constant 30 Ry cutoff then converges. Live runs
    took 51 optimizer steps to ``fmax`` 0.01, and every energy and force
    on that trajectory was finite, so the step cap is 51.
    """
    if not opencpmd:
        pytest.skip("the reference well has no retained OpenCPMD orbitals")
    park_orbitals(params)
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=17)
    energies: list[float] = []
    forces: list[np.ndarray] = []

    def record() -> None:
        energies.append(float(atoms.get_potential_energy()))
        forces.append(np.array(atoms.get_forces(), dtype=float, copy=True))

    # Two identical live OpenCPMD runs, one BLAS thread, seed 17.
    step_cap = 51
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        with BFGS(FrechetCellFilter(atoms), logfile=None) as opt:
            opt.attach(record)
            converged = opt.run(fmax=0.01, steps=step_cap)
            nsteps = opt.nsteps
    assert energies and forces
    flat_forces = np.concatenate([sample.ravel() for sample in forces])
    assert converged, f"nsteps={nsteps}"
    assert nsteps <= step_cap
    assert np.all(np.isfinite(energies))
    assert np.all(np.isfinite(flat_forces))


def test_a_fixed_basis_cell_filter_collapses_the_cell(periodic_calc, opencpmd):
    """Without a constant cutoff the OpenCPMD stress is Pulay stress.

    The 40 Ry stress is about 16 eV over the volume. The optimizer shrinks
    the cell along it while the energy at the fixed cutoff rises.
    """
    if not opencpmd:
        pytest.skip("the reference well has no basis and no Pulay stress")
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=17)
    atoms.calc = periodic_calc
    assert atoms.get_stress()[:3].mean() > 1.0e-2
    _, energies = relax_cell(atoms, steps=4)
    assert atoms.get_volume() < 0.75 * np.linalg.det(CUBE)
    assert energies[-1] > energies[1]


@pytest.mark.parametrize("constant_cutoff", [None, 30.0], ids=["fixed-basis", "constant-cutoff"])
def test_stress_follows_the_energy_over_a_finite_strain(
    periodic_params, opencpmd, constant_cutoff
):
    """The mean stress against dE/dV over a 1.5% isotropic strain.

    With a fixed basis the OpenCPMD stress is the derivative at a fixed
    number of plane waves, about 1.6e-2 eV/Angstrom^3 above the slope of the
    energy at the fixed cutoff. The constant-cutoff stress follows it.
    """
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(periodic_params, constant_cutoff=constant_cutoff) as calc:
        atoms.calc = calc
        pressure = atoms.get_stress()[:3].mean()
        points = []
        for scale in (0.985, 1.015):
            strained = atoms.copy()
            strained.set_cell(atoms.cell.array * scale, scale_atoms=True)
            strained.calc = calc
            points.append((strained.get_volume(), strained.get_potential_energy()))
    (v_low, e_low), (v_high, e_high) = points
    slope = (e_high - e_low) / (v_high - v_low)
    if not opencpmd:
        assert pressure == pytest.approx(slope, rel=1.0e-3, abs=1.0e-12)
    elif constant_cutoff is None:
        assert pressure - slope > 1.0e-2
    else:
        assert abs(pressure - slope) < 5.0e-4


def strained_steps(calc, scales, cell=CUBE):
    """Energy, forces, and stress of the water with the cell scaled by each factor, in turn."""
    steps = []
    for scale in scales:
        atoms = water(cell=cell, pbc=True)
        atoms.set_cell(np.multiply(cell, scale), scale_atoms=True)
        atoms.calc = calc
        steps.append((atoms.get_volume(), atoms.get_potential_energy(), atoms.get_forces(), atoms.get_stress()))
    return steps


def warm_starts(err: str) -> list[int]:
    return [int(s) for s in re.findall(r"^cpmdc_warm_start\s+(-?\d+)", err, re.M)]


def test_the_constant_cutoff_stress_is_the_volume_derivative_of_the_energy(
    periodic_params, opencpmd, monkeypatch, capfd
):
    """-integral P dV over a 2% stretch against the energy change.

    CPMDC_CELL_IN_PLACE keeps one set of 40 Ry plane waves through the
    stretched cells, so the three energies lie on one smooth surface and
    Simpson's rule on the three pressures integrates its slope. The
    fixed-basis stress of this deck is about 1.6e-2 eV/Angstrom^3 off that
    slope, close to 1 eV over the stretch.
    """
    if not opencpmd:
        pytest.skip("the reference well keeps no plane waves between calls")
    monkeypatch.setenv("CPMDC_CELL_IN_PLACE", "0.05")
    monkeypatch.setenv("CPMDC_SCF_STEPS", "1")
    scales = (1.0, 1.01, 1.02)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        steps = strained_steps(calc, scales)
    assert warm_starts(capfd.readouterr().err)[1:] == [0, 0]
    volume0 = steps[0][0]
    pressures = [-np.mean(stress[:3]) for _, _, _, stress in steps]
    dvds = [3.0 * scale**2 * volume0 for scale in scales]
    work = (scales[2] - scales[0]) / 6.0 * (
        pressures[0] * dvds[0] + 4.0 * pressures[1] * dvds[1] + pressures[2] * dvds[2]
    )
    assert steps[2][1] - steps[0][1] == pytest.approx(-work, rel=1.0e-2, abs=1.0e-3)


def test_a_constant_cutoff_cell_moves_in_place_within_the_strain_limit(
    periodic_params, opencpmd, monkeypatch, capfd
):
    """CPMDC_CELL_IN_PLACE against a new setup for every cell.

    The 1% stretch and the return to the setup cell move in place and start
    from the stored orbitals; the 3% stretch is 3% beyond the setup cell, so
    it sets CPMD up again. In place, the plane waves between the 36 Ry
    sphere and the 40 Ry cutoff are those of the setup cell, which on this
    deck moves the 1% step by about 3 meV, 1.5 meV/Angstrom, and
    3e-6 eV/Angstrom^3 from a new setup. The process keeps one CPMD, so the
    pass without in-place moves runs first and leaves the 3% cell set up.
    """
    if not opencpmd:
        pytest.skip("the reference well keeps no plane waves between calls")
    monkeypatch.setenv("CPMDC_SCF_STEPS", "1")
    scales = (1.0, 1.01, 1.0, 1.03)
    runs = {}
    for limit in ("0", "0.02"):
        monkeypatch.setenv("CPMDC_CELL_IN_PLACE", limit)
        with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
            steps = strained_steps(calc, scales)
        err = capfd.readouterr().err
        iterations = [int(n) for n in re.findall(r"^cpmdc_scf_steps\s+(\d+)", err, re.M)]
        runs[limit] = (steps, warm_starts(err), iterations)
    assert runs["0"][1][1:] == [-2, -2, -2]
    assert runs["0.02"][1] == [-2, 0, 0, -2]
    assert runs["0.02"][2][1] < runs["0"][2][1]
    back, start = runs["0.02"][0][2], runs["0.02"][0][0]
    assert back[1] == pytest.approx(start[1], rel=0, abs=1.0e-4)
    np.testing.assert_allclose(back[3], start[3], rtol=0, atol=1.0e-5)
    for moved, fresh in zip(runs["0.02"][0], runs["0"][0]):
        assert moved[1] == pytest.approx(fresh[1], rel=0, abs=1.0e-2)
        np.testing.assert_allclose(moved[2], fresh[2], rtol=0, atol=5.0e-3)
        np.testing.assert_allclose(moved[3], fresh[3], rtol=0, atol=2.0e-5)


def test_a_triclinic_cell_moves_in_place_after_an_orthorhombic_one(
    periodic_params, opencpmd, monkeypatch, capfd
):
    """A CELL VECTORS setup folds through its own cell after a SYMMETRY 8 one.

    setsc writes the minimum-image half lengths of a CELL deck and none for
    CELL VECTORS, so the lengths of the orthorhombic box must not outlive
    its setup. With them, pbc3 folds the triclinic cell as a monoclinic one
    of the old lengths and the in-place check refuses the 1% stretch.
    """
    if not opencpmd:
        pytest.skip("the reference well keeps no plane waves between calls")
    monkeypatch.setenv("CPMDC_SCF_STEPS", "1")
    monkeypatch.setenv("CPMDC_CELL_IN_PLACE", "0.02")
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms = water(cell=ORTHORHOMBIC, pbc=True)
        atoms.calc = calc
        atoms.get_potential_energy()
        strained_steps(calc, (1.0, 1.01), cell=TRICLINIC)
    assert warm_starts(capfd.readouterr().err)[1:] == [-2, 0]


def test_constant_cutoff_is_a_calculator_parameter(periodic_params, calls):
    with CPMDCCalculator(periodic_params, constant_cutoff=30) as calc:
        assert calc.parameters["constant_cutoff"] == 30.0
        assert calc.todict() == {"params": str(periodic_params), "constant_cutoff": 30.0}
        json.dumps(calc.todict())
        atoms = water(cell=CUBE, pbc=True)
        atoms.calc = calc
        atoms.get_stress()
        first = calc.potential
        explicit = {"cutoff_ry": 30.0, "height_ry": 2000.0, "width_ry": 3.0}
        assert calc.set(constant_cutoff=ConstantCutoff(**explicit)) == {
            "constant_cutoff": explicit
        }
        assert calc.potential is not first
        assert calc.results == {}
        json.dumps(calc.todict())
        with CPMDCCalculator(**calc.todict()) as rebuilt:
            assert rebuilt.parameters == calc.parameters
        assert calc.set(constant_cutoff=None) == {"constant_cutoff": None}
        assert calc.todict() == {"params": str(periodic_params)}
    assert len(calls) == 1


def test_a_constant_cutoff_the_deck_cannot_take_raises(periodic_params):
    with pytest.raises(CPMDCError, match="exceeds the plane-wave cutoff 40 Ry"):
        CPMDCCalculator(periodic_params, constant_cutoff=38.0)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        first = calc.potential
        with pytest.raises(CPMDCError, match="exceeds the plane-wave cutoff"):
            calc.set(constant_cutoff={"cutoff_ry": 30.0, "width_ry": 6.0})
        assert calc.parameters["constant_cutoff"] == 30.0
        assert calc.potential is first


def test_one_periodic_flag_covers_every_axis():
    from cpmdc.ase import _pbc_flags

    assert _pbc_flags(True) == [True, True, True]
    assert _pbc_flags([False]) == [False, False, False]
    with pytest.raises(ValueError, match="one flag or three"):
        _pbc_flags([True, False])


def test_library_and_bytearray_params_round_trip(params):
    blob = Path(params).read_bytes()
    library = str(cpmdc.find_library())
    with CPMDCCalculator(bytearray(blob), library=library) as calc:
        assert calc.parameters["params"] == blob
        assert calc.parameters["library"] == library
        atoms = water()
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        calc.calculate(atoms, properties=[])
        assert calc.results["energy"] == energy
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            calc.set(schema=None)
        calc.parameters["params"] = None
        with pytest.raises(TypeError, match="params is required"):
            calc.set(library=str(Path(library).parent / "missing.so"))
        assert calc.parameters["library"] == library
        assert atoms.get_potential_energy() == energy


@pytest.mark.parametrize(
    ("value", "error", "message"),
    [
        ({"base64": "QQ=="}, TypeError, "cpmdparams_base64"),
        ({"cpmdparams_base64": 1}, TypeError, "cpmdparams_base64"),
        ({"cpmdparams_base64": "***"}, ValueError, "valid standard base64"),
    ],
)
def test_todict_params_mapping_is_validated(params, value, error, message):
    with CPMDCCalculator(params) as calc:
        with pytest.raises(error, match=message):
            calc.set(params=value)


@pytest.mark.parametrize(
    "value, error",
    [(-1.0, ValueError), (float("nan"), ValueError), (True, TypeError), ("30", TypeError)],
)
def test_constant_cutoff_values_are_checked(value, error):
    with pytest.raises(error):
        ConstantCutoff(value)
    with pytest.raises(error):
        ConstantCutoff(30.0, width_ry=value)
    assert ConstantCutoff.coerce(30) == ConstantCutoff(30.0)
    assert ConstantCutoff.coerce({"cutoff_ry": 30, "width_ry": 2}).todict() == {
        "cutoff_ry": 30.0,
        "width_ry": 2.0,
    }
    assert ConstantCutoff.coerce(None) is None


def test_get_calculator_class_returns_cpmdc_calculator():
    assert get_calculator_class("cpmdc") is CPMDCCalculator
    assert "cpmdc" in calculator_names


def _import_env(*extra: Path) -> dict[str, str]:
    """Environment for a child process that imports this checkout's ``cpmdc``."""
    env = os.environ.copy()
    paths = [str(path) for path in extra]
    if not INSTALLED:
        paths.append(str(ROOT / "python"))
    if env.get("PYTHONPATH"):
        paths.append(env["PYTHONPATH"])
    if paths:
        env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def test_repeated_registration_is_idempotent():
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import importlib, cpmdc.ase\n"
            "from ase.calculators.calculator import get_calculator_class, names\n"
            "for _ in range(2):\n"
            "    importlib.reload(cpmdc.ase)\n"
            "assert get_calculator_class('cpmdc') is cpmdc.ase.CPMDCCalculator\n"
            "assert names.count('cpmdc') == 1\n",
        ],
        env=_import_env(),
        check=True,
    )


def test_global_plugin_collection_lists_cpmdc_after_cpmdc_ase_import(tmp_path):
    pytest.importorskip("ase._4.plugins", reason="ASE v4 plugin API is not present")
    extra = []
    if not INSTALLED:
        dist = tmp_path / "cpmdc-0.0.0.dist-info"
        dist.mkdir()
        (dist / "METADATA").write_text("Metadata-Version: 2.1\nName: cpmdc\nVersion: 0.0.0\n")
        (dist / "entry_points.txt").write_text("[ase.plugins]\ncpmdc = cpmdc.ase\n")
        extra.append(tmp_path)
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import cpmdc.ase\n"
            "from ase._4.plugins import plugins\n"
            "by_name = {plugin.name: plugin for plugin in plugins.calculators}\n"
            "assert by_name['cpmdc']() is cpmdc.ase.CPMDCCalculator, sorted(by_name)\n",
        ],
        env=_import_env(*extra),
        check=True,
    )


@pytest.mark.parametrize("periodic", [False, True], ids=["cluster", "periodic"])
def test_calculator_built_via_lookup_matches_direct(params, periodic_params, opencpmd, periodic):
    cls = get_calculator_class("cpmdc")
    deck = periodic_params if periodic else params
    atoms = water(cell=TRICLINIC, pbc=True) if periodic else water()
    other = atoms.copy()
    with cls(deck) as looked_up, CPMDCCalculator(deck) as direct:
        atoms.calc = looked_up
        other.calc = direct
        energy, forces = atoms.get_potential_energy(), atoms.get_forces()
        ref_energy, ref_forces = other.get_potential_energy(), other.get_forces()
        scf = dict(rel=0, abs=1e-5) if opencpmd else dict(rel=1e-14)
        force_atol = 1e-3 if opencpmd else 0.0
        assert energy == pytest.approx(ref_energy, **scf)
        np.testing.assert_allclose(forces, ref_forces, rtol=1e-14, atol=force_atol)
        if periodic:
            stress_atol = 1e-4 if opencpmd else 0.0
            np.testing.assert_allclose(
                atoms.get_stress(), other.get_stress(), rtol=1e-12, atol=stress_atol
            )


def test_ase_plugin_resolves_cpmdc_calculator(params, opencpmd):
    pytest.importorskip("ase._4.plugins", reason="ASE v4 plugin API is not present")
    from ase._4.plugins import PluginCollection

    import cpmdc.ase as ase_mod

    plugins = PluginCollection(ase_mod.__ase_plugins__).calculators
    by_name = {plugin.name: plugin for plugin in plugins}
    assert "cpmdc" in by_name
    cls = by_name["cpmdc"]()
    assert cls is ase_mod.CPMDCCalculator
    atoms = water()
    other = water()
    with cls(params) as via_plugin, ase_mod.CPMDCCalculator(params) as direct:
        atoms.calc = via_plugin
        other.calc = direct
        scf = dict(rel=0, abs=1e-5) if opencpmd else dict(rel=1e-14)
        force_atol = 1e-3 if opencpmd else 0.0
        assert atoms.get_potential_energy() == pytest.approx(
            other.get_potential_energy(), **scf
        )
        np.testing.assert_allclose(
            atoms.get_forces(), other.get_forces(), rtol=1e-14, atol=force_atol
        )


CHARGE_TWO = ("charge = 0,", "charge = 2,")
MULTIPLICITY_THREE = ("multiplicity = 1,", "multiplicity = 3,")


def tagged_water(charges=None, magmoms=None) -> Atoms:
    atoms = water()
    if charges is not None:
        atoms.set_initial_charges(charges)
    if magmoms is not None:
        atoms.set_initial_magnetic_moments(magmoms)
    return atoms


def energy_of(blob: Path, atoms: Atoms) -> float:
    with CPMDCCalculator(blob) as calculator:
        atoms.calc = calculator
        return atoms.get_potential_energy()


def test_all_zero_atoms_leave_the_method_to_cpmdparams(cluster_variant, calls):
    plain = cluster_variant()
    with CPMDCCalculator(plain) as calc:
        bare = water()
        bare.calc = calc
        energy = bare.get_potential_energy()
        zeros = tagged_water([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        zeros.calc = calc
        assert zeros.get_potential_energy() == energy
        assert calc._override == (None, None)
    assert energy_of(plain, tagged_water(np.zeros(3), np.zeros(3))) == energy
    # A cancelling set of partial charges is a total charge of zero.
    assert energy_of(plain, tagged_water([0.5, -0.25, -0.25])) == energy
    # Zero moments of the non-collinear shape are not a moment either.
    assert energy_of(plain, tagged_water(None, np.zeros((3, 3)))) == energy
    charged = cluster_variant(*CHARGE_TWO)
    assert energy_of(charged, tagged_water(np.zeros(3), np.zeros(3))) == energy_of(
        charged, water()
    )


@pytest.mark.parametrize(
    "variant, charges, magmoms",
    [
        (CHARGE_TWO, [2.0, 0.0, 0.0], None),
        (CHARGE_TWO, [0.5, 1.0, 0.5], None),
        (MULTIPLICITY_THREE, None, [2.0, 0.0, 0.0]),
        (MULTIPLICITY_THREE, None, [-1.0, -1.0, 0.0]),
    ],
    ids=["charge", "charge-spread", "multiplicity", "multiplicity-negative"],
)
def test_atoms_set_the_charge_and_multiplicity(
    cluster_variant, opencpmd, variant, charges, magmoms
):
    if opencpmd:
        pytest.skip("a charged or open-shell water is a separate OpenCPMD SCF")
    plain = cluster_variant()
    expected = energy_of(cluster_variant(*variant), water())
    assert expected != energy_of(plain, water())
    assert energy_of(plain, tagged_water(charges, magmoms)) == expected


@pytest.mark.parametrize("variant", [CHARGE_TWO, MULTIPLICITY_THREE], ids=["charge", "multiplicity"])
def test_the_edited_message_is_the_message_capnp_encodes(cluster_variant, variant):
    from cpmdc.ase import _params_with_method

    plain = cluster_variant().read_bytes()
    encoded = cluster_variant(*variant).read_bytes()
    charge, multiplicity = (2, None) if variant is CHARGE_TWO else (None, 3)
    assert _params_with_method(plain, charge, multiplicity) == encoded
    assert _params_with_method(plain, None, None) is plain


def test_a_message_without_the_scalar_words_is_refused():
    from cpmdc.ase import _params_with_method

    root = (1 << 32) | 0
    short = struct.pack("<II", 0, 2) + struct.pack("<Q", root) + bytes(8)
    with pytest.raises(CPMDCError, match="data bytes"):
        _params_with_method(short, 1, None)
    with pytest.raises(CPMDCError, match="segment"):
        _params_with_method(b"\0\0\0\0", 1, None)
    far = struct.pack("<II", 0, 1) + struct.pack("<Q", 2)
    with pytest.raises(CPMDCError, match="struct pointer"):
        _params_with_method(far, 1, None)
    with pytest.raises(ValueError, match="Int32"):
        _params_with_method(short, 2**31, None)


def test_agreeing_atoms_and_cpmdparams_are_silent(cluster_variant, opencpmd):
    if opencpmd:
        pytest.skip("a charged or open-shell water is a separate OpenCPMD SCF")
    both = cluster_variant(
        "charge = 0,\n  multiplicity = 1,", "charge = 2,\n  multiplicity = 3,"
    )
    expected = energy_of(both, water())
    assert energy_of(both, tagged_water([2.0, 0.0, 0.0], [2.0, 0.0, 0.0])) == expected
    assert energy_of(both, tagged_water([2.0, 0.0, 0.0])) == expected
    assert energy_of(both, tagged_water(None, [-2.0, 0.0, 0.0])) == expected


def test_a_charge_that_disagrees_with_cpmdparams_raises(cluster_variant, calls):
    charged = cluster_variant(*CHARGE_TWO)
    with CPMDCCalculator(charged) as calc:
        atoms = tagged_water([1.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError) as caught:
            atoms.get_potential_energy()
        message = str(caught.value)
        assert "charge 1" in message and "Atoms.initial_charges" in message
        assert "charge 2" in message and "CPMDParams" in message
        assert calc.results == {}
        with pytest.raises(ValueError, match="disagrees"):
            atoms.get_potential_energy()
        assert calls == []


def test_a_multiplicity_that_disagrees_with_cpmdparams_raises(cluster_variant, calls):
    triplet = cluster_variant(*MULTIPLICITY_THREE)
    with CPMDCCalculator(triplet) as calc:
        atoms = tagged_water(None, [4.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError) as caught:
            atoms.get_potential_energy()
        message = str(caught.value)
        assert "multiplicity 5" in message and "Atoms.initial_magmoms" in message
        assert "multiplicity 3" in message and "CPMDParams" in message
        assert calls == []


def test_explicit_zero_charge_with_nonzero_initial_charges_keeps_zero(
    cluster_variant, calls, opencpmd
):
    explicit = cluster_variant("charge = 0,", "charge = 0,\n  hasCharge = true,")
    from cpmdc.ase import _top_level_has_charge

    assert _top_level_has_charge(explicit.read_bytes()) is True
    assert cpmdc.params_method(explicit).charge == 0
    neutral = energy_of(explicit, water())
    tol = dict(rel=0, abs=1e-5 if opencpmd else 0)
    assert energy_of(explicit, tagged_water([0.0, 0.0, 0.0])) == pytest.approx(neutral, **tol)
    before = len(calls)
    with CPMDCCalculator(explicit) as calc:
        atoms = tagged_water([1.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError, match="disagrees") as caught:
            atoms.get_potential_energy()
        message = str(caught.value)
        assert "charge 1" in message and "Atoms.initial_charges" in message
        assert "charge 0" in message and "CPMDParams" in message
        assert len(calls) == before
        assert calc.results == {}
        assert calc._override == (None, None)
    assert energy_of(explicit, water()) == pytest.approx(neutral, **tol)


def test_an_explicit_nonzero_charge_is_kept(cluster_variant, calls, opencpmd):
    charged = cluster_variant(*CHARGE_TWO)
    with CPMDCCalculator(charged) as calc:
        atoms = tagged_water([1.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError, match="disagrees") as caught:
            atoms.get_potential_energy()
        message = str(caught.value)
        assert "charge 1" in message and "charge 2" in message
        assert calls == []
        assert calc._override == (None, None)
    if opencpmd:
        pytest.skip("a charged water is a separate OpenCPMD SCF")
    expected = energy_of(charged, water())
    assert energy_of(charged, tagged_water([0.0, 0.0, 0.0])) == expected
    assert energy_of(charged, tagged_water([2.0, 0.0, 0.0])) == expected
    with CPMDCCalculator(charged) as calc:
        assert calc._method_override(tagged_water([2.0, 0.0, 0.0])) == (None, None)
        assert calc._method_override(tagged_water([0.0, 0.0, 0.0])) == (None, None)


def test_unset_charge_takes_the_sum_of_initial_charges(cluster_variant, opencpmd):
    from cpmdc.ase import _top_level_has_charge

    plain = cluster_variant()
    assert _top_level_has_charge(plain.read_bytes()) is False
    assert cpmdc.params_method(plain).charge == 0
    with CPMDCCalculator(plain) as calc:
        assert calc._method_override(tagged_water([2.0, 0.0, 0.0])) == (2, None)
        assert calc._method_override(tagged_water([0.5, 1.0, 0.5])) == (2, None)
    if opencpmd:
        pytest.skip("a charged water is a separate OpenCPMD SCF")
    expected = energy_of(cluster_variant(*CHARGE_TWO), water())
    assert energy_of(plain, tagged_water([2.0, 0.0, 0.0])) == expected


def test_all_zero_initial_charges_leave_the_charge_neutral(cluster_variant, opencpmd):
    plain = cluster_variant()
    explicit = cluster_variant("charge = 0,", "charge = 0,\n  hasCharge = true,")
    tol = dict(rel=0, abs=1e-5 if opencpmd else 0)
    neutral = energy_of(plain, water())
    assert energy_of(plain, tagged_water([0.0, 0.0, 0.0])) == pytest.approx(neutral, **tol)
    assert energy_of(explicit, water()) == pytest.approx(neutral, **tol)
    assert energy_of(explicit, tagged_water([0.0, 0.0, 0.0])) == pytest.approx(neutral, **tol)
    # Partial charges that cancel are still a total of zero.
    assert energy_of(explicit, tagged_water([0.5, -0.25, -0.25])) == pytest.approx(
        neutral, **tol
    )
    with CPMDCCalculator(plain) as calc, CPMDCCalculator(explicit) as marked:
        assert calc._method_override(tagged_water(np.zeros(3))) == (None, None)
        assert marked._method_override(tagged_water(np.zeros(3))) == (None, None)


def test_a_typed_system_explicit_zero_is_kept(cluster_variant, calls, opencpmd):
    pinned = cluster_variant(
        "symmetry = 0,",
        "symmetry = 0, charge = 0, hasCharge = true,",
    )
    assert cpmdc.params_method(pinned).charge == 0
    tol = dict(rel=0, abs=1e-5 if opencpmd else 0)
    neutral = energy_of(cluster_variant(), water())
    assert energy_of(pinned, water()) == pytest.approx(neutral, **tol)
    before = len(calls)
    with CPMDCCalculator(pinned) as calc:
        atoms = tagged_water([2.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError, match="fixes it") as caught:
            atoms.get_potential_energy()
        message = str(caught.value)
        assert "charge 2" in message and "charge 0" in message
        assert len(calls) == before
        assert calc.results == {}
    assert energy_of(pinned, tagged_water([0.0, 0.0, 0.0])) == pytest.approx(neutral, **tol)


def test_a_charge_that_a_typed_section_pins_is_refused(tmp_path, calls):
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("put capnp on PATH to encode a CPMDParams variant")
    text = (ROOT / "tests" / "params_system_has_charge_zero.capnp.txt").read_bytes()
    blob = tmp_path / "pinned.bin"
    blob.write_bytes(
        subprocess.run(
            [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
            input=text,
            stdout=subprocess.PIPE,
            check=True,
        ).stdout
    )
    assert cpmdc.params_method(blob).charge == 0
    with CPMDCCalculator(blob) as calc:
        atoms = tagged_water([2.0, 0.0, 0.0])
        atoms.calc = calc
        with pytest.raises(ValueError, match="fixes it"):
            atoms.get_potential_energy()
        assert calls == []


@pytest.mark.parametrize(
    "charges, magmoms, message",
    [
        ([0.5, 0.0, 0.0], None, "Atoms.initial_charges sums to 0.5"),
        ([0.001, 0.0, 0.0], None, "not an integer"),
        ([0.1] * 3, None, "not an integer"),
        (None, [0.5, 0.0, 0.0], "Atoms.initial_magmoms sums to 0.5"),
        (None, [1.0, 0.25, 0.0], "not an integer"),
        ([float("nan"), 0.0, 0.0], None, "non-finite"),
        (None, [float("inf"), 0.0, 0.0], "non-finite"),
    ],
)
def test_a_non_integral_total_raises(cluster_variant, calls, charges, magmoms, message):
    with CPMDCCalculator(cluster_variant()) as calc:
        atoms = tagged_water(charges, magmoms)
        atoms.calc = calc
        with pytest.raises(ValueError, match=message):
            atoms.get_potential_energy()
        assert calls == []
        assert calc.results == {}


def test_a_sum_within_float64_rounding_of_an_integer_is_accepted(cluster_variant, opencpmd):
    if opencpmd:
        pytest.skip("a charged water is a separate OpenCPMD SCF")
    plain = cluster_variant()
    expected = energy_of(cluster_variant(*CHARGE_TWO), water())
    third = [2.0 / 3.0] * 3
    assert energy_of(plain, tagged_water(third)) == expected
    tenths = Atoms("H10", positions=np.arange(30.0).reshape(10, 3))
    tenths.set_initial_charges([0.1] * 10)
    with CPMDCCalculator(plain) as calc:
        assert calc._method_override(tenths) == (1, None)


def test_non_collinear_moments_raise(cluster_variant, calls):
    with CPMDCCalculator(cluster_variant()) as calc:
        atoms = tagged_water(None, [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        atoms.calc = calc
        with pytest.raises(ValueError, match="non-collinear"):
            atoms.get_potential_energy()
        assert calls == []


def test_a_change_of_initial_charges_evaluates_again(cluster_variant, calls, opencpmd):
    with CPMDCCalculator(cluster_variant()) as calc:
        atoms = water()
        atoms.calc = calc
        neutral = atoms.get_potential_energy()
        session = calc.potential
        assert len(calls) == 1
        atoms.set_initial_charges([2.0, 0.0, 0.0])
        charged = atoms.get_potential_energy()
        assert len(calls) == 2
        assert calc.potential is not session
        if not opencpmd:
            assert charged != neutral
        reopened = calc.potential
        atoms.set_initial_charges([0.0, 2.0, 0.0])
        assert calc.check_state(atoms) == ["initial_charges"]
        assert atoms.get_potential_energy() == charged
        assert len(calls) == 3
        assert calc.potential is reopened
        atoms.set_initial_charges([0.0, 0.0, 0.0])
        assert atoms.get_potential_energy() == pytest.approx(
            neutral, rel=0, abs=1e-5 if opencpmd else 0
        )
        assert len(calls) == 4
        assert calc.potential is not reopened


def test_a_change_of_params_resets_the_atoms_override(cluster_variant, calls, opencpmd):
    if opencpmd:
        pytest.skip("a charged water is a separate OpenCPMD SCF")
    plain = cluster_variant()
    charged = cluster_variant(*CHARGE_TWO)
    with CPMDCCalculator(plain) as calc:
        atoms = tagged_water([2.0, 0.0, 0.0])
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        assert calc.set(params=charged) != {}
        assert calc._override == (None, None)
        assert atoms.get_potential_energy() == energy
        assert calc._override == (None, None)
        bad = tagged_water([1.0, 0.0, 0.0])
        bad.calc = calc
        with pytest.raises(ValueError, match="disagrees"):
            bad.get_potential_energy()


def test_the_atoms_override_survives_todict_and_pickle(cluster_variant, opencpmd):
    if opencpmd:
        pytest.skip("a charged water is a separate OpenCPMD SCF")
    plain = cluster_variant()
    expected = energy_of(cluster_variant(*CHARGE_TWO), water())
    with CPMDCCalculator(plain) as calc:
        atoms = tagged_water([2.0, 0.0, 0.0])
        atoms.calc = calc
        assert atoms.get_potential_energy() == expected
        saved = calc.todict()
        assert set(saved) == {"params"}
        restored = pickle.loads(pickle.dumps(calc))
        copied = copy.deepcopy(calc)
    with restored, copied, CPMDCCalculator(**saved) as rebuilt:
        assert restored.todict() == saved
        assert restored.results["energy"] == expected
        for other in (restored, copied, rebuilt):
            again = tagged_water([2.0, 0.0, 0.0])
            again.calc = other
            assert again.get_potential_energy() == expected
            neutral = water()
            neutral.calc = other
            assert neutral.get_potential_energy() != expected
