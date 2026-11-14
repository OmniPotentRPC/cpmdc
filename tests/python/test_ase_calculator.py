"""CPMDCCalculator against ASE's Calculator contract.

The default build evaluates the reference well: a harmonic well about the
origin, ``E = sum_i 0.5 k Z_i |r_i|^2`` plus small deck and cell terms.
It is exact to differentiate but not translation invariant.
"""

from __future__ import annotations

import copy
import json
import pickle

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, PropertyNotImplementedError
from ase.calculators.fd import calculate_numerical_forces, calculate_numerical_stress
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from conftest import WATER

from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Potential, Session
from cpmdc.ase import CPMDCCalculator

CUBE = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
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


def test_is_an_ase_calculator(calc, params):
    assert isinstance(calc, Calculator)
    assert calc.implemented_properties == ["energy", "free_energy", "forces", "stress"]
    assert calc.default_parameters == {"params": None, "library": None}
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


@pytest.mark.parametrize("periodic", [False, True])
def test_forces_match_finite_differences(params, periodic_params, periodic):
    if periodic:
        atoms = water(cell=TRICLINIC, pbc=True)
        atoms.calc = CPMDCCalculator(periodic_params)
    else:
        atoms = water()
        atoms.calc = CPMDCCalculator(params)
    atoms.rattle(0.05, seed=3)
    forces = atoms.get_forces()
    numeric = calculate_numerical_forces(atoms, eps=1.0e-4)
    np.testing.assert_allclose(forces, numeric, rtol=1.0e-6, atol=1.0e-9)
    atoms.calc.close()


@pytest.mark.parametrize("cell", [CUBE, TRICLINIC], ids=["cubic", "triclinic"])
def test_stress_matches_numerical_stress(periodic_calc, cell):
    atoms = water(cell=cell, pbc=True)
    atoms.rattle(0.05, seed=5)
    atoms.calc = periodic_calc
    stress = atoms.get_stress()
    assert stress.shape == (6,)
    numeric = calculate_numerical_stress(atoms, eps=1.0e-5)
    np.testing.assert_allclose(stress, numeric, rtol=1.0e-5, atol=1.0e-10)
    full = atoms.get_stress(voigt=False)
    np.testing.assert_allclose(full, full.T)


def test_stress_has_the_ase_sign(periodic_calc):
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
    assert sxx == pytest.approx(expected, rel=1e-6)


def test_translation_moves_along_the_forces(calc):
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
    assert e1 - e0 == pytest.approx(-0.5 * np.sum((f0 + f1) * shift), rel=1e-9, abs=1e-14)


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
    np.testing.assert_allclose(f0.sum(axis=0), 0.0, atol=1e-3)
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
    atoms.set_initial_charges([0.5, 0.0, 0.0])
    atoms.set_initial_magnetic_moments([1.0, 0.0, 0.0])
    atoms.get_forces()
    assert len(calls) == 1


@pytest.mark.parametrize(
    "change",
    ["positions", "numbers", "cell", "pbc"],
)
def test_a_system_change_recomputes(periodic_calc, calls, change):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    first = atoms.get_potential_energy()
    if change == "positions":
        atoms.positions[1, 0] += 0.01
    elif change == "numbers":
        atoms.numbers = [8, 1, 3]
    elif change == "cell":
        atoms.set_cell(np.array(CUBE) * 1.01)
    else:
        atoms.pbc = False
    assert periodic_calc.check_state(atoms) == [change]
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


@pytest.mark.parametrize("name", ["energies", "magmom", "magmoms", "charges", "dipole"])
def test_unsupported_properties_raise(calc, name):
    atoms = water()
    atoms.calc = calc
    with pytest.raises(PropertyNotImplementedError):
        calc.get_property(name, atoms)


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
    assert calc.todict() == {}
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


def test_unknown_constructor_arguments_raise(params):
    with pytest.raises(TypeError, match="unknown calculator arguments"):
        CPMDCCalculator(params, xc="PBE")


@pytest.mark.parametrize("how", ["copy", "deepcopy", "pickle"])
def test_copies_keep_results_and_own_a_session(periodic_calc, how):
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
    np.testing.assert_allclose(clone.get_stress(atoms), stress)
    assert periodic_calc.results.get("energy") == energy
    clone.close()


def test_close_then_reuse(calc):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    calc.close()
    calc.close()
    atoms.positions[0, 0] += 0.01
    atoms.get_potential_energy()
    atoms.positions[0, 0] -= 0.01
    assert atoms.get_potential_energy() == energy


def test_a_new_composition_opens_a_new_session(calc):
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
    assert atoms.get_potential_energy() == first


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
