"""CPMDCCalculator against ASE's own tools, on the default build.

``test_ase_calculator.py`` covers what the calculator adds to ASE. This
file asks the questions an ASE user asks of any calculator, and answers
them with ASE's machinery: ``check_state``, ``system_changes``,
``get_property``, ``ase.io``, ``ase.units``, the finite-difference tools,
``BFGS`` with ``FixAtoms``, and ``copy`` and ``pickle``. The default build
evaluates the reference well, so the numbers are exact; with OpenCPMD
linked the comparisons take the SCF tolerance.
"""

from __future__ import annotations

import copy
import os
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.calculator import (
    CalculationFailed,
    CalculatorError,
    PropertyNotImplementedError,
    all_changes,
    all_properties,
)
from ase.calculators.fd import calculate_numerical_forces, calculate_numerical_stress
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.stress import voigt_6_to_full_3x3_stress
from ase.units import create_units
from conftest import WATER

import cpmdc
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, CPMDCError, Potential, Session
from cpmdc.ase import CPMDCCalculator
from cpmdc.ase import CalculationFailed as CPMDCCalculationFailed

CUBE = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]
STEP = {"energy", "free_energy", "forces"}


def water(**kwargs: Any) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


@pytest.fixture
def calc(params):
    with CPMDCCalculator(params) as calculator:
        yield calculator


@pytest.fixture
def periodic_calc(periodic_params):
    with CPMDCCalculator(periodic_params) as calculator:
        yield calculator


@pytest.fixture
def opened(monkeypatch: pytest.MonkeyPatch) -> list[Session]:
    """Every Session that any calculator opens."""
    seen: list[Session] = []
    init = Session.__init__

    def counted(self, *args, **kwargs):
        init(self, *args, **kwargs)
        seen.append(self)

    monkeypatch.setattr(Session, "__init__", counted)
    return seen


@pytest.fixture
def steps(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """One entry per evaluation by the library."""
    seen: list[int] = []
    evaluate = Potential.evaluate

    def counted(self, *args, **kwargs):
        seen.append(1)
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", counted)
    return seen


@pytest.fixture
def tolerance(opencpmd) -> float:
    return 1e-5 if opencpmd else 0.0


# --- implemented_properties -------------------------------------------------


def test_implemented_properties_are_known_to_ase(calc):
    assert set(calc.implemented_properties) <= set(all_properties)
    assert len(set(calc.implemented_properties)) == len(calc.implemented_properties)


def test_a_cluster_returns_exactly_energy_free_energy_and_forces(calc):
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    assert set(calc.results) == STEP
    for name in STEP:
        assert name in calc.implemented_properties


def test_a_periodic_step_returns_exactly_what_it_lists(periodic_calc):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    for name in STEP | {"stress"}:
        assert name in periodic_calc.implemented_properties
        assert periodic_calc.get_property(name, atoms) is not None
    assert set(periodic_calc.results) <= set(periodic_calc.implemented_properties)
    assert set(periodic_calc.results) == STEP | {"stress"}


def test_every_returned_value_has_the_ase_shape(periodic_calc):
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.calc = periodic_calc
    atoms.get_stress()
    results = periodic_calc.results
    assert isinstance(results["energy"], float)
    assert isinstance(results["free_energy"], float)
    assert results["forces"].shape == (3, 3)
    assert results["stress"].shape == (6,)
    assert all(np.asarray(value).dtype == np.float64 for value in results.values())
    assert all(np.all(np.isfinite(value)) for value in results.values())


@pytest.mark.parametrize("name", ["energies", "stresses", "charges", "magmom", "magmoms"])
def test_get_property_refuses_what_the_backend_does_not_return(calc, steps, name):
    atoms = water()
    atoms.calc = calc
    with pytest.raises(PropertyNotImplementedError):
        calc.get_property(name, atoms)
    assert steps == []
    assert name not in calc.implemented_properties


def test_the_atoms_api_refuses_unsupported_properties(calc, steps):
    atoms = water()
    atoms.calc = calc
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_charges()
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_magnetic_moment()
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_potential_energies()
    with pytest.raises(PropertyNotImplementedError):
        atoms.get_stresses()
    assert steps == []


def test_stress_on_free_boundaries_is_refused_not_invented(calc):
    atoms = water()
    atoms.calc = calc
    with pytest.raises(PropertyNotImplementedError, match="fully periodic"):
        atoms.get_stress()
    assert "stress" not in calc.results


def test_extxyz_stores_only_what_the_calculator_returned(
    calc, periodic_calc, tmp_path: Path, tolerance
):
    cluster = water()
    cluster.calc = calc
    energy = cluster.get_potential_energy()
    write(tmp_path / "cluster.extxyz", cluster)
    back = read(tmp_path / "cluster.extxyz")
    assert set(back.calc.results) == STEP
    assert back.get_potential_energy() == pytest.approx(energy, rel=1e-8)

    periodic = water(cell=TRICLINIC, pbc=True)
    periodic.calc = periodic_calc
    periodic.get_stress()
    write(tmp_path / "periodic.extxyz", periodic)
    back = read(tmp_path / "periodic.extxyz")
    assert set(back.calc.results) == STEP | {"stress"}
    np.testing.assert_allclose(back.get_stress(), periodic.get_stress(), rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(back.get_forces(), periodic.get_forces(), rtol=1e-6, atol=1e-8)


# --- caching and invalidation ----------------------------------------------


def changed(atoms: Atoms, what: str, opencpmd: bool = False) -> Atoms:
    other = atoms.copy()
    if what == "positions":
        other.positions[1, 0] += 0.01
    elif what == "numbers":
        other.numbers = [1, 8, 1] if opencpmd else [8, 1, 3]
    elif what == "cell":
        other.set_cell(np.array(CUBE) * 1.01)
    elif what == "pbc":
        other.pbc = False
    else:
        raise AssertionError(what)
    return other


@pytest.mark.parametrize("what", ["positions", "numbers", "cell", "pbc"])
def test_check_state_names_the_change(periodic_calc, opencpmd, what):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    atoms.get_potential_energy()
    assert periodic_calc.check_state(atoms) == []
    assert periodic_calc.check_state(changed(atoms, what, opencpmd)) == [what]


@pytest.mark.parametrize("what", ["positions", "numbers", "cell", "pbc"])
def test_calculate_receives_the_system_changes(periodic_calc, monkeypatch, opencpmd, what):
    seen: list[list[str]] = []
    calculate = CPMDCCalculator.calculate

    def spy(self, atoms=None, properties=None, system_changes=all_changes):
        seen.append(list(system_changes))
        return calculate(self, atoms, properties, system_changes)

    monkeypatch.setattr(CPMDCCalculator, "calculate", spy)
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    atoms.get_potential_energy()
    other = changed(atoms, what, opencpmd)
    other.calc = periodic_calc
    other.get_potential_energy()
    assert seen == [list(all_changes), [what]]


@pytest.mark.parametrize("what", ["positions", "numbers", "cell", "pbc"])
def test_a_change_discards_the_stored_step(periodic_calc, steps, opencpmd, what):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    atoms.get_stress()
    before = {key: np.copy(value) for key, value in periodic_calc.results.items()}
    other = changed(atoms, what, opencpmd)
    assert periodic_calc.get_property("forces", other, allow_calculation=False) is None
    assert periodic_calc.results == {}
    assert periodic_calc.calculation_required(other, ["energy"])
    other.calc = periodic_calc
    other.get_forces()
    assert len(steps) == 2
    if what == "pbc":
        # The periodic deck holds the cube that the structure had. A step
        # without a cell is that same system, so the energy can repeat and
        # the stress of the replaced step is what must be gone.
        assert "stress" not in periodic_calc.results
    else:
        assert periodic_calc.results["energy"] != before["energy"]


def test_a_change_of_one_part_of_a_step_is_never_served_stale(periodic_calc, steps):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    forces, stress = atoms.get_forces(), atoms.get_stress()
    atoms.set_cell(np.array(CUBE) * 1.05, scale_atoms=True)
    assert not np.allclose(atoms.get_stress(), stress, rtol=0, atol=1e-12)
    assert not np.allclose(atoms.get_forces(), forces, rtol=0, atol=1e-12)
    assert len(steps) == 2


def test_a_tiny_change_is_still_a_change(calc, steps):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    atoms.positions[0, 2] += 1e-10
    assert calc.check_state(atoms) == ["positions"]
    atoms.get_forces()
    assert len(steps) == 2


def test_an_equal_structure_reuses_the_step(calc, steps):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    rebuilt = Atoms("OHH", positions=np.array(WATER).copy())
    assert calc.check_state(rebuilt) == []
    assert calc.get_property("energy", rebuilt) == energy
    assert calc.get_property("forces", rebuilt, allow_calculation=False) is not None
    assert len(steps) == 1


def test_initial_charges_and_moments_do_not_invalidate(calc, steps):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    atoms.set_initial_charges([1.0, -0.5, -0.5])
    atoms.set_initial_magnetic_moments([2.0, 0.0, 0.0])
    assert calc.check_state(atoms) == []
    atoms.get_forces()
    assert len(steps) == 1


def test_returned_arrays_do_not_alias_the_stored_step(calc):
    atoms = water()
    atoms.calc = calc
    forces = atoms.get_forces()
    stored = calc.results["forces"].copy()
    forces[:] = 123.0
    np.testing.assert_array_equal(atoms.get_forces(), stored)
    np.testing.assert_array_equal(calc.results["forces"], stored)


def test_the_calculator_keeps_its_own_copy_of_the_atoms(calc, steps):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    assert calc.atoms is not atoms
    atoms.positions[0, 0] += 0.1
    assert calc.check_state(atoms) == ["positions"]
    calc.atoms.positions[0, 0] += 0.1
    assert calc.check_state(atoms) == []


def test_a_direct_step_replaces_every_stored_result(periodic_calc):
    periodic = water(cell=CUBE, pbc=True)
    periodic.calc = periodic_calc
    periodic.get_stress()
    assert "stress" in periodic_calc.results
    cluster = water()
    periodic_calc.calculate(cluster, ["energy"], ["cell", "pbc"])
    assert set(periodic_calc.results) == STEP
    cluster.calc = periodic_calc
    with pytest.raises(PropertyNotImplementedError):
        cluster.get_stress()


def test_a_stress_returned_for_a_cluster_is_not_stored(calc, monkeypatch):
    evaluate = Potential.evaluate

    def with_stress(self, *args, **kwargs):
        done = evaluate(self, *args, **kwargs)
        done["stress"] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        return done

    monkeypatch.setattr(Potential, "evaluate", with_stress)
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    assert set(calc.results) == STEP
    with pytest.raises(PropertyNotImplementedError, match="fully periodic"):
        atoms.get_stress()


def test_a_failed_direct_step_leaves_no_stale_results(calc, monkeypatch):
    atoms = water()
    atoms.calc = calc
    first = atoms.get_potential_energy()
    moved = water()
    moved.positions[0, 2] += 0.2
    evaluate = Potential.evaluate
    broken = [True]

    def flaky(self, *args, **kwargs):
        if broken[0]:
            raise CPMDCError("SCF did not converge")
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", flaky)
    with pytest.raises(CalculationFailed):
        calc.calculate(moved, ["energy"], ["positions"])
    assert calc.results == {}
    assert calc.get_property("energy", moved, allow_calculation=False) is None
    broken[0] = False
    moved.calc = calc
    assert moved.get_potential_energy() != first


# --- one session for a run of steps ----------------------------------------


def test_repeated_calculate_calls_reuse_the_session(opened, steps, calc):
    atoms = water()
    for shift in range(5):
        atoms.positions[0, 2] = WATER[0][2] + 0.01 * shift
        calc.calculate(atoms, ["energy"], ["positions"])
    session = calc.potential
    calc.calculate(atoms, ["forces"], all_changes)
    assert calc.potential is session
    assert len(opened) == 1
    assert len(steps) == 6


def test_repeated_property_requests_reuse_the_session(opened, steps, periodic_calc):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    for shift in range(4):
        atoms.positions[1, 1] = WATER[1][1] + 0.02 * shift
        atoms.get_forces()
        atoms.get_stress()
        atoms.get_potential_energy()
    assert len(opened) == 1
    assert len(steps) == 4


def test_the_session_opens_with_the_parameters_and_survives_reset(opened, params):
    with CPMDCCalculator(params) as calc:
        assert len(opened) == 1
        atoms = water()
        atoms.calc = calc
        atoms.get_forces()
        calc.reset()
        atoms.get_forces()
        assert len(opened) == 1


def test_setting_the_same_parameters_keeps_the_session(opened, calc, params):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    session = calc.potential
    before = len(opened)
    assert calc.set(params=str(params)) == {}
    assert calc.set(library=None, constant_cutoff=None) == {}
    assert calc.potential is session
    assert len(opened) == before
    atoms.get_potential_energy(force_consistent=True)
    assert len(opened) == before


def test_a_composition_change_opens_one_session_and_stays_on_it(opened, calc):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    other = Atoms("OHHH", positions=WATER + [[0.0, 0.0, 2.0]])
    other.calc = calc
    for shift in range(3):
        other.positions[3, 2] = 2.0 + 0.1 * shift
        other.get_forces()
    assert len(opened) == 2


# --- failures ----------------------------------------------------------------


def test_a_real_error_of_the_library_is_a_calculation_failure(calc, tmp_path, monkeypatch):
    """The session already holds two atoms, so the library refuses three."""
    calc.potential.evaluate(WATER[:2], [8, 1])
    monkeypatch.chdir(tmp_path)
    here = Path.cwd()
    atoms = water()
    atoms.calc = calc
    with pytest.raises(CalculationFailed, match="topology change requires a new session") as err:
        atoms.get_potential_energy()
    assert isinstance(err.value, CPMDCCalculationFailed)
    assert isinstance(err.value, CalculatorError)
    assert isinstance(err.value, CPMDCError)
    assert isinstance(err.value.__cause__, CPMDCError)
    assert str(here) in str(err.value)
    assert Path.cwd() == here
    assert calc.results == {}


def test_a_failed_step_is_not_cached(calc, steps):
    calc.potential.evaluate(WATER[:2], [8, 1])
    steps.clear()
    atoms = water()
    atoms.calc = calc
    for _ in range(2):
        with pytest.raises(CalculationFailed):
            atoms.get_forces()
    assert len(steps) == 2
    assert calc.results == {}


def test_the_calculator_recovers_after_a_failure(calc, params, tolerance):
    calc.potential.evaluate(WATER[:2], [8, 1])
    atoms = water()
    atoms.calc = calc
    with pytest.raises(CalculationFailed):
        atoms.get_potential_energy()
    calc.close()
    with Session(params) as session:
        native = session.energy_forces(WATER, [8, 1, 1])
    assert atoms.get_potential_energy() == pytest.approx(
        native.energy_h * HARTREE_TO_EV, rel=0, abs=tolerance
    )


@pytest.mark.parametrize("scale", [1e200, -1e200])
def test_a_non_finite_result_is_a_calculation_failure(calc, scale, opencpmd, cold_start):
    """The reference well returns an infinite energy for this geometry.

    OpenCPMD returns no result. A cold step writes the coordinates into
    the deck, CPMD reads asterisks and stops, and the library gives its
    own reason, which the ``cold_start`` fixture makes the one to expect.
    ASE fixes the exception, not the text.
    """
    atoms = water()
    atoms.positions *= scale
    atoms.calc = calc
    reason = "CPMD stopgm during embed SCF" if opencpmd else "non-finite"
    where = re.escape(os.path.abspath(calc.directory))
    pattern = f"^cpmdc could not evaluate the step in {where}: .*{reason}"
    with pytest.raises(CalculationFailed, match=pattern):
        atoms.get_potential_energy()
    assert calc.results == {}


@pytest.mark.parametrize("name", ["energy", "forces", "stress"])
def test_a_non_finite_value_from_the_library_is_a_calculation_failure(
    periodic_params, monkeypatch, name
):
    def poisoned(self, *args, **kwargs):
        done = {"energy": 0.0, "forces": [0.0] * 9, "stress": [0.0] * 9}
        done[name] = float("nan") if name == "energy" else [float("inf")] * 9
        return done

    monkeypatch.setattr(Potential, "evaluate", poisoned)
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        with pytest.raises(CalculationFailed, match=f"non-finite {name}"):
            atoms.get_stress()
        assert calc.results == {}


def test_a_failure_in_a_relaxation_reaches_the_caller(calc, monkeypatch):
    evaluate = Potential.evaluate
    done = [0]

    def dies_on_the_third_step(self, *args, **kwargs):
        done[0] += 1
        if done[0] == 3:
            raise CPMDCError("SCF did not converge")
        return evaluate(self, *args, **kwargs)

    monkeypatch.setattr(Potential, "evaluate", dies_on_the_third_step)
    atoms = water()
    atoms.rattle(0.1, seed=3)
    atoms.calc = calc
    with pytest.raises(CalculationFailed, match="SCF did not converge"):
        BFGS(atoms, logfile=None).run(fmax=1e-6, steps=50)


def test_the_error_class_is_ases_and_cpmdcs(calc):
    assert issubclass(CPMDCCalculationFailed, CalculationFailed)
    assert issubclass(CPMDCCalculationFailed, CPMDCError)


# --- relaxation with a constraint --------------------------------------------


def test_bfgs_with_fixatoms_holds_the_fixed_atom(opened, steps, calc):
    atoms = water()
    atoms.rattle(0.1, seed=7)
    atoms.set_constraint(FixAtoms([0]))
    fixed = atoms.positions[0].copy()
    free = atoms.positions[1:].copy()
    atoms.calc = calc
    start = atoms.get_potential_energy()
    with BFGS(atoms, logfile=None) as opt:
        opt.run(fmax=1e-12, steps=5)
    assert opt.nsteps == 5
    np.testing.assert_array_equal(atoms.positions[0], fixed)
    assert not np.allclose(atoms.positions[1:], free)
    assert atoms.get_potential_energy() < start
    assert len(opened) == 1


def test_the_constraint_is_applied_by_ase_and_not_by_the_calculator(calc):
    atoms = water()
    atoms.rattle(0.1, seed=9)
    atoms.set_constraint(FixAtoms([0]))
    atoms.calc = calc
    constrained = atoms.get_forces()
    raw = atoms.get_forces(apply_constraint=False)
    np.testing.assert_array_equal(constrained[0], 0.0)
    np.testing.assert_array_equal(constrained[1:], raw[1:])
    assert np.abs(raw[0]).max() > 0.0
    np.testing.assert_array_equal(calc.results["forces"], raw)


def test_bfgs_with_fixatoms_converges_on_the_free_atoms(calc):
    atoms = water()
    atoms.rattle(0.1, seed=11)
    atoms.set_constraint(FixAtoms([0]))
    fixed = atoms.positions[0].copy()
    atoms.calc = calc
    with BFGS(atoms, logfile=None) as opt:
        assert opt.run(fmax=0.01, steps=200)
    np.testing.assert_array_equal(atoms.positions[0], fixed)
    assert np.linalg.norm(atoms.get_forces(), axis=1).max() < 0.01


def test_a_cell_relaxation_keeps_a_constraint(periodic_params, steps):
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=5)
    atoms.set_constraint(FixAtoms([1]))
    fixed = atoms.positions[1].copy()
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        with BFGS(atoms, logfile=None) as opt:
            opt.run(fmax=1e-12, steps=3)
    np.testing.assert_array_equal(atoms.positions[1], fixed)
    assert len(steps) >= 3


# --- copies and pickles --------------------------------------------------------


def test_atoms_with_a_calculator_copy_and_pickle(periodic_calc, tolerance):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    energy = atoms.get_potential_energy()
    for clone in (copy.deepcopy(atoms), pickle.loads(pickle.dumps(atoms))):
        assert clone.calc is not periodic_calc
        assert clone.get_potential_energy() == energy
        clone.positions[0, 2] += 0.01
        assert clone.get_potential_energy() != energy
        clone.calc.close()
    assert atoms.get_potential_energy() == energy


def test_a_copy_does_not_need_the_params_file(params, tmp_path, tolerance):
    path = tmp_path / "water.params.bin"
    path.write_bytes(Path(params).read_bytes())
    atoms = water()
    with CPMDCCalculator(path) as calc:
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        blob = pickle.dumps(calc)
    path.unlink()
    clone = pickle.loads(blob)
    moved = water()
    moved.calc = clone
    assert moved.get_potential_energy() == pytest.approx(energy, rel=0, abs=tolerance)
    clone.close()


def test_pickling_leaves_the_original_session_working(opened, calc):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    session = calc.potential
    pickle.dumps(calc)
    copy.copy(calc)
    atoms.positions[0, 2] += 0.01
    atoms.get_forces()
    assert calc.potential is session
    assert len(opened) == 1


def test_a_copy_keeps_the_constant_cutoff(periodic_params):
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        for clone in (copy.deepcopy(calc), pickle.loads(pickle.dumps(calc))):
            assert clone.parameters["constant_cutoff"] == 30.0
            assert clone.todict() == calc.todict()
            clone.close()


def test_a_pickle_holds_no_live_handle(calc):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    state = calc.__getstate__()
    assert state["_potential"] is None
    assert state["_composition"] is None
    assert b"ctypes" not in pickle.dumps(calc)


# --- units ----------------------------------------------------------------------


def test_the_conversion_constants_are_codata_2018():
    exact = create_units("2018")
    assert HARTREE_TO_EV == pytest.approx(exact.Hartree, rel=1e-10)
    assert BOHR_TO_ANGSTROM == pytest.approx(exact.Bohr, rel=1e-10)
    assert HARTREE_TO_EV == pytest.approx(units.Hartree, rel=1e-7)
    assert BOHR_TO_ANGSTROM == pytest.approx(units.Bohr, rel=1e-7)
    assert cpmdc.forces_to_ev_per_angstrom([[1.0, 0.0, 0.0]])[0][0] == pytest.approx(
        exact.Hartree / exact.Bohr, rel=1e-10
    )
    assert cpmdc.stress_to_ev_per_angstrom3([1.0] * 9)[0] == pytest.approx(
        exact.Hartree / exact.Bohr**3, rel=1e-10
    )


def test_energy_and_forces_cross_the_boundary_in_ev_and_ev_per_angstrom(params):
    exact = create_units("2018")
    atoms = water()
    with Session(params) as session:
        native = session.energy_forces(WATER, [8, 1, 1])
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        energy, forces = atoms.get_potential_energy(), atoms.get_forces()
    assert energy == pytest.approx(native.energy_h * exact.Hartree, rel=1e-10)
    assert energy == pytest.approx(native.energy_h * units.Hartree, rel=1e-7)
    expected = np.array(native.forces_h_bohr) * (exact.Hartree / exact.Bohr)
    np.testing.assert_allclose(forces, expected, rtol=1e-10)
    np.testing.assert_allclose(forces, np.array(native.forces_h_bohr) * units.Hartree / units.Bohr, rtol=1e-7)


def test_stress_crosses_the_boundary_in_ev_per_angstrom_cubed(periodic_params):
    exact = create_units("2018")
    atoms = water(cell=np.array(TRICLINIC) * 0.7, pbc=True)
    with Session(periodic_params) as session:
        native = session.energy_forces(atoms.positions, [8, 1, 1], atoms.cell.array)
    assert native.stress_h_bohr3 is not None
    tensor = np.array(native.stress_h_bohr3).reshape(3, 3)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        stress = atoms.get_stress()
        full = atoms.get_stress(voigt=False)
    expected = -0.5 * (tensor + tensor.T) * (exact.Hartree / exact.Bohr**3)
    assert np.abs(expected).max() > 0.0
    np.testing.assert_allclose(full, expected, rtol=1e-9, atol=1e-18)
    np.testing.assert_allclose(voigt_6_to_full_3x3_stress(stress), full, rtol=0, atol=1e-18)


def test_forces_are_the_negative_gradient_of_the_energy_in_the_same_units(periodic_calc, tolerance):
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=1)
    atoms.calc = periodic_calc
    numerical = calculate_numerical_forces(atoms, eps=1e-3)
    margin = 1e-9 if tolerance == 0.0 else 2e-3
    assert np.abs(atoms.get_forces() - numerical).max() < margin
    numerical = calculate_numerical_stress(atoms, eps=1e-5)
    margin = 1e-9 if tolerance == 0.0 else 1e-4
    assert np.abs(atoms.get_stress() - numerical).max() < margin


def test_an_asymmetric_library_tensor_is_averaged_into_voigt(periodic_params, monkeypatch):
    raw = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    evaluate = Potential.evaluate

    def asymmetric(self, *args, **kwargs):
        done = evaluate(self, *args, **kwargs)
        done["stress"] = raw
        return done

    monkeypatch.setattr(Potential, "evaluate", asymmetric)
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        stress = atoms.get_stress()
    np.testing.assert_allclose(stress, [1.0, 5.0, 9.0, 7.0, 5.0, 3.0])
    np.testing.assert_allclose(atoms.get_stress(voigt=False), atoms.get_stress(voigt=False).T)


# --- the ASE calculator object ----------------------------------------------


def test_the_directory_is_created_like_ase_does_it(params, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = water()
    with CPMDCCalculator(params, directory=Path("deep") / "run") as calc:
        atoms.calc = calc
        atoms.get_forces()
    assert (tmp_path / "deep" / "run").is_dir()
    assert Path.cwd() == tmp_path


def test_a_missing_atoms_object_is_not_a_cache_hit(calc, steps):
    assert calc.get_property("energy", allow_calculation=False) is None
    assert steps == []
