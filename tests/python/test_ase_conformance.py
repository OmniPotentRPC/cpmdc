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
import json
import os
import pickle
import re
from pathlib import Path
from typing import Any

import ase.db
import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.calculator import (
    CalculationFailed,
    CalculatorError,
    PropertyNotImplementedError,
    ReadError,
    all_changes,
    all_properties,
)
from ase.calculators.mixing import SumCalculator
from ase_codata import derived_matches_recommended
from ase_fd import calculate_numerical_forces, calculate_numerical_stress
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.md.verlet import VelocityVerlet
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
    # The calculator multiplies by cpmdc's constants, the CODATA 2018
    # recommended values. ASE's derived Bohr and Hartree match them when
    # the 2018 table stores the CODATA 2018 mu0, which it does from ASE 3.24.
    # ASE 3.23 still stores mu0 = 4*pi*1e-7, so the derived values miss by
    # about 1e-9. The charge, Planck constant, Boltzmann constant, and
    # Avogadro number identify the 2018 table on every version.
    assert exact["_e"] == 1.602176634e-19
    assert exact["_hplanck"] == 6.62607015e-34
    assert exact["_k"] == 1.380649e-23
    assert exact["_Nav"] == 6.02214076e23
    if derived_matches_recommended(exact):
        assert HARTREE_TO_EV == pytest.approx(exact.Hartree, rel=1e-10)
        assert BOHR_TO_ANGSTROM == pytest.approx(exact.Bohr, rel=1e-10)
        assert cpmdc.forces_to_ev_per_angstrom([[1.0, 0.0, 0.0]])[0][0] == pytest.approx(
            exact.Hartree / exact.Bohr, rel=1e-10
        )
        assert cpmdc.stress_to_ev_per_angstrom3([1.0] * 9)[0] == pytest.approx(
            exact.Hartree / exact.Bohr**3, rel=1e-10
        )
    else:
        assert cpmdc.forces_to_ev_per_angstrom([[1.0, 0.0, 0.0]])[0][0] == pytest.approx(
            HARTREE_TO_EV / BOHR_TO_ANGSTROM
        )
        assert cpmdc.stress_to_ev_per_angstrom3([1.0] * 9)[0] == pytest.approx(
            HARTREE_TO_EV / BOHR_TO_ANGSTROM**3
        )
    assert HARTREE_TO_EV == pytest.approx(units.Hartree, rel=1e-7)
    assert BOHR_TO_ANGSTROM == pytest.approx(units.Bohr, rel=1e-7)


def test_energy_and_forces_cross_the_boundary_in_ev_and_ev_per_angstrom(params):
    atoms = water()
    with Session(params) as session:
        native = session.energy_forces(WATER, [8, 1, 1])
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        energy, forces = atoms.get_potential_energy(), atoms.get_forces()
    assert energy == pytest.approx(native.energy_h * HARTREE_TO_EV, rel=1e-10)
    assert energy == pytest.approx(native.energy_h * units.Hartree, rel=1e-7)
    expected = np.array(native.forces_h_bohr) * (HARTREE_TO_EV / BOHR_TO_ANGSTROM)
    np.testing.assert_allclose(forces, expected, rtol=1e-10)
    np.testing.assert_allclose(forces, np.array(native.forces_h_bohr) * units.Hartree / units.Bohr, rtol=1e-7)


def test_stress_crosses_the_boundary_in_ev_per_angstrom_cubed(periodic_params):
    atoms = water(cell=np.array(TRICLINIC) * 0.7, pbc=True)
    with Session(periodic_params) as session:
        native = session.energy_forces(atoms.positions, [8, 1, 1], atoms.cell.array)
    assert native.stress_h_bohr3 is not None
    tensor = np.array(native.stress_h_bohr3).reshape(3, 3)
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        stress = atoms.get_stress()
        full = atoms.get_stress(voigt=False)
    expected = -0.5 * (tensor + tensor.T) * (HARTREE_TO_EV / BOHR_TO_ANGSTROM**3)
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


# --- adapted from ASE's own calculator tests ----------------------------------
#
# ASE has no single conformance suite. These cases follow what its tests ask
# of every calculator: ase/test/calculator/test_traj.py (results survive a
# trajectory, calculation_required after a step), test_calculator.py
# (get_spin_polarized), test_singlepointcalc.py and fio/test_extxyz.py
# (results through SinglePointCalculator and extended XYZ), and
# calculator/vasp/test_vasp_check_state.py (a changed parameter discards the
# step, get_atoms keeps the results).
#
# ASE rebuilds a calculator from its name and its constructor keywords: a
# trajectory and an ase.db row store calc.name and calc.todict(), and
# PythonSubProcessCalculator sends the keywords, not the calculator, because
# it assumes a calculator may not pickle. ASE asks no pickling of Calculator.
# copy.deepcopy(atoms) and the pickle that ase.visualize.view sends to its
# viewer process copy a live calculator; the copies-and-pickles section
# above covers those. What ASE needs is that todict() is JSON data and that
# CPMDCCalculator(**calc.todict()) evaluates the same method wherever it is
# called.

VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
# ase.io.extxyz writes float columns with "%16.8f".
EXTXYZ_COLUMN_DECIMALS = 8


def test_unchanged_atoms_make_no_second_library_call(periodic_calc, steps):
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.calc = periodic_calc
    energy = atoms.get_potential_energy()
    assert atoms.get_potential_energy() == energy
    atoms.get_forces()
    atoms.get_stress()
    atoms.get_potential_energy(force_consistent=True)
    # Dipole is implemented, and this deck does not produce one. Asking for
    # it reports that without evaluating the step again.
    stored = [name for name in periodic_calc.implemented_properties if name != "dipole"]
    for name in stored:
        assert periodic_calc.get_property(name, atoms) is not None
    assert not periodic_calc.calculation_required(atoms, stored)
    with pytest.raises(PropertyNotImplementedError):
        periodic_calc.get_property("dipole", atoms)
    assert len(steps) == 1


def changed_parameter(what: str, periodic_params: Path, tmp_path: Path) -> dict[str, Any]:
    """A new value of one parameter that evaluates the same deck, or the constant cutoff."""
    if what == "params":
        other = tmp_path / "copy.periodic.bin"
        other.write_bytes(Path(periodic_params).read_bytes())
        return {"params": str(other)}
    if what == "library":
        # The dynamic linker maps a link to the loaded library onto the same image.
        link = tmp_path / "libcpmdc-link.so"
        link.symlink_to(cpmdc.find_library().resolve())
        return {"library": str(link)}
    if what == "constant_cutoff":
        return {"constant_cutoff": 30.0}
    raise AssertionError(what)


@pytest.mark.parametrize("what", ["params", "library", "constant_cutoff"])
def test_a_changed_parameter_discards_the_stored_step(
    periodic_calc, periodic_params, tmp_path, steps, tolerance, what
):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    energy = atoms.get_potential_energy()
    change = changed_parameter(what, periodic_params, tmp_path)
    assert periodic_calc.set(**change) == change
    assert periodic_calc.results == {}
    assert periodic_calc.get_property("energy", atoms, allow_calculation=False) is None
    assert periodic_calc.calculation_required(atoms, ["energy"])
    assert len(steps) == 1
    again = atoms.get_potential_energy()
    assert len(steps) == 2
    assert periodic_calc.set(**change) == {}
    assert atoms.get_potential_energy() == again
    assert len(steps) == 2
    if what != "constant_cutoff":
        assert again == pytest.approx(energy, rel=0, abs=tolerance)


def test_get_atoms_keeps_the_stored_step(periodic_calc, steps):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    atoms.get_stress()
    before = dict(periodic_calc.results)
    copied = periodic_calc.get_atoms()
    assert copied.calc is periodic_calc
    assert periodic_calc.results == before
    assert copied.get_potential_energy() == before["energy"]
    assert len(steps) == 1


@pytest.mark.parametrize("with_atoms", [True, False], ids=["atoms", "calculator-atoms"])
def test_get_property_without_calculation_returns_none_until_a_step(calc, steps, with_atoms):
    atoms = water()
    atoms.calc = calc
    target = atoms if with_atoms else None
    for name in calc.implemented_properties:
        assert calc.get_property(name, target, allow_calculation=False) is None
    atoms.get_potential_energy()
    for name in STEP:
        assert calc.get_property(name, target, allow_calculation=False) is not None
    # A cluster step stores no stress, so the stored step has none to give.
    assert calc.get_property("stress", target, allow_calculation=False) is None
    calc.reset()
    assert calc.get_property("energy", atoms, allow_calculation=False) is None
    assert len(steps) == 1


@pytest.mark.parametrize("name", ["energies", "stresses", "charges", "magmom", "magmoms"])
def test_an_unimplemented_property_raises_even_without_calculation(calc, steps, name):
    atoms = water()
    atoms.calc = calc
    atoms.get_potential_energy()
    with pytest.raises(PropertyNotImplementedError):
        calc.get_property(name, atoms, allow_calculation=False)
    assert len(steps) == 1


def test_the_spin_question_has_the_ase_default_answer(calc):
    for candidate in (calc, pickle.loads(pickle.dumps(calc)), copy.deepcopy(calc)):
        with pytest.warns(FutureWarning):
            assert candidate.get_spin_polarized() is False


def snapshot(atoms: Atoms) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        atoms.positions.copy(),
        atoms.numbers.copy(),
        atoms.cell.array.copy(),
        atoms.pbc.copy(),
    )


@pytest.mark.parametrize("route", ["property", "calculate", "constructor"])
def test_mutating_the_callers_atoms_leaves_the_calculator_copy(periodic_params, steps, route):
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(
        periodic_params, atoms=atoms if route == "constructor" else None
    ) as calc:
        if route == "calculate":
            calc.calculate(atoms, ["energy"], all_changes)
        else:
            atoms.calc = calc
            atoms.get_forces()
        assert calc.atoms is not atoms
        kept = snapshot(calc.atoms)
        atoms.positions[0, 2] += 0.1
        atoms.numbers[1] = 3
        atoms.cell.array[0, 0] *= 1.1
        atoms.pbc[2] = False
        for stored, now in zip(snapshot(calc.atoms), kept):
            np.testing.assert_array_equal(stored, now)
        assert set(calc.check_state(atoms)) == {"positions", "numbers", "cell", "pbc"}
        returned = calc.get_atoms()
        returned.positions[:] = 0.0
        returned.cell.array[:] = 0.0
        np.testing.assert_array_equal(calc.atoms.positions, kept[0])
        np.testing.assert_array_equal(calc.atoms.cell.array, kept[2])
        assert calc.get_property("energy", allow_calculation=False) is not None
        assert len(steps) == 1


def test_the_ase_trajectory_contract(periodic_calc, tmp_path):
    """ase/test/calculator/test_traj.py, with the periodic water for H2."""
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.calc = periodic_calc
    energy = atoms.get_potential_energy()
    assert not periodic_calc.calculation_required(atoms, ["energy"])
    forces = atoms.get_forces()
    assert not periodic_calc.calculation_required(atoms, ["energy", "forces"])
    stress = atoms.get_stress()
    write(tmp_path / "water.traj", atoms)
    back = read(tmp_path / "water.traj")
    assert abs(energy - back.get_potential_energy()) < 1e-12
    assert abs(forces - back.get_forces()).max() < 1e-12
    assert abs(stress - back.get_stress()).max() < 1e-12
    assert back.calc.name == periodic_calc.name
    assert back.calc.parameters == periodic_calc.todict()


def test_extxyz_round_trips_a_single_point_step(periodic_calc, tmp_path):
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.rattle(0.05, seed=3)
    atoms.calc = periodic_calc
    energy, forces, stress = atoms.get_potential_energy(), atoms.get_forces(), atoms.get_stress()
    frozen = atoms.copy()
    frozen.calc = SinglePointCalculator(frozen, energy=energy, forces=forces, stress=stress)
    write(tmp_path / "step.extxyz", frozen)
    back = read(tmp_path / "step.extxyz")
    assert isinstance(back.calc, SinglePointCalculator)
    assert set(back.calc.results) == {"energy", "forces", "stress"}
    # Info values keep their repr; float columns keep the printed decimals.
    assert back.get_potential_energy() == energy
    np.testing.assert_array_equal(back.get_stress(), stress)
    np.testing.assert_allclose(
        back.get_forces(), forces, rtol=0, atol=0.5 * 10.0**-EXTXYZ_COLUMN_DECIMALS
    )
    np.testing.assert_array_equal(back.cell.array, atoms.cell.array)
    assert back.pbc.all()


def test_an_ase_database_row_rebuilds_the_calculator(periodic_params, tmp_path, tolerance):
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        stress = atoms.get_stress()
        with ase.db.connect(tmp_path / "rows.json") as db:
            db.write(atoms)
        row = ase.db.connect(tmp_path / "rows.json").get(id=1)
        assert row.calculator == calc.name
        assert row.calculator_parameters == calc.todict()
        stored = row.toatoms()
        np.testing.assert_array_equal(stored.get_stress(), stress)
        assert stored.get_potential_energy() == atoms.get_potential_energy()
        with CPMDCCalculator(**row.calculator_parameters) as rebuilt:
            assert rebuilt.parameters == calc.parameters
            fresh = water(cell=CUBE, pbc=True)
            fresh.calc = rebuilt
            np.testing.assert_allclose(fresh.get_stress(), stress, rtol=0, atol=tolerance)


@pytest.mark.parametrize("what", ["params", "library"])
def test_todict_rebuilds_the_calculator_from_another_directory(
    params, periodic_params, tmp_path, monkeypatch, tolerance, what
):
    """A relative path is read where it is set, and todict() names that file.

    The other directory holds a different ``deck.bin``, the cluster message,
    so a stored relative path would rebuild another method without an error.
    """
    first, other = tmp_path / "first", tmp_path / "other"
    for where in (first, other):
        (where / "lib").mkdir(parents=True)
    (first / "deck.bin").write_bytes(Path(periodic_params).read_bytes())
    (other / "deck.bin").write_bytes(Path(params).read_bytes())
    (first / "lib" / "libcpmdc.so").symlink_to(cpmdc.find_library().resolve())
    keywords = {"params": "deck.bin"}
    if what == "library":
        keywords = {"params": str(first / "deck.bin"), "library": "lib/libcpmdc.so"}
    monkeypatch.chdir(first)
    atoms = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(**keywords) as calc:
        atoms.calc = calc
        stress = atoms.get_stress()
        saved = json.loads(json.dumps(calc.todict()))
        assert Path(saved[what]).is_absolute()
        assert Path(saved[what]).samefile(first / keywords[what])
    monkeypatch.chdir(other)
    with CPMDCCalculator(**saved) as rebuilt:
        assert rebuilt.todict() == saved
        moved = water(cell=CUBE, pbc=True)
        moved.calc = rebuilt
        np.testing.assert_allclose(moved.get_stress(), stress, rtol=0, atol=tolerance)
        clone = pickle.loads(pickle.dumps(rebuilt))
        clone.close()
        clone.reset()
        assert clone.get_stress(moved) == pytest.approx(stress, rel=0, abs=tolerance)
        clone.close()


def test_a_bytes_message_round_trips_through_todict(params, tolerance):
    blob = Path(params).read_bytes()
    with CPMDCCalculator(blob) as calc:
        saved = json.loads(json.dumps(calc.todict()))
    assert saved["params"].keys() == {"cpmdparams_base64"}
    rebuilt = Atoms("OHH", positions=np.array(WATER).copy())
    with CPMDCCalculator(**saved) as calc:
        rebuilt.calc = calc
        through = rebuilt.get_potential_energy()
    direct = Atoms("OHH", positions=np.array(WATER).copy())
    with CPMDCCalculator(blob) as calc:
        direct.calc = calc
        expected = direct.get_potential_energy()
    assert through == pytest.approx(expected, rel=0, abs=tolerance)


def strain_energy(atoms: Atoms, calc: CPMDCCalculator, pair: tuple[int, int], delta: float) -> float:
    """Energy of ``atoms`` under the symmetric strain with ``eps_ij + eps_ji = 2 delta``."""
    i, j = pair
    deformation = np.eye(3)
    deformation[i, j] += delta
    if i != j:
        deformation[j, i] += delta
    strained = atoms.copy()
    strained.set_cell(atoms.cell.array @ deformation, scale_atoms=True)
    strained.calc = calc
    return float(strained.get_potential_energy())


def test_stress_is_the_central_difference_of_the_energy_over_strain(periodic_calc, opencpmd):
    """Each Voigt component against a central difference of the energy.

    The difference is taken here, in eV and Angstrom^3, so the comparison
    checks the sign, the units, the Voigt order, and the factor of the shear
    components at once. ``eps_ii = d`` changes the energy by
    ``d V sigma_ii``; ``eps_ij = eps_ji = d`` changes it by
    ``2 d V sigma_ij``. On OpenCPMD the margin is that of
    ``calculate_numerical_stress`` above.
    """
    atoms = water(cell=TRICLINIC, pbc=True)
    atoms.rattle(0.05, seed=1)
    atoms.calc = periodic_calc
    stress = atoms.get_stress()
    volume = atoms.get_volume()
    delta = 1e-5
    numerical = np.array(
        [
            (strain_energy(atoms, periodic_calc, pair, delta) - strain_energy(atoms, periodic_calc, pair, -delta))
            / ((2.0 if pair[0] == pair[1] else 4.0) * delta * volume)
            for pair in VOIGT_PAIRS
        ]
    )
    margin = 1e-4 if opencpmd else 1e-9
    np.testing.assert_allclose(stress, numerical, rtol=0, atol=margin)
    signed = np.abs(numerical) > 10.0 * margin
    assert signed.any()
    np.testing.assert_array_equal(np.sign(stress[signed]), np.sign(numerical[signed]))


def test_a_cell_filter_relaxation_lowers_the_energy_and_the_stress(periodic_params, steps):
    """A few BFGS steps of FrechetCellFilter on the smallest periodic deck.

    The constant cutoff gives OpenCPMD a stress that follows the energy, as
    ``docs/orgmode/howto/ase.org`` explains; on the reference well it changes
    nothing that this test reads.
    """
    atoms = water(cell=CUBE, pbc=True)
    atoms.rattle(0.05, seed=17)
    with CPMDCCalculator(periodic_params, constant_cutoff=30.0) as calc:
        atoms.calc = calc
        start_energy = atoms.get_potential_energy()
        start_stress = np.abs(atoms.get_stress()).max()
        with BFGS(FrechetCellFilter(atoms), logfile=None) as opt:
            converged = opt.run(fmax=0.01, steps=8)
        end_stress = np.abs(atoms.get_stress()).max()
        assert atoms.get_potential_energy() < start_energy
    assert converged or end_stress < start_stress
    assert not np.allclose(atoms.cell.array, CUBE)
    assert len(steps) == opt.nsteps + 1


# --- read and restart ------------------------------------------------------


def test_read_refuses_to_restore_what_the_calculator_never_stored(calc):
    """ASE's default ``read`` sets the label and returns, as if it had restored."""
    with pytest.raises(ReadError, match="keeps no calculator file"):
        calc.read("old/run")
    assert calc.label is None
    assert calc.directory == "."
    assert calc.results == {}
    assert calc.atoms is None


def test_a_read_failure_is_an_ase_calculator_error(calc):
    with pytest.raises(CalculatorError):
        calc.read("old/run")


def test_read_names_the_label_it_refused(calc):
    with pytest.raises(ReadError, match="'old/run'"):
        calc.read("old/run")


def test_a_refused_read_keeps_the_calculator_usable(calc, steps):
    atoms = water()
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    with pytest.raises(ReadError):
        calc.read("old/run")
    assert atoms.get_potential_energy() == energy
    assert len(steps) == 1


# --- state that is not a system change --------------------------------------


def _tag(atoms: Atoms) -> None:
    atoms.set_tags([1, 2, 3])


def _momenta(atoms: Atoms) -> None:
    atoms.set_momenta(np.ones((3, 3)))


def _constraint(atoms: Atoms) -> None:
    atoms.set_constraint(FixAtoms(indices=[0]))


def _info(atoms: Atoms) -> None:
    atoms.info["label"] = "kept"


@pytest.mark.parametrize("touch", [_tag, _momenta, _constraint, _info])
def test_state_outside_the_system_changes_keeps_the_step(calc, steps, touch):
    atoms = water()
    atoms.calc = calc
    atoms.get_forces()
    touch(atoms)
    assert calc.check_state(atoms) == []
    atoms.get_potential_energy()
    assert len(steps) == 1


def test_an_equal_copy_of_the_atoms_reuses_the_step(periodic_calc, steps):
    atoms = water(cell=CUBE, pbc=True)
    atoms.calc = periodic_calc
    stress = atoms.get_stress()
    twin = atoms.copy()
    assert periodic_calc.get_property("stress", twin) is not None
    np.testing.assert_array_equal(periodic_calc.get_stress(twin), stress)
    assert len(steps) == 1


# --- parameters --------------------------------------------------------------


@pytest.mark.parametrize(
    "cutoff",
    [None, 30.0, {"cutoff_ry": 30.0, "height_ry": 2.0, "width_ry": 5.0}],
    ids=["none", "e0", "mapping"],
)
def test_parameters_round_trip_through_todict_and_a_parameters_file(
    periodic_params, tmp_path, cutoff
):
    with CPMDCCalculator(periodic_params, constant_cutoff=cutoff) as calc:
        saved = tmp_path / "calc.json"
        calc.parameters.write(saved)
        with CPMDCCalculator(**calc.todict()) as rebuilt:
            assert rebuilt.parameters == calc.parameters
            assert rebuilt.todict() == calc.todict()
        with CPMDCCalculator(parameters=str(saved)) as loaded:
            assert loaded.parameters == calc.parameters
        json.dumps(calc.todict())


# --- other ASE tools ----------------------------------------------------------


def _summed(params_path: Path, atoms: Atoms) -> tuple[SumCalculator, CPMDCCalculator, CPMDCCalculator]:
    first, second = CPMDCCalculator(params_path), CPMDCCalculator(params_path)
    atoms.calc = SumCalculator([first, second])
    return atoms.calc, first, second


def test_sum_calculator_adds_two_calculators_and_calls_each_once(params, steps, tolerance):
    atoms = water()
    reference = water()
    with CPMDCCalculator(params) as alone:
        reference.calc = alone
        energy, forces = reference.get_potential_energy(), reference.get_forces()
    del steps[:]
    summed, first, second = _summed(params, atoms)
    try:
        np.testing.assert_allclose(atoms.get_potential_energy(), 2.0 * energy, rtol=0, atol=tolerance)
        np.testing.assert_allclose(atoms.get_forces(), 2.0 * forces, rtol=0, atol=tolerance)
        assert len(steps) == 2
        assert first.results.keys() == second.results.keys() == STEP
    finally:
        first.close()
        second.close()


def test_sum_calculator_adds_the_stress_of_a_periodic_cell(periodic_params, tolerance):
    atoms = water(cell=CUBE, pbc=True)
    reference = water(cell=CUBE, pbc=True)
    with CPMDCCalculator(periodic_params) as alone:
        reference.calc = alone
        stress = reference.get_stress()
    summed, first, second = _summed(periodic_params, atoms)
    try:
        np.testing.assert_allclose(atoms.get_stress(), 2.0 * stress, rtol=0, atol=tolerance)
    finally:
        first.close()
        second.close()


def test_sum_calculator_passes_on_the_refusal_of_stress_for_a_cluster(params):
    atoms = water()
    summed, first, second = _summed(params, atoms)
    try:
        with pytest.raises(PropertyNotImplementedError):
            atoms.get_stress()
    finally:
        first.close()
        second.close()


def test_velocity_verlet_conserves_the_energy_with_one_step_per_move(calc, steps, opencpmd):
    atoms = water()
    atoms.calc = calc
    atoms.set_momenta(np.array([[0.0, 0.0, 0.02], [0.0, 0.03, -0.01], [0.0, -0.03, -0.01]]))
    with VelocityVerlet(atoms, 0.2 * units.fs, logfile=None) as dyn:
        start = atoms.get_total_energy()
        dyn.run(10)
    drift = abs(atoms.get_total_energy() - start)
    assert drift < (1e-3 if opencpmd else 1e-6)
    assert len(steps) == 1 + 10
