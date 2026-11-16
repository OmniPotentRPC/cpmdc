"""bench/constant_cutoff_npt.py: the thermostats and barostats of run, their states, and the scaling of pressure."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.md.langevinbaoab import LangevinBAOAB
from ase.md.velocitydistribution import thermalize_momenta

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("constant_cutoff_npt", ROOT / "bench" / "constant_cutoff_npt.py")
npt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(npt)

SETUP = argparse.Namespace(dt=2.0, temperature=300.0, pressure=1.0, t_tau=20.0, p_tau=200.0, seed=3)


def copper() -> Atoms:
    atoms = bulk("Cu", cubic=True).repeat(2)
    atoms.calc = EMT()
    thermalize_momenta(atoms, 300.0, rng=np.random.default_rng(1))
    return atoms


def run(dyn, atoms: Atoms, steps: int) -> list[float]:
    temperatures = []
    dyn.attach(lambda: temperatures.append(atoms.get_temperature()), interval=1)
    dyn.run(steps)
    return temperatures


def rate_of(dyn) -> float:
    if isinstance(dyn, npt.MelchionnaNPT):
        return np.trace(dyn.get_strain_rate()) / 3.0
    return dyn.p_eps / dyn.barostat_mass


def test_rescaling_replaces_the_langevin_atoms_and_keeps_the_piston():
    atoms = copper()
    dyn = npt.integrator(SETUP, {"thermostat": "svr"}, atoms)
    assert isinstance(dyn, npt.RescalingBAOAB)
    assert dyn.gamma == 0.0
    assert dyn.barostat_gamma > 0.0
    temperatures = run(dyn, atoms, 400)
    assert abs(np.mean(temperatures[100:]) - 300.0) < 40.0
    langevin = npt.integrator(SETUP, {}, copper())
    assert type(langevin) is LangevinBAOAB and langevin.gamma > 0.0


def test_parrinello_rahman_moves_only_the_volume():
    atoms = copper()
    side = atoms.cell[0, 0]
    dyn = npt.integrator(SETUP, {"barostat": "pr", "bulk_modulus": 140.0}, atoms)
    assert isinstance(dyn, npt.MelchionnaNPT)
    run(dyn, atoms, 50)
    cell = atoms.cell.array
    assert abs(cell[0, 0] - side) > 1e-6
    assert np.allclose(cell, cell[0, 0] * np.identity(3), rtol=0.0, atol=1e-12)


def test_each_state_keeps_the_strain_rate(tmp_path):
    out = tmp_path / "run.json"
    for meta in ({}, {"thermostat": "svr"}, {"barostat": "pr", "bulk_modulus": 140.0}):
        atoms = copper()
        dyn = npt.integrator(SETUP, meta, atoms)
        run(dyn, atoms, 5)
        rate = rate_of(dyn)
        assert rate != 0.0
        npt.save(out, {"args": meta}, atoms, dyn)
        saved = np.load(out.with_suffix(".state.npz"))
        assert npt.strain_rate(saved) == rate
        for other in ({}, {"barostat": "pr", "bulk_modulus": 140.0}):
            again = npt.integrator(SETUP, other, copper())
            npt.set_strain_rate(again, npt.strain_rate(saved))
            assert np.isclose(rate_of(again), rate, rtol=1e-14, atol=0.0)
            if isinstance(again, npt.MelchionnaNPT):
                assert np.allclose(again.get_strain_rate(), rate * np.identity(3), rtol=1e-14, atol=0.0)
                run(again, again.atoms, 2)


def test_the_record_keeps_the_decimals_of_each_field(tmp_path):
    out = tmp_path / "run.json"
    atoms = copper()
    record = {"args": {}, "time_ps": [0.35000000000000003], "Volume": [1.1252695252636975], "Density": [806.41234567],
              "Pressure": [-1566.123456], "Temperature": [352.123456], "Energy": [-15003.3266123456],
              "step_wall_s": [7.32091234], "scf_steps": [9], "setups": [0]}
    npt.save(out, record, atoms, npt.integrator(SETUP, {}, atoms))
    assert json.loads(out.read_text()) == {
        "args": {}, "time_ps": [0.35], "Volume": [1.12527], "Density": [806.412], "Pressure": [-1566.12],
        "Temperature": [352.12], "Energy": [-15003.326612], "step_wall_s": [7.3209], "scf_steps": [9], "setups": [0]}
    assert record["Volume"] == [1.1252695252636975]


def test_a_state_without_a_strain_rate_takes_it_from_the_piston():
    assert npt.strain_rate({"p_eps": np.asarray(3.0), "barostat_mass": np.asarray(4.0)}) == 0.75


def test_scaling_keeps_each_molecule_rigid_across_the_cell_face():
    side = 5.0
    atoms = Atoms("OH2OH2", cell=[side, side, side], pbc=True,
                  positions=[[1.0, 1.0, 1.0], [1.757, 1.586, 1.0], [0.243, 1.586, 1.0],
                             [0.2, 3.0, 3.0], [0.957, 3.586, 3.0], [4.443, 3.586, 3.0]])

    def bonds(a: Atoms) -> np.ndarray:
        r = a.positions.reshape(-1, 3, 3)
        d = r - r[:, :1]
        return d - np.round(d / a.cell[0, 0]) * a.cell[0, 0]

    def centres(a: Atoms) -> np.ndarray:
        m = a.get_masses().reshape(-1, 3, 1)
        return a.positions.reshape(-1, 3, 3)[:, 0] + np.sum(m * bonds(a), axis=1) / np.sum(m, axis=1)

    shape, middle = bonds(atoms), centres(atoms)
    npt.scale_molecules(atoms, 1.1)
    assert np.isclose(atoms.get_volume(), (1.1 * side) ** 3)
    assert np.allclose(bonds(atoms), shape, rtol=0.0, atol=1e-12)
    assert np.allclose(centres(atoms), 1.1 * middle, rtol=0.0, atol=1e-12)
