"""Hirshfeld charges on Session, Potential, and the ASE calculator.

The reference evaluator leaves them unset. A request then raises
``PropertyNotImplementedError`` once, and a second request does not
evaluate again. OpenCPMD stores q = Z_v - n. Neutral water sums to 0
within 2e-3 e, and a permutation of the atoms carries the same charges.
"""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError

from conftest import WATER
from cpmdc import Potential, Session
from cpmdc.ase import CPMDCCalculator

# This 50 Ry Hockney deck leaves a real-space residual of about 1.2e-3 e.
SUM_TOL = 2.0e-3


def _molecule() -> Atoms:
    return Atoms("OHH", positions=WATER)


def test_reference_charges_stay_unset(params, opencpmd):
    if opencpmd:
        pytest.skip("OpenCPMD stores Hirshfeld charges for this water deck")
    with Session(params) as session:
        assert session.charges() is None
        assert session.charges_possible() is False
        done = session.energy_forces(WATER, [8, 1, 1])
        assert done.charges_e is None
        assert session.charges() is None
        assert session.charges_possible() is False
    with Potential(params) as pot:
        result = pot.evaluate(WATER, [8, 1, 1])
        assert "charges" not in result
        assert pot.caps()["charges"] is False


def test_ase_refuses_missing_charges_without_a_second_evaluation(params, opencpmd):
    if opencpmd:
        pytest.skip("OpenCPMD stores Hirshfeld charges for this water deck")
    atoms = _molecule()
    calls = {"n": 0}
    real = Potential.evaluate

    def counting(self, *args, **kwargs):
        calls["n"] += 1
        return real(self, *args, **kwargs)

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(Potential, "evaluate", counting)
    try:
        with CPMDCCalculator(params) as calc:
            atoms.calc = calc
            with pytest.raises(PropertyNotImplementedError, match="per-atom charges"):
                atoms.get_charges()
            with pytest.raises(PropertyNotImplementedError, match="per-atom charges"):
                atoms.get_charges()
            assert calls["n"] == 1
            assert "charges" not in calc.results
            assert "charges" in calc.implemented_properties
    finally:
        monkeypatch.undo()


def test_opencpmd_hirshfeld_charges_sum_and_follow_atoms(params, opencpmd):
    if not opencpmd:
        pytest.skip("the reference evaluator leaves Hirshfeld charges unset")
    perm_pos = [WATER[1], WATER[0], WATER[2]]
    with Session(params) as session:
        assert session.charges() is None
        assert session.charges_possible() is True
        first = session.energy_forces(WATER, [8, 1, 1])
        assert first.charges_e is not None
        charges = np.asarray(first.charges_e, dtype=float)
        assert charges.shape == (3,)
        assert abs(float(charges.sum())) < SUM_TOL
        assert session.charges_possible() is True
        np.testing.assert_allclose(session.charges(), charges)
        cached = session.energy_forces(perm_pos, [1, 8, 1])
        assert cached.charges_e is not None
        np.testing.assert_allclose(
            cached.charges_e,
            [charges[1], charges[0], charges[2]],
            rtol=0,
            atol=1e-12,
        )
    with Session(params) as session:
        cold = session.energy_forces(perm_pos, [1, 8, 1])
    assert cold.charges_e is not None
    np.testing.assert_allclose(
        cold.charges_e,
        [charges[1], charges[0], charges[2]],
        rtol=0,
        atol=1e-8,
    )
    with Potential(params) as pot:
        provided = pot.evaluate(WATER, [8, 1, 1])
        assert pot.caps()["charges"] is True
        np.testing.assert_allclose(provided["charges"], charges, rtol=0, atol=1e-8)
    atoms = _molecule()
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        got = atoms.get_charges()
        np.testing.assert_allclose(got, charges, rtol=0, atol=1e-8)
        assert got.dtype == np.float64
        assert got.shape == (3,)


def test_charge_snapshot_can_be_turned_off(params, opencpmd, monkeypatch):
    if not opencpmd:
        pytest.skip("the reference evaluator leaves Hirshfeld charges unset")
    monkeypatch.setenv("CPMDC_CHARGE_SNAPSHOT", "0")
    with Session(params) as session:
        done = session.energy_forces(WATER, [8, 1, 1])
        assert done.charges_e is None
        assert session.charges() is None
        assert session.charges_possible() is False
    atoms = _molecule()
    with CPMDCCalculator(params) as calc:
        atoms.calc = calc
        with pytest.raises(PropertyNotImplementedError, match="per-atom charges"):
            atoms.get_charges()
        with pytest.raises(PropertyNotImplementedError, match="per-atom charges"):
            atoms.get_charges()
