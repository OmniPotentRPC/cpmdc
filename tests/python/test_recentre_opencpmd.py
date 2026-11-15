"""Rigid translations of a cluster in the isolated cell of an OpenCPMD deck.

The cluster water deck has a 10 Angstrom Hockney cell with its origin at a
corner. ``CPMDC_RECENTER`` (default ``mass``) moves the atoms to the middle
of that cell before each evaluation and projects the net force out, so a
host that lets the QM region drift sees the same energy, forces and
``sum_i x_i (x) f_i`` wherever the region sits. These cases need a library
linked to OpenCPMD and ``CPMDC_PSEUDO_DIR``; they skip otherwise.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
from conftest import WATER

from cpmdc import CPMDCError, Session

NUMBERS = [8, 1, 1]
BOX = 10.0
# Largest |y| of the water is 0.7572 Angstrom, so this puts a hydrogen
# 0.5 Angstrom from the y = BOX face.
NEAR_FACE_Y = BOX - 0.5 - 0.7572
SHIFTS = {
    "centred": (5.0, 5.0, 5.0),
    "two_angstrom": (5.0, BOX - 2.0 - 0.7572, 5.0),
    "near_face": (5.0, NEAR_FACE_Y, 5.0),
    "straddling": (0.0, 0.0, 0.0),
    "next_image": (15.0, 5.0, -5.0),
}


@pytest.fixture(autouse=True)
def _needs_opencpmd(opencpmd, monkeypatch):
    if not opencpmd:
        pytest.skip("the reference well is not translation invariant by design")
    if not os.environ.get("CPMDC_PSEUDO_DIR"):
        pytest.skip("CPMDC_PSEUDO_DIR is not set")
    monkeypatch.delenv("CPMDC_RECENTER", raising=False)
    monkeypatch.delenv("CPMDC_RECENTER_MARGIN", raising=False)


def _water(shift):
    return np.asarray(WATER) + np.asarray(shift)


def _evaluate(session, positions):
    result = session.energy_forces(positions.tolist(), NUMBERS)
    forces = np.asarray(result.forces_h_bohr)
    # The virial the GROMACS shim reports for CELL NONE, in Angstrom Hartree/Bohr.
    virial = positions.T @ forces
    return result.energy_h, forces, 0.5 * (virial + virial.T)


def test_translation_leaves_energy_forces_and_virial(params):
    with Session(params) as session:
        results = {name: _evaluate(session, _water(shift)) for name, shift in SHIFTS.items()}
    energy0, forces0, virial0 = results["centred"]
    for name, (energy, forces, virial) in results.items():
        assert abs(energy - energy0) < 1.0e-7, name
        np.testing.assert_allclose(forces, forces0, atol=1.0e-5, err_msg=name)
        np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1.0e-10, err_msg=name)
        np.testing.assert_allclose(virial, virial0, atol=1.0e-4, err_msg=name)


def test_geometry_centre_is_also_invariant(params, monkeypatch):
    monkeypatch.setenv("CPMDC_RECENTER", "geometry")
    with Session(params) as session:
        energy0, forces0, _ = _evaluate(session, _water(SHIFTS["centred"]))
        energy, forces, _ = _evaluate(session, _water(SHIFTS["straddling"]))
    assert abs(energy - energy0) < 1.0e-7
    np.testing.assert_allclose(forces, forces0, atol=1.0e-5)
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1.0e-10)


def test_off_shows_the_fixed_cell(params, monkeypatch):
    monkeypatch.setenv("CPMDC_RECENTER", "off")
    monkeypatch.setenv("CPMDC_RECENTER_MARGIN", "-1")
    with Session(params) as session:
        energy0, forces0, _ = _evaluate(session, _water(SHIFTS["centred"]))
        near, _, _ = _evaluate(session, _water(SHIFTS["near_face"]))
        across, forces, _ = _evaluate(session, _water(SHIFTS["straddling"]))
    assert abs(near - energy0) > 1.0e-4
    assert abs(across - energy0) > 1.0e-3
    assert np.linalg.norm(forces.sum(axis=0)) > 0.02
    assert np.linalg.norm(forces0.sum(axis=0)) < 1.0e-2


def test_margin_refuses_an_atom_near_a_face(params, monkeypatch):
    monkeypatch.setenv("CPMDC_RECENTER", "off")
    monkeypatch.setenv("CPMDC_RECENTER_MARGIN", "1.0")
    with Session(params) as session, pytest.raises(CPMDCError, match="CPMDC_RECENTER_MARGIN"):
        session.energy_forces(_water(SHIFTS["near_face"]).tolist(), NUMBERS)


def test_off_refuses_an_atom_outside_the_cell(params, monkeypatch):
    monkeypatch.setenv("CPMDC_RECENTER", "off")
    with Session(params) as session, pytest.raises(CPMDCError, match=r"-0\.757 Angstrom"):
        session.energy_forces(_water(SHIFTS["straddling"]).tolist(), NUMBERS)


def test_margin_applies_after_the_shift(params, monkeypatch):
    monkeypatch.setenv("CPMDC_RECENTER_MARGIN", "4.5")
    with Session(params) as session, pytest.raises(CPMDCError, match="enlarge the cell"):
        session.energy_forces(_water(SHIFTS["near_face"]).tolist(), NUMBERS)


@pytest.mark.parametrize(("name", "value"), [("CPMDC_RECENTER", "middle"), ("CPMDC_RECENTER_MARGIN", "wide")])
def test_unparsed_setting_fails_the_call(params, monkeypatch, name, value):
    monkeypatch.setenv(name, value)
    with Session(params) as session, pytest.raises(CPMDCError, match="CPMDC_RECENTER"):
        session.energy_forces(_water(SHIFTS["centred"]).tolist(), NUMBERS)
