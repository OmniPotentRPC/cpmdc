"""Opt-in periodic Gaussian embedding on the reference evaluator.

The oracle is the conducting-boundary Ewald field written in the C ABI
reference. A charge on a face centre of a cubic cell puts every sine at zero,
so the induction energy and the forces are zero.
"""

from __future__ import annotations

import math
import os

import pytest

import cpmdc
from cpmdc import Session

BOHR = cpmdc.BOHR_TO_ANGSTROM
HALF_DECAY = 36.0
INDEX_CAP = 48


def _induction(centres, numbers, charges, cell_bohr):
    """Return induction energy, QM forces (-dE/dR), and MM forces (-dE/dR)."""
    rows = [cell_bohr[3 * r : 3 * r + 3] for r in range(3)]
    det = (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )
    volume = abs(det)
    scale = 2.0 * math.pi / det

    def cross(a, b):
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    recip = [
        tuple(scale * value for value in cross(rows[1], rows[2])),
        tuple(scale * value for value in cross(rows[2], rows[0])),
        tuple(scale * value for value in cross(rows[0], rows[1])),
    ]
    sigma_min = min(charge[2] for charge in charges)
    gcut = math.sqrt(2.0 * HALF_DECAY) / sigma_min
    bounds = []
    for axis in rows:
        length = math.sqrt(sum(value * value for value in axis))
        ratio = gcut * length / (2.0 * math.pi)
        if ratio > INDEX_CAP:
            ratio = float(INDEX_CAP)
        bound = math.ceil(ratio) if ratio > 0.0 else 0
        bounds.append(min(INDEX_CAP, max(0, bound)))
    four_pi = 4.0 * math.pi
    energy = 0.0
    qm_forces = []
    mm_forces = [[0.0, 0.0, 0.0] for _ in charges]

    def vectors():
        for h in range(-bounds[0], bounds[0] + 1):
            for k in range(-bounds[1], bounds[1] + 1):
                for ell in range(-bounds[2], bounds[2] + 1):
                    if h == 0 and k == 0 and ell == 0:
                        continue
                    gvec = tuple(
                        h * recip[0][c] + k * recip[1][c] + ell * recip[2][c]
                        for c in range(3)
                    )
                    g2 = sum(value * value for value in gvec)
                    if g2 <= 0.0:
                        continue
                    if 0.5 * sigma_min * sigma_min * g2 > HALF_DECAY:
                        continue
                    yield gvec, g2

    shell = list(vectors())
    for centre, atomic_number in zip(centres, numbers):
        alpha = 0.20 * max(1, atomic_number)
        field = [0.0, 0.0, 0.0]
        for position, q, sigma in charges:
            delta = tuple(centre[c] - position[c] for c in range(3))
            for gvec, g2 in shell:
                phase = sum(gvec[c] * delta[c] for c in range(3))
                pref = (four_pi * q / volume) * math.exp(-0.5 * sigma * sigma * g2) / g2
                sine = math.sin(phase)
                for c in range(3):
                    field[c] += pref * gvec[c] * sine
        energy += -0.5 * alpha * sum(value * value for value in field)
        qm = [0.0, 0.0, 0.0]
        for index, (position, q, sigma) in enumerate(charges):
            delta = tuple(centre[c] - position[c] for c in range(3))
            dfield = [[0.0, 0.0, 0.0] for _ in range(3)]
            for gvec, g2 in shell:
                phase = sum(gvec[c] * delta[c] for c in range(3))
                pref = (four_pi * q / volume) * math.exp(-0.5 * sigma * sigma * g2) / g2
                cosine = math.cos(phase)
                for b in range(3):
                    for a in range(3):
                        dfield[b][a] += pref * gvec[b] * cosine * gvec[a]
            for a in range(3):
                contrib = alpha * sum(field[b] * dfield[b][a] for b in range(3))
                qm[a] += contrib
                mm_forces[index][a] -= contrib
        qm_forces.extend(qm)
    return energy, qm_forces, [value for row in mm_forces for value in row]


def _angstrom(bohr):
    return [value * BOHR for value in bohr]


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_face_centre_charge_has_zero_induction(periodic_params, monkeypatch, cpmd_backend):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    side = 10.0
    cell = [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side]
    charge = _angstrom([side / 2.0, 0.0, 0.0])
    width = [1.0 * BOHR]
    with Session(periodic_params) as session:
        bare, bare_forces, _ = session.calculate_au([0.0, 0.0, 0.0], [1], cell)
        session.set_external_charges(charge, [0.5], width)
        energy, forces, _ = session.calculate_au([0.0, 0.0, 0.0], [1], cell)
        mm = session.external_charge_forces_au(1)
    assert energy - bare == pytest.approx(0.0, abs=1e-12)
    for got, base in zip(forces, bare_forces):
        assert got - base == pytest.approx(0.0, abs=1e-12)
    for component in mm:
        assert component == pytest.approx(0.0, abs=1e-12)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_ewald_induction_matches_the_reciprocal_sum(periodic_params, monkeypatch, cpmd_backend):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    side = 8.0
    cell = [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side]
    centres = [[0.2, -0.3, 0.4], [1.1, 0.5, -0.7]]
    numbers = [8, 1]
    charges = [
        ([2.0, 0.4, -0.2], 0.40, 1.20),
        ([1.0, 2.0, -0.5], -0.25, 0.80),
    ]
    flat = [value for centre in centres for value in centre]
    expected_energy, expected_qm, expected_mm = _induction(centres, numbers, charges, cell)
    with Session(periodic_params) as session:
        bare, bare_forces, virial_bare = session.calculate_au(flat, numbers, cell)
        session.set_external_charges(
            [value * BOHR for charge in charges for value in charge[0]],
            [charge[1] for charge in charges],
            [charge[2] * BOHR for charge in charges],
        )
        energy, forces, virial = session.calculate_au(flat, numbers, cell)
        mm = list(session.external_charge_forces_au(len(charges)))
    assert energy - bare == pytest.approx(expected_energy, rel=0.0, abs=1e-9)
    for got, base, want in zip(forces, bare_forces, expected_qm):
        assert got - base == pytest.approx(want, rel=0.0, abs=1e-9)
    for got, want in zip(mm, expected_mm):
        assert got == pytest.approx(want, rel=0.0, abs=1e-9)
    assert list(virial) == [0.0] * 9
    assert any(value != 0.0 for value in virial_bare)

    step = 1.0e-5
    moved = [flat[0] + step, *flat[1:]]
    with Session(periodic_params) as session:
        session.set_external_charges(
            [value * BOHR for charge in charges for value in charge[0]],
            [charge[1] for charge in charges],
            [charge[2] * BOHR for charge in charges],
        )
        plus, _, _ = session.calculate_au(moved, numbers, cell)
        moved[0] -= 2.0 * step
        minus, _, _ = session.calculate_au(moved, numbers, cell)
    finite = -(plus - minus) / (2.0 * step)
    assert forces[0] == pytest.approx(finite, rel=0.0, abs=1e-7)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_periodic_embedding_stays_off_without_the_variable(periodic_params, monkeypatch, cpmd_backend):
    monkeypatch.delenv("CPMDC_PERIODIC_EMBEDDING", raising=False)
    cell = [8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0]
    with Session(periodic_params) as session:
        session.set_external_charges([1.0, 0.0, 0.0], [0.2], [0.3])
        with pytest.raises(cpmdc.CPMDCUnsupportedPeriodicEmbedding):
            session.calculate_au([0.0, 0.0, 0.0], [1], cell)
    assert os.environ.get("CPMDC_PERIODIC_EMBEDDING") is None


def _rows_volume(cell_bohr):
    rows = [cell_bohr[3 * r : 3 * r + 3] for r in range(3)]
    det = (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )
    return rows, abs(det)


def _direct_field(delta, charge, sigma, rows, volume, nshell):
    """Spherical sum of the isolated Gaussian field, minus the Lorentz field.

    That surface term is what makes the spherical sample agree with the
    conducting boundary of the reciprocal Ewald sum.
    """
    alpha = 1.0 / (sigma * math.sqrt(2.0))
    field = [0.0, 0.0, 0.0]
    limit = nshell * nshell
    for n1 in range(-nshell, nshell + 1):
        for n2 in range(-nshell, nshell + 1):
            for n3 in range(-nshell, nshell + 1):
                if n1 * n1 + n2 * n2 + n3 * n3 > limit:
                    continue
                disp = [
                    delta[c] - (n1 * rows[0][c] + n2 * rows[1][c] + n3 * rows[2][c])
                    for c in range(3)
                ]
                radius = math.sqrt(sum(value * value for value in disp))
                if radius == 0.0:
                    continue
                scaled = alpha * radius
                screened = 0.0 if scaled > 20.0 else math.exp(-scaled * scaled)
                lam = math.erf(scaled) / radius**3 - (
                    (2.0 / math.sqrt(math.pi)) * alpha * screened / radius**2
                )
                for component in range(3):
                    field[component] += charge * lam * disp[component]
    lorentz = 4.0 * math.pi * charge / (3.0 * volume)
    for component in range(3):
        field[component] -= lorentz * delta[component]
    return field


def _direct_induction(centres, numbers, charges, cell_bohr, nshell=8):
    rows, volume = _rows_volume(cell_bohr)
    energy = 0.0
    for centre, atomic_number in zip(centres, numbers):
        polar = 0.20 * max(1, atomic_number)
        field = [0.0, 0.0, 0.0]
        for position, charge, sigma in charges:
            delta = tuple(centre[c] - position[c] for c in range(3))
            piece = _direct_field(delta, charge, sigma, rows, volume, nshell)
            for component in range(3):
                field[component] += piece[component]
        energy += -0.5 * polar * sum(value * value for value in field)
    return energy


def _charged_induction(session_params, centre, numbers, cell, position, charge, width):
    with Session(session_params) as session:
        bare, _, _ = session.calculate_au(centre, numbers, cell)
        session.set_external_charges(position, [charge], width)
        energy, forces, _ = session.calculate_au(centre, numbers, cell)
    return energy - bare, forces


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_sheared_cell_matches_the_reduced_basis(periodic_params, monkeypatch, cpmd_backend):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    sheared = [10.0, 0.0, 0.0, 10.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    reduced = [10.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    centre = [0.0, 0.0, 0.0]
    numbers = [1]
    charges = [([1.0, 0.0, 0.0], 0.5, 1.0)]
    position = _angstrom([1.0, 0.0, 0.0])
    width = [1.0 * BOHR]
    sheared_energy, forces = _charged_induction(
        periodic_params, centre, numbers, sheared, position, 0.5, width
    )
    reduced_energy, _ = _charged_induction(
        periodic_params, centre, numbers, reduced, position, 0.5, width
    )
    expected, _, _ = _induction([centre], numbers, charges, sheared)
    assert sheared_energy == pytest.approx(reduced_energy, rel=0.0, abs=1e-9)
    assert sheared_energy == pytest.approx(expected, rel=0.0, abs=1e-9)

    step = 1.0e-5
    with Session(periodic_params) as session:
        session.set_external_charges(position, [0.5], width)
        plus, _, _ = session.calculate_au([step, 0.0, 0.0], numbers, sheared)
        minus, _, _ = session.calculate_au([-step, 0.0, 0.0], numbers, sheared)
    finite = -(plus - minus) / (2.0 * step)
    assert forces[0] == pytest.approx(finite, rel=0.0, abs=1e-7)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_narrow_width_matches_a_direct_sum(periodic_params, monkeypatch, cpmd_backend):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    side = 8.0
    cell = [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side]
    centre = [0.2, -0.1, 0.3]
    numbers = [1]
    sigma = 1.0e-12 / BOHR
    charges = [([1.5, 0.4, -0.2], 0.5, sigma)]
    expected = _direct_induction([centre], numbers, charges, cell)
    energy, forces = _charged_induction(
        periodic_params,
        centre,
        numbers,
        cell,
        _angstrom(charges[0][0]),
        0.5,
        [1.0e-12],
    )
    assert energy == pytest.approx(expected, rel=0.0, abs=1e-8)

    step = 1.0e-5
    moved = [centre[0] + step, centre[1], centre[2]]
    with Session(periodic_params) as session:
        session.set_external_charges(_angstrom(charges[0][0]), [0.5], [1.0e-12])
        plus, _, _ = session.calculate_au(moved, numbers, cell)
        moved[0] -= 2.0 * step
        minus, _, _ = session.calculate_au(moved, numbers, cell)
    finite = -(plus - minus) / (2.0 * step)
    assert forces[0] == pytest.approx(finite, rel=0.0, abs=1e-6)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_default_width_in_a_long_cell_matches_a_direct_sum(
    periodic_params, monkeypatch, cpmd_backend
):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    side = 20.0 / BOHR
    cell = [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side]
    centre = [0.0, 0.0, 0.0]
    numbers = [1]
    sigma = 0.30 / BOHR
    charges = [([1.0, 0.4, -0.3], 1.0, sigma)]
    expected = _direct_induction([centre], numbers, charges, cell)
    energy, _ = _charged_induction(
        periodic_params,
        centre,
        numbers,
        cell,
        _angstrom(charges[0][0]),
        1.0,
        [0.30],
    )
    assert energy == pytest.approx(expected, rel=0.0, abs=1e-8)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_anisotropic_cell_reports_an_unconverged_embedding(
    periodic_params, monkeypatch, cpmd_backend
):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    cell = [400.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    with Session(periodic_params) as session:
        session.set_external_charges(_angstrom([1.0, 0.0, 0.0]), [0.5], [1.0 * BOHR])
        with pytest.raises(cpmdc.CPMDCEvaluationFailed, match="did not converge"):
            session.calculate_au([0.0, 0.0, 0.0], [1], cell)


@pytest.mark.cpmd_backend(
    "reference",
    reason="the Ewald oracle describes the reference backend, not an OpenCPMD SCF",
)
def test_periodic_deck_without_a_cell_stays_unsupported(
    periodic_params, monkeypatch, cpmd_backend
):
    monkeypatch.setenv("CPMDC_PERIODIC_EMBEDDING", "1")
    with Session(periodic_params) as session:
        session.set_external_charges([1.0, 0.0, 0.0], [0.2], [0.3])
        with pytest.raises(cpmdc.CPMDCUnsupportedPeriodicEmbedding):
            session.calculate_au([0.0, 0.0, 0.0], [1])
