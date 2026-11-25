"""Electrostatic QM/MM external-charge embedding and typed error tests."""

from __future__ import annotations

import ctypes
import math
import re
import shutil
import subprocess

import mpmath
import numpy as np
import pytest
from conftest import ROOT, WATER

import cpmdc
from cpmdc import CPMDCError, Session

# Relative part of the finite-difference comparison. pytest.approx
# allows max(rel * |expected|, abs).
_FORCE_FD_REL = 2.0e-6
# convergenceOrbitals in tests/water_blyp_embedding_params.capnp.txt.
# The analytic MM force is the Hellmann-Feynman derivative of the grid
# energy. A re-solved SCF difference still carries a remainder of this
# size: the orbital gradient is stopped here, and that remainder is
# first order in the gradient. The direct erf/r coefficient matches a
# fixed-density central difference, so this is not a kernel error.
_EMBED_ORBITAL_GRADIENT = 1.5e-7
_MACHINE_EPS = math.ulp(1.0)


def _finite_difference_allowance(plus, minus, step, opencpmd):
    """Roundoff of a central difference, plus the SCF residual."""
    roundoff = (
        _MACHINE_EPS
        * (abs(plus) + abs(minus))
        / (2.0 * step)
        * cpmdc.BOHR_TO_ANGSTROM
    )
    if opencpmd:
        return max(roundoff, _EMBED_ORBITAL_GRADIENT)
    return roundoff


def _encode_params_text(tmp_path, text):
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is not on PATH")
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode("utf-8"),
        stdout=subprocess.PIPE,
        check=True,
    )
    path = tmp_path / "params.bin"
    path.write_bytes(done.stdout)
    return path


def test_external_charge_qm_and_mm_forces_are_energy_derivatives(
    embedding_params, monkeypatch, capfd
):
    width = [0.30]
    charge = [0.7]
    mm = [3.0, 0.2, -0.1]
    step = 1.0e-5
    monkeypatch.setenv("CPMDC_SCF_STEPS", "1")
    with Session(embedding_params) as session:
        opencpmd = session.available()
        session.set_external_charges(mm, charge, width)
        done = session.energy_forces(WATER, [8, 1, 1])
        moved = [list(row) for row in WATER]
        moved[0][0] += step
        plus = session.energy_forces(moved, [8, 1, 1]).energy_h
        moved[0][0] -= 2.0 * step
        minus = session.energy_forces(moved, [8, 1, 1]).energy_h
        qm_fd = -(plus - minus) / (2.0 * step) * cpmdc.BOHR_TO_ANGSTROM
        qm_abs = _finite_difference_allowance(plus, minus, step, opencpmd)
        assert done.forces_h_bohr[0][0] == pytest.approx(
            qm_fd, rel=_FORCE_FD_REL, abs=qm_abs
        )

        shifted = list(mm)
        shifted[0] += step
        session.set_external_charges(shifted, charge, width)
        plus = session.energy_forces(WATER, [8, 1, 1]).energy_h
        shifted[0] -= 2.0 * step
        session.set_external_charges(shifted, charge, width)
        minus = session.energy_forces(WATER, [8, 1, 1]).energy_h
        mm_fd = -(plus - minus) / (2.0 * step) * cpmdc.BOHR_TO_ANGSTROM
        session.set_external_charges(mm, charge, width)
        session.energy_forces(WATER, [8, 1, 1])
        mm_force = session.external_charge_forces_au(1)
        mm_abs = _finite_difference_allowance(plus, minus, step, opencpmd)
        # The allowance is smaller than the force, so a sign flip
        # (gap 2*|F|) or a factor-of-two error (gap |F|) still fails.
        assert mm_abs < abs(mm_fd)
        assert mm_force[0] == pytest.approx(mm_fd, rel=_FORCE_FD_REL, abs=mm_abs)

    if opencpmd:
        scf_steps = [
            int(value)
            for value in re.findall(r"^cpmdc_scf_steps (\d+)$", capfd.readouterr().err, re.M)
        ]
        assert len(scf_steps) >= 5
        assert all(0 < value < 200 for value in scf_steps)


def test_external_charge_limits_are_finite_and_recover_bare_energy(embedding_params):
    with Session(embedding_params) as session:
        bare = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([1.0e6, 0.0, 0.0], [1.0], [0.30])
        far = session.energy_forces(WATER, [8, 1, 1]).energy_h
        assert far == pytest.approx(bare, abs=1e-12)
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [1.0e-6])
        tiny = session.energy_forces(WATER, [8, 1, 1])
        assert math.isfinite(tiny.energy_h)
        assert all(math.isfinite(value) for row in tiny.forces_h_bohr for value in row)


def test_clearing_external_charges_restores_bare_energy(embedding_params):
    with Session(embedding_params) as session:
        bare = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        embedded = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([], [])
        cleared = session.energy_forces(WATER, [8, 1, 1]).energy_h
        assert abs(embedded - bare) > 1.0e-8
        assert cleared == pytest.approx(bare, abs=5.0e-7)


@pytest.mark.cpmd_backend(
    "opencpmd", reason="the reference well is anchored at the origin"
)
def test_external_charge_translation_invariance(embedding_params, cpmd_backend):
    shift = [0.23, -0.17, 0.31]
    mm = [3.0, 0.2, -0.1]
    with Session(embedding_params) as session:
        session.set_external_charges(mm, [0.7], [0.30])
        original = session.energy_forces(WATER, [8, 1, 1])

        moved_water = [
            [value + delta for value, delta in zip(row, shift)] for row in WATER
        ]
        moved_mm = [value + delta for value, delta in zip(mm, shift)]
        session.set_external_charges(moved_mm, [0.7], [0.30])
        moved = session.energy_forces(moved_water, [8, 1, 1])

    assert moved.energy_h == pytest.approx(original.energy_h, abs=1.0e-8)
    for moved_row, original_row in zip(moved.forces_h_bohr, original.forces_h_bohr):
        assert moved_row == pytest.approx(original_row, rel=1.0e-4, abs=2.0e-8)


class _PropertySnapshot(ctypes.Structure):
    _fields_ = [
        ("valid", ctypes.c_int),
        ("hessian_count", ctypes.c_size_t),
        ("hessian", ctypes.c_double * 4096),
        ("dipole_count", ctypes.c_size_t),
        ("dipole", ctypes.c_double * 3),
        ("polarizability_count", ctypes.c_size_t),
        ("polarizability", ctypes.c_double * 9),
    ]


@pytest.mark.cpmd_backend(
    "reference", reason="OpenCPMD does not export the induced dipole snapshot"
)
def test_reference_embedding_step_has_no_fabricated_dipole(embedding_params, cpmd_backend):
    """The reference evaluator exports no dipole, with or without charges."""
    with Session(embedding_params) as session:
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        session.energy_forces(WATER, [8, 1, 1])
        snapshot = _PropertySnapshot()
        session._lib.cpmdc_session_last_property_snapshot.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(_PropertySnapshot),
        ]
        session._lib.cpmdc_session_last_property_snapshot.restype = ctypes.c_int
        assert (
            session._lib.cpmdc_session_last_property_snapshot(
                session._handle, ctypes.byref(snapshot)
            )
            == 0
        )
        assert snapshot.valid
        assert snapshot.dipole_count == 0
        assert list(snapshot.dipole) == [0.0, 0.0, 0.0]


def test_periodic_external_charges_have_typed_error(periodic_params):
    with Session(periodic_params) as session:
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(
            cpmdc.CPMDCUnsupportedPeriodicEmbedding, match="periodic electrostatic embedding"
        ) as caught:
            session.calculate_au(
                [value / cpmdc.BOHR_TO_ANGSTROM for row in WATER for value in row],
                [8, 1, 1],
                [20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0],
            )
        assert session._lib.cpmdc_last_error_code.argtypes == []
        assert session._lib.cpmdc_last_error_code.restype is ctypes.c_int
        assert caught.value.code == cpmdc.CPMDCUnsupportedPeriodicEmbedding.code
        assert session._lib.cpmdc_last_error_code() == caught.value.code


def test_typed_failure_replaces_previous_error_code(periodic_params):
    with Session(periodic_params) as session:
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(cpmdc.CPMDCUnsupportedPeriodicEmbedding):
            session.calculate_au(
                [value / cpmdc.BOHR_TO_ANGSTROM for row in WATER for value in row],
                [8, 1, 1],
                [20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0],
            )
        assert session._lib.cpmdc_last_error_code() == 3

        assert not session._lib.cpmdc_session_create(None, 0)
        failure = session._fail("cpmdc_session_create failed")
        assert type(failure) is cpmdc.CPMDCInvalidArgument
        assert failure.code == cpmdc.CPMDCInvalidArgument.code
        assert session._lib.cpmdc_last_error_code() == failure.code


def test_isolated_cell_external_charges_evaluate(embedding_params):
    mm = [3.0, 0.2, -0.1]
    charge = [0.7]
    width = [0.30]
    cell_bohr = [20.0, 0.0, 0.0, 0.0, 22.0, 0.0, 0.0, 0.0, 24.0]
    cell_ang = [value * cpmdc.BOHR_TO_ANGSTROM for value in cell_bohr]
    water_bohr = [value / cpmdc.BOHR_TO_ANGSTROM for row in WATER for value in row]
    step = 1.0e-5
    with Session(embedding_params) as session:
        opencpmd = session.available()
        session.set_external_charges(mm, charge, width)
        energy_au, qm_forces_au, virial = session.calculate_au(
            water_bohr, [8, 1, 1], cell_bohr
        )
        assert session._lib.cpmdc_last_error_code() == 0
        assert math.isfinite(energy_au)
        assert all(math.isfinite(value) for value in qm_forces_au)
        assert list(virial) == [0.0] * 9
        mm_forces_au = session.external_charge_forces_au(1)
        assert all(math.isfinite(value) for value in mm_forces_au)

        done = session.energy_forces(WATER, [8, 1, 1], cell=cell_ang)
        assert done.stress_h_bohr3 is None
        assert done.energy_h == pytest.approx(energy_au, rel=1e-12, abs=1e-12)
        mm_forces_msg = session.external_charge_forces_au(1)
        assert list(mm_forces_msg) == pytest.approx(list(mm_forces_au), rel=1e-12, abs=1e-12)

        moved = [list(row) for row in WATER]
        moved[0][0] += step
        plus = session.energy_forces(moved, [8, 1, 1], cell=cell_ang).energy_h
        moved[0][0] -= 2.0 * step
        minus = session.energy_forces(moved, [8, 1, 1], cell=cell_ang).energy_h
        qm_fd = -(plus - minus) / (2.0 * step) * cpmdc.BOHR_TO_ANGSTROM
        qm_abs = _finite_difference_allowance(plus, minus, step, opencpmd)
        assert done.forces_h_bohr[0][0] == pytest.approx(
            qm_fd, rel=_FORCE_FD_REL, abs=qm_abs
        )

        shifted = list(mm)
        shifted[0] += step
        session.set_external_charges(shifted, charge, width)
        plus = session.energy_forces(WATER, [8, 1, 1], cell=cell_ang).energy_h
        shifted[0] -= 2.0 * step
        session.set_external_charges(shifted, charge, width)
        minus = session.energy_forces(WATER, [8, 1, 1], cell=cell_ang).energy_h
        mm_fd = -(plus - minus) / (2.0 * step) * cpmdc.BOHR_TO_ANGSTROM
        mm_abs = _finite_difference_allowance(plus, minus, step, opencpmd)
        assert mm_abs < abs(mm_fd)
        assert mm_forces_msg[0] == pytest.approx(mm_fd, rel=_FORCE_FD_REL, abs=mm_abs)

        session.set_external_charges(mm, charge, width)
        zero_energy, zero_forces, _ = session.calculate_au(
            water_bohr, [8, 1, 1], [0.0] * 9
        )
        assert math.isfinite(zero_energy)
        assert all(math.isfinite(value) for value in zero_forces)
        assert all(math.isfinite(value) for value in session.external_charge_forces_au(1))


def test_smearing_width_changes_embedding_and_origin_is_finite(embedding_params):
    origin = [value for row in WATER[:1] for value in row]
    with Session(embedding_params) as session:
        session.set_external_charges([3.0, 0.2, -0.1], [0.7], [0.20])
        narrow = session.energy_forces(WATER, [8, 1, 1])
        narrow_mm = session.external_charge_forces_au(1)
        session.set_external_charges([3.0, 0.2, -0.1], [0.7], [0.60])
        broad = session.energy_forces(WATER, [8, 1, 1])
        broad_mm = session.external_charge_forces_au(1)
        session.set_external_charges(origin, [1.0e-6], [0.30])
        coincident = session.energy_forces(WATER, [8, 1, 1])
        coincident_mm = session.external_charge_forces_au(1)

    assert math.isfinite(narrow.energy_h)
    assert math.isfinite(broad.energy_h)
    assert all(math.isfinite(value) for row in narrow.forces_h_bohr for value in row)
    assert all(math.isfinite(value) for row in broad.forces_h_bohr for value in row)
    assert all(math.isfinite(value) for value in narrow_mm + broad_mm)
    assert abs(narrow.energy_h - broad.energy_h) > 1.0e-10
    assert math.isfinite(coincident.energy_h)
    assert all(math.isfinite(value) for row in coincident.forces_h_bohr for value in row)
    assert all(math.isfinite(value) for value in coincident_mm)


def test_keyed_session_with_external_charges(embedding_params):
    water_second = [[row[0] + 0.05, row[1], row[2]] for row in WATER]
    mm_first = [3.0, 0.2, -0.1]
    mm_second = [-2.8, 0.1, 0.2]
    with Session(embedding_params) as session:
        session.select_orbitals(1)
        bare_first = session.energy_forces(WATER, [8, 1, 1])
        session.set_external_charges(mm_first, [0.7], [0.30])
        emb_first = session.energy_forces(WATER, [8, 1, 1])
        mm_forces_first = session.external_charge_forces_au(1)

        session.select_orbitals(2)
        session.set_external_charges(mm_second, [-0.5], [0.30])
        emb_second = session.energy_forces(water_second, [8, 1, 1])
        mm_forces_second = session.external_charge_forces_au(1)

        session.select_orbitals(1)
        session.set_external_charges([], [])
        cleared_first = session.energy_forces(WATER, [8, 1, 1])

    assert abs(emb_first.energy_h - bare_first.energy_h) > 1.0e-8
    assert cleared_first.energy_h == pytest.approx(bare_first.energy_h, abs=5.0e-7)
    assert math.isfinite(emb_second.energy_h)
    assert all(math.isfinite(value) for value in mm_forces_first + mm_forces_second)


def _gaussian_reference(radius, sigma, charge):
    with mpmath.workdps(80):
        r = mpmath.mpf(radius)
        width = mpmath.mpf(sigma)
        q = mpmath.mpf(charge)
        root = mpmath.sqrt(2 / mpmath.pi)
        erfv = mpmath.erf(r / (mpmath.sqrt(2) * width))
        expo = mpmath.exp(-mpmath.mpf("0.5") * (r / width) ** 2)
        coefficient = erfv / r**3 - root * expo / (width * r**2)
        derivative = (
            3 * root * expo / (width * r**3)
            + root * expo / (width**3 * r)
            - 3 * erfv / r**4
        )
        field = -q * coefficient * r
        alpha = mpmath.mpf("0.2")
        energy = -mpmath.mpf("0.5") * alpha * field**2
        force = -alpha * q * (coefficient + derivative * r) * field
        return float(energy), float(force)


@pytest.mark.cpmd_backend(
    "reference",
    reason=(
        "the closed-form oracle describes the reference backend, "
        "not an OpenCPMD SCF"
    ),
)
def test_gaussian_small_r_series_is_continuous_and_conservative(embedding_params, cpmd_backend):
    sigma_bohr = 8.0e-12
    charge = sigma_bohr**2
    series_switch = sigma_bohr / 8.0
    old_width_guard = sigma_bohr * math.sqrt(math.sqrt(math.ulp(1.0)))

    with Session(embedding_params) as session:
        bare = session.energy_forces([[0.0, 0.0, 0.0]], [1]).energy_h

        def observed(radius):
            session.set_external_charges(
                [radius * cpmdc.BOHR_TO_ANGSTROM, 0.0, 0.0],
                [charge],
                [sigma_bohr * cpmdc.BOHR_TO_ANGSTROM],
            )
            result = session.energy_forces([[0.0, 0.0, 0.0]], [1])
            return result.energy_h - bare, session.external_charge_forces_au(1)[0]

        assert observed(0.0) == pytest.approx((0.0, 0.0), abs=1e-18)
        for radius in (
            series_switch * (1.0 - 2.0**-10),
            series_switch * (1.0 + 2.0**-10),
            old_width_guard * 0.5,
            old_width_guard * 2.0,
        ):
            energy, force = observed(radius)
            expected_energy, expected_force = _gaussian_reference(
                radius, sigma_bohr, charge
            )
            assert energy == pytest.approx(expected_energy, rel=3e-13, abs=2e-15)
            assert force == pytest.approx(expected_force, rel=3e-13)

            step = radius * 2.0**-12
            plus, _ = observed(radius + step)
            minus, _ = observed(radius - step)
            finite_difference = -(plus - minus) / (2.0 * step)
            assert force == pytest.approx(finite_difference, rel=2e-7)

        below = observed(series_switch * (1.0 - 2.0**-10))
        above = observed(series_switch * (1.0 + 2.0**-10))
        assert below[0] == pytest.approx(above[0], rel=5e-3)
        assert below[1] == pytest.approx(above[1], rel=5e-3)


@pytest.mark.cpmd_backend(
    "reference", reason="the closed-form oracle describes the reference backend"
)
def test_reference_smearing_width_is_converted_from_angstrom(embedding_params, cpmd_backend):
    mm = [3.0, 0.2, -0.1]
    charge = 0.7

    def polarization(width):
        total = 0.0
        sigma = width / cpmdc.BOHR_TO_ANGSTROM
        for xyz, atomic_number in zip(WATER, [8, 1, 1]):
            delta = [
                (value - centre) / cpmdc.BOHR_TO_ANGSTROM
                for value, centre in zip(xyz, mm)
            ]
            radius = math.sqrt(sum(value * value for value in delta))
            kernel = (
                math.erf(radius / (math.sqrt(2.0) * sigma)) / radius**3
                - math.sqrt(2.0 / math.pi)
                * math.exp(-0.5 * (radius / sigma) ** 2)
                / (sigma * radius**2)
            )
            field2 = sum((charge * kernel * value) ** 2 for value in delta)
            total -= 0.5 * (0.20 * atomic_number) * field2
        return total

    with Session(embedding_params) as session:
        session.set_external_charges(mm, [charge], [0.30])
        narrow = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges(mm, [charge], [0.60])
        broad = session.energy_forces(WATER, [8, 1, 1]).energy_h

    expected = polarization(0.30) - polarization(0.60)
    assert narrow - broad == pytest.approx(expected, rel=2e-12, abs=1e-14)


def test_external_charge_force_reader_reports_each_failure(embedding_params):
    with Session(embedding_params) as session:
        with pytest.raises(cpmdc.CPMDCInvalidArgument, match="not available") as missing:
            session.external_charge_forces_au(1)
        assert missing.value.code == cpmdc.CPMDCInvalidArgument.code
        assert session._lib.cpmdc_last_error_code() == missing.value.code

        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(cpmdc.CPMDCInvalidArgument, match="not available"):
            session.external_charge_forces_au(1)

        session.energy_forces(WATER, [8, 1, 1])
        with pytest.raises(cpmdc.CPMDCInvalidArgument, match="wrong length") as wrong:
            session.external_charge_forces_au(2)
        assert wrong.value.code == cpmdc.CPMDCInvalidArgument.code

        forces = session.external_charge_forces_au(1)
        assert session._lib.cpmdc_last_error_code() == cpmdc.CPMDCError.code
        assert session._error() == ""
        assert len(forces) == 3
        assert all(math.isfinite(value) for value in forces)


def test_closed_session_force_reader_stays_generic(params):
    session = Session(params)
    session.close()
    with pytest.raises(CPMDCError, match="closed") as caught:
        session.external_charge_forces_au(0)
    assert type(caught.value) is CPMDCError


def _deck_with_flag(flag, directory):
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is needed to encode a deck with EXTERNAL POTENTIAL")
    source = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text(encoding="utf-8")
    needle = "optimizeWavefunction = true,"
    edited = source.replace(needle, f"{needle}\n        {flag} = true,", 1)
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=edited.encode(),
        stdout=subprocess.PIPE,
        check=True,
    )
    path = directory / f"water_{flag}.bin"
    path.write_bytes(done.stdout)
    return path


@pytest.mark.parametrize("flag", ["externalPotential", "externalPotentialAdd"])
def test_deck_external_potential_is_its_own_error(flag, tmp_path):
    params = _deck_with_flag(flag, tmp_path)
    with Session(params) as session:
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(
            cpmdc.CPMDCDeckExternalPotential, match="EXTERNAL POTENTIAL"
        ) as caught:
            session.energy_forces(WATER, [8, 1, 1])
        assert caught.value.code == cpmdc.CPMDCDeckExternalPotential.code
        assert caught.value.code != cpmdc.CPMDCEvaluationFailed.code
        assert session._lib.cpmdc_last_error_code() == caught.value.code


class _AbiZeroLibrary:
    """A libcpmdc that predates electrostatic embedding.

    It forwards every symbol of the real library except the three that
    embedding added, so ``hasattr`` on those fails as it does on an older
    shared object.
    """

    HIDDEN = frozenset(
        {
            "cpmdc_session_set_external_charges",
            "cpmdc_session_last_external_charge_forces",
            "cpmdc_last_error_code",
        }
    )

    load = staticmethod(ctypes.CDLL)

    def __init__(self, path: str) -> None:
        self._real = self.load(path)

    def __getattr__(self, name: str):
        if name in self.HIDDEN:
            raise AttributeError(name)
        return getattr(self._real, name)


@pytest.fixture
def abi_zero_library(monkeypatch):
    monkeypatch.setattr(ctypes, "CDLL", _AbiZeroLibrary)


def test_older_library_still_evaluates_and_fails_with_plain_errors(
    params, abi_zero_library
):
    with Session(params) as session:
        for name in _AbiZeroLibrary.HIDDEN:
            assert not hasattr(session._lib, name)
        done = session.energy_forces(WATER, [8, 1, 1])
        assert np.isfinite(done.energy_h)

    with pytest.raises(CPMDCError) as caught:
        Session(b"not a CPMDParams message")
    assert type(caught.value) is CPMDCError
    assert caught.value.code == 0
    assert str(caught.value)


def test_older_library_refuses_external_charges_by_name(params, abi_zero_library):
    with Session(params) as session:
        with pytest.raises(CPMDCError, match="cpmdc_session_set_external_charges") as caught:
            session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        assert type(caught.value) is CPMDCError
        with pytest.raises(CPMDCError, match="cpmdc_session_set_external_charges"):
            session.set_external_charges([], [])
        with pytest.raises(CPMDCError, match="cpmdc_session_last_external_charge_forces"):
            session.external_charge_forces_au(1)


@pytest.mark.parametrize("width", [float("nan"), float("inf")])
def test_non_finite_gaussian_widths_are_refused_before_the_library(embedding_params, width):
    with Session(embedding_params) as session:
        with pytest.raises(ValueError, match="widths_angstrom must be finite"):
            session.set_external_charges([3.0, 0.2, -0.1], [0.7], [width])


@pytest.mark.parametrize("width", [0.0, -0.1, 1.0e-13])
def test_unusable_gaussian_widths_are_refused_and_keep_the_charges(embedding_params, width):
    with Session(embedding_params) as session:
        session.set_external_charges([3.0, 0.2, -0.1], [0.7], [0.30])
        before = session.energy_forces(WATER, [8, 1, 1])
        with pytest.raises(cpmdc.CPMDCInvalidArgument) as caught:
            session.set_external_charges([3.0, 0.2, -0.1], [0.7], [width])
        assert caught.value.code == cpmdc.CPMDCInvalidArgument.code
        assert "width" in str(caught.value)
        after = session.energy_forces(WATER, [8, 1, 1])
        assert after.energy_h == before.energy_h
        assert np.isfinite(after.energy_h)


def test_the_smallest_allowed_width_gives_a_finite_energy(embedding_params):
    with Session(embedding_params) as session:
        session.set_external_charges([3.0, 0.2, -0.1], [0.7], [1.0e-12])
        done = session.energy_forces(WATER, [8, 1, 1])
        forces = np.asarray(session.external_charge_forces_au(1))
    assert np.isfinite(done.energy_h)
    assert np.all(np.isfinite(forces))


def test_periodic_deck_is_refused_without_a_cell_argument(periodic_params):
    """A periodic deck is periodic whether or not the call passes a cell."""
    with Session(periodic_params) as session:
        bare = session.energy_forces(WATER, [8, 1, 1])
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(
            cpmdc.CPMDCUnsupportedPeriodicEmbedding, match="periodic electrostatic embedding"
        ) as caught:
            session.energy_forces(WATER, [8, 1, 1])
        assert caught.value.code == cpmdc.CPMDCUnsupportedPeriodicEmbedding.code
        session.set_external_charges([], [])
        assert session.energy_forces(WATER, [8, 1, 1]).energy_h == bare.energy_h


def test_cluster_deck_and_unset_periodicity_cell_rules(tmp_path):
    cluster_blob = _encode_params_text(
        tmp_path,
        """(
  functional = "BLYP",
  cutOffRy = 50.0,
  charge = 0,
  multiplicity = 1,
  task = "gradient",
  inputSections = [
    ( cpmd = ( optimizeWavefunction = true, convergenceOrbitals = 1.0e-5, maxStep = 120, maxIter = 120 ) ),
    ( dft = ( functional = "BLYP", lsd = false ) ),
    ( raw = "&SYSTEM\\n CLUSTER\\n POISSON SOLVER HOCKNEY\\n ANGSTROM\\n CELL\\n  10.0 1.0 1.0 0.0 0.0 0.0\\n CUTOFF\\n  50.0\\n&END\\n" ),
    ( atoms = ( pseudopotentials = [
        ( element = "O", path = "O_MT_BLYP.psp", lmax = 1 ),
        ( element = "H", path = "H_CVB_BLYP.psp", lmax = 0 )
      ] ) )
  ]
)""",
    ).read_bytes()
    unset_blob = _encode_params_text(
        tmp_path,
        """(
  functional = "BLYP",
  cutOffRy = 50.0,
  charge = 0,
  multiplicity = 1,
  task = "gradient",
  inputSections = [
    ( cpmd = ( optimizeWavefunction = true, convergenceOrbitals = 1.0e-5, maxStep = 120, maxIter = 120 ) ),
    ( dft = ( functional = "BLYP", lsd = false ) ),
    ( atoms = ( pseudopotentials = [
        ( element = "O", path = "O_MT_BLYP.psp", lmax = 1 ),
        ( element = "H", path = "H_CVB_BLYP.psp", lmax = 0 )
      ] ) )
  ]
)""",
    ).read_bytes()
    water_bohr = [value / cpmdc.BOHR_TO_ANGSTROM for row in WATER for value in row]
    cell_bohr = [20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0]

    with Session(cluster_blob) as cluster:
        cluster.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        energy_h, qm_forces, _ = cluster.calculate_au(water_bohr, [8, 1, 1], cell_bohr)
        assert math.isfinite(energy_h)
        assert all(math.isfinite(value) for value in qm_forces)
        assert all(math.isfinite(value) for value in cluster.external_charge_forces_au(1))

    with Session(unset_blob) as unset:
        unset.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        energy_zero, _, _ = unset.calculate_au(water_bohr, [8, 1, 1], [0.0] * 9)
        assert math.isfinite(energy_zero)
        with pytest.raises(
            cpmdc.CPMDCUnsupportedPeriodicEmbedding, match="periodic electrostatic embedding"
        ) as caught:
            unset.calculate_au(water_bohr, [8, 1, 1], cell_bohr)
        assert caught.value.code == cpmdc.CPMDCUnsupportedPeriodicEmbedding.code
