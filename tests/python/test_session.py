"""Session and Potential below the ASE layer."""

from __future__ import annotations

import shutil
import subprocess
from array import array
import ctypes
import math
import re

import mpmath
import pytest
from conftest import ROOT, WATER

import cpmdc
from cpmdc import CPMDCError, Potential, Session


def _force_input_text(flat, numbers, cell) -> str:
    lines = [
        "(",
        "  pos = [%s]," % ", ".join(f"{value:.17g}" for value in flat),
        "  atmnrs = [%s]," % ", ".join(str(value) for value in numbers),
    ]
    if cell is not None:
        lines.append("  box = [%s]," % ", ".join(f"{value:.17g}" for value in cell))
    lines += ['  lengthUnit = "angstrom",', '  energyUnit = "hartree",', ")"]
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize(
    ("atoms", "cell"),
    [
        (3, None),
        (3, [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0]),
        (1, [9.0, 0.0, 0.0, 1.2, 8.5, 0.0, -0.7, 0.9, 9.4]),
        (4, None),
    ],
)
def test_force_input_matches_capnp_encode(atoms, cell):
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is not on PATH")
    flat = array("d", [0.1 * i - 0.37 for i in range(3 * atoms)])
    numbers = array("i", [8, 1, 1, 6][:atoms])
    expected = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "ForceInput"],
        input=_force_input_text(flat, numbers, cell).encode(),
        stdout=subprocess.PIPE,
        check=True,
    ).stdout
    assert cpmdc._force_input_message(flat, numbers, cell) == expected


def test_periodic_step_needs_no_capnp(periodic_params, monkeypatch):
    # Open MPI's singleton startup otherwise looks up "ssh" on PATH before
    # entering the library, which would defeat this no-executables check.
    monkeypatch.setenv("OMPI_MCA_plm", "isolated")
    monkeypatch.setenv("PATH", "")
    with Potential(periodic_params) as pot:
        done = pot.evaluate(WATER, [8, 1, 1], [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]])
    assert len(done["stress"]) == 9


def test_closed_session_raises(params):
    session = Session(params)
    session.close()
    with pytest.raises(CPMDCError, match="closed"):
        session.energy_forces(WATER, [8, 1, 1])


def test_old_error_name_is_the_same_class():
    with pytest.warns(DeprecationWarning, match="CPMDCError"):
        old = cpmdc.CpmdcError
    assert old is CPMDCError
    assert "CpmdcError" not in cpmdc.__all__
    with pytest.raises(AttributeError):
        cpmdc.NoSuchName


def test_repeated_geometry_is_deterministic(params):
    with Session(params) as session:
        first = session.energy_forces(WATER, [8, 1, 1])
        again = session.energy_forces(WATER, [8, 1, 1])
    assert first.energy_h == again.energy_h
    assert first.forces_h_bohr == again.forces_h_bohr
    assert first.stress_h_bohr3 is None


# Relative part of the finite-difference comparison. pytest.approx
# allows max(rel * |expected|, abs).
_FORCE_FD_REL = 2.0e-6
# convergenceOrbitals in tests/water_blyp_embedding_params.capnp.txt.
# The analytic MM force is the Hellmann-Feynman derivative of the grid
# energy. A re-solved SCF difference still carries a remainder of this
# size: the orbital gradient is stopped here, and that remainder is
# first order in the gradient. The direct erf/r coefficient matches a
# fixed-density central difference, so this is not a kernel error.
_EMBED_ORBITAL_GRADIENT = 1.0e-7
_MACHINE_EPS = math.ulp(1.0)


def _finite_difference_allowance(plus, minus, step, opencpmd):
    """Roundoff of a central difference, plus the SCF residual.

    The subtracted energies each carry about eps*|E|. Dividing by
    twice the step (Angstrom) and converting to Ha/Bohr is the
    roundoff eps*|E|/h. Truncation is O(step^2 * d^3E). A step ten
    times larger moves this derivative by far less than the roundoff
    at `step`, so the h^2 term is not added again. On the linked
    build the Hellmann-Feynman remainder above dominates both.
    """
    roundoff = (
        _MACHINE_EPS
        * (abs(plus) + abs(minus))
        / (2.0 * step)
        * cpmdc.BOHR_TO_ANGSTROM
    )
    if opencpmd:
        return max(roundoff, _EMBED_ORBITAL_GRADIENT)
    return roundoff


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


def test_external_charge_translation_invariance(embedding_params):
    shift = [0.23, -0.17, 0.31]
    mm = [3.0, 0.2, -0.1]
    with Session(embedding_params) as session:
        if not session.available():
            pytest.skip("the reference well is anchored at the origin")
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


def test_density_dipole_responds_toward_external_field(embedding_params):
    with Session(embedding_params) as session:
        if session.available():
            pytest.skip("OpenCPMD does not export the induced dipole snapshot")
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
        assert snapshot.dipole[0] < 0.0
        assert abs(snapshot.dipole[0]) > abs(snapshot.dipole[1])


def test_periodic_external_charges_have_typed_error(embedding_params):
    with Session(embedding_params) as session:
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


def test_typed_failure_replaces_previous_error_code(embedding_params):
    with Session(embedding_params) as session:
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


def test_gaussian_small_r_series_is_continuous_and_conservative(embedding_params):
    sigma_bohr = 8.0e-12
    charge = sigma_bohr**2
    series_switch = sigma_bohr / 8.0
    old_width_guard = sigma_bohr * math.sqrt(math.sqrt(math.ulp(1.0)))

    with Session(embedding_params) as session:
        # _gaussian_reference is the reference well's induced dipole,
        # -1/2 alpha |field|^2. OpenCPMD returns a BLYP total energy, so
        # the same closed form does not describe that build and the
        # comparison cannot be moved onto it.
        if session.available():
            pytest.skip(
                "the closed-form oracle describes the reference backend, "
                "not an OpenCPMD SCF"
            )
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


def test_reference_smearing_width_is_converted_from_angstrom(embedding_params):
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
        if session.available():
            pytest.skip("the closed-form oracle describes the reference backend")
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
