"""Session and Potential below the ASE layer."""

from __future__ import annotations

import shutil
import subprocess
from array import array
import ctypes
import math

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


def _reference_session(params):
    session = Session(params)
    if session._lib.cpmdc_available():
        session.close()
        pytest.skip("analytic embedding oracle is a reference-backend test")
    return session


def test_external_charge_qm_and_mm_forces_are_energy_derivatives(params):
    width = [0.30]
    charge = [0.7]
    mm = [3.0, 0.2, -0.1]
    step = 1.0e-5
    with _reference_session(params) as session:
        session.set_external_charges(mm, charge, width)
        done = session.energy_forces(WATER, [8, 1, 1])
        moved = [list(row) for row in WATER]
        moved[0][0] += step
        plus = session.energy_forces(moved, [8, 1, 1]).energy_h
        moved[0][0] -= 2.0 * step
        minus = session.energy_forces(moved, [8, 1, 1]).energy_h
        qm_fd = -(plus - minus) / (2.0 * step) * cpmdc.BOHR_TO_ANGSTROM
        assert done.forces_h_bohr[0][0] == pytest.approx(qm_fd, rel=2e-5, abs=2e-8)

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
        assert mm_force[0] == pytest.approx(mm_fd, rel=2e-5, abs=2e-8)


def test_external_charge_limits_are_finite_and_recover_bare_energy(params):
    with _reference_session(params) as session:
        bare = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([1.0e6, 0.0, 0.0], [1.0], [0.30])
        far = session.energy_forces(WATER, [8, 1, 1]).energy_h
        assert far == pytest.approx(bare, abs=1e-12)
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [1.0e-6])
        tiny = session.energy_forces(WATER, [8, 1, 1])
        assert math.isfinite(tiny.energy_h)
        assert all(math.isfinite(value) for row in tiny.forces_h_bohr for value in row)


def test_clearing_external_charges_restores_bare_energy(params):
    with Session(params) as session:
        bare = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        embedded = session.energy_forces(WATER, [8, 1, 1]).energy_h
        session.set_external_charges([], [])
        cleared = session.energy_forces(WATER, [8, 1, 1]).energy_h
        assert abs(embedded - bare) > 1.0e-8
        assert cleared == pytest.approx(bare, abs=1.0e-8)


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


def test_density_dipole_responds_toward_external_field(params):
    with _reference_session(params) as session:
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


def test_periodic_external_charges_have_typed_error(params):
    with Session(params) as session:
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


def test_external_charge_force_reader_reports_each_failure(params):
    with Session(params) as session:
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
        bare = session.energy_forces(WATER, [8, 1, 1])
        assert math.isfinite(bare.energy_h)
        session.set_external_charges([3.0, 0.0, 0.0], [1.0], [0.30])
        with pytest.raises(
            cpmdc.CPMDCDeckExternalPotential, match="EXTERNAL POTENTIAL"
        ) as caught:
            session.energy_forces(WATER, [8, 1, 1])
        assert caught.value.code == cpmdc.CPMDCDeckExternalPotential.code
        assert caught.value.code != cpmdc.CPMDCEvaluationFailed.code
        assert session._lib.cpmdc_last_error_code() == caught.value.code
        session.set_external_charges([], [])
        cleared = session.energy_forces(WATER, [8, 1, 1])
        assert cleared.energy_h == pytest.approx(bare.energy_h)
