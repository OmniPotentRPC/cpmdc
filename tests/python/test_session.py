"""Session and Potential below the ASE layer."""

from __future__ import annotations

import ctypes
import shutil
import subprocess
from array import array

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


def test_orbital_key_methods(params):
    session = Session(params)
    session.select_orbitals(11)
    session.energy_forces(WATER, [8, 1, 1])
    session.copy_orbitals(11, 12)
    session.close()
    with pytest.raises(CPMDCError, match="closed"):
        session.select_orbitals(1)
    with pytest.raises(CPMDCError, match="closed"):
        session.copy_orbitals(1, 2)


def test_evaluation_has_a_readable_repr_and_keeps_identity_semantics():
    first = cpmdc.Evaluation(1.5, [[0.0, 0.0, 0.0]], None)
    second = cpmdc.Evaluation(1.5, [[0.0, 0.0, 0.0]], None)
    assert repr(first) == (
        "Evaluation(energy_h=1.5, forces_h_bohr=[[0.0, 0.0, 0.0]], stress_h_bohr3=None)"
    )
    assert first != second
    assert len({first, second}) == 2
    assert cpmdc.Evaluation(energy_h=2.0, forces_h_bohr=[], stress_h_bohr3=[0.0] * 9).energy_h == 2.0


def test_session_rejects_bad_inputs(params, tmp_path, monkeypatch):
    import numpy as np

    with pytest.raises(CPMDCError, match="empty"):
        Session(b"")
    with pytest.raises(CPMDCError, match="parse failed"):
        Session(b"not a CPMDParams message")
    bogus = tmp_path / "notlib.so"
    bogus.write_bytes(b"not a shared library")
    with pytest.raises(CPMDCError, match="failed to load"):
        Session(params, library=bogus)
    with pytest.raises(CPMDCError, match="library not found"):
        cpmdc.find_library(tmp_path / "missing.so")
    with pytest.raises(CPMDCError, match="schema not found"):
        cpmdc.find_schema(tmp_path / "missing.capnp")

    real = cpmdc.find_library()
    monkeypatch.setenv("CPMDC_LIBRARY", str(tmp_path / "missing.so"))
    monkeypatch.setenv("RGPOT_CPMDC_ENGINE", str(real))
    assert cpmdc.find_library() == real
    monkeypatch.delenv("CPMDC_LIBRARY")
    monkeypatch.delenv("RGPOT_CPMDC_ENGINE")
    monkeypatch.delenv("RGPOT_CPMD_ENGINE", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(CPMDCError, match="library not found"):
        cpmdc.find_library()
    monkeypatch.undo()

    with pytest.raises(ValueError, match="9 components"):
        cpmdc.stress_to_ev_per_angstrom3([0.0] * 8)

    session = Session(params)
    assert isinstance(session.abi_version(), int)

    with pytest.raises(ValueError, match="finite"):
        session.energy_forces([[float("nan"), 0.0, 0.0]], [1])
    with pytest.raises(ValueError, match="at least one atom"):
        session.energy_forces([], [])
    with pytest.raises(ValueError, match="three coordinates"):
        session.energy_forces([[0.0, 0.0]], [1])
    with pytest.raises(ValueError, match="one entry per atom"):
        session.energy_forces(WATER, [8, 1])
    with pytest.raises(ValueError, match="3x3 or 9"):
        session.energy_forces(WATER, [8, 1, 1], [1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="3x3 or 9"):
        session.energy_forces(
            WATER, [8, 1, 1], [[1.0, 0.0, 0.0], [0.0, 1.0], [0.0, 0.0, 1.0]]
        )
    with pytest.raises(ValueError, match="finite"):
        session.energy_forces(
            WATER,
            [8, 1, 1],
            [[1.0, 0.0, 0.0], [0.0, float("nan"), 0.0], [0.0, 0.0, 1.0]],
        )

    bohr = [[coord / cpmdc.BOHR_TO_ANGSTROM for coord in row] for row in WATER]
    rows = [tuple(row) for row in bohr]
    energy, forces, virial = session.calculate_au(rows, [8, 1, 1])
    done = session.energy_forces(WATER, [8, 1, 1])
    assert energy == done.energy_h
    assert len(list(forces)) == 9
    assert len(list(virial)) == 9
    readonly = np.zeros(9, dtype=np.float64)
    readonly.setflags(write=False)
    with pytest.raises(ValueError, match="writable"):
        session.calculate_au(rows, [8, 1, 1], forces=readonly)
    with pytest.raises(ValueError, match="positions_bohr"):
        session.calculate_au([0.0, 0.0], [8, 1, 1])
    with pytest.raises(ValueError, match="at least one atom"):
        session.calculate_au([], [])
    session.close()
    with pytest.raises(CPMDCError, match="closed"):
        session.calculate_au(rows, [8, 1, 1])


def test_available_matches_the_build(params):
    library = ctypes.CDLL(str(cpmdc.find_library()))
    library.cpmdc_available.restype = ctypes.c_int
    with Session(params) as session:
        assert session.available() is bool(library.cpmdc_available())
    with Potential(params) as pot:
        assert pot.caps() == {"stress": True, "energies": False}
