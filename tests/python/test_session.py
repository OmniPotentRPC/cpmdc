"""Session and Potential below the ASE layer."""

from __future__ import annotations

import ctypes
import shutil
import subprocess
import sys
from array import array

import pytest
from conftest import ROOT, WATER

import cpmdc
from cpmdc import BOHR_TO_ANGSTROM, CPMDCError, Potential, Session, dipole_to_e_angstrom
from cpmdc._rank import _cell_volume


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


def test_session_preserves_an_accepted_small_triclinic_cell(periodic_params):
    """Tilt components far below the lengths ride along with a physical cell.

    The cell passes the rank guard with a tilt of 2.58e-12 and 1e-13 Angstrom
    beside 6 Angstrom edges. It must give the energy, forces and stress of
    the same cell with the tilts set to zero, on the reference evaluator and
    on OpenCPMD alike.
    """
    tilted = [[6.0, 0.0, 0.0], [2.58e-12, 6.0, 0.0], [1.0e-13, 0.0, 6.0]]
    straight = [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]
    with Session(periodic_params) as session:
        done = session.energy_forces(WATER, [8, 1, 1], tilted)
    with Session(periodic_params) as session:
        base = session.energy_forces(WATER, [8, 1, 1], straight)
    assert done.stress_h_bohr3 is not None
    assert done.energy_h == pytest.approx(base.energy_h, abs=1e-9)
    flat = [value for row in done.forces_h_bohr for value in row]
    flat_base = [value for row in base.forces_h_bohr for value in row]
    assert flat == pytest.approx(flat_base, abs=1e-9)
    assert done.stress_h_bohr3 == pytest.approx(base.stress_h_bohr3, abs=1e-9)


def test_session_accepts_a_cell_thinner_than_any_electronic_structure(periodic_params, opencpmd):
    """The rank guard takes this cell, so only the backend can judge it.

    The third vector is 1e-13 Angstrom long. The reference evaluator
    returns a stress for it. OpenCPMD is handed the cell as written, with
    the lattice vectors, reciprocal vectors and the 8 by 8 by 4 mesh the
    cell implies, and its plane wave basis shrinks to 5 waves for 4 orbitals
    at a density near 1e13 electrons per cubic bohr. The first electronic
    gradient is already 1e13 on the first step of a cold start and the SCF
    cannot converge. The error comes from CPMD, not from the cell guard.
    """
    cell = [[1.0, 0.0, 0.0], [2.58e-12, 1.0, 0.0], [1.0e-13, 0.0, 1.0e-13]]
    with Session(periodic_params) as session:
        if not opencpmd:
            done = session.energy_forces(WATER, [8, 1, 1], cell)
            assert done.stress_h_bohr3 is not None
            return
        with pytest.raises(CPMDCError, match=r"CPMD (SCF did not converge|stopgm)") as failure:
            session.energy_forces(WATER, [8, 1, 1], cell)
    assert "independent vectors" not in str(failure.value)
    assert "too small" not in str(failure.value)


def _water_bohr() -> list[float]:
    return [coord / BOHR_TO_ANGSTROM for row in WATER for coord in row]


def test_virial_uses_the_scaled_lu_volume(periodic_params):
    """The virial is that volume times the stress, bit for bit.

    For this cell the unscaled three-term determinant is one unit in the
    last place away from the LU volume, so the two are told apart.
    """
    cell = [18.7, 0.2, 0.0, 0.3, 19.3, 0.4, 0.7, -1.1, 18.1]
    rows = [cell[3 * k : 3 * k + 3] for k in range(3)]
    volume = _cell_volume(rows)
    three_term = abs(
        cell[0] * (cell[4] * cell[8] - cell[5] * cell[7])
        - cell[1] * (cell[3] * cell[8] - cell[5] * cell[6])
        + cell[2] * (cell[3] * cell[7] - cell[4] * cell[6])
    )
    assert three_term != volume
    positions = _water_bohr()
    with Session(periodic_params) as session:
        _, _, virial = session.calculate_au(positions, [8, 1, 1], cell)
        done = session.energy_forces(
            [[value * BOHR_TO_ANGSTROM for value in positions[3 * i : 3 * i + 3]] for i in range(3)],
            [8, 1, 1],
            [[value * BOHR_TO_ANGSTROM for value in row] for row in rows],
        )
    assert done.stress_h_bohr3 is not None and any(done.stress_h_bohr3)
    assert list(virial) == [volume * component for component in done.stress_h_bohr3]


def test_near_singular_cancelling_cell_gets_stress_and_virial(periodic_params, opencpmd):
    if opencpmd:
        pytest.skip("uses the reference evaluator's stress formula")
    near = 1.0 + 28.0 * sys.float_info.epsilon
    cell_bohr = [1.0, 1.0, 1.0, 1.0, near, 1.0, 1.0, 1.0, near]
    rows_bohr = [cell_bohr[3 * k : 3 * k + 3] for k in range(3)]
    three_term = abs(
        cell_bohr[0] * (cell_bohr[4] * cell_bohr[8] - cell_bohr[5] * cell_bohr[7])
        - cell_bohr[1] * (cell_bohr[3] * cell_bohr[8] - cell_bohr[5] * cell_bohr[6])
        + cell_bohr[2] * (cell_bohr[3] * cell_bohr[7] - cell_bohr[4] * cell_bohr[6])
    )
    assert three_term == 0.0
    volume = _cell_volume(rows_bohr)
    assert volume > 0.0
    rows_ang = [[value * BOHR_TO_ANGSTROM for value in row] for row in rows_bohr]
    with Session(periodic_params) as session:
        _, _, virial = session.calculate_au(_water_bohr(), [8, 1, 1], cell_bohr)
        done = session.energy_forces(WATER, [8, 1, 1], rows_ang)
    assert done.stress_h_bohr3 is not None and any(done.stress_h_bohr3)
    assert any(virial)
    assert list(virial) == [volume * component for component in done.stress_h_bohr3]


def test_overflowing_periodic_cell_fails_instead_of_dropping_stress(periodic_params, opencpmd):
    if opencpmd:
        pytest.skip("uses the reference evaluator's stress formula")
    side = sys.float_info.max / 2.0
    cell_rows = [[side, 0.0, 0.0], [0.0, side, 0.0], [0.0, 0.0, side]]
    with Session(periodic_params) as session:
        with pytest.raises(CPMDCError, match="cell volume is too large for binary64"):
            session.energy_forces(WATER, [8, 1, 1], cell_rows)
        with pytest.raises(CPMDCError, match="cell volume is too large for binary64"):
            session.calculate_au(
                _water_bohr(),
                [8, 1, 1],
                [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side],
                virial=None,
            )


def test_virial_volume_above_the_largest_float_is_refused(periodic_params):
    side = sys.float_info.max / 2.0
    with Session(periodic_params) as session:
        with pytest.raises(CPMDCError, match="cell volume is too large for binary64"):
            session.calculate_au(
                _water_bohr(), [8, 1, 1], [side, 0.0, 0.0, 0.0, side, 0.0, 0.0, 0.0, side]
            )


def test_virial_that_overflows_is_refused(periodic_params, opencpmd):
    if opencpmd:
        pytest.skip("uses the reference evaluator's stress formula")
    positions = list(_water_bohr())
    positions[0] = 2.0**520
    with Session(periodic_params) as session:
        with pytest.raises(CPMDCError, match="virial is not finite"):
            session.calculate_au(
                positions, [8, 1, 1], [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0]
            )


def test_reference_step_has_no_fabricated_dipole(params):
    with Session(params) as session:
        done = session.energy_forces(WATER, [8, 1, 1])
    assert done.dipole_e_bohr is None
    with Potential(params) as pot:
        result = pot.evaluate(WATER, [8, 1, 1])
    assert "dipole" not in result
    assert pot.caps() == {"stress": True, "dipole": False, "energies": False}


def test_dipole_unit_conversion():
    assert dipole_to_e_angstrom([1.0, -2.0, 0.5]) == pytest.approx(
        [BOHR_TO_ANGSTROM, -2.0 * BOHR_TO_ANGSTROM, 0.5 * BOHR_TO_ANGSTROM]
    )
    with pytest.raises(ValueError, match="3 components"):
        dipole_to_e_angstrom([1.0, 2.0])


def test_property_gradient_limit_fails_without_truncation(params):
    positions = [[0.0, 0.0, 0.0]] * 1366
    with Session(params) as session:
        with pytest.raises(CPMDCError, match="at most 1365 atoms"):
            session.energy_forces(positions, [1] * 1366)


def _encode_params(text: str) -> bytes:
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is not on PATH")
    return subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode("utf-8"),
        stdout=subprocess.PIPE,
        check=True,
    ).stdout


def _periodic_source() -> str:
    return (ROOT / "tests" / "water_blyp_periodic_params.capnp.txt").read_text(encoding="utf-8")


GAMMA_KPOINT = "kpoints = [( coordinates = [0.0, 0.0, 0.0], weight = 1.0 )],"


def test_dipole_dynamics_with_kpoints_is_rejected_up_front():
    source = _periodic_source().replace("maxIter = 40,", "maxIter = 40, dipoleDynamics = true,")
    source = source.replace("symmetry = 1,", f"symmetry = 1, {GAMMA_KPOINT}")
    assert "dipoleDynamics" in source and "kpoints" in source
    with pytest.raises(CPMDCError, match=r"DIPOLE DYNAMICS.*KPOINTS.*Gamma-point"):
        Session(_encode_params(source))
    with pytest.raises(CPMDCError, match=r"DIPOLE DYNAMICS.*KPOINTS"):
        Potential(_encode_params(source))


def test_dipole_dynamics_with_kpoints_in_a_raw_block_is_rejected():
    raw = (
        "&CPMD\\n OPTIMIZE WAVEFUNCTION\\n DIPOLE DYNAMICS\\n&END\\n"
        "&SYSTEM\\n SYMMETRY\\n  1\\n KPOINTS\\n 1\\n 0.0 0.0 0.0 1.0\\n&END\\n"
    )
    source = f'(functional = "BLYP", cutOffRy = 40.0, task = "gradient", inputBlocks = ["{raw}"])'
    with pytest.raises(CPMDCError, match=r"DIPOLE DYNAMICS.*KPOINTS"):
        Session(_encode_params(source))


def _raw_params(block: str) -> bytes:
    source = (
        f'(functional = "BLYP", cutOffRy = 40.0, task = "gradient", '
        f'inputBlocks = ["{block}"])'
    )
    return _encode_params(source)


def test_a_commented_dipole_dynamics_line_is_not_a_kpoint_conflict():
    commented = (
        "&CPMD\\n OPTIMIZE WAVEFUNCTION\\n ! DIPOLE DYNAMICS\\n&END\\n"
        "&SYSTEM\\n SYMMETRY\\n  1\\n KPOINTS\\n 1\\n 0.0 0.0 0.0 1.0\\n&END\\n"
    )
    trailing = (
        "&CPMD\\n DIPOLE DYNAMICS ! berry\\n&END\\n"
        "&SYSTEM\\n SYMMETRY\\n  1\\n KPOINTS ! mesh\\n 1\\n"
        " 0.0 0.0 0.0 1.0\\n&END\\n"
    )
    with Session(_raw_params(commented)):
        pass
    with pytest.raises(CPMDCError, match=r"DIPOLE DYNAMICS.*KPOINTS"):
        Session(_raw_params(trailing))


def test_kpoints_and_dipole_dynamics_alone_are_accepted():
    kpoints = _periodic_source().replace("symmetry = 1,", f"symmetry = 1, {GAMMA_KPOINT}")
    dipole = _periodic_source().replace("maxIter = 40,", "maxIter = 40, dipoleDynamics = true,")
    for source in (kpoints, dipole):
        with Session(_encode_params(source)):
            pass


def test_topology_change_is_not_reported_as_the_gradient_limit(params):
    big = [[0.0, 0.0, 0.0]] * 1366
    with Session(params) as session:
        session.energy_forces(WATER, [8, 1, 1])
        with pytest.raises(CPMDCError, match="topology change"):
            session.energy_forces(big, [1] * 1366)
        again = session.energy_forces(WATER, [8, 1, 1])
    assert again.forces_h_bohr is not None


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


def test_params_method_reads_effective_scalars(params):
    method = cpmdc.params_method(params)
    assert method == cpmdc.MethodInfo(
        functional="BLYP", charge=0, multiplicity=1
    )


def _encode_params_text(tmp_path, text):
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip("capnp is not on PATH")
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode(),
        stdout=subprocess.PIPE,
        check=True,
    )
    path = tmp_path / "params.bin"
    path.write_bytes(done.stdout)
    return path


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (
            """(
  functional = "BLYP",
  inputSections = [ ( raw = "&DFT\\n FUNCTIONAL\\n  PBE\\n&END\\n" ) ]
)""",
            "raw section opens &DFT",
        ),
        (
            """(
  functional = "BLYP",
  inputSections = [ ( generic = ( name = "SYSTEM", directives = [ ( keyword = "CHARGE", args = ["2"] ) ] ) ) ]
)""",
            "generic section opens &SYSTEM",
        ),
        (
            """(
  functional = "BLYP",
  inputSections = [ ( set = ( key = "SYSTEM.CHARGE", value = "2" ) ) ]
)""",
            "SYSTEM scalar set directive",
        ),
    ],
    ids=["raw-dft", "generic-system", "set-charge"],
)
def test_params_method_rejects_opaque_scalar_sections(tmp_path, text, message):
    params = _encode_params_text(tmp_path, text)
    with pytest.raises(CPMDCError, match=message):
        cpmdc.params_method(params)


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
        "Evaluation(energy_h=1.5, forces_h_bohr=[[0.0, 0.0, 0.0]], "
        "stress_h_bohr3=None, dipole_e_bohr=None)"
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
    # The assertion below is the empty search. A wheel still ships
    # libcpmdc under _WHEEL_LIBS, and ctypes.util.find_library can return
    # that soname. Neither location is part of this check.
    empty_libs = tmp_path / "wheel-libs"
    empty_libs.mkdir()
    monkeypatch.setattr(cpmdc, "_WHEEL_LIBS", empty_libs)
    monkeypatch.setattr(ctypes.util, "find_library", lambda _name: None)
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
        assert pot.caps() == {
            "stress": True,
            "dipole": bool(library.cpmdc_available()),
            "energies": False,
        }


def test_load_rejects_library_without_available_symbol(monkeypatch):
    class NoAvailableProxy:
        def __init__(self, inner):
            self._inner = inner

        def __getattr__(self, name):
            if name == "cpmdc_available":
                raise AttributeError(name)
            return getattr(self._inner, name)

    original = ctypes.CDLL(str(cpmdc.find_library()))
    monkeypatch.setattr(cpmdc.ctypes, "CDLL", lambda _path: NoAvailableProxy(original))
    with pytest.raises(CPMDCError, match="has no cpmdc_available symbol"):
        cpmdc._load(cpmdc.find_library())


def test_older_library_without_params_method_loads_and_reports_unavailable(
    params, monkeypatch
):
    class NoParamsMethodProxy:
        def __init__(self, inner):
            self._inner = inner

        def __getattr__(self, name):
            if name == "cpmdc_params_method":
                raise AttributeError(name)
            return getattr(self._inner, name)

    original = ctypes.CDLL(str(cpmdc.find_library()))
    monkeypatch.setattr(
        cpmdc.ctypes, "CDLL", lambda _path: NoParamsMethodProxy(original)
    )
    loaded = cpmdc._load(cpmdc.find_library())
    assert not hasattr(loaded, "cpmdc_params_method")
    with Session(params) as session:
        assert session.abi_version() == 0
    with pytest.raises(CPMDCError, match="no cpmdc_params_method"):
        cpmdc.params_method(params)
