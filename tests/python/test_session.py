"""Session and Potential below the ASE layer."""

from __future__ import annotations

import contextlib
import ctypes
import math
import re
import shutil
import struct
import subprocess
import sys
from array import array
from collections.abc import Iterator

import mpmath
import numpy as np
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
    # Open MPI's singleton startup otherwise looks up "ssh" on PATH before
    # entering the library, which would defeat this no-executables check.
    monkeypatch.setenv("OMPI_MCA_plm", "isolated")
    monkeypatch.setenv("PATH", "")
    with Potential(periodic_params) as pot:
        done = pot.evaluate(WATER, [8, 1, 1], [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]])
    assert len(done["stress"]) == 9


def _decode_force_input_box(raw: bytes) -> list[float]:
    """Decode ForceInput.box directly from Cap'n Proto single-segment bytes."""
    seg_count_minus_one, seg_words = struct.unpack_from("<II", raw, 0)
    assert seg_count_minus_one == 0
    assert len(raw) == 8 + 8 * seg_words
    root_word = struct.unpack_from("<Q", raw, 8)[0]
    assert (root_word & 0x3) == 0
    struct_offset = (root_word >> 2) & 0x3FFFFFFF
    if struct_offset & 0x20000000:
        struct_offset -= 0x40000000
    data_words = (root_word >> 32) & 0xFFFF
    ptr_words = (root_word >> 48) & 0xFFFF
    assert ptr_words >= 3
    box_ptr_byte = 8 + 8 * (1 + struct_offset + data_words + 2)
    box_ptr = struct.unpack_from("<Q", raw, box_ptr_byte)[0]
    assert (box_ptr & 0x3) == 1
    list_offset = (box_ptr >> 2) & 0x3FFFFFFF
    if list_offset & 0x20000000:
        list_offset -= 0x40000000
    elem_size = (box_ptr >> 32) & 0x7
    elem_count = (box_ptr >> 35) & 0x1FFFFFFF
    assert elem_size == 5
    assert elem_count == 9
    target_byte = box_ptr_byte + 8 + 8 * list_offset
    return list(struct.unpack_from("<9d", raw, target_byte))


@contextlib.contextmanager
def _captured_force_input_boxes(session: Session) -> Iterator[list[list[float]]]:
    boxes: list[list[float]] = []
    original = session._lib.cpmdc_session_calculate_forces

    def capture(handle, payload, length, forces, capacity):
        boxes.append(_decode_force_input_box(ctypes.string_at(payload, length)))
        return original(handle, payload, length, forces, capacity)

    session._lib.cpmdc_session_calculate_forces = capture
    try:
        yield boxes
    finally:
        session._lib.cpmdc_session_calculate_forces = original


def test_session_preserves_an_accepted_small_triclinic_cell(periodic_params):
    """Tilt components far below the lengths ride along with a physical cell.

    The cell passes the rank guard with a tilt of 2.58e-12 and 1e-13 Angstrom
    beside 6 Angstrom edges. It must give the energy, forces and stress of
    the same cell with the tilts set to zero, on the reference evaluator and
    on OpenCPMD alike, while still keeping the off-diagonal tilts in ForceInput.
    """
    tilted = [[6.0, 0.0, 0.0], [2.58e-12, 6.0, 0.0], [1.0e-13, 0.0, 6.0]]
    straight = [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]
    assert tilted != straight
    assert tilted[1][0] != 0.0 and tilted[2][0] != 0.0
    with Session(periodic_params) as session:
        with _captured_force_input_boxes(session) as boxes:
            done = session.energy_forces(WATER, [8, 1, 1], tilted)
    assert boxes == [[6.0, 0.0, 0.0, 2.58e-12, 6.0, 0.0, 1.0e-13, 0.0, 6.0]]
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
    assert cell[1][0] != 0.0 and cell[2][0] != 0.0
    with Session(periodic_params) as session:
        with _captured_force_input_boxes(session) as boxes:
            if not opencpmd:
                done = session.energy_forces(WATER, [8, 1, 1], cell)
                assert boxes == [[1.0, 0.0, 0.0, 2.58e-12, 1.0, 0.0, 1.0e-13, 0.0, 1.0e-13]]
                assert done.stress_h_bohr3 is not None
                return
            with pytest.raises(CPMDCError, match=r"CPMD (SCF did not converge|stopgm)") as failure:
                session.energy_forces(WATER, [8, 1, 1], cell)
    assert boxes == [[1.0, 0.0, 0.0, 2.58e-12, 1.0, 0.0, 1.0e-13, 0.0, 1.0e-13]]
    assert "independent vectors" not in str(failure.value)
    assert "too small" not in str(failure.value)


def test_session_triclinic_tilt_oracle_and_rejects_transpose(periodic_params, opencpmd):
    """A 1 Angstrom triclinic tilt matches an independent numpy oracle and fails when transposed or zeroed."""
    # Expected lattice vectors (rows in Angstrom) computed independently of cpmdc.
    oracle_vectors = np.array(
        [
            [6.5, 0.0, 0.0],
            [1.2, 6.2, 0.0],
            [-0.9, 0.8, 6.4],
        ],
        dtype=np.float64,
    )
    tilted = [[6.5, 0.0, 0.0], [1.2, 6.2, 0.0], [-0.9, 0.8, 6.4]]
    straight = [[6.5, 0.0, 0.0], [0.0, 6.2, 0.0], [0.0, 0.0, 6.4]]
    transposed_vectors = oracle_vectors.T
    transposed = transposed_vectors.tolist()

    # Hand-computed volume of the lower-triangular cell: 6.5 * 6.2 * 6.4 = 257.92 A^3.
    hand_volume_ang3 = 6.5 * 6.2 * 6.4
    oracle_volume_ang3 = abs(
        float(np.dot(oracle_vectors[0], np.cross(oracle_vectors[1], oracle_vectors[2])))
    )
    assert oracle_volume_ang3 == pytest.approx(hand_volume_ang3, rel=1e-14, abs=1e-14)
    assert abs(float(np.linalg.det(oracle_vectors))) == pytest.approx(
        hand_volume_ang3, rel=1e-14, abs=1e-14
    )

    # Independent fractional <-> Cartesian mapping in plain numpy.
    water_np = np.array(WATER, dtype=np.float64)
    inv_vectors = np.linalg.inv(oracle_vectors)
    frac_water = water_np @ inv_vectors
    np.testing.assert_allclose(frac_water @ oracle_vectors, water_np, rtol=1e-14, atol=1e-15)

    shifted_frac = frac_water.copy()
    shifted_frac[1] += np.array([0.0, 1.0, 0.0])
    shifted_cart = shifted_frac @ oracle_vectors
    np.testing.assert_allclose(
        shifted_cart[1] - water_np[1], np.array([1.2, 6.2, 0.0]), rtol=1e-14, atol=1e-15
    )
    assert not np.allclose(shifted_frac @ transposed_vectors, shifted_cart, atol=1e-3)

    with Session(periodic_params) as session:
        with _captured_force_input_boxes(session) as boxes:
            done = session.energy_forces(WATER, [8, 1, 1], tilted)
            shifted_done = session.energy_forces(shifted_cart.tolist(), [8, 1, 1], tilted)
            straight_done = session.energy_forces(WATER, [8, 1, 1], straight)
            trans_done = session.energy_forces(WATER, [8, 1, 1], transposed)
            trans_shifted = session.energy_forces(
                shifted_cart.tolist(), [8, 1, 1], transposed
            )
            no_cell_done = session.energy_forces(WATER, [8, 1, 1], None)

        assert done.stress_h_bohr3 is not None
        assert trans_done.stress_h_bohr3 is not None
        assert straight_done.stress_h_bohr3 is not None

        if opencpmd:
            # Periodic-image invariance holds only under the un-transposed, tilted cell.
            assert shifted_done.energy_h == pytest.approx(done.energy_h, abs=1e-9)
            np.testing.assert_allclose(
                shifted_done.forces_h_bohr, done.forces_h_bohr, rtol=1e-9, atol=1e-11
            )
            np.testing.assert_allclose(
                shifted_done.stress_h_bohr3, done.stress_h_bohr3, rtol=1e-9, atol=1e-11
            )

            # Zeroing the off-diagonal tilts changes both energy and stress significantly.
            assert abs(done.energy_h - straight_done.energy_h) > 1e-4
            assert (
                np.max(
                    np.abs(np.array(done.stress_h_bohr3) - np.array(straight_done.stress_h_bohr3))
                )
                > 1e-6
            )

            # A transposed cell reaching the session fails cell-orientation and lattice-shift checks.
            assert abs(trans_done.energy_h - done.energy_h) > 1e-4
            assert (
                np.max(
                    np.abs(np.array(trans_done.stress_h_bohr3) - np.array(done.stress_h_bohr3))
                )
                > 1e-6
            )
            assert abs(trans_shifted.energy_h - trans_done.energy_h) > 1e-2
        else:
            # Independent plain-numpy oracle for the reference evaluator (CODATA 2018 Bohr).
            bohr_ang = 0.529177210903
            vol_bohr3 = hand_volume_ang3 / bohr_ang**3
            pos_bohr = water_np / bohr_ang
            atomic_numbers = np.array([8, 1, 1], dtype=np.float64)
            k = 1.0e-3 * max(0.1, 40.0 / 70.0)
            oracle_cell_energy = 1.0e-10 * float(np.sum(oracle_vectors**2))
            oracle_forces = -(k * atomic_numbers[:, None] * pos_bohr)
            ded_atoms = k * (pos_bohr.T @ (atomic_numbers[:, None] * pos_bohr))
            ded_cell = 2.0e-10 * (oracle_vectors.T @ oracle_vectors)
            oracle_stress = -(ded_atoms + ded_cell) / vol_bohr3
            oracle_virial = vol_bohr3 * oracle_stress

            assert (done.energy_h - no_cell_done.energy_h) == pytest.approx(
                oracle_cell_energy, rel=1e-13, abs=1e-15
            )
            np.testing.assert_allclose(done.forces_h_bohr, oracle_forces, rtol=1e-13, atol=1e-15)
            np.testing.assert_allclose(
                done.stress_h_bohr3, oracle_stress.ravel(), rtol=1e-12, atol=1e-15
            )

            au_energy, au_forces, au_virial = session.calculate_au(
                pos_bohr.ravel().tolist(), [8, 1, 1], (oracle_vectors / bohr_ang).ravel().tolist()
            )
            assert (au_energy - no_cell_done.energy_h) == pytest.approx(
                oracle_cell_energy, rel=1e-13, abs=1e-15
            )
            np.testing.assert_allclose(
                list(au_forces), oracle_forces.ravel(), rtol=1e-13, atol=1e-15
            )
            np.testing.assert_allclose(
                list(au_virial), oracle_virial.ravel(), rtol=1e-12, atol=1e-15
            )

            # Zeroing off-diagonals changes both reference energy and stress beyond oracle tolerance.
            assert abs(done.energy_h - straight_done.energy_h) > 1e-11
            assert (
                np.max(
                    np.abs(np.array(done.stress_h_bohr3) - np.array(straight_done.stress_h_bohr3))
                )
                > 1e-13
            )

            # A transposed cell reaching the session fails the independent stress oracle.
            with pytest.raises(AssertionError):
                np.testing.assert_allclose(
                    trans_done.stress_h_bohr3, oracle_stress.ravel(), rtol=1e-12, atol=1e-15
                )

        assert boxes[0] == oracle_vectors.ravel().tolist()
        assert boxes[0] != transposed_vectors.ravel().tolist()
        assert boxes[2] == [6.5, 0.0, 0.0, 0.0, 6.2, 0.0, 0.0, 0.0, 6.4]
        assert boxes[3] == transposed_vectors.ravel().tolist()

    with Potential(periodic_params) as pot:
        pot_res = pot.evaluate(WATER, [8, 1, 1], tilted)
    assert pot_res["energy"] == pytest.approx(done.energy_h * cpmdc.HARTREE_TO_EV, rel=1e-12)
    np.testing.assert_allclose(
        pot_res["stress"],
        [-value for value in cpmdc.stress_to_ev_per_angstrom3(done.stress_h_bohr3)],
        rtol=1e-12,
        atol=1e-15,
    )


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
    with pytest.raises(CPMDCError, match=message) as caught:
        cpmdc.params_method(params)
    assert caught.value.code == cpmdc.CPMDCInvalidArgument.code


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


def test_reference_embedding_step_has_no_fabricated_dipole(embedding_params):
    """The reference evaluator exports no dipole, with or without charges."""
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
        assert snapshot.dipole_count == 0
        assert list(snapshot.dipole) == [0.0, 0.0, 0.0]


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
    params, abi_zero_library, tmp_path
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
