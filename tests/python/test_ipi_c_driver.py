"""cpmdc-ipi-c, the C i-PI client, against ASE, i-PI, a minimal server, and cpmdc-ipi.

Meson builds the client from ``tools/cpmdc_ipi.c`` on
``cpmdc_session_calculate_au()``. ``CPMDC_IPI_C`` names it; otherwise
``build/cpmdc-ipi-c`` of this checkout is used, and the tests skip when
neither exists. Every step is compared with the Python driver, cpmdc-ipi,
which ``test_ipi_driver.py`` checks against the direct calculator.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.socketio import SocketIOCalculator
from conftest import ROOT, WATER
from test_ipi_driver import (  # noqa: F401  (fixtures are used by name)
    CUBE,
    IPI_INPUT,
    TIMEOUT,
    TRICLINIC,
    Client,
    MiniServer,
    codata2018,
    free_port,
    mini,
    native,
    socket_dir,
    water,
)

from cpmdc import BOHR_TO_ANGSTROM
from cpmdc.ipi import Driver


def c_client_path() -> Path | None:
    named = os.environ.get("CPMDC_IPI_C")
    candidate = Path(named) if named else ROOT / "build" / "cpmdc-ipi-c"
    return candidate if candidate.is_file() else None


C_CLIENT = c_client_path()
pytestmark = pytest.mark.skipif(C_CLIENT is None, reason="cpmdc-ipi-c is not built")


class CClient(Client):
    """``cpmdc-ipi-c`` in a subprocess."""

    def __init__(self, *args: str) -> None:
        self.proc = subprocess.Popen(
            [str(C_CLIENT), *args],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )


@pytest.fixture
def drivers():
    """Start a client: ``drivers("c", ...)`` or ``drivers("python", ...)``."""
    started: list[Client] = []

    def start(kind: str, *args: str) -> Client:
        client = CClient(*args) if kind == "c" else Client(*args)
        started.append(client)
        return client

    yield start
    for client in started:
        client.kill()


def socket_pair(kind: str, drivers, params: Path):
    """One SocketIOCalculator per driver, C first, over UNIX or TCP."""
    calcs = []
    for driver in ("c", "python"):
        if kind == "unix":
            name = "cpmdc-" + uuid.uuid4().hex[:12]
            calc = SocketIOCalculator(unixsocket=name, timeout=TIMEOUT)
            drivers(driver, "--unix", name, "--params", str(params), "--species", "O H H")
        else:
            port = free_port()
            calc = SocketIOCalculator(port=port, timeout=TIMEOUT)
            drivers(driver, "--host", "localhost", "--port", str(port), "--params",
                    str(params), "--species", "OH2")
        calcs.append(calc)
    return calcs


def assert_close(left, right) -> None:
    np.testing.assert_allclose(left, right, rtol=1e-10, atol=1e-10)


# ASE SocketIOCalculator, both drivers


@pytest.mark.parametrize("kind", ["unix", "inet"])
def test_ase_cluster_matches_python_driver(kind, params, drivers, codata2018):
    c_calc, py_calc = socket_pair(kind, drivers, params)
    with c_calc, py_calc:
        for seed in (None, 1, 2):
            atoms = water()
            if seed is not None:
                atoms.rattle(0.05, seed=seed)
            through_c, through_py = atoms.copy(), atoms.copy()
            through_c.calc, through_py.calc = c_calc, py_calc
            assert_close(through_c.get_potential_energy(), through_py.get_potential_energy())
            assert_close(through_c.get_forces(), through_py.get_forces())
            with pytest.raises(Exception):
                through_c.get_stress()


@pytest.mark.parametrize("kind", ["unix", "inet"])
def test_ase_periodic_matches_python_driver(kind, periodic_params, drivers, codata2018):
    c_calc, py_calc = socket_pair(kind, drivers, periodic_params)
    with c_calc, py_calc:
        for cell, seed in ((TRICLINIC, None), (TRICLINIC, 4), (CUBE, 4), (np.array(CUBE) * 1.07, 5)):
            atoms = water(cell=cell, pbc=True)
            if seed is not None:
                atoms.rattle(0.05, seed=seed)
            through_c, through_py = atoms.copy(), atoms.copy()
            through_c.calc, through_py.calc = c_calc, py_calc
            assert_close(through_c.get_potential_energy(), through_py.get_potential_energy())
            assert_close(through_c.get_forces(), through_py.get_forces())
            assert_close(through_c.get_stress(), through_py.get_stress())
            assert np.abs(through_c.get_stress()).max() > 0.0


def test_exits_cleanly_when_ase_closes(params, drivers, codata2018):
    name = "cpmdc-" + uuid.uuid4().hex[:12]
    with SocketIOCalculator(unixsocket=name, timeout=TIMEOUT) as calc:
        client = drivers("c", "--unix", name, "--params", str(params), "--species", "OHH",
                         "--verbose")
        atoms = water()
        atoms.calc = calc
        atoms.get_forces()
        atoms.positions[0, 2] += 0.01
        atoms.get_forces()
    code, err = client.finish()
    assert code == 0, err
    assert "the server closed the connection" in err
    assert "2 steps" in err


# the minimal i-PI-style server, against the Python driver in process


def wire(cell_rows_ang, positions_ang) -> tuple[np.ndarray, np.ndarray]:
    h = (np.asarray(cell_rows_ang, dtype=float).T / BOHR_TO_ANGSTROM).ravel()
    return h, (np.asarray(positions_ang, dtype=float) / BOHR_TO_ANGSTROM).ravel()


def python_step(params: Path, cell_rows_ang, positions_ang, cluster: bool = False):
    h, pos = wire(cell_rows_ang, positions_ang)
    with Driver(params, [8, 1, 1], cluster=cluster) as driver:
        energy, forces, virial = driver.compute(h, pos)
    # The driver sends the virial laid out as the cell, which MiniServer
    # transposes back.
    return energy, np.reshape(forces, (3, 3)), np.reshape(virial, (3, 3))


def test_full_handshake_and_exit(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    assert mini.status() == "NEEDINIT"
    mini.init(b"cpmdc test", bead=3)
    energy, forces, virial, extra = mini.step(np.zeros((3, 3)), WATER)
    expected = python_step(params, np.zeros((3, 3)), WATER)
    assert energy == expected[0]
    np.testing.assert_array_equal(forces, expected[1])
    assert not virial.any()
    assert extra == b""
    assert mini.status() == "READY"
    mini.send(mini.header("EXIT"))
    code, err = client.finish()
    assert code == 0, err


def test_cell_change_matches_python_driver(mini, socket_dir, periodic_params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(periodic_params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    for cell in (TRICLINIC, CUBE, np.array(CUBE) * 0.95):
        energy, forces, virial, _ = mini.step(cell, WATER)
        expected = python_step(periodic_params, cell, WATER)
        assert energy == expected[0]
        np.testing.assert_array_equal(forces, expected[1])
        np.testing.assert_allclose(virial, expected[2], rtol=1e-10, atol=1e-16)
        assert np.abs(virial).max() > 0.0
    mini.send(mini.header("EXIT"))
    assert client.finish()[0] == 0


def test_pipelined_posdata_and_getforce(mini, socket_dir, params, drivers):
    drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    assert mini.status() == "NEEDINIT"
    mini.init()
    assert mini.status() == "READY"
    for shift in (0.0, 0.02, 0.04):
        moved = np.array(WATER) + [0.0, 0.0, shift]
        energy, forces, _, _ = mini.getforce(pipelined=mini.posdata(np.zeros((3, 3)), moved))
        expected = python_step(params, np.zeros((3, 3)), moved)
        assert energy == pytest.approx(expected[0], rel=1e-12)
        np.testing.assert_allclose(forces, expected[1], rtol=1e-11, atol=1e-14)


def test_cluster_flag_ignores_the_cell(mini, socket_dir, params, drivers):
    drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH",
            "--cluster")
    mini.accept()
    mini.status()
    mini.init()
    energy, _, virial, _ = mini.step(CUBE, WATER)
    assert energy == python_step(params, CUBE, WATER, cluster=True)[0]
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)
    assert not virial.any()


def test_periodic_cell_without_stress_warns_once(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    for scale in (1.0, 1.02):
        energy, _, virial, _ = mini.step(np.array(CUBE) * scale, WATER)
        assert energy == python_step(params, np.array(CUBE) * scale, WATER)[0]
        assert not virial.any()
    mini.send(mini.header("EXIT"))
    code, err = client.finish()
    assert code == 0, err
    assert err.count("the library computed no stress for this periodic cell") == 1


def test_socket_prefix(socket_dir, params, drivers):
    server = MiniServer(socket_dir / "pre_s")
    try:
        drivers("c", "--unix", "s", "--socket-prefix", str(socket_dir / "pre_"),
                "--params", str(params), "--species", "OHH")
        server.accept()
        assert server.status() == "NEEDINIT"
    finally:
        server.close()


def test_server_dropped_between_messages(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.step(np.zeros((3, 3)), WATER)
    mini.conn.close()
    code, err = client.finish()
    assert code == 0, err


def test_server_dropped_mid_message(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.send(mini.posdata(np.zeros((3, 3)), WATER)[:40])
    mini.conn.close()
    code, err = client.finish()
    assert code == 1
    assert "closed the connection in the middle of POSDATA" in err


def test_changed_atom_count_is_an_error(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.step(np.zeros((3, 3)), WATER)
    mini.send(mini.posdata(np.zeros((3, 3)), WATER[:2]))
    code, err = client.finish()
    assert code == 1
    assert "sent 2 atoms" in err and "started for 3" in err
    assert "carries no species" in err


def test_batched_init_is_refused(mini, socket_dir, params, drivers):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init(b"batch_size:4")
    code, err = client.finish()
    assert code == 1
    assert "batch_size 4" in err


@pytest.mark.parametrize(
    "message, expected",
    [
        (MiniServer.header("GETFORCE"), "GETFORCE arrived before POSDATA"),
        (MiniServer.header("BOGUS"), "unknown message 'BOGUS'"),
        (
            MiniServer.posdata([[5.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 0.0, 5.0]], WATER),
            "without three independent vectors",
        ),
        (
            MiniServer.posdata(np.zeros((3, 3)), WATER)
            + MiniServer.posdata(np.zeros((3, 3)), WATER),
            "POSDATA arrived before GETFORCE collected the last step",
        ),
    ],
    ids=["getforce-first", "unknown", "flat-cell", "posdata-twice"],
)
def test_bad_messages_stop_the_driver(mini, socket_dir, params, drivers, message, expected):
    client = drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params),
                     "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.send(message)
    code, err = client.finish()
    assert code == 1
    assert expected in err


def test_posdata_without_init(mini, socket_dir, params, drivers):
    drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.send(mini.posdata(np.zeros((3, 3)), WATER))
    assert mini.status() == "HAVEDATA"
    energy, *_ = mini.getforce()
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)


@pytest.mark.parametrize("species", ["O H H", "O,H,H", "OH2", "8 1 1"])
def test_species_spellings(mini, socket_dir, params, drivers, species):
    drivers("c", "--unix", str(socket_dir / "s"), "--params", str(params), "--species", species)
    mini.accept()
    mini.status()
    mini.init()
    energy, *_ = mini.step(np.zeros((3, 3)), WATER)
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)


def test_cli_errors_before_connecting(params, drivers, socket_dir):
    for species, message in (("Qq", "unknown element"), ("o h h", "cannot read species"),
                             ("119", "outside 1..118"), ("", "at least one atom")):
        code, err = drivers("c", "--unix", "x", "--params", str(params),
                            "--species", species).finish()
        assert code == 1 and message in err, err
    code, err = drivers(
        "c", "--unix", str(socket_dir / "missing"), "--params", str(params), "--species", "OHH"
    ).finish()
    assert code == 1 and "cannot connect" in err
    code, err = drivers("c", "--port", "70000", "--params", str(params),
                        "--species", "OHH").finish()
    assert code == 1 and "outside 1..65535" in err
    code, err = drivers(
        "c", "--unix", "x", "--params", str(socket_dir / "none.bin"), "--species", "OHH"
    ).finish()
    assert code == 1 and "cannot read" in err
    code, err = drivers("c", "--params", str(params), "--species", "OHH").finish()
    assert code == 2 and "exactly one of --unix and --port" in err


# i-PI itself, both drivers


def ipi_run(kind: str, params: Path, drivers, socket_dir: Path, work: Path):
    ipi = Path(sys.executable).parent / "i-pi"
    if not ipi.is_file():
        found = shutil.which("i-pi")
        if found is None:
            pytest.skip("the i-pi command is not installed")
        ipi = Path(found)
    work.mkdir()
    bohr = np.array(WATER) / BOHR_TO_ANGSTROM
    lines = ["3", "# CELL(abcABC): 20.0 20.0 20.0 90.0 90.0 90.0 positions{atomic_unit} "
             "cell{atomic_unit}"]
    lines += [f"{sym} {x:.17g} {y:.17g} {z:.17g}" for sym, (x, y, z) in zip("OHH", bohr)]
    (work / "init.xyz").write_text("\n".join(lines) + "\n")
    prefix = str(socket_dir) + os.sep + kind + "_"
    (work / "input.xml").write_text(IPI_INPUT.format(prefix=prefix))
    server = subprocess.Popen(
        [sys.executable, str(ipi), "input.xml"],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        client = drivers(kind, "--unix", "water", "--socket-prefix", prefix, "--params",
                         str(params), "--species", "OHH", "--cluster", "--wait", str(TIMEOUT))
        code, err = client.finish()
        out, _ = server.communicate(timeout=TIMEOUT)
    finally:
        if server.poll() is None:
            server.kill()
            server.communicate()
    assert code == 0, err
    assert server.returncode == 0, out
    return (work / "sim.out").read_text(), (work / "sim.for_0.xyz").read_text()


def test_ipi_nve_run_matches_python_driver(params, drivers, socket_dir, tmp_path):
    pytest.importorskip("ipi")
    c_out, c_forces = ipi_run("c", params, drivers, socket_dir, tmp_path / "c")
    py_out, py_forces = ipi_run("python", params, drivers, socket_dir, tmp_path / "python")
    rows = np.loadtxt(tmp_path / "c" / "sim.out", ndmin=2)
    assert rows.shape[0] == 4
    assert np.all(np.diff(rows[:, 1]) != 0.0)
    # i-PI prints the potential with nine significant digits and forces with six.
    assert rows[0, 1] == pytest.approx(native(params, WATER).energy_h, rel=1e-8)
    assert c_out == py_out
    assert c_forces == py_forces
