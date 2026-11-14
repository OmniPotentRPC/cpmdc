"""cpmdc-ipi against ASE's socket calculator, i-PI, and a minimal i-PI server.

The minimal server follows the i-PI documentation, "Communication
protocol" (https://docs.ipi-code.org/distributed.html), and i-PI's
``ipi/interfaces/sockets.py``; ASE's side is ``ase/calculators/socketio.py``.

ASE converts with CODATA 2014 constants and cpmdc with CODATA 2018 ones.
The ``codata2018`` fixture gives ASE cpmdc's Bohr and Hartree, so a step
through the socket and a direct step do the same arithmetic.
"""

from __future__ import annotations

import os
import shutil
import socket
import struct
import subprocess
import sys
import tempfile
import threading
import uuid
from pathlib import Path

import ase.units
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.socketio import SocketIOCalculator, actualunixsocketname
from ase.optimize import BFGS
from conftest import INSTALLED, ROOT, WATER

import cpmdc.ipi
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Session
from cpmdc.ase import CPMDCCalculator
from cpmdc.ipi import Driver, ProtocolError, species_numbers, unix_socket_path

TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]
CUBE = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
TIMEOUT = 60.0


def water(**kwargs) -> Atoms:
    return Atoms("OHH", positions=WATER, **kwargs)


def driver_env() -> dict[str, str]:
    env = os.environ.copy()
    if not INSTALLED:
        env["PYTHONPATH"] = os.pathsep.join(
            [str(ROOT / "python")] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
        )
    return env


class Client:
    """``python -m cpmdc.ipi`` in a subprocess."""

    def __init__(self, *args: str) -> None:
        self.proc = subprocess.Popen(
            [sys.executable, "-m", "cpmdc.ipi", *args],
            env=driver_env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def finish(self) -> tuple[int, str]:
        try:
            _, err = self.proc.communicate(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            _, err = self.proc.communicate()
            raise
        return self.proc.returncode, err

    def kill(self) -> None:
        if self.proc.poll() is None:
            self.proc.kill()
            self.proc.communicate()


@pytest.fixture
def clients():
    started: list[Client] = []

    def start(*args: str) -> Client:
        client = Client(*args)
        started.append(client)
        return client

    yield start
    for client in started:
        client.kill()


@pytest.fixture
def codata2018(monkeypatch: pytest.MonkeyPatch) -> None:
    # ASE derives its 2018 Bohr from other constants; it agrees to 1e-11.
    units = ase.units.create_units("2018")
    assert units["Bohr"] == pytest.approx(BOHR_TO_ANGSTROM, rel=1e-11)
    assert units["Hartree"] == pytest.approx(HARTREE_TO_EV, rel=1e-11)
    monkeypatch.setattr(ase.units, "Bohr", BOHR_TO_ANGSTROM)
    monkeypatch.setattr(ase.units, "Ha", HARTREE_TO_EV)


@pytest.fixture
def socket_dir():
    # AF_UNIX paths are limited to about 100 bytes, so not tmp_path.
    path = Path(tempfile.mkdtemp(prefix="ipi"))
    yield path
    shutil.rmtree(path, ignore_errors=True)


def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("localhost", 0))
        return probe.getsockname()[1]


class MiniServer:
    """The server half of the i-PI protocol, one client at a time."""

    def __init__(self, path: Path | None = None, conn: socket.socket | None = None) -> None:
        self.listener: socket.socket | None = None
        if path is not None:
            self.listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            self.listener.bind(str(path))
            self.listener.listen(1)
            self.listener.settimeout(TIMEOUT)
        self.conn = conn
        if conn is not None:
            conn.settimeout(TIMEOUT)

    def accept(self) -> None:
        assert self.listener is not None
        self.conn, _ = self.listener.accept()
        self.conn.settimeout(TIMEOUT)

    def send(self, data: bytes) -> None:
        assert self.conn is not None
        self.conn.sendall(data)

    def recv(self, nbytes: int) -> bytes:
        assert self.conn is not None
        out = b""
        while len(out) < nbytes:
            chunk = self.conn.recv(nbytes - len(out))
            assert chunk, "client closed the connection"
            out += chunk
        return out

    @staticmethod
    def header(name: str) -> bytes:
        return name.encode("ascii").ljust(12)

    def status(self) -> str:
        self.send(self.header("STATUS"))
        return self.recv(12).decode("ascii").strip()

    def init(self, text: bytes = b"", bead: int = 0) -> None:
        self.send(self.header("INIT") + struct.pack("=ii", bead, len(text)) + text)

    @staticmethod
    def posdata(cell_rows_ang, positions_ang) -> bytes:
        h = np.asarray(cell_rows_ang, dtype=float).T / BOHR_TO_ANGSTROM
        ih = np.linalg.pinv(h) if h.any() else np.zeros((3, 3))
        pos = np.asarray(positions_ang, dtype=float) / BOHR_TO_ANGSTROM
        return (
            MiniServer.header("POSDATA")
            + h.tobytes()
            + ih.tobytes()
            + struct.pack("=i", len(pos))
            + pos.tobytes()
        )

    def getforce(self, pipelined: bytes = b"") -> tuple[float, np.ndarray, np.ndarray, bytes]:
        self.send(pipelined + self.header("GETFORCE"))
        assert self.recv(12).decode("ascii").strip() == "FORCEREADY"
        (energy,) = struct.unpack("=d", self.recv(8))
        (n_atoms,) = struct.unpack("=i", self.recv(4))
        forces = np.frombuffer(self.recv(24 * n_atoms), dtype=float).reshape(n_atoms, 3)
        virial = np.frombuffer(self.recv(72), dtype=float).reshape(3, 3).T
        (extra,) = struct.unpack("=i", self.recv(4))
        return energy, forces, virial, self.recv(extra)

    def step(self, cell_rows_ang, positions_ang):
        assert self.status() == "READY"
        self.send(self.posdata(cell_rows_ang, positions_ang))
        assert self.status() == "HAVEDATA"
        return self.getforce()

    def close(self) -> None:
        if self.conn is not None:
            self.conn.close()
        if self.listener is not None:
            self.listener.close()


@pytest.fixture
def mini(socket_dir: Path):
    server = MiniServer(socket_dir / "s")
    yield server
    server.close()


def native(params: Path, positions_ang, cell=None):
    with Session(params) as session:
        return session.energy_forces(positions_ang, [8, 1, 1], cell)


# ASE SocketIOCalculator against the direct calculator


def socket_calculator(kind: str, clients, params: Path, *extra: str) -> SocketIOCalculator:
    if kind == "unix":
        name = "cpmdc-" + uuid.uuid4().hex[:12]
        calc = SocketIOCalculator(unixsocket=name, timeout=TIMEOUT)
        clients("--unix", name, "--params", str(params), "--species", "O H H", *extra)
    else:
        port = free_port()
        calc = SocketIOCalculator(port=port, timeout=TIMEOUT)
        clients("--host", "localhost", "--port", str(port), "--params", str(params),
                "--species", "OH2", *extra)
    return calc


def assert_same(left: Atoms, right: Atoms, stress: bool) -> None:
    assert left.get_potential_energy() == pytest.approx(
        right.get_potential_energy(), rel=1e-10, abs=1e-10
    )
    np.testing.assert_allclose(left.get_forces(), right.get_forces(), rtol=1e-10, atol=1e-10)
    if stress:
        np.testing.assert_allclose(left.get_stress(), right.get_stress(), rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("kind", ["unix", "inet"])
def test_cluster_round_trip_matches_direct(kind, params, clients, codata2018):
    calc = socket_calculator(kind, clients, params)
    direct = CPMDCCalculator(params)
    with calc, direct:
        for seed in (None, 1, 2):
            atoms = water()
            if seed is not None:
                atoms.rattle(0.05, seed=seed)
            through = atoms.copy()
            through.calc = calc
            atoms.calc = direct
            assert_same(through, atoms, stress=False)
            with pytest.raises(Exception):
                through.get_stress()


@pytest.mark.parametrize("kind", ["unix", "inet"])
def test_periodic_round_trip_matches_direct(kind, periodic_params, clients, codata2018):
    calc = socket_calculator(kind, clients, periodic_params)
    direct = CPMDCCalculator(periodic_params)
    with calc, direct:
        for cell, seed in ((TRICLINIC, None), (TRICLINIC, 4), (CUBE, 4), (np.array(CUBE) * 1.07, 5)):
            atoms = water(cell=cell, pbc=True)
            if seed is not None:
                atoms.rattle(0.05, seed=seed)
            through = atoms.copy()
            through.calc = calc
            atoms.calc = direct
            assert_same(through, atoms, stress=True)
            assert np.abs(through.get_stress()).max() > 0.0


def test_driver_exits_cleanly_when_ase_closes(params, clients, codata2018):
    name = "cpmdc-" + uuid.uuid4().hex[:12]
    with SocketIOCalculator(unixsocket=name, timeout=TIMEOUT) as calc:
        client = clients("--unix", name, "--params", str(params), "--species", "OHH", "--verbose")
        atoms = water()
        atoms.calc = calc
        atoms.get_forces()
        atoms.positions[0, 2] += 0.01
        atoms.get_forces()
    code, err = client.finish()
    assert code == 0, err
    assert "the server closed the connection" in err
    assert "2 steps" in err


def test_bfgs_through_socket_matches_direct(params, clients, codata2018):
    through = water()
    through.rattle(0.02, seed=7)
    direct = through.copy()
    with socket_calculator("unix", clients, params) as calc:
        through.calc = calc
        with BFGS(through, logfile=None) as opt:
            seen_socket = [through.get_potential_energy()]
            opt.attach(lambda: seen_socket.append(through.get_potential_energy()))
            opt.run(fmax=0.05, steps=6)
    with CPMDCCalculator(params) as calc:
        direct.calc = calc
        with BFGS(direct, logfile=None) as opt:
            seen_direct = [direct.get_potential_energy()]
            opt.attach(lambda: seen_direct.append(direct.get_potential_energy()))
            opt.run(fmax=0.05, steps=6)
    assert len(seen_socket) == len(seen_direct) > 3
    np.testing.assert_allclose(seen_socket, seen_direct, rtol=1e-9)
    np.testing.assert_allclose(through.positions, direct.positions, atol=1e-9)
    assert seen_socket[-1] < seen_socket[0]


def test_stock_ase_units_differ_by_codata_only(periodic_params, clients):
    atoms = water(cell=TRICLINIC, pbc=True)
    with socket_calculator("unix", clients, periodic_params) as calc:
        through = atoms.copy()
        through.calc = calc
        energy = through.get_potential_energy()
    with CPMDCCalculator(periodic_params) as calc:
        atoms.calc = calc
        reference = atoms.get_potential_energy()
    assert energy != reference
    assert energy == pytest.approx(reference, rel=1e-7)


# the minimal i-PI-style server


def test_full_handshake_and_exit(mini, socket_dir, params, clients):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    assert mini.status() == "NEEDINIT"
    mini.init(b"cpmdc test", bead=3)
    energy, forces, virial, extra = mini.step(np.zeros((3, 3)), WATER)
    expected = native(params, WATER)
    assert energy == pytest.approx(expected.energy_h, rel=1e-13)
    np.testing.assert_allclose(forces, expected.forces_h_bohr, rtol=1e-12, atol=1e-15)
    assert not virial.any()
    assert extra == b""
    assert mini.status() == "READY"
    mini.send(mini.header("EXIT"))
    code, err = client.finish()
    assert code == 0, err


def test_pipelined_posdata_and_getforce(mini, socket_dir, params, clients):
    """i-PI 3 sends POSDATA and GETFORCE in one write, with no STATUS between."""
    clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    assert mini.status() == "NEEDINIT"
    mini.init()
    assert mini.status() == "READY"
    for shift in (0.0, 0.02, 0.04):
        moved = np.array(WATER) + [0.0, 0.0, shift]
        energy, forces, _, _ = mini.getforce(pipelined=mini.posdata(np.zeros((3, 3)), moved))
        expected = native(params, moved)
        assert energy == pytest.approx(expected.energy_h, rel=1e-12)
        np.testing.assert_allclose(forces, expected.forces_h_bohr, rtol=1e-11, atol=1e-14)


def test_cell_change_between_steps(mini, socket_dir, periodic_params, clients):
    client = clients(
        "--unix", str(socket_dir / "s"), "--params", str(periodic_params), "--species", "OHH"
    )
    mini.accept()
    mini.status()
    mini.init()
    with CPMDCCalculator(periodic_params) as direct:
        for cell in (TRICLINIC, CUBE, np.array(CUBE) * 0.95):
            energy, forces, virial, _ = mini.step(cell, WATER)
            atoms = water(cell=cell, pbc=True)
            atoms.calc = direct
            assert energy * HARTREE_TO_EV == pytest.approx(
                atoms.get_potential_energy(), rel=1e-12
            )
            np.testing.assert_allclose(
                forces * HARTREE_TO_EV / BOHR_TO_ANGSTROM, atoms.get_forces(), rtol=1e-11,
                atol=1e-13,
            )
            # The virial is -volume * stress, as ASE's SocketClient sends it.
            expected = -atoms.get_volume() * atoms.get_stress(voigt=False) / HARTREE_TO_EV
            np.testing.assert_allclose(virial, expected, rtol=1e-10, atol=1e-16)
            np.testing.assert_allclose(virial, virial.T)
    mini.send(mini.header("EXIT"))
    assert client.finish()[0] == 0


def test_cluster_flag_ignores_the_cell(mini, socket_dir, params, clients):
    clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH",
            "--cluster")
    mini.accept()
    mini.status()
    mini.init()
    energy, _, virial, _ = mini.step(CUBE, WATER)
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)
    assert not virial.any()


def test_server_dropped_between_messages(mini, socket_dir, params, clients):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.step(np.zeros((3, 3)), WATER)
    mini.conn.close()
    code, err = client.finish()
    assert code == 0, err


def test_server_dropped_mid_message(mini, socket_dir, params, clients):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.send(mini.posdata(np.zeros((3, 3)), WATER)[:40])
    mini.conn.close()
    code, err = client.finish()
    assert code == 1
    assert "closed the connection in the middle of POSDATA" in err


def test_changed_atom_count_is_an_error(mini, socket_dir, params, clients):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.step(np.zeros((3, 3)), WATER)
    mini.send(mini.posdata(np.zeros((3, 3)), WATER[:2]))
    code, err = client.finish()
    assert code == 1
    assert "sent 2 atoms" in err and "started for 3" in err
    assert "carries no species" in err


def test_batched_init_is_refused(mini, socket_dir, params, clients):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
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
    ],
    ids=["getforce-first", "unknown", "flat-cell"],
)
def test_bad_messages_stop_the_driver(mini, socket_dir, params, clients, message, expected):
    client = clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.status()
    mini.init()
    mini.send(message)
    code, err = client.finish()
    assert code == 1
    assert expected in err


def test_posdata_without_init(mini, socket_dir, params, clients):
    """ASE's client starts READY and some servers skip INIT; that works."""
    clients("--unix", str(socket_dir / "s"), "--params", str(params), "--species", "OHH")
    mini.accept()
    mini.send(mini.posdata(np.zeros((3, 3)), WATER))
    assert mini.status() == "HAVEDATA"
    energy, *_ = mini.getforce()
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)


# in process


def test_one_session_serves_every_step(params, monkeypatch):
    opened: list[int] = []
    evaluated: list[int] = []

    class Counting(Session):
        def __init__(self, *args, **kwargs):
            opened.append(1)
            super().__init__(*args, **kwargs)

        def energy_forces(self, *args, **kwargs):
            evaluated.append(1)
            return super().energy_forces(*args, **kwargs)

    monkeypatch.setattr(cpmdc.ipi, "Session", Counting)
    server_end, client_end = socket.socketpair()
    with Driver(params, [8, 1, 1]) as driver:
        worker = threading.Thread(target=driver.serve, args=(client_end,))
        worker.start()
        server = MiniServer(conn=server_end)
        server.status()
        server.init()
        for shift in range(4):
            server.step(np.zeros((3, 3)), np.array(WATER) + 0.01 * shift)
        server.send(server.header("EXIT"))
        worker.join(TIMEOUT)
        assert not worker.is_alive()
    server.close()
    client_end.close()
    assert opened == [1]
    assert len(evaluated) == 4 == driver.steps


def test_compute_virial_sign(periodic_params):
    h = (np.array(TRICLINIC).T / BOHR_TO_ANGSTROM).ravel()
    pos = (np.array(WATER) / BOHR_TO_ANGSTROM).ravel()
    with Driver(periodic_params, [8, 1, 1]) as driver:
        energy, forces, virial = driver.compute(h, pos)
    with Session(periodic_params) as session:
        done = session.energy_forces(WATER, [8, 1, 1], TRICLINIC)
    volume = abs(np.linalg.det(np.array(TRICLINIC))) / BOHR_TO_ANGSTROM**3
    assert energy == pytest.approx(done.energy_h, rel=1e-13)
    np.testing.assert_allclose(virial, volume * np.array(done.stress_h_bohr3), rtol=1e-12)
    assert len(forces) == 9


def test_compute_rejects_another_atom_count(params):
    with Driver(params, [8, 1, 1]) as driver:
        with pytest.raises(ProtocolError, match="carries no species"):
            driver.compute([0.0] * 9, [0.0] * 6)


def test_species_parsing():
    assert species_numbers("O H H") == [8, 1, 1]
    assert species_numbers("O,H,H") == [8, 1, 1]
    assert species_numbers("OH2") == [8, 1, 1]
    assert species_numbers("Si3N4") == [14] * 3 + [7] * 4
    assert species_numbers("8 1 1") == [8, 1, 1]
    for bad in ("", "Xx", "o h h", "0", "119", "O-H"):
        with pytest.raises(ValueError):
            species_numbers(bad)


def test_unix_socket_path_matches_ase():
    assert unix_socket_path("water") == actualunixsocketname("water")
    assert unix_socket_path("water", "pre_") == "pre_water"
    assert unix_socket_path(os.path.join("run", "sock")) == os.path.join("run", "sock")


def test_cli_errors_before_connecting(params, clients, socket_dir):
    code, err = clients("--unix", "x", "--params", str(params), "--species", "Qq").finish()
    assert code == 1 and "unknown element" in err
    code, err = clients(
        "--unix", str(socket_dir / "missing"), "--params", str(params), "--species", "OHH"
    ).finish()
    assert code == 1 and "cannot connect" in err
    code, err = clients("--port", "70000", "--params", str(params), "--species", "OHH").finish()
    assert code == 1 and "outside 1..65535" in err
    code, err = clients(
        "--unix", "x", "--params", str(socket_dir / "none.bin"), "--species", "OHH"
    ).finish()
    assert code == 1


def test_structure_file_gives_species(mini, socket_dir, params, clients, tmp_path):
    path = tmp_path / "water.xyz"
    water().write(path)
    clients("--unix", str(socket_dir / "s"), "--params", str(params), "--structure", str(path))
    mini.accept()
    mini.status()
    mini.init()
    energy, *_ = mini.step(np.zeros((3, 3)), WATER)
    assert energy == pytest.approx(native(params, WATER).energy_h, rel=1e-13)


# i-PI itself

IPI_INPUT = """<simulation verbosity="low" safe_stride="1000" sockets_prefix="{prefix}">
  <output prefix="sim">
    <properties stride="1" filename="out"> [ step, potential{{atomic_unit}} ] </properties>
    <trajectory filename="for" stride="1" format="xyz"> forces </trajectory>
  </output>
  <total_steps> 3 </total_steps>
  <prng><seed>31415</seed></prng>
  <ffsocket name="cpmdc" mode="unix" pbc="false">
    <address> water </address>
  </ffsocket>
  <system>
    <initialize nbeads="1">
      <file mode="xyz"> init.xyz </file>
      <velocities mode="thermal" units="kelvin"> 0 </velocities>
    </initialize>
    <forces><force forcefield="cpmdc"/></forces>
    <ensemble><temperature units="kelvin"> 0 </temperature></ensemble>
    <motion mode="dynamics">
      <dynamics mode="nve"><timestep units="femtosecond"> 0.25 </timestep></dynamics>
    </motion>
  </system>
</simulation>
"""


def test_ipi_nve_run(params, clients, socket_dir, tmp_path):
    pytest.importorskip("ipi")
    ipi = Path(sys.executable).parent / "i-pi"
    if not ipi.is_file():
        found = shutil.which("i-pi")
        if found is None:
            pytest.skip("the i-pi command is not installed")
        ipi = Path(found)
    bohr = np.array(WATER) / BOHR_TO_ANGSTROM
    lines = ["3", "# CELL(abcABC): 20.0 20.0 20.0 90.0 90.0 90.0 positions{atomic_unit} "
             "cell{atomic_unit}"]
    lines += [f"{sym} {x:.17g} {y:.17g} {z:.17g}" for sym, (x, y, z) in zip("OHH", bohr)]
    (tmp_path / "init.xyz").write_text("\n".join(lines) + "\n")
    prefix = str(socket_dir) + os.sep + "ipi_"
    (tmp_path / "input.xml").write_text(IPI_INPUT.format(prefix=prefix))
    server = subprocess.Popen(
        [sys.executable, str(ipi), "input.xml"],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        client = clients("--unix", "water", "--socket-prefix", prefix, "--params", str(params),
                         "--species", "OHH", "--cluster", "--wait", str(TIMEOUT))
        code, err = client.finish()
        out, _ = server.communicate(timeout=TIMEOUT)
    finally:
        if server.poll() is None:
            server.kill()
            server.communicate()
    assert code == 0, err
    assert server.returncode == 0, out
    rows = np.loadtxt(tmp_path / "sim.out", ndmin=2)
    assert rows.shape[0] == 4
    expected = native(params, WATER)
    # i-PI prints the potential with nine significant digits and forces with six.
    assert rows[0, 1] == pytest.approx(expected.energy_h, rel=1e-8)
    frame = (tmp_path / "sim.for_0.xyz").read_text().splitlines()[2:5]
    forces = np.array([[float(v) for v in line.split()[1:4]] for line in frame])
    np.testing.assert_allclose(forces, expected.forces_h_bohr, rtol=1e-5, atol=1e-9)
    assert np.all(np.diff(rows[:, 1]) != 0.0)
