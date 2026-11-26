"""A recomputed molecular dipole is the i-PI extra string.

i-PI's communication protocol appends a length-prefixed extra string to
``FORCEREADY`` (https://docs.ipi-code.org/distributed.html).
``ipi.interfaces.utils.parse_extra`` JSON-decodes that string, and the
dummy driver sends ``{"dipole": [...]}``
(https://github.com/i-pi/i-pi/blob/master/ipi/pes/dummy.py).
ASE's ``SocketClient`` reads the announced number of extra bytes and
ignores their contents
(https://gitlab.com/ase/ase/-/blob/master/ase/calculators/socketio.py).

cpmdc already returns the dipole in electron Bohr from
``cpmdc_session_last_dipole`` / ``Evaluation.dipole_e_bohr``. These tests
require both socket clients to forward a finite dipole and to send no
extra bytes when the library did not recompute one.
"""

from __future__ import annotations

import json
import math
import os
import socket
import threading
from pathlib import Path

import numpy as np
import pytest

from test_ipi_driver import MiniServer, WATER, mini, socket_dir  # noqa: F401

from cpmdc import BOHR_TO_ANGSTROM, CPMDCError, Evaluation
from cpmdc._rank import _cell_volume
from cpmdc.ipi import Driver

DIPOLE = (1.0, -2.5, 0.0)
DIPOLE_EXTRA = (
    b'{"dipole":[1.0000000000000000e+00,-2.5000000000000000e+00,'
    b"0.0000000000000000e+00]}"
)
FORCES = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
TRICLINIC = [[9.0, 0.0, 0.0], [1.2, 8.5, 0.0], [-0.7, 0.9, 9.4]]


def _evaluation(dipole: list[float] | None, stress: list[float] | None = None) -> Evaluation:
    return Evaluation(-1.5, FORCES, stress, dipole)


def _serve(driver: Driver):
    server_end, client_end = socket.socketpair()
    worker = threading.Thread(target=driver.serve, args=(client_end,))
    worker.start()
    server = MiniServer(conn=server_end)
    return server, worker, client_end


def _finish(driver: Driver, server: MiniServer, worker: threading.Thread, client_end: socket.socket) -> None:
    server.send(server.header("EXIT"))
    worker.join(60.0)
    assert not worker.is_alive()
    server.close()
    client_end.close()
    driver.close()


def test_recomputed_dipole_is_the_extra_string(params: Path) -> None:
    """Finite dipole JSON follows the virial; a later step without one is empty."""
    parse_extra = None
    try:
        from ipi.interfaces.utils import parse_extra as _parse_extra
    except ImportError:
        _parse_extra = None
    else:
        parse_extra = _parse_extra
    dipoles: list[list[float] | None] = [list(DIPOLE), None]
    with Driver(params, [8, 1, 1]) as driver:
        driver._session.energy_forces = lambda *_a, **_k: _evaluation(dipoles.pop(0))
        server, worker, client_end = _serve(driver)
        assert server.status() == "NEEDINIT"
        server.init(b"cpmdc dipole")
        energy, forces, virial, extra = server.step(np.zeros((3, 3)), WATER)
        assert energy == pytest.approx(-1.5)
        np.testing.assert_allclose(forces, FORCES)
        assert not virial.any()
        assert extra == DIPOLE_EXTRA
        parsed = json.loads(extra.decode("utf-8"))
        assert parsed == {"dipole": [1.0, -2.5, 0.0]}
        if parse_extra is not None:
            decoded = parse_extra(extra.decode("utf-8"))
            assert decoded["dipole"] == [1.0, -2.5, 0.0]
        _energy, _forces, _virial, again = server.step(np.zeros((3, 3)), WATER)
        assert again == b""
        _finish(driver, server, worker, client_end)
    assert driver.steps == 2


def test_dipole_extra_does_not_shift_the_virial(params: Path) -> None:
    stress = [0.5] * 9
    with Driver(params, [8, 1, 1]) as driver:
        driver._session.energy_forces = lambda *_a, **_k: _evaluation(list(DIPOLE), stress)
        server, worker, client_end = _serve(driver)
        server.status()
        server.init()
        _energy, _forces, virial, extra = server.step(TRICLINIC, WATER)
        rows = [[TRICLINIC[k][i] / BOHR_TO_ANGSTROM for i in range(3)] for k in range(3)]
        volume = _cell_volume(rows)
        np.testing.assert_allclose(virial, volume * np.array(stress).reshape(3, 3), rtol=1e-12)
        assert extra == DIPOLE_EXTRA
        _finish(driver, server, worker, client_end)


def test_negative_zero_dipole_keeps_its_sign(params: Path) -> None:
    with Driver(params, [8, 1, 1]) as driver:
        driver._session.energy_forces = lambda *_a, **_k: _evaluation([-0.0, 0.0, 1.0])
        server, worker, client_end = _serve(driver)
        server.status()
        server.init()
        _energy, _forces, _virial, extra = server.step(np.zeros((3, 3)), WATER)
        assert extra == (
            b'{"dipole":[-0.0000000000000000e+00,0.0000000000000000e+00,'
            b"1.0000000000000000e+00]}"
        )
        loaded = json.loads(extra.decode("utf-8"))["dipole"][0]
        assert loaded == 0.0
        assert math.copysign(1.0, loaded) < 0.0
        _finish(driver, server, worker, client_end)


def test_nonfinite_dipole_is_not_counted(params: Path) -> None:
    positions = [0.0] * 9
    cell = [0.0] * 9
    with Driver(params, [8, 1, 1]) as driver:
        driver._session.energy_forces = lambda *_a, **_k: _evaluation([float("nan"), 0.0, 0.0])
        with pytest.raises(CPMDCError, match="dipole is not finite"):
            driver.compute(positions, cell)
        assert driver.steps == 0
        driver._session.energy_forces = lambda *_a, **_k: _evaluation([0.0, float("inf"), 0.0])
        with pytest.raises(CPMDCError, match="dipole is not finite"):
            driver.compute(positions, cell)
        assert driver.steps == 0
        driver._session.energy_forces = lambda *_a, **_k: _evaluation([0.0, 0.0])
        with pytest.raises(CPMDCError, match="dipole must have three components"):
            driver.compute(positions, cell)
        assert driver.steps == 0


def _c_client(params: Path, socket_path: Path, spec: str):
    import subprocess

    client = os.environ.get("CPMDC_IPI_C")
    interpose = os.environ.get("CPMDC_IPI_DIPOLE_INTERPOSE")
    if not client or not interpose:
        pytest.skip("cpmdc-ipi-c dipole interposer is not built")
    env = os.environ.copy()
    # ASan refuses to start unless its runtime is the first preloaded library.
    # The interposer still supplies cpmdc_session_last_dipole, which libasan
    # does not define.
    existing = env.get("LD_PRELOAD", "")
    env["LD_PRELOAD"] = f"{existing} {interpose}".strip() if existing else interpose
    env["CPMDC_TEST_DIPOLE"] = spec
    return subprocess.Popen(
        [client, "--unix", str(socket_path), "--params", str(params), "--species", "OHH"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def test_c_driver_sends_the_same_dipole_extra(params: Path, mini: MiniServer, socket_dir: Path) -> None:
    client = _c_client(params, socket_dir / "s", "once")
    try:
        mini.accept()
        assert mini.status() == "NEEDINIT"
        mini.init()
        _energy, forces, virial, extra = mini.step(np.zeros((3, 3)), WATER)
        assert extra == DIPOLE_EXTRA
        _energy, _forces, _virial, again = mini.step(np.zeros((3, 3)), WATER)
        assert again == b""
        assert forces.shape == (3, 3)
        assert virial.shape == (3, 3)
        mini.send(mini.header("EXIT"))
        _out, err = client.communicate(timeout=60.0)
        assert client.returncode == 0, err
    finally:
        if client.poll() is None:
            client.kill()
            client.communicate()


def test_c_driver_rejects_a_nonfinite_dipole(params: Path, mini: MiniServer, socket_dir: Path) -> None:
    client = _c_client(params, socket_dir / "s", "nan")
    try:
        mini.accept()
        mini.status()
        mini.init()
        mini.send(mini.posdata(np.zeros((3, 3)), WATER))
        _out, err = client.communicate(timeout=60.0)
        assert client.returncode == 1, err
        assert "dipole is not finite" in err
    finally:
        if client.poll() is None:
            client.kill()
            client.communicate()
