"""TCP round trip of the MDI engine against the stub library.

The driver and the engine are subprocesses, because ``MDI_Init`` is once
per process. ``pymdi`` is optional: the round trip skips when it is missing,
and the command-line tests do not import it.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import INSTALLED, ROOT, WATER

from cpmdc import BOHR_TO_ANGSTROM, Session
from cpmdc.mdi import COMMANDS

TIMEOUT = 60.0


def driver_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if not INSTALLED:
        env["PYTHONPATH"] = os.pathsep.join(
            [str(ROOT / "python")] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
        )
    return env


def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def bohr_of(rows: list[list[float]]) -> list[float]:
    return [value / BOHR_TO_ANGSTROM for row in rows for value in row]


def det3(m: list[float]) -> float:
    return (
        m[0] * (m[4] * m[8] - m[5] * m[7])
        - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
    )


def direct(params: Path, coords: list[float], numbers: list[int], cell: list[float] | None):
    with Session(params) as session:
        energy, forces, virial = session.calculate_au(coords, numbers, cell)
    forces = [float(value) for value in forces]
    virial = [float(value) for value in virial]
    if cell is None:
        stress = [0.0] * 9
    else:
        volume = abs(det3(cell))
        stress = [value / volume for value in virial] if volume else [0.0] * 9
    return float(energy), forces, stress


def require_pymdi() -> None:
    try:
        import mdi  # noqa: F401
    except ImportError:
        pytest.skip("pymdi is not installed")


class Pair:
    def __init__(
        self, port: int, script: list, params: Path, work: Path, *engine_args: str
    ) -> None:
        self.port = port
        self.params = params
        self.engine_args = engine_args
        script_path = work / f"mdi-script-{port}.json"
        self.result_path = work / f"mdi-result-{port}.json"
        script_path.write_text(json.dumps(script), encoding="utf-8")
        self.driver = subprocess.Popen(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--driver",
                str(port),
                str(script_path),
                str(self.result_path),
            ],
            env=driver_env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.engine: subprocess.Popen[str] | None = None

    def start_engine(self) -> None:
        self.engine = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "cpmdc.mdi",
                "-mdi",
                (
                    "-role ENGINE -name CPMDC -method TCP "
                    f"-hostname 127.0.0.1 -port {self.port}"
                ),
                "--params",
                str(self.params),
                "--verbose",
                *self.engine_args,
            ],
            env=driver_env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def finish(self) -> tuple[int, str, int, str, dict]:
        assert self.engine is not None
        try:
            _, engine_err = self.engine.communicate(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            self.kill()
            raise
        try:
            _, driver_err = self.driver.communicate(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            self.kill()
            raise
        result = {}
        if self.result_path.is_file():
            result = json.loads(self.result_path.read_text(encoding="utf-8"))
        return (
            self.engine.returncode if self.engine is not None else 1,
            engine_err,
            self.driver.returncode,
            driver_err,
            result,
        )

    def kill(self) -> None:
        for proc in (self.engine, self.driver):
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.communicate()


@pytest.fixture
def pairs():
    started: list[Pair] = []

    def start(port: int, script: list, params: Path, work: Path, *engine_args: str) -> Pair:
        pair = Pair(port, script, params, work, *engine_args)
        started.append(pair)
        assert pair.driver.stdout is not None
        line = pair.driver.stdout.readline()
        if line.strip() != "LISTENING":
            pair.kill()
            raise AssertionError(f"driver did not listen: {line!r} {pair.driver.stderr.read()}")
        pair.start_engine()
        return pair

    yield start
    for pair in started:
        pair.kill()


def _exchange(commands: list) -> list:
    script = []
    for command in commands:
        script.append(["check", command])
    script.append(["check", "<KE_ELEC"])
    script.append(["check", ">TYPES"])
    return script


def test_cli_help() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "cpmdc.mdi", "--help"],
        env=driver_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert proc.returncode == 0
    assert "cpmdc-mdi" in proc.stdout
    assert "-mdi" in proc.stdout


def test_cli_rejects_missing_options(params: Path) -> None:
    env = driver_env()
    env.pop("CPMDC_PARAMS", None)
    missing = subprocess.run(
        [sys.executable, "-m", "cpmdc.mdi", "--params", str(params)],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert missing.returncode == 1
    assert "-mdi" in missing.stderr
    bare = subprocess.run(
        [
            sys.executable,
            "-m",
            "cpmdc.mdi",
            "-mdi",
            "-role ENGINE -name CPMDC -method TCP -hostname 127.0.0.1 -port 1",
        ],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert bare.returncode == 1
    assert "CPMDC_PARAMS" in bare.stderr


def test_cli_rejects_mpi_and_link(params: Path) -> None:
    env = driver_env()
    mpi = subprocess.run(
        [
            sys.executable,
            "-m",
            "cpmdc.mdi",
            "-mdi",
            "-role ENGINE -name CPMDC -method MPI",
            "--params",
            str(params),
        ],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert mpi.returncode == 1
    assert "MPI" in mpi.stderr
    link = subprocess.run(
        [
            sys.executable,
            "-m",
            "cpmdc.mdi",
            "-mdi",
            "-role ENGINE -name CPMDC -method LINK",
            "--params",
            str(params),
        ],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert link.returncode == 1
    assert "MDI_Plugin_init_cpmdc" in link.stderr


def test_import_does_not_load_pymdi() -> None:
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, cpmdc; assert 'mdi' not in sys.modules; "
            "assert 'cpmdc.mdi' not in sys.modules; "
            "from cpmdc import MDI_Plugin_init_cpmdc; "
            "assert 'cpmdc.mdi' not in sys.modules",
        ],
        env=driver_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_tcp_commands_units_and_cell(
    params: Path, periodic_params: Path, pairs, tmp_path: Path
) -> None:
    require_pymdi()
    coords = bohr_of(WATER)
    nudged = [coords[0] + 0.05, *coords[1:]]
    # H, O, H: the same elements in another order, with the coordinates moved to match.
    reordered_numbers = [1, 8, 1]
    reordered = coords[3:6] + coords[0:3] + coords[6:9]
    cube = [8.0 / BOHR_TO_ANGSTROM, 0.0, 0.0, 0.0, 9.0 / BOHR_TO_ANGSTROM, 0.0, 0.0, 0.0, 10.0 / BOHR_TO_ANGSTROM]
    script = [
        *_exchange(COMMANDS),
        ["command", "<NAME"],
        ["recv_name", "name"],
        ["command", "<@"],
        ["recv_name", "node"],
        ["command", ">NATOMS"],
        ["send_int", 2],
        ["command", ">NATOMS"],
        ["send_int", 3],
        ["command", "<NATOMS"],
        ["recv_int", "natoms"],
        ["command", ">ELEMENTS"],
        ["send_ints", [8, 1, 1]],
        ["command", ">COORDS"],
        ["send_doubles", coords],
        ["command", "<COORDS"],
        ["recv_doubles", "coords", 9],
        ["command", ">CELL"],
        ["send_doubles", [0.0] * 9],
        ["command", "<CELL"],
        ["recv_doubles", "cell", 9],
        ["command", ">CELL_DISPL"],
        ["send_doubles", [0.4, -0.2, 0.7]],
        ["command", ">CHARGES"],
        ["send_doubles", [0.0, 0.0, 0.0]],
        ["command", ">CHARGES"],
        ["send_doubles", [8.0, 1.0, 1.0]],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy", 1],
        ["command", "<PE"],
        ["recv_doubles", "pe", 1],
        ["command", "<FORCES"],
        ["recv_doubles", "forces", 9],
        ["command", "<STRESS"],
        ["recv_doubles", "stress", 9],
        ["command", ">NATOMS"],
        ["send_int", 3],
        ["command", ">COORDS"],
        ["send_doubles", nudged],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy2", 1],
        ["command", "<FORCES"],
        ["recv_doubles", "forces2", 9],
        ["command", ">ELEMENTS"],
        ["send_ints", reordered_numbers],
        ["command", ">COORDS"],
        ["send_doubles", reordered],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy3", 1],
        ["command", "<FORCES"],
        ["recv_doubles", "forces3", 9],
        ["exit"],
    ]
    pair = pairs(free_port(), script, params, tmp_path)
    engine_code, engine_err, driver_code, driver_err, result = pair.finish()
    assert engine_code == 0, engine_err
    assert driver_code == 0, driver_err
    assert result["error"] is None, result
    out = result["out"]
    checks = {name: flag for name, flag in out["checks"]}
    for command in COMMANDS:
        assert checks[command] == 1
    assert checks["<KE_ELEC"] == 0
    assert checks[">TYPES"] == 0
    assert out["name"] == ["CPMDC"]
    assert out["node"] == ["@DEFAULT"]
    assert out["natoms"] == [3]
    assert out["coords"] == [coords]
    assert out["cell"] == [[0.0] * 9]
    energy, forces, stress = direct(params, coords, [8, 1, 1], None)
    assert out["energy"] == [[energy]]
    assert out["pe"] == [[energy]]
    assert out["forces"] == [forces]
    assert out["stress"] == [stress]
    assert stress == [0.0] * 9
    energy2, forces2, _stress2 = direct(params, nudged, [8, 1, 1], None)
    assert energy2 != energy
    assert out["energy2"] == [[energy2]]
    assert out["forces2"] == [forces2]
    energy3, forces3, _stress3 = direct(params, reordered, reordered_numbers, None)
    assert out["energy3"] == [[energy3]]
    assert out["forces3"] == [forces3]
    assert "2 steps" not in engine_err
    assert "3 steps" in engine_err
    # A periodic cell on the isolated deck is a separate connection.
    _periodic_cell(periodic_params, pairs, coords, cube, tmp_path)


def _periodic_cell(
    params: Path, pairs, coords: list[float], cell: list[float], work: Path
) -> None:
    taller = [cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8] + 1.0]
    script = [
        ["command", ">NATOMS"],
        ["send_int", 3],
        ["command", ">ELEMENTS"],
        ["send_ints", [8, 1, 1]],
        ["command", ">COORDS"],
        ["send_doubles", coords],
        ["command", ">CELL"],
        ["send_doubles", cell],
        ["command", "<CELL"],
        ["recv_doubles", "cell", 9],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy", 1],
        ["command", "<FORCES"],
        ["recv_doubles", "forces", 9],
        ["command", "<STRESS"],
        ["recv_doubles", "stress", 9],
        ["command", ">CELL"],
        ["send_doubles", taller],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy2", 1],
        ["command", "<STRESS"],
        ["recv_doubles", "stress2", 9],
        ["exit"],
    ]
    pair = pairs(free_port(), script, params, work)
    engine_code, engine_err, driver_code, driver_err, result = pair.finish()
    assert engine_code == 0, engine_err
    assert driver_code == 0, driver_err
    assert result["error"] is None, result
    out = result["out"]
    energy, forces, stress = direct(params, coords, [8, 1, 1], cell)
    assert any(stress)
    assert out["cell"] == [cell]
    assert out["energy"] == [[energy]]
    assert out["forces"] == [forces]
    assert out["stress"] == [stress]
    energy2, _forces2, stress2 = direct(params, coords, [8, 1, 1], taller)
    assert energy2 != energy
    assert out["energy2"] == [[energy2]]
    assert out["stress2"] == [stress2]
    assert "warning:" not in engine_err


def test_cluster_flag_ignores_the_cell(params: Path, pairs, tmp_path: Path) -> None:
    require_pymdi()
    coords = bohr_of(WATER)
    cell = [6.0 / BOHR_TO_ANGSTROM, 0.0, 0.0, 0.2, 7.0 / BOHR_TO_ANGSTROM, 0.0, 0.0, 0.1, 8.0 / BOHR_TO_ANGSTROM]
    script = [
        ["command", ">NATOMS"],
        ["send_int", 3],
        ["command", ">ELEMENTS"],
        ["send_ints", [8, 1, 1]],
        ["command", ">COORDS"],
        ["send_doubles", coords],
        ["command", ">CELL"],
        ["send_doubles", cell],
        ["command", "<CELL"],
        ["recv_doubles", "cell", 9],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy", 1],
        ["command", "<STRESS"],
        ["recv_doubles", "stress", 9],
        ["exit"],
    ]
    pair = pairs(free_port(), script, params, tmp_path, "--cluster")
    engine_code, engine_err, driver_code, driver_err, result = pair.finish()
    assert engine_code == 0, engine_err
    assert driver_code == 0, driver_err
    assert result["error"] is None
    energy, _forces, stress = direct(params, coords, [8, 1, 1], None)
    assert result["out"]["cell"] == [cell]
    assert result["out"]["energy"] == [[energy]]
    assert result["out"]["stress"] == [stress]
    assert stress == [0.0] * 9


def _fail(script: list, params: Path, pairs, work: Path, needle: str) -> None:
    require_pymdi()
    pair = pairs(free_port(), script, params, work)
    engine_code, engine_err, driver_code, _driver_err, result = pair.finish()
    assert engine_code == 1, engine_err
    assert needle in engine_err
    assert driver_code == 0
    assert result["error"]


def _system(coords: list[float]) -> list:
    return [
        ["command", ">NATOMS"],
        ["send_int", 3],
        ["command", ">ELEMENTS"],
        ["send_ints", [8, 1, 1]],
        ["command", ">COORDS"],
        ["send_doubles", coords],
        ["command", "<ENERGY"],
        ["recv_doubles", "energy", 1],
    ]


def test_error_paths(params: Path, pairs, tmp_path: Path) -> None:
    require_pymdi()
    coords = bohr_of(WATER)
    _fail(
        [["command", "<ENERGY"], ["recv_doubles", "energy", 1], ["exit"]],
        params,
        pairs,
        tmp_path,
        "cannot evaluate yet",
    )
    _fail(
        [["command", "<KE_ELEC"], ["recv_doubles", "energy", 1], ["exit"]],
        params,
        pairs,
        tmp_path,
        "unsupported MDI command",
    )
    # A > command has no acknowledgement. The following < command is what
    # fails in the driver once the engine has closed the connection.
    _fail(
        [
            *_system(coords),
            ["command", ">NATOMS"],
            ["send_int", 4],
            ["command", "<NATOMS"],
            ["recv_int", "natoms"],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "topology change requires a new session",
    )
    _fail(
        [
            *_system(coords),
            ["command", ">ELEMENTS"],
            ["send_ints", [8, 6, 1]],
            ["command", "<ENERGY"],
            ["recv_doubles", "energy2", 1],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "topology change requires a new session",
    )
    _fail(
        [
            ["command", ">NATOMS"],
            ["send_int", 3],
            ["command", ">CHARGES"],
            ["send_doubles", [0.2, 0.0, -0.2]],
            ["command", "<NATOMS"],
            ["recv_int", "natoms"],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "CPMDParams",
    )
    _fail(
        [
            ["command", ">CELL"],
            ["send_doubles", [1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0]],
            ["command", "<CELL"],
            ["recv_doubles", "cell", 9],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "independent",
    )
    _fail(
        [
            ["command", ">NATOMS"],
            ["send_int", 0],
            ["command", "<NATOMS"],
            ["recv_int", "natoms"],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "at least one atom",
    )
    _fail(
        [
            ["command", ">NATOMS"],
            ["send_int", 1],
            ["command", ">COORDS"],
            ["send_nan"],
            ["command", "<COORDS"],
            ["recv_doubles", "coords", 3],
            ["exit"],
        ],
        params,
        pairs,
        tmp_path,
        "finite",
    )


def test_plugin_skips_without_mpi() -> None:
    require_pymdi()
    try:
        import mpi4py  # noqa: F401
    except ImportError:
        pytest.skip(
            "MDI_Launch_plugin requires an MPI communicator; mpi4py is not installed"
        )


def test_lammps_mdi_qm_smoke(params: Path, tmp_path: Path) -> None:
    require_pymdi()
    import shutil

    lmp = shutil.which("lmp") or shutil.which("lmp_serial") or shutil.which("lmp_mpi")
    if lmp is None:
        pytest.skip("LAMMPS with the MDI package is not installed")
    probe = subprocess.run(
        [lmp, "-mdi", "probe", "-log", "none"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=30,
        check=False,
    )
    probe_text = probe.stdout + probe.stderr
    if "Invalid command-line argument: -mdi" in probe_text:
        pytest.skip("this LAMMPS was built without the MDI package")
    port = free_port()
    script = tmp_path / "in.mdi"
    script.write_text(
        """
units           metal
atom_style      atomic
atom_modify     map array
boundary        f f f
region          box block 0 10 0 10 0 10
create_box      2 box
mass            1 16.0
mass            2 1.0
create_atoms    1 single 0 0 0.1173
create_atoms    2 single 0 0.7572 -0.4692
create_atoms    2 single 0 -0.7572 -0.4692
fix             qm all mdi/qm elements O H
run             0
""",
        encoding="utf-8",
    )
    engine = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "cpmdc.mdi",
            "-mdi",
            f"-role ENGINE -name CPMDC -method TCP -hostname 127.0.0.1 -port {port}",
            "--params",
            str(params),
        ],
        env=driver_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        proc = subprocess.run(
            [lmp, "-mdi", f"-role DRIVER -name LMP -method TCP -port {port}", "-in", str(script)],
            cwd=tmp_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=TIMEOUT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        engine.kill()
        engine.communicate()
        raise
    text = proc.stdout + proc.stderr
    if proc.returncode != 0 and (
        "Unrecognized fix style" in text or "Invalid command-line argument: -mdi" in text
    ):
        engine.kill()
        engine.communicate()
        pytest.skip("this LAMMPS was built without the MDI package")
    if proc.returncode != 0:
        engine.kill()
        engine.communicate()
        raise AssertionError(text)
    engine_code, engine_err = _wait(engine)
    assert engine_code == 0, engine_err


def _wait(proc: subprocess.Popen[str]) -> tuple[int, str]:
    try:
        _, err = proc.communicate(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        proc.kill()
        _, err = proc.communicate()
        raise
    return proc.returncode, err


def _driver_main() -> int:
    import mdi

    port = int(sys.argv[2])
    script = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
    result_path = Path(sys.argv[4])
    mdi.MDI_Init(f"-role DRIVER -name driver -method TCP -hostname 127.0.0.1 -port {port}")
    print("LISTENING", flush=True)
    comm = mdi.MDI_Accept_Communicator()
    out: dict = {}
    error = None
    try:
        for step in script:
            op = step[0]
            if op == "check":
                flag = mdi.MDI_Check_Command_Exists("@DEFAULT", step[1], comm)
                out.setdefault("checks", []).append([step[1], int(flag)])
            elif op == "command":
                mdi.MDI_Send_Command(step[1], comm)
            elif op == "send_int":
                mdi.MDI_Send(step[1], 1, mdi.MDI_INT, comm)
            elif op == "send_ints":
                mdi.MDI_Send(step[1], len(step[1]), mdi.MDI_INT, comm)
            elif op == "send_doubles":
                values = step[1]
                mdi.MDI_Send(values, len(values), mdi.MDI_DOUBLE, comm)
            elif op == "send_nan":
                mdi.MDI_Send([float("nan"), 0.0, 0.0], 3, mdi.MDI_DOUBLE, comm)
            elif op == "recv_int":
                out.setdefault(step[1], []).append(int(mdi.MDI_Recv(1, mdi.MDI_INT, comm)))
            elif op == "recv_doubles":
                count = step[2]
                raw = mdi.MDI_Recv(count, mdi.MDI_DOUBLE, comm)
                values = [float(raw)] if count == 1 else [float(value) for value in raw]
                out.setdefault(step[1], []).append(values)
            elif op == "recv_name":
                out.setdefault(step[1], []).append(mdi.MDI_Recv(mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm))
            elif op == "exit":
                mdi.MDI_Send_Command("EXIT", comm)
            else:
                raise AssertionError(op)
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
    result_path.write_text(json.dumps({"error": error, "out": out}), encoding="utf-8")
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--driver":
        try:
            raise SystemExit(_driver_main())
        except Exception as exc:
            print(f"driver failed: {exc}", file=sys.stderr)
            raise
