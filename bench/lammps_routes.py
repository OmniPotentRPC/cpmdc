#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Per-step overhead of coupling LAMMPS to cpmdc, in process and over i-PI sockets.

Every route integrates the same system: two water molecules in a 9 Angstrom
cubic box, metal units, NVE at 0.25 fs from 300 K velocities. The routes:

  direct        Session.calculate_au in a Python loop, no integrator
  lammps-zero   LAMMPS with pair_style zero (LAMMPS alone)
  plugin        LAMMPS with pair_style cpmdc
  fix-external  LAMMPS with cpmdc.lammps.FixExternal
  fix-socket    FixExternal whose evaluation goes through ASE's SocketServer
                to cpmdc-ipi-c over a UNIX socket; only the transport differs
                from fix-external
  ipi           i-PI integrates; cpmdc-ipi-c computes over a UNIX socket
  ipi-lammps    i-PI integrates; cpmdc-ipi-c and LAMMPS fix ipi with
                pair_style zero are both force fields, as in the socket how-to

In-process routes time ``run STEPS`` after a ``run 0`` warm-up. The i-PI
routes read i-PI's own per-step wall time (verbosity medium prints it every
100 steps) and average the windows after step 0, which holds the client
handshake; process startup and teardown stay out of every number. Each
repeat starts fresh; the summary is the median per step with every raw
value. With the reference evaluator of a build without OpenCPMD, the
numbers are coupling overhead.

usage: lammps_routes.py [--build DIR] [--steps N] [--repeats N]
                        [--routes a,b,...] [--json FILE]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))

from cpmdc import BOHR_TO_ANGSTROM, Session
from cpmdc.lammps import FixExternal

EDGE = 9.0
WATERS = [
    [1.0, 1.0, 1.1173],
    [1.0, 1.7572, 0.5308],
    [1.0, 0.2428, 0.5308],
    [5.0, 4.0, 4.1],
    [5.6, 4.7, 4.2],
    [4.3, 4.4, 3.6],
]
SYMBOLS = ["O", "H", "H", "O", "H", "H"]
NUMBERS = [8, 1, 1, 8, 1, 1]
TYPES = [1, 2, 2, 1, 2, 2]
CMDARGS = ["-log", "none", "-screen", "none", "-nocite"]
ROUTES = (
    "direct",
    "lammps-zero",
    "plugin",
    "fix-external",
    "fix-socket",
    "ipi",
    "ipi-lammps",
)


def lammps_system():
    from lammps import lammps

    lmp = lammps(cmdargs=CMDARGS)
    lmp.commands_string(
        "\n".join(
            [
                "units metal",
                "atom_modify map array sort 0 0",
                "boundary p p p",
                f"region box block 0 {EDGE} 0 {EDGE} 0 {EDGE}",
                "create_box 2 box",
                "mass 1 15.999",
                "mass 2 1.008",
            ]
        )
    )
    flat = np.asarray(WATERS).ravel().tolist()
    lmp.create_atoms(len(WATERS), list(range(1, len(WATERS) + 1)), TYPES, flat)
    return lmp


DYNAMICS = [
    "velocity all create 300 4928 loop geom",
    "fix nve all nve",
    "timestep 0.00025",
    "thermo 0",
]


class SocketFix(FixExternal):
    """FixExternal that evaluates through an i-PI socket server instead of a Session."""

    def __init__(self, lmp, params, elements, client, library=None) -> None:
        from ase.calculators.socketio import SocketServer

        name = "cpmdc-bench-" + uuid.uuid4().hex[:8]
        self._server = SocketServer(unixsocket=name, timeout=60)
        self._client = subprocess.Popen(
            [
                str(client),
                "--unix",
                name,
                "--params",
                str(params),
                "--species",
                " ".join(SYMBOLS),
                "--wait",
                "60",
            ]
        )
        super().__init__(lmp, params, elements, library=library)

    def _evaluate(self, tags, xu, numbers, rows):
        index = tags - 1
        positions = np.empty_like(xu)
        positions[index] = xu
        if self._server.protocol is None:
            self._server._accept()
        out = self._server.protocol.calculate(positions, rows)
        return (
            out["energy"],
            np.asarray(out["forces"])[index],
            np.asarray(out["virial"]),
        )

    def close(self) -> None:
        super().close()
        self._server.close()
        self._client.wait(timeout=60)


def time_lammps(args, route: str) -> float:
    lmp = lammps_system()
    driver = None
    if route == "lammps-zero":
        lmp.commands_string("pair_style zero 3.0\npair_coeff * *")
    elif route == "plugin":
        lmp.command(f"plugin load {args.plugin}")
        lmp.command(f"pair_style cpmdc {args.params}")
        lmp.command("pair_coeff * * O H")
    elif route == "fix-external":
        driver = FixExternal(lmp, args.params, "O H", library=args.library)
    elif route == "fix-socket":
        driver = SocketFix(lmp, args.params, "O H", args.client, library=args.library)
    lmp.commands_string("\n".join(DYNAMICS))
    run = driver.run if driver is not None else (lambda n: lmp.command(f"run {n}"))
    run(0)
    start = time.perf_counter()
    run(args.steps)
    elapsed = time.perf_counter() - start
    if driver is not None:
        driver.close()
    lmp.close()
    return elapsed / args.steps


def time_direct(args) -> float:
    session = Session(args.params, library=args.library)
    positions = np.asarray(WATERS) / BOHR_TO_ANGSTROM
    cell = np.eye(3) * EDGE / BOHR_TO_ANGSTROM
    forces = np.zeros(positions.size)
    virial = np.zeros(9)
    numbers = np.asarray(NUMBERS, dtype=np.int32)
    session.calculate_au(positions, numbers, cell, forces=forces, virial=virial)
    start = time.perf_counter()
    for _ in range(args.steps):
        session.calculate_au(positions, numbers, cell, forces=forces, virial=virial)
    elapsed = time.perf_counter() - start
    session.close()
    return elapsed / args.steps


IPI_INPUT = """<simulation verbosity="medium" safe_stride="100000000">
  <output prefix="sim">
    <properties stride="100000000" filename="out"> [ step ] </properties>
  </output>
  <total_steps> {steps} </total_steps>
  <prng><seed>31415</seed></prng>
  <ffsocket name="cpmdc" mode="unix"><address> {cpmdc} </address></ffsocket>
  {lammps_socket}
  <system>
    <initialize nbeads="1">
      <file mode="xyz"> init.xyz </file>
      <velocities mode="thermal" units="kelvin"> 300 </velocities>
    </initialize>
    <forces><force forcefield="cpmdc"/>{lammps_force}</forces>
    <ensemble><temperature units="kelvin"> 300 </temperature></ensemble>
    <motion mode="dynamics">
      <dynamics mode="nve"><timestep units="femtosecond"> 0.25 </timestep></dynamics>
    </motion>
  </system>
</simulation>
"""

LAMMPS_IPI_INPUT = """units metal
atom_style atomic
atom_modify map array sort 0 0
boundary p p p
region box block 0 {edge} 0 {edge} 0 {edge}
create_box 2 box
mass 1 15.999
mass 2 1.008
{atoms}
pair_style zero 3.0
pair_coeff * *
thermo 0
fix ipi all ipi {address} 0 unix
run 100000000
"""


def wait_for_sockets(
    server: subprocess.Popen, names: list[str], log: list[str]
) -> None:
    """Read i-PI's log until every forcefield polls, then keep reading it into ``log``."""
    pending = {f"@ForceField ({name}): Starting the polling thread" for name in names}
    for line in server.stdout:
        log.append(line)
        pending = {mark for mark in pending if mark not in line}
        if not pending:
            break
    else:
        raise RuntimeError("i-PI stopped before its sockets listened:\n" + "".join(log))
    threading.Thread(target=lambda: log.extend(server.stdout), daemon=True).start()


TIMING = re.compile(r"Average timings at MD step\s+(\d+)\. t/step:\s+(\S+)")


def time_ipi(args, route: str) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        return ipi_run(args, route, args.steps, Path(tmp))


def ipi_run(args, route: str, steps: int, work: Path) -> float:
    """i-PI's own wall time per step, over every 100-step window after the first."""
    tag = uuid.uuid4().hex[:8]
    cpmdc_name = f"cpmdc-{tag}"
    lammps_name = f"lmp-{tag}"
    with_lammps = route == "ipi-lammps"
    lines = [str(len(WATERS))]
    lines.append(
        f"# CELL(abcABC): {EDGE} {EDGE} {EDGE} 90.0 90.0 90.0 positions{{angstrom}} cell{{angstrom}}"
    )
    lines += [f"{s} {x!r} {y!r} {z!r}" for s, (x, y, z) in zip(SYMBOLS, WATERS)]
    (work / "init.xyz").write_text("\n".join(lines) + "\n")
    (work / "input.xml").write_text(
        IPI_INPUT.format(
            steps=steps,
            cpmdc=cpmdc_name,
            lammps_socket=(
                f'<ffsocket name="lammps" mode="unix"><address> {lammps_name} </address></ffsocket>'
                if with_lammps
                else ""
            ),
            lammps_force='<force forcefield="lammps"/>' if with_lammps else "",
        )
    )
    procs = []
    log: list[str] = []
    server = subprocess.Popen(
        [args.ipi, "input.xml"],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )
    wait_for_sockets(server, ["cpmdc", "lammps"] if with_lammps else ["cpmdc"], log)
    command = [str(args.client), "--unix", cpmdc_name, "--params", str(args.params)]
    command += ["--species", " ".join(SYMBOLS), "--wait", "60"]
    procs.append(subprocess.Popen(command, cwd=work))
    if with_lammps:
        atoms = "\n".join(
            f"create_atoms {t} single {x!r} {y!r} {z!r}"
            for t, (x, y, z) in zip(TYPES, WATERS)
        )
        (work / "in.lammps").write_text(
            LAMMPS_IPI_INPUT.format(edge=EDGE, atoms=atoms, address=lammps_name)
        )
        procs.append(
            subprocess.Popen(
                [args.lmp, "-in", "in.lammps", *CMDARGS],
                cwd=work,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        )
    try:
        if server.wait(timeout=1800) != 0:
            raise RuntimeError(f"i-PI exited with {server.returncode}")
        for proc in procs:
            proc.wait(timeout=60)
    finally:
        for proc in [server, *procs]:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
    windows = [
        float(m.group(2)) for m in map(TIMING.search, log) if m and int(m.group(1)) > 0
    ]
    if not windows:
        raise RuntimeError("i-PI printed no step timings:\n" + "".join(log[-20:]))
    return statistics.fmean(windows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--params", type=Path)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--plugin", type=Path)
    parser.add_argument("--client", type=Path)
    parser.add_argument("--ipi", default=shutil.which("i-pi"))
    parser.add_argument("--lmp", default=shutil.which("lmp"))
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--routes", default=",".join(ROUTES))
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    args.params = args.params or args.build / "water_blyp_periodic_params.bin"
    args.library = args.library or args.build / "libcpmdc.so"
    args.plugin = args.plugin or args.build / "cpmdcplugin.so"
    args.client = args.client or args.build / "cpmdc-ipi-c"
    os.environ.setdefault("CPMDC_LIBRARY", str(args.library))

    results = {}
    for route in args.routes.split(","):
        if route not in ROUTES:
            parser.error(f"unknown route {route!r}; choose from {', '.join(ROUTES)}")
        raw = []
        for _ in range(args.repeats):
            if route == "direct":
                raw.append(time_direct(args))
            elif route in ("ipi", "ipi-lammps"):
                raw.append(time_ipi(args, route))
            else:
                raw.append(time_lammps(args, route))
        results[route] = {
            "median_us": statistics.median(raw) * 1e6,
            "raw_us": [v * 1e6 for v in raw],
        }
        print(
            f"{route:13s} median {results[route]['median_us']:9.2f} us/step  raw "
            + " ".join(f"{v:.2f}" for v in results[route]["raw_us"]),
            flush=True,
        )
    if args.json:
        meta = {"steps": args.steps, "repeats": args.repeats, "atoms": len(WATERS)}
        args.json.write_text(
            json.dumps({"setup": meta, "routes": results}, indent=2) + "\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
