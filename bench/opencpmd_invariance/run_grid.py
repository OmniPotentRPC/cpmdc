#!/usr/bin/env python3
"""Run the OpenCPMD rank/thread grid and write one CSV row set per run.

Evidence only. Each run uses a fresh directory so RESTART is not reused.
"""

from __future__ import annotations

import argparse
import os
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
AU_KB = 294210.1080
ENERGIES = [
    "etot", "ekin", "epseu", "enl", "eht", "ehep", "ehee", "ehii", "exc", "vxc",
    "egc", "esr", "eself", "eeig", "eband", "entropy", "eext", "ehsic", "ecnstr",
    "erestr", "eefield", "ebogo", "amu", "etddft",
]


def parse_steps(text: str) -> tuple[int, float, str]:
    """Return SCF iterations, median step seconds, and a short status."""
    times = []
    last_nfi = 0
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != 6:
            continue
        try:
            nfi = int(parts[0])
            tcpu = float(parts[5])
        except ValueError:
            continue
        if nfi < last_nfi:
            continue
        last_nfi = nfi
        times.append(tcpu)
    deltas = []
    prev = 0.0
    for t in times:
        deltas.append(t - prev)
        prev = t
    wall = statistics.median(deltas) if deltas else float("nan")
    status = "ok"
    if "NO CONVERGENCE" in text:
        status = "maxstep"
    if "PROGRAM STOPS" in text or "STOPGM" in text:
        status = "stop"
    return last_nfi, wall, status


def step_wall(path: Path, fallback: float) -> float:
    """Median SCF step seconds.

    The normal wavefunction loop resets time1 at each iteration, so each
    logged TCPU value is already that step. The printed column is 0.01 s.
    """
    import struct

    if not path.exists():
        return fallback
    data = path.read_bytes()
    n = len(data) // 8
    if n == 0:
        return fallback
    times = struct.unpack("<" + "d" * n, data[: n * 8])
    return statistics.median(times) if times else fallback


def read_evidence(path: Path) -> dict:
    import struct

    data = path.read_bytes()
    off = 0

    def i32() -> int:
        nonlocal off
        val = struct.unpack_from("<i", data, off)[0]
        off += 4
        return val

    def f64() -> float:
        nonlocal off
        val = struct.unpack_from("<d", data, off)[0]
        off += 8
        return val

    def raw64() -> bytes:
        nonlocal off
        blob = data[off : off + 8]
        off += 8
        return blob

    magic = i32()
    version = i32()
    nat = i32()
    nener = i32()
    if magic != 0x45564431 or version != 1:
        raise RuntimeError(f"bad evidence header {magic:#x} {version}")
    energies = []
    for _ in range(nener):
        name = data[off : off + 16].decode("ascii", "replace").strip()
        off += 16
        blob = raw64()
        val = struct.unpack("<d", blob)[0]
        energies.append((name, val, blob))
    nforce = i32()
    forces = []
    for _ in range(nforce):
        blob = raw64()
        forces.append((struct.unpack("<d", blob)[0], blob))
    nstress = i32()
    paiu = []
    for _ in range(nstress):
        blob = raw64()
        paiu.append((struct.unpack("<d", blob)[0], blob))
    omega = f64()
    gemax = f64()
    cnorm = f64()
    nproc = i32()
    planes = []
    for _ in range(nproc):
        planes.append((i32(), i32(), i32()))
    return {
        "nat": nat,
        "energies": energies,
        "forces": forces,
        "paiu": paiu,
        "omega": omega,
        "gemax": gemax,
        "cnorm": cnorm,
        "planes": planes,
    }


def hex_of(blob: bytes) -> str:
    return blob[::-1].hex().upper()


def csv_rows(system: str, ranks: int, threads: int, rep: int, variant: str, iters: int, wall: float, ev: dict) -> list[str]:
    rows = []

    def add(quantity: str, index: int, value: float, blob: bytes) -> None:
        rows.append(
            ",".join(
                [
                    system,
                    str(ranks),
                    str(threads),
                    str(rep),
                    quantity,
                    str(index),
                    hex_of(blob),
                    f"{value:.16e}",
                    str(iters),
                    f"{wall:.8f}",
                    variant,
                ]
            )
        )

    import struct

    for i, (name, value, blob) in enumerate(ev["energies"]):
        add(name if name else ENERGIES[i], i, value, blob)
    for i, (value, blob) in enumerate(ev["forces"]):
        add("force", i, value, blob)
    omega = ev["omega"]
    for i, (value, blob) in enumerate(ev["paiu"]):
        add("stress_paiu", i, value, blob)
        gpa = 0.0 if omega == 0.0 else value / omega * AU_KB * 0.1
        add("stress_gpa", i, gpa, struct.pack("<d", gpa))
    for r, (ngw, lo, hi) in enumerate(ev["planes"]):
        add("ngw", r, float(ngw), struct.pack("<d", float(ngw)))
        add("nrxpl_lo", r, float(lo), struct.pack("<d", float(lo)))
        add("nrxpl_hi", r, float(hi), struct.pack("<d", float(hi)))
    return rows


def launch(bin_path: Path, inp: Path, pp: Path, ranks: int, threads: int, rundir: Path, ncpu: int) -> subprocess.CompletedProcess:
    rundir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_WAIT_POLICY"] = "PASSIVE"
    env["CPMDC_EVIDENCE_BIN"] = str(rundir / "evidence.bin")
    env["CPMDC_EVIDENCE_STEPS"] = str(rundir / "steps.bin")
    over = ranks * threads > ncpu
    cmd = ["mpirun", "--bind-to", "none" if over else "core"]
    if over:
        cmd.append("--oversubscribe")
    else:
        cmd.extend(["--map-by", f"slot:PE={threads}"])
    cmd.extend(["-np", str(ranks), str(bin_path), str(inp), str(pp)])
    (rundir / "cmd.txt").write_text(" ".join(cmd) + f"\noversubscribe={int(over)}\n", encoding="utf-8")
    return subprocess.run(cmd, cwd=rundir, env=env, capture_output=True, text=True, timeout=3600)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", default="stock")
    parser.add_argument("--systems", default="si8,water")
    parser.add_argument("--ranks", default="1,2,3,4")
    parser.add_argument("--threads", default="1,2,8")
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    variants = args.variants.split(",")
    systems = args.systems.split(",")
    ranks = [int(x) for x in args.ranks.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    ncpu = os.cpu_count() or 1
    pp = ROOT / "pp"
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for rep in range(args.reps):
        for system in systems:
            for rank in ranks:
                for thread in threads:
                    for variant in variants:
                        jobs.append((system, rank, thread, rep, variant))
    header = "system,ranks,threads,rep,quantity,index,value_hex,value,iters,wall_step_s,variant\n"
    for system, rank, thread, rep, variant in jobs:
        shard = out / f"{system}_{variant}.csv"
        if not shard.exists():
            shard.write_text(header, encoding="utf-8")
        tag = f"{system} r{rank} t{thread} rep{rep} {variant}"
        rundir = out / "runs" / variant / system / f"r{rank}_t{thread}_rep{rep}"
        done = rundir / "evidence.bin"
        if done.exists() and (rundir / "logged.txt").exists():
            print("skip", tag, flush=True)
            continue
        binary = ROOT / "bin" / f"cpmd.x.{variant}"
        if not binary.exists():
            print("missing", binary, flush=True)
            sys.exit(1)
        inp = ROOT / "inputs" / f"{system}.inp"
        print("run", tag, flush=True)
        try:
            proc = launch(binary, inp, pp, rank, thread, rundir, ncpu)
        except subprocess.TimeoutExpired as exc:
            (rundir / "stdout.txt").write_text(exc.stdout or "", encoding="utf-8")
            (rundir / "stderr.txt").write_text((exc.stderr or "") + "\nTIMEOUT\n", encoding="utf-8")
            print("timeout", tag, flush=True)
            continue
        (rundir / "stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
        (rundir / "stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
        iters, wall, status = parse_steps(proc.stdout or "")
        wall = step_wall(rundir / "steps.bin", wall)
        (rundir / "status.txt").write_text(f"{status} exit={proc.returncode} iters={iters} wall={wall}\n", encoding="utf-8")
        if proc.returncode != 0 or not done.exists():
            print("fail", tag, "exit", proc.returncode, status, flush=True)
            tail = "\n".join((proc.stdout or "").splitlines()[-30:])
            (rundir / "tail.txt").write_text(tail, encoding="utf-8")
            continue
        ev = read_evidence(done)
        rows = csv_rows(system, rank, thread, rep, variant, iters, wall, ev)
        with shard.open("a", encoding="utf-8") as fh:
            fh.write("\n".join(rows) + "\n")
        (rundir / "logged.txt").write_text("ok\n", encoding="utf-8")
        print("ok", tag, "iters", iters, "wall", f"{wall:.4f}", status, flush=True)


if __name__ == "__main__":
    main()
