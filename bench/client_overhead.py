#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Client side cost of one energy and force evaluation, per call.

cpmdc runs in process: a Python call packs positions, crosses ctypes into
``libcpmdc``, and unpacks forces, stress and dipole. This harness times that
path for N atoms with the SCF taken out, in two backends:

``mock``
    ``bench/client_overhead_mock.c``, built here into a temporary directory.
    It does one pass over the positions and one over the forces, so what
    remains is the cost of the Python client. It has no atom count limit.
``reference``
    The ``libcpmdc`` of a build without OpenCPMD, whose deterministic
    reference evaluator is what the test suite runs. It adds the real C
    session layer (topology guard, unit conversion, cache bookkeeping) and
    accepts at most 1365 atoms; larger sizes are reported as refused.

Each path is one callable per call. ``floor`` calls
``cpmdc_session_calculate_au`` through bare ctypes with prebuilt buffers, so
``path - floor`` is what the Python layer adds. A ring of four distinct
geometries is cycled so the session never answers from its last-geometry
cache. Timing is round-robin over the paths: every run times every path once
before the next run starts, so drift in machine load lands on all paths.
One run is a batch of calls; its figure is the mean time per call. The
tables give the median over runs and the spread (quartiles, minimum,
maximum).

A change is called real only when the interquartile ranges of the two sides,
pooled over the files of each side, do not overlap and the medians differ by
more than MATERIAL_RATIO. The ``floor`` row runs identical code on both sides,
so it shows how far a row can move with no change behind it.

usage: client_overhead.py [--library LIBCPMDC --params CPMDPARAMS
                           --periodic-params CPMDPARAMS]
                          [--atoms 8,64,512,4096] [--runs 15]
                          [--batch-ms 50] [--pin CPU] [--json FILE]
                          [--profile PATH]
       client_overhead.py --base-json FILE... --tip-json FILE...
"""

from __future__ import annotations

import argparse
import ctypes
import cProfile
import gc
import json
import os
import platform
import pstats
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("CPMDC_BENCH_PYTHON", str(ROOT / "python")))

import cpmdc  # noqa: E402

RING = 4
SEED = 20260101
BOHR_PER_ANGSTROM = 1.0 / cpmdc.BOHR_TO_ANGSTROM
REFERENCE_ATOM_LIMIT = 1365
ATOMS_PER_VOLUME = 0.1
MIN_EDGE_ANGSTROM = 6.0
PROBE_CALLS = 3
MATERIAL_RATIO = 1.05


class Geometry:
    """N atoms in a cubic cell, as the ring of distinct frames a run cycles."""

    def __init__(self, n_atoms: int) -> None:
        rng = np.random.default_rng(SEED + n_atoms)
        self.n_atoms = n_atoms
        edge = max(MIN_EDGE_ANGSTROM, (n_atoms / ATOMS_PER_VOLUME) ** (1.0 / 3.0))
        self.cell = np.diag([edge, edge, edge])
        base = rng.uniform(0.0, edge, size=(n_atoms, 3))
        self.frames = [base + rng.normal(0.0, 0.05, size=(n_atoms, 3)) for _ in range(RING)]
        self.numbers = np.array([8 if i % 3 == 0 else 1 for i in range(n_atoms)], dtype=np.int32)
        self.numbers_list = [int(z) for z in self.numbers]
        self.frame_lists = [frame.tolist() for frame in self.frames]
        self.flat_bohr = [np.ascontiguousarray(f * BOHR_PER_ANGSTROM).reshape(-1) for f in self.frames]
        self.flat_bohr_lists = [flat.tolist() for flat in self.flat_bohr]
        self.cell_bohr = np.ascontiguousarray(self.cell * BOHR_PER_ANGSTROM).reshape(-1)
        self.cell_list = self.cell.tolist()


class RawLibrary:
    """The three entry points the floor path calls, bound with bare ctypes."""

    def __init__(self, library: Path, params: bytes) -> None:
        self.lib = ctypes.CDLL(str(library))
        dp = ctypes.POINTER(ctypes.c_double)
        self.lib.cpmdc_session_create.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        self.lib.cpmdc_session_create.restype = ctypes.c_void_p
        self.lib.cpmdc_session_destroy.argtypes = [ctypes.c_void_p]
        self.lib.cpmdc_session_destroy.restype = None
        self.lib.cpmdc_session_calculate_au.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            dp,
            ctypes.POINTER(ctypes.c_int),
            dp,
            dp,
            dp,
            dp,
        ]
        self.lib.cpmdc_session_calculate_au.restype = ctypes.c_int
        self.buffer = ctypes.create_string_buffer(params, len(params))
        self.handle = self.lib.cpmdc_session_create(self.buffer, len(params))
        if not self.handle:
            raise RuntimeError("cpmdc_session_create failed")

    def close(self) -> None:
        self.lib.cpmdc_session_destroy(self.handle)


def floor_path(raw: RawLibrary, geo: Geometry) -> Callable[[int], object]:
    n = geo.n_atoms
    numbers = (ctypes.c_int * n)(*geo.numbers_list)
    positions = [(ctypes.c_double * (3 * n)).from_buffer_copy(f) for f in geo.flat_bohr]
    cell = (ctypes.c_double * 9).from_buffer_copy(geo.cell_bohr)
    forces = (ctypes.c_double * (3 * n))()
    virial = (ctypes.c_double * 9)()
    energy = ctypes.c_double()
    call = raw.lib.cpmdc_session_calculate_au

    def run(i: int) -> object:
        status = call(raw.handle, n, positions[i % RING], numbers, cell,
                      ctypes.byref(energy), forces, virial)
        if status != 0:
            raise RuntimeError("floor call failed")
        return energy.value

    return run


def build_paths(
    library: Path, params: bytes, periodic_params: bytes, geo: Geometry
) -> tuple[dict[str, Callable[[int], object]], list[Callable[[], None]]]:
    """The timed callables for one atom count, and the closers to run after."""
    closers: list[Callable[[], None]] = []
    paths: dict[str, Callable[[int], object]] = {}

    raw = RawLibrary(library, periodic_params)
    closers.append(raw.close)
    paths["floor"] = floor_path(raw, geo)

    cluster = cpmdc.Session(params, library=library)
    periodic = cpmdc.Session(periodic_params, library=library)
    potential = cpmdc.Potential(periodic_params, library=library)
    closers += [cluster.close, periodic.close, potential.close]
    forces_out = np.zeros(3 * geo.n_atoms)
    virial_out = np.zeros(9)

    paths["Session.calculate_au (numpy)"] = lambda i: periodic.calculate_au(
        geo.flat_bohr[i % RING], geo.numbers, geo.cell_bohr,
        forces=forces_out, virial=virial_out)
    paths["Session.calculate_au (lists)"] = lambda i: periodic.calculate_au(
        geo.flat_bohr_lists[i % RING], geo.numbers_list, geo.cell_list)
    paths["Session.energy_forces (cluster, numpy)"] = lambda i: cluster.energy_forces(
        geo.frames[i % RING], geo.numbers_list)
    paths["Session.energy_forces (periodic, numpy)"] = lambda i: periodic.energy_forces(
        geo.frames[i % RING], geo.numbers_list, geo.cell)
    paths["Session.energy_forces (periodic, lists)"] = lambda i: periodic.energy_forces(
        geo.frame_lists[i % RING], geo.numbers_list, geo.cell_list)
    paths["Potential.evaluate (periodic, numpy)"] = lambda i: potential.evaluate(
        geo.frames[i % RING], geo.numbers_list, geo.cell)

    try:
        from cpmdc.lammps import _SessionPotential
    except ImportError:
        pass
    else:
        driver = _SessionPotential(periodic_params, library, True)
        closers.append(driver.close)
        paths["lammps _SessionPotential.evaluate (periodic)"] = lambda i: driver.evaluate(
            geo.frames[i % RING], geo.numbers, geo.cell)

    try:
        from ase import Atoms

        from cpmdc.ase import CPMDCCalculator
    except ImportError:
        return paths, closers
    calc = CPMDCCalculator(params=periodic_params, library=str(library))
    closers.append(calc.close)
    atoms = Atoms(numbers=geo.numbers_list, positions=geo.frames[0], cell=geo.cell, pbc=True)
    atoms.calc = calc

    def ase_step(i: int) -> object:
        atoms.set_positions(geo.frames[i % RING])
        return (atoms.get_potential_energy(), atoms.get_forces(), atoms.get_stress())

    paths["CPMDCCalculator (energy, forces, stress)"] = ase_step
    return paths, closers


def probe(paths: dict[str, Callable[[int], object]]) -> str | None:
    """None when every path runs; otherwise the first library refusal."""
    try:
        for run in paths.values():
            for i in range(PROBE_CALLS):
                run(i)
    except (cpmdc.CPMDCError, RuntimeError) as exc:
        return str(exc)
    return None


def calibrate(run: Callable[[int], object], batch_ms: float) -> int:
    start = time.perf_counter_ns()
    count = 0
    while time.perf_counter_ns() - start < 5_000_000:
        run(count)
        count += 1
    each = (time.perf_counter_ns() - start) / count
    return max(3, int(batch_ms * 1e6 / each))


def timed_batch(run: Callable[[int], object], reps: int) -> float:
    gc.collect()
    start = time.perf_counter_ns()
    for i in range(reps):
        run(i)
    return (time.perf_counter_ns() - start) / reps / 1e3


def summarize(samples: list[float]) -> dict[str, float]:
    q = np.percentile(samples, [25, 50, 75])
    return {
        "median_us": float(q[1]),
        "q1_us": float(q[0]),
        "q3_us": float(q[2]),
        "min_us": float(min(samples)),
        "max_us": float(max(samples)),
        "iqr_rel": float((q[2] - q[0]) / q[1]),
        "runs": len(samples),
    }


def measure(
    library: Path, params: bytes, periodic_params: bytes, n_atoms: int, runs: int, batch_ms: float
) -> dict[str, object]:
    geo = Geometry(n_atoms)
    paths, closers = build_paths(library, params, periodic_params, geo)
    try:
        refusal = probe(paths)
        if refusal is not None:
            return {"refused": refusal}
        reps = {name: calibrate(run, batch_ms) for name, run in paths.items()}
        samples: dict[str, list[float]] = {name: [] for name in paths}
        for _ in range(runs):
            for name, run in paths.items():
                samples[name].append(timed_batch(run, reps[name]))
        out = {name: summarize(values) | {"batch": reps[name]} for name, values in samples.items()}
        floor = out["floor"]["median_us"]
        for stats in out.values():
            stats["over_floor_us"] = stats["median_us"] - floor
        return out
    finally:
        for close in closers:
            close()


def build_mock(directory: Path) -> Path:
    target = directory / "libcpmdc_bench_mock.so"
    compiler = os.environ.get("CC", "cc")
    subprocess.run(
        [compiler, "-O2", "-shared", "-fPIC", "-fvisibility=hidden", "-std=c11",
         str(ROOT / "bench" / "client_overhead_mock.c"), "-o", str(target)],
        check=True,
    )
    return target


def compiler_banner() -> str:
    compiler = os.environ.get("CC", "cc")
    done = subprocess.run([compiler, "--version"], check=True, capture_output=True, text=True)
    return done.stdout.splitlines()[0]


def machine_info() -> dict[str, str]:
    model = "unknown"
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {
        "cpu": model,
        "cores": str(os.cpu_count()),
        "kernel": platform.release(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "cc": compiler_banner(),
    }


def org_table(title: str, results: dict[str, dict[str, object]]) -> str:
    lines = [f"* {title}", ""]
    for n_atoms, per_path in results.items():
        lines.append(f"** N = {n_atoms}")
        lines.append("")
        if "refused" in per_path:
            lines += [f"Refused by the library: {per_path['refused']}", ""]
            continue
        lines.append("| path | median us/call | q1 | q3 | min | max | IQR/median | over floor us |")
        lines.append("|-")
        for name, st in per_path.items():
            lines.append(
                f"| {name} | {st['median_us']:.2f} | {st['q1_us']:.2f} | {st['q3_us']:.2f} "
                f"| {st['min_us']:.2f} | {st['max_us']:.2f} | {st['iqr_rel']:.3f} "
                f"| {st['over_floor_us']:.2f} |"
            )
        lines.append("")
    return "\n".join(lines)


def pooled(files: list[Path], backend: str, size: str, name: str) -> dict[str, float] | None:
    """Median of the file medians, and the widest interquartile range over the files."""
    found = []
    for path in files:
        stats = json.loads(path.read_text())["results"].get(backend, {}).get(size, {}).get(name)
        if stats is not None:
            found.append(stats)
    if not found:
        return None
    return {
        "median_us": statistics.median(item["median_us"] for item in found),
        "q1_us": min(item["q1_us"] for item in found),
        "q3_us": max(item["q3_us"] for item in found),
    }


def compare(base_files: list[Path], tip_files: list[Path]) -> str:
    lines = ["#+TITLE: client overhead, base against tip", ""]
    first = json.loads(base_files[0].read_text())
    for backend, per_size in first["results"].items():
        lines += [f"* {backend}", ""]
        for size, per_path in per_size.items():
            if "refused" in per_path:
                continue
            lines += [f"** N = {size}", "", "| path | base us | tip us | base/tip | beyond noise |", "|-"]
            for name in per_path:
                before = pooled(base_files, backend, size, name)
                after = pooled(tip_files, backend, size, name)
                if before is None or after is None:
                    continue
                ratio = before["median_us"] / after["median_us"]
                apart = (after["q3_us"] < before["q1_us"] or before["q3_us"] < after["q1_us"]) and (
                    ratio > MATERIAL_RATIO or ratio < 1.0 / MATERIAL_RATIO
                )
                lines.append(
                    f"| {name} | {before['median_us']:.2f} | {after['median_us']:.2f} "
                    f"| {ratio:.2f} | {'yes' if apart else 'no'} |"
                )
            lines.append("")
    return "\n".join(lines)


def profile_one(library: Path, params: bytes, periodic_params: bytes, n_atoms: int, name: str) -> str:
    geo = Geometry(n_atoms)
    paths, closers = build_paths(library, params, periodic_params, geo)
    try:
        if name not in paths:
            raise SystemExit(f"unknown path {name!r}; choose from {sorted(paths)}")
        run = paths[name]
        reps = calibrate(run, 200.0)
        profiler = cProfile.Profile()
        profiler.enable()
        for i in range(reps):
            run(i)
        profiler.disable()
        out = tempfile.TemporaryFile("w+")
        pstats.Stats(profiler, stream=out).sort_stats("tottime").print_stats(15)
        out.seek(0)
        return f"{reps} calls\n{out.read()}"
    finally:
        for close in closers:
            close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", type=Path, help="libcpmdc of a reference build; omit for mock only")
    parser.add_argument("--params", type=Path, help="CPMDParams message for the cluster session")
    parser.add_argument("--periodic-params", type=Path, help="CPMDParams message for the periodic session")
    parser.add_argument("--atoms", default="8,64,512,4096")
    parser.add_argument("--runs", type=int, default=15)
    parser.add_argument("--batch-ms", type=float, default=50.0)
    parser.add_argument("--pin", type=int, help="pin the process to this CPU")
    parser.add_argument("--json", type=Path)
    parser.add_argument("--profile", metavar="PATH", help="cProfile one path of the mock backend")
    parser.add_argument("--profile-atoms", type=int, default=512)
    parser.add_argument("--base-json", type=Path, nargs="+", help="result files of the base tree")
    parser.add_argument("--tip-json", type=Path, nargs="+", help="result files of the tip tree")
    args = parser.parse_args()

    if args.base_json or args.tip_json:
        if not (args.base_json and args.tip_json):
            parser.error("--base-json and --tip-json go together")
        print(compare(args.base_json, args.tip_json))
        return 0

    if args.pin is not None:
        os.sched_setaffinity(0, {args.pin})
    sizes = [int(item) for item in args.atoms.split(",")]
    if bool(args.library) != bool(args.params) or bool(args.params) != bool(args.periodic_params):
        parser.error("--library, --params and --periodic-params go together")

    with tempfile.TemporaryDirectory() as scratch:
        mock = build_mock(Path(scratch))
        blob = b"\0" * 16
        if args.profile:
            print(profile_one(mock, blob, blob, args.profile_atoms, args.profile))
            return 0
        backends: dict[str, tuple[Path, bytes, bytes]] = {"mock": (mock, blob, blob)}
        if args.library:
            backends["reference"] = (
                args.library, args.params.read_bytes(), args.periodic_params.read_bytes())
        results: dict[str, dict[str, dict[str, object]]] = {}
        for backend, (library, params, periodic_params) in backends.items():
            results[backend] = {}
            for n_atoms in sizes:
                if backend == "reference" and n_atoms > REFERENCE_ATOM_LIMIT:
                    results[backend][str(n_atoms)] = {
                        "refused": f"the reference evaluator takes at most {REFERENCE_ATOM_LIMIT} atoms"
                    }
                    continue
                print(f"{backend} N={n_atoms}", file=sys.stderr, flush=True)
                results[backend][str(n_atoms)] = measure(
                    library, params, periodic_params, n_atoms, args.runs, args.batch_ms)

    info = machine_info()
    print("#+TITLE: cpmdc client overhead per call")
    print("")
    print("| machine | " + "; ".join(f"{k}={v}" for k, v in info.items()) + " |")
    print(f"| runs per path | {args.runs}, batch about {args.batch_ms} ms |")
    print("")
    titles = {"mock": "Mock library (Python client cost only)", "reference": "Reference libcpmdc (client plus C session layer)"}
    for backend, table in results.items():
        print(org_table(titles[backend], table))
    if args.json:
        args.json.write_text(json.dumps({"machine": info, "results": results}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
