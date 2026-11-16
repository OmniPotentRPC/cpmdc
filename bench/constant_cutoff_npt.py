#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Density of liquid water from NPT dynamics through CPMDCCalculator.

start writes the box a run begins from: WATERS molecules at DENSITY g/cm^3,
oxygens on the sites of a simple cubic lattice taken in checkerboard order,
seeded random orientations, then PS ps of Langevin dynamics (0.5 fs steps,
friction 0.01/fs) at TEMPERATURE in a classical flexible model: the TIP3P
charges and Lennard-Jones of ASE's TIP3P, with harmonic O-H bonds and H-O-H
angle that take the SPC/Fw force constants about the TIP3P geometry.

run integrates hydrostatic Langevin NPT (ASE's LangevinBAOAB) at PRESSURE
bar and TEMPERATURE K through CPMDCCalculator: PBE, the Troullier-Martins
O_MT_PBE and H_MT_PBE pseudopotentials of the CPMD library, basis cutoff
CUTOFF, orbital convergence CONV, and orbitals predicted by ASPC of order 4.
--e0 adds the constant-cutoff functional at E0 (A = 50 x CUTOFF and sigma =
E0 / 10 unless --height and --width are given); without it each new cell is
a new setup at CUTOFF, whose stress is that of a fixed number of plane
waves. Hydrogen carries --hydrogen-mass, deuterium by default, so that 1 fs
steps resolve the O-H stretch; configurational averages such as the density
do not depend on the masses. The density is that of H2O, the water count
times 18.015 u over the volume, and the pressure is the instantaneous one
that the barostat holds, with the kinetic term. A run stops when the density
leaves --density-range, as it does when the cell collapses. --start takes the
box of start, whose momenta the run draws at TEMPERATURE, or the state file of
another run, whose positions, momenta, cell and barostat momentum it
continues from with its own settings; the record then holds the state's step
count times DT as start_ps, where the figure draws the run.

Each step appends the time, volume, density, pressure, temperature,
potential energy, wall time, SCF iterations and new setups to the record
OUT. The last two come from the CPMDC_SCF_STEPS lines of the error stream,
which rank 0 sends to OUT with the suffix .err. Every --checkpoint steps
rank 0 rewrites OUT and the state file (suffix .state.npz), and --resume
continues the run from them with the settings recorded in OUT.
CPMDC_CELL_IN_PLACE applies as set.

pressure evaluates one configuration, a box or the state file of a run, at
each basis it is given, CUTOFF for a fixed cutoff or CUTOFF/E0[/WIDTH] for a
constant cutoff with the default A, and writes an Org table of the energy and
of two virial pressures without the kinetic term: that of the atoms, and that
of the molecular centres of mass (with the masses of the state file), which
leaves out the intramolecular forces.

table writes Org tables of the block-averaged volume, density and pressure of
each run (bench/npt_block_average.py: the detected equilibration discarded,
at least --min-discard-ps of a run that starts from a box and
--branch-discard-ps of one that starts from a state, Flyvbjerg-Petersen
standard error); of the runs that stopped or hold less than
--min-production-ps after that, with the density at their start and end; of
the cost per step; and of the density every --branch-step-ps after the start
of each run that starts from a state. --json also writes the means and the
window means that docs/render_constant_cutoff_npt.py draws.

capnp must be on PATH and CPMDC_PSEUDO_DIR must name the pseudopotential
directory. Under mpirun with CPMDC_RANK_CONSENSUS=1 every rank integrates the
same trajectory and rank 0 writes.

usage: constant_cutoff_npt.py start --out BOX.xyz [--waters 32] [--density 1.0]
                                    [--ps 10] [--temperature 350] [--seed 7]
       constant_cutoff_npt.py run --library LIB --out RUN.json
                                  (--start BOX.xyz|RUN.state.npz --cutoff RY [--e0 RY]
                                   [--height RY] [--width RY] [--conv 1e-6]
                                   [--dt 1] [--temperature 350] [--pressure 1]
                                   [--t-tau 100] [--p-tau 500] [--seed 11]
                                   | --resume) [--steps 3000] [--checkpoint 10]
       constant_cutoff_npt.py pressure --library LIB --start BOX.xyz|RUN.state.npz
                                       [--conv 1e-6] [--output FILE] BASIS [BASIS ...]
       constant_cutoff_npt.py table RUN.json [RUN.json ...] [--json OUT]
                                    [--discard-ps PS] [--min-discard-ps 0.5]
                                    [--branch-discard-ps 0.05] [--min-production-ps 0.1]
                                    [--window-ps 0.1] [--plot-window-ps 0.02]
                                    [--branch-step-ps 0.01] [--output FILE]
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import re
import subprocess
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "bench"))

import ase  # noqa: E402
from ase import Atoms, units  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402
from ase.calculators.mixing import SumCalculator  # noqa: E402
from ase.calculators.tip3p import TIP3P, angleHOH, rOH  # noqa: E402
from ase.data import atomic_masses  # noqa: E402
from ase.io import read, write  # noqa: E402
from ase.md.langevin import Langevin  # noqa: E402
from ase.md.langevinbaoab import LangevinBAOAB  # noqa: E402
from ase.md.velocitydistribution import Stationary, thermalize_momenta  # noqa: E402

from constant_cutoff_scan import cpu_model  # noqa: E402
from cpmdc import find_schema  # noqa: E402
from cpmdc.ase import CPMDCCalculator  # noqa: E402
from npt_block_average import analyse, curve_rows, detail_rows, summary_rows  # noqa: E402

WATER_MASS = float(atomic_masses[8] + 2.0 * atomic_masses[1])
BOND_K = 1059.162 * units.kcal / units.mol
ANGLE_K = 75.90 * units.kcal / units.mol
THETA0 = np.radians(angleHOH)
FIELDS = ("time_ps", "Volume", "Density", "Pressure", "Temperature", "Energy", "step_wall_s", "scf_steps",
          "setups")
PHYSICS = ("cutoff", "e0", "height", "width", "conv", "dt", "temperature", "pressure", "t_tau", "p_tau",
           "hydrogen_mass", "seed", "density_range")


class Intramolecular(Calculator):
    """Harmonic O-H bonds and H-O-H angle of molecules stored as O, H, H."""

    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        r = self.atoms.positions.reshape(-1, 3, 3)
        cell = self.atoms.cell.array
        inverse = np.linalg.inv(cell)

        def nearest(d: np.ndarray) -> np.ndarray:
            s = d @ inverse
            return (s - np.round(s)) @ cell

        forces = np.zeros_like(r)
        energy = 0.0
        bonds = [nearest(r[:, j] - r[:, 0]) for j in (1, 2)]
        for j, bond in zip((1, 2), bonds):
            length = np.linalg.norm(bond, axis=1)
            energy += 0.5 * BOND_K * float(np.sum((length - rOH) ** 2))
            g = (BOND_K * (length - rOH) / length)[:, None] * bond
            forces[:, j] -= g
            forces[:, 0] += g
        n1 = np.linalg.norm(bonds[0], axis=1)
        n2 = np.linalg.norm(bonds[1], axis=1)
        cos = np.sum(bonds[0] * bonds[1], axis=1) / (n1 * n2)
        theta = np.arccos(np.clip(cos, -1.0, 1.0))
        energy += 0.5 * ANGLE_K * float(np.sum((theta - THETA0) ** 2))
        scale = -ANGLE_K * (theta - THETA0) / np.sqrt(np.maximum(1.0 - cos * cos, 1e-12))
        g1 = scale[:, None] * (bonds[1] / (n1 * n2)[:, None] - (cos / n1**2)[:, None] * bonds[0])
        g2 = scale[:, None] * (bonds[0] / (n1 * n2)[:, None] - (cos / n2**2)[:, None] * bonds[1])
        forces[:, 1] -= g1
        forces[:, 2] -= g2
        forces[:, 0] += g1 + g2
        self.results = {"energy": energy, "forces": forces.reshape(-1, 3)}


def water_box(waters: int, density: float, rng: np.random.Generator) -> Atoms:
    side = (waters * WATER_MASS / units.mol / density * 1e24) ** (1.0 / 3.0)
    m = int(np.ceil(waters ** (1.0 / 3.0) - 1e-9))
    sites = sorted(product(range(m), repeat=3), key=lambda s: (sum(s) % 2, s))[:waters]
    x, y = rOH * np.sin(THETA0 / 2.0), rOH * np.cos(THETA0 / 2.0)
    molecule = np.array([[0.0, 0.0, 0.0], [x, y, 0.0], [-x, y, 0.0]])
    positions = []
    for site in sites:
        q = rng.normal(size=4)
        a, b, c, d = q / np.linalg.norm(q)
        rotation = np.array([
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],
        ])
        positions.extend((np.asarray(site) + 0.5) * side / m + molecule @ rotation.T)
    return Atoms("OH2" * waters, positions=positions, cell=[side] * 3, pbc=True)


def start(args: argparse.Namespace) -> int:
    rng = np.random.default_rng(args.seed)
    atoms = water_box(args.waters, args.density, rng)
    atoms.calc = SumCalculator([TIP3P(rc=atoms.cell[0, 0] / 2.0 - 0.3), Intramolecular()])
    thermalize_momenta(atoms, args.temperature, rng=rng)
    Stationary(atoms)
    dyn = Langevin(atoms, 0.5 * units.fs, temperature_K=args.temperature, friction=0.01 / units.fs, fixcm=False,
                   rng=rng)
    dyn.run(int(round(args.ps * 2000)))
    print(f"{args.waters} waters in a {atoms.cell[0, 0]:.4f} Angstrom cube, {atoms.get_temperature():.0f} K "
          f"after {args.ps:g} ps", flush=True)
    atoms.calc = None
    atoms.wrap()
    write(args.out, atoms, format="extxyz")
    return 0


def deck(cutoff: float, side: float, conv: float) -> bytes:
    text = f"""(
  functional = "PBE", cutOffRy = {cutoff:.2f}, charge = 0, multiplicity = 1, task = "gradient",
  inputSections = [
    ( cpmd = ( optimizeWavefunction = true, convergenceOrbitals = {conv:g}, maxStep = 200, maxIter = 200,
        extrapolateWfnOptions = "ASPC", extrapolateWfnOrder = 4,
        directives = [ ( keyword = "CENTER MOLECULE OFF", args = [] ) ] ) ),
    ( dft = ( functional = "PBE", lsd = false, directives = [ ( keyword = "OLDCODE", args = [] ) ] ) ),
    ( system = ( symmetry = 1, angstrom = true, cutOffRy = {cutoff:.2f},
        cell = [{side:.6f}, 1.0, 1.0, 0.0, 0.0, 0.0] ) ),
    ( atoms = ( pseudopotentials = [
        ( element = "O", path = "O_MT_PBE.psp", lmax = 1 ),
        ( element = "H", path = "H_MT_PBE.psp", lmax = 0 ) ] ) ) ] )
"""
    schema = find_schema(ROOT / "schema" / "Potentials.capnp")
    done = subprocess.run(["capnp", "encode", str(schema), "CPMDParams"], input=text.encode("utf-8"),
                          stdout=subprocess.PIPE, check=True)
    return done.stdout


def mpi_rank() -> int:
    for name in ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK"):
        if name in os.environ:
            return int(os.environ[name])
    return 0


def mpi_size() -> int:
    for name in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE"):
        if name in os.environ:
            return int(os.environ[name])
    return 1


def finalize_mpi() -> None:
    """MPI_Finalize through the libmpi that CPMD loaded; cpmdc does not call it."""
    try:
        maps = Path("/proc/self/maps").read_text().split()
    except OSError:
        return
    for path in sorted({word for word in maps if re.search(r"/libmpi\.so[.0-9]*$", word)}):
        mpi = ctypes.CDLL(path)
        started, finished = ctypes.c_int(0), ctypes.c_int(0)
        mpi.MPI_Initialized(ctypes.byref(started))
        mpi.MPI_Finalized(ctypes.byref(finished))
        if started.value and not finished.value:
            mpi.MPI_Finalize()


class ScfLines:
    """The CPMDC_SCF_STEPS lines of each call, from the error stream sent to a file."""

    def __init__(self, path: Path) -> None:
        sys.stderr.flush()
        console = os.dup(2)
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
        os.dup2(fd, 2)
        os.close(fd)
        sys.stderr = os.fdopen(console, "w", buffering=1)
        self.reader = path.open(encoding="utf-8", errors="replace")
        self.reader.seek(0, os.SEEK_END)

    def take(self) -> tuple[int, int]:
        """SCF iterations, and evaluations that started with nothing stored, since the last take."""
        text = self.reader.read()
        steps = sum(int(n) for n in re.findall(r"^cpmdc_scf_steps\s+(\d+)", text, re.M))
        setups = len(re.findall(r"^cpmdc_warm_start\s+-2\b", text, re.M))
        return steps, setups


class LeftRange(Exception):
    """The density left --density-range."""


def write_atomic(path: Path, write_to) -> None:
    tmp = path.with_name(path.name + ".tmp")
    write_to(tmp)
    os.replace(tmp, path)


def save(out: Path, record: dict, atoms: Atoms, dyn: LangevinBAOAB) -> None:
    write_atomic(out, lambda p: p.write_text(json.dumps(record) + "\n"))

    def state(p: Path) -> None:
        with p.open("wb") as handle:
            np.savez(handle, numbers=atoms.numbers, positions=atoms.positions, momenta=atoms.get_momenta(),
                     masses=atoms.get_masses(), cell=atoms.cell.array, p_eps=np.asarray(dyn.p_eps),
                     nsteps=dyn.nsteps, rng=json.dumps(dyn.rng.bit_generator.state))

    write_atomic(out.with_suffix(".state.npz"), state)


def run(args: argparse.Namespace) -> int:
    rank = mpi_rank()
    out = args.out
    if rank == 0:
        out.parent.mkdir(parents=True, exist_ok=True)
    os.environ["CPMDC_SCF_STEPS"] = "1"
    scf = ScfLines(out.with_suffix(".err")) if rank == 0 else None
    if args.resume:
        record = json.loads(out.read_text())
        setup = argparse.Namespace(**record["args"])
        saved = np.load(out.with_suffix(".state.npz"))
        atoms = Atoms(numbers=saved["numbers"], positions=saved["positions"], cell=saved["cell"], pbc=True)
        atoms.set_masses(saved["masses"])
        atoms.set_momenta(saved["momenta"])
    else:
        if args.start is None or args.cutoff is None:
            raise SystemExit("run needs --start and --cutoff, or --resume")
        setup = args
        start_ps = 0.0
        if args.start.suffix == ".npz":
            saved = np.load(args.start)
            atoms = Atoms(numbers=saved["numbers"], positions=saved["positions"], cell=saved["cell"], pbc=True)
            atoms.set_masses(saved["masses"])
            atoms.set_momenta(saved["momenta"])
            start_ps = int(saved["nsteps"]) * args.dt * 1e-3
        else:
            atoms = read(args.start)
            masses = atoms.get_masses()
            masses[atoms.numbers == 1] = args.hydrogen_mass
            atoms.set_masses(masses)
            thermalize_momenta(atoms, args.temperature, rng=np.random.default_rng(args.seed))
            Stationary(atoms)
        waters = int(np.sum(atoms.numbers == 8))
        record = {
            "args": {**{name: getattr(args, name) for name in PHYSICS}, "start": args.start.name,
                     "start_ps": start_ps, "waters": waters, "ranks": mpi_size(), "cpu": cpu_model(),
                     "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
                     "cell_in_place": os.environ.get("CPMDC_CELL_IN_PLACE", ""), "ase": ase.__version__},
            **{name: [] for name in FIELDS},
            "stopped": None,
        }
    waters = record["args"]["waters"]
    cc = None if setup.e0 is None else {"cutoff_ry": setup.e0, "height_ry": setup.height,
                                        "width_ry": setup.width}
    calc = CPMDCCalculator(deck(setup.cutoff, float(atoms.cell[0, 0]), setup.conv), library=str(args.library),
                           constant_cutoff=cc)
    atoms.calc = calc
    clock = [time.perf_counter()]
    dyn = LangevinBAOAB(atoms, setup.dt * units.fs, temperature_K=setup.temperature,
                        externalstress=-setup.pressure * units.bar, hydrostatic=True,
                        T_tau=setup.t_tau * units.fs, P_tau=setup.p_tau * units.fs,
                        rng=np.random.default_rng(setup.seed + 1))
    if args.resume:
        dyn.p_eps = float(saved["p_eps"])
        dyn.nsteps = int(saved["nsteps"])
        dyn.rng.bit_generator.state = json.loads(str(saved["rng"]))
    elif args.start.suffix == ".npz":
        dyn.p_eps = float(saved["p_eps"])
    lo, hi = setup.density_range

    def frame() -> None:
        now = time.perf_counter()
        steps, setups = scf.take() if scf else (0, 0)
        volume = atoms.get_volume()
        density = waters * WATER_MASS / units.mol / volume * 1e27
        pressure = -np.trace(atoms.get_stress(voigt=False, include_ideal_gas=True)) / 3.0 / units.bar
        row = (dyn.nsteps * setup.dt * 1e-3, volume * 1e-3, density, pressure, atoms.get_temperature(),
               atoms.get_potential_energy(), now - clock[0], steps, setups)
        clock[0] = now
        for name, value in zip(FIELDS, row):
            record[name].append(float(value) if name not in ("scf_steps", "setups") else int(value))
        if rank == 0 and (dyn.nsteps % args.checkpoint == 0 or not lo <= density <= hi):
            save(out, record, atoms, dyn)
        if not lo <= density <= hi:
            raise LeftRange(f"density {density:.0f} kg/m^3 at {row[0]:g} ps")

    dyn.attach(frame, interval=1)
    try:
        dyn.run(max(0, args.steps - dyn.nsteps))
    except LeftRange as exc:
        record["stopped"] = str(exc)
    if rank == 0:
        save(out, record, atoms, dyn)
    calc.close()
    finalize_mpi()
    return 0


def load_configuration(path: Path) -> Atoms:
    if path.suffix != ".npz":
        return read(path)
    saved = np.load(path)
    atoms = Atoms(numbers=saved["numbers"], positions=saved["positions"], cell=saved["cell"], pbc=True)
    atoms.set_masses(saved["masses"])
    return atoms


def molecular_virial(atoms: Atoms, forces: np.ndarray) -> float:
    """sum over atoms of (r - R) . F, eV, with R the centre of mass of the atom's O, H, H."""
    cell = atoms.cell.array
    inverse = np.linalg.inv(cell)
    r = atoms.positions.reshape(-1, 3, 3)
    d = r - r[:, :1]
    d = d - np.round(d @ inverse) @ cell
    m = atoms.get_masses().reshape(-1, 3, 1)
    rel = d - np.sum(m * d, axis=1, keepdims=True) / np.sum(m, axis=1, keepdims=True)
    return float(np.sum(rel * forces.reshape(-1, 3, 3)))


def pressure(args: argparse.Namespace) -> int:
    rank = mpi_rank()
    atoms = load_configuration(args.start)
    if np.any(atoms.numbers.reshape(-1, 3) != [8, 1, 1]):
        raise SystemExit("pressure needs the molecules stored as O, H, H")
    volume = atoms.get_volume()
    waters = len(atoms) // 3
    density = waters * WATER_MASS / units.mol / volume * 1e27
    lines = [
        "# Generated by bench/constant_cutoff_npt.py pressure; do not edit by hand.",
        f"# {waters} waters from {args.start.name}, {volume:.2f} Angstrom^3, {density:.1f} kg/m^3; PBE, "
        f"O_MT_PBE and H_MT_PBE, orbital convergence {args.conv:g}; {mpi_size()} rank(s) on {cpu_model()}",
        "",
        "| Cutoff, Ry | E0, Ry | A, Ry | sigma, Ry | Energy, eV | Atomic pressure, kbar | "
        "Molecular pressure, kbar | Wall, s |",
        "|------------+--------+-------+-----------+------------+-----------------------+"
        "--------------------------+---------|",
    ]
    for basis in args.bases:
        fields = [float(x) for x in basis.split("/")]
        cutoff, e0 = fields[0], fields[1] if len(fields) > 1 else None
        width = fields[2] if len(fields) > 2 else None
        cc = None if e0 is None else {"cutoff_ry": e0, "width_ry": width}
        probe = atoms.copy()
        clock = time.perf_counter()
        with CPMDCCalculator(deck(cutoff, float(atoms.cell[0, 0]), args.conv), library=str(args.library),
                             constant_cutoff=cc) as calc:
            probe.calc = calc
            energy = probe.get_potential_energy()
            stress = probe.get_stress(voigt=False)
            forces = probe.get_forces()
        wall = time.perf_counter() - clock
        atomic = -np.trace(stress) / 3.0 / units.GPa * 10.0
        molecular = atomic - molecular_virial(probe, forces) / (3.0 * volume) / units.GPa * 10.0
        height = "" if e0 is None else f"{50.0 * cutoff:g}"
        sigma = "" if e0 is None else f"{width if width is not None else 0.1 * e0:g}"
        lines.append(f"| {cutoff:g} | {'' if e0 is None else f'{e0:g}'} | {height} | {sigma} | {energy:.4f} | "
                     f"{atomic:.2f} | {molecular:.2f} | {wall:.1f} |")
        if rank == 0:
            sys.stderr.write(lines[-1] + "\n")
    if rank == 0:
        text = "\n".join(lines) + "\n"
        if args.output:
            args.output.write_text(text)
        sys.stdout.write(text)
    finalize_mpi()
    return 0


def label(meta: dict) -> str:
    if meta["e0"] is None:
        return f"{meta['cutoff']:g} Ry, no constant cutoff"
    return f"{meta['cutoff']:g} Ry, E0 {meta['e0']:g} Ry"


def windows(run: dict, width_ps: float) -> list[tuple[float, float]]:
    t = np.asarray(run["time_ps"])
    rho = np.asarray(run["Density"])
    rows = []
    for lo in np.arange(0.0, t[-1], width_ps):
        sel = (t >= lo) & (t < lo + width_ps)
        if sel.any():
            rows.append((float(lo + 0.5 * width_ps), float(rho[sel].mean())))
    return rows


def density_window_rows(runs: list[tuple[str, dict]], width_ps: float) -> list[str]:
    columns = [dict(windows(data, width_ps)) for _, data in runs]
    centres = sorted(set().union(*columns))
    head = " | ".join(f"{name}, kg/m^3" for name, _ in runs)
    lines = [f"| Window, ps | {head} |", "|-" + "-+-" * len(runs) + "-|"]
    for c in centres:
        cells = [f"{column[c]:.1f}" if c in column else "" for column in columns]
        lines.append(f"| {c - 0.5 * width_ps:g}-{c + 0.5 * width_ps:g} | " + " | ".join(cells) + " |")
    return lines


def start_of(meta: dict) -> str:
    if meta.get("start_ps"):
        return f"{meta['start']} at {meta['start_ps']:g} ps"
    return meta["start"]


def cost_rows(runs: list[tuple[str, dict]]) -> list[str]:
    lines = ["| Run | Start | Ranks | Steps | Wall per step, s | SCF iterations per step | New setups | Wall, h |",
             "|-----+-------+-------+-------+------------------+-------------------------+------------+---------|"]
    for name, data in runs:
        wall = np.asarray(data["step_wall_s"][1:])
        scf = np.asarray(data["scf_steps"][1:])
        lines.append(f"| {name} | {start_of(data['args'])} | {data['args']['ranks']} | {len(wall)} | "
                     f"{wall.mean():.1f} | {scf.mean():.1f} | {int(np.sum(data['setups'][1:]))} | "
                     f"{np.sum(data['step_wall_s']) / 3600.0:.2f} |")
    return lines


def branch_rows(runs: list[tuple[str, dict]], step_ps: float) -> list[str]:
    """The density every step_ps after the start of each run that starts from the state of another."""
    branches = [(name, data) for name, data in runs if data["args"].get("start_ps")]
    if not branches:
        return []
    columns = [{round(float(t), 6): float(rho) for t, rho in zip(data["time_ps"], data["Density"])}
               for _, data in branches]
    end = max(float(data["time_ps"][-1]) for _, data in branches)
    marks = [round(k * step_ps, 6) for k in range(int(np.floor(end / step_ps + 1e-9)) + 1)]
    head = " | ".join(f"{name}, kg/m^3" for name, _ in branches)
    lines = [f"| After the start, ps | {head} |", "|-" + "-+-" * len(branches) + "-|"]
    for mark in marks:
        cells = [f"{column[mark]:.1f}" if mark in column else "" for column in columns]
        lines.append(f"| {mark:g} | " + " | ".join(cells) + " |")
    return lines


def table(args: argparse.Namespace) -> int:
    loaded = [json.loads(path.read_text()) for path in args.runs]
    meta = loaded[0]["args"]
    finished, unaveraged, plotted = [], [], []
    for data in loaded:
        name = label(data["args"])
        data["wall_s"] = float(np.sum(data["step_wall_s"]))
        start_ps = data["args"].get("start_ps", 0.0)
        floor = args.branch_discard_ps if start_ps else args.min_discard_ps
        length = float(data["time_ps"][-1])
        entry = {"label": name, "cutoff": data["args"]["cutoff"], "e0": data["args"]["e0"], "start_ps": start_ps,
                 "length_ps": length, "windows": windows(data, args.plot_window_ps), "stopped": data["stopped"]}
        if data["stopped"] is None and length - floor >= args.min_production_ps - 1e-9:
            result = analyse(data, args.discard_ps, floor, args.window_ps, args.min_blocks)
            finished.append((name, result))
            entry.update(mean=result["terms"]["Density"]["mean"], sem=result["terms"]["Density"]["sem"],
                         discard_ps=result["discard_ps"])
        else:
            why = data["stopped"] or (f"{length:g} ps, less than {floor:g} ps to discard and "
                                      f"{args.min_production_ps:g} ps to average")
            unaveraged.append((name, data, why))
        plotted.append(entry)
    lines = [
        "# Generated by bench/constant_cutoff_npt.py table; do not edit by hand.",
        f"# {meta['waters']} waters, PBE, O_MT_PBE and H_MT_PBE, orbital convergence {meta['conv']:g}, "
        f"ASPC order 4; LangevinBAOAB hydrostatic, {meta['temperature']:g} K, {meta['pressure']:g} bar, "
        f"T_tau {meta['t_tau']:g} fs, P_tau {meta['p_tau']:g} fs, {meta['dt']:g} fs steps, "
        f"hydrogen mass {meta['hydrogen_mass']:g} u; CPMDC_CELL_IN_PLACE={meta['cell_in_place'] or 'unset'}",
        f"# {meta['ranks']} rank(s) on {meta['cpu']}, ASE {meta['ase']}",
    ]
    if finished:
        lines += ["", *summary_rows(finished)]
    if unaveraged:
        lines += ["", "| Run | Start | Steps | No mean | Density at the start, kg/m^3 | Density at the end, kg/m^3 | "
                  "Mean pressure, bar |",
                  "|-----+-------+-------+---------+------------------------------+----------------------------+"
                  "--------------------|"]
        for name, data, why in unaveraged:
            lines.append(f"| {name} | {start_of(data['args'])} | {len(data['time_ps']) - 1} | {why} | "
                         f"{data['Density'][0]:.1f} | {data['Density'][-1]:.1f} | "
                         f"{np.mean(data['Pressure'][1:]):.0f} |")
    named = [(label(d["args"]), d) for d in loaded]
    lines += ["", *cost_rows(named)]
    if finished:
        lines += ["", *detail_rows(finished), "", "Standard error of the mean against block length:", "",
                  *curve_rows(finished, args.min_blocks)]
    lines += ["", "Mean density over consecutive windows:", "", *density_window_rows(named, args.window_ps)]
    rows = branch_rows(named, args.branch_step_ps)
    if rows:
        lines += ["", "Density after the start of each run that starts from the state of another:", "", *rows]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
    if args.json:
        args.json.write_text(json.dumps({"settings": lines[1:3], "runs": plotted}, indent=1) + "\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("start")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--waters", type=int, default=32)
    p.add_argument("--density", type=float, default=1.0)
    p.add_argument("--ps", type=float, default=10.0)
    p.add_argument("--temperature", type=float, default=350.0)
    p.add_argument("--seed", type=int, default=7)
    p = sub.add_parser("run")
    p.add_argument("--library", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--start", type=Path)
    p.add_argument("--cutoff", type=float)
    p.add_argument("--e0", type=float)
    p.add_argument("--height", type=float)
    p.add_argument("--width", type=float)
    p.add_argument("--conv", type=float, default=1e-6)
    p.add_argument("--dt", type=float, default=1.0, help="fs")
    p.add_argument("--temperature", type=float, default=350.0)
    p.add_argument("--pressure", type=float, default=1.0, help="bar")
    p.add_argument("--t-tau", type=float, default=100.0, help="fs")
    p.add_argument("--p-tau", type=float, default=500.0, help="fs")
    p.add_argument("--hydrogen-mass", type=float, default=2.014)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--density-range", type=float, nargs=2, default=(500.0, 1500.0), help="kg/m^3")
    p.add_argument("--steps", type=int, default=3000)
    p.add_argument("--checkpoint", type=int, default=10)
    p = sub.add_parser("pressure")
    p.add_argument("--library", required=True)
    p.add_argument("--start", type=Path, required=True)
    p.add_argument("--conv", type=float, default=1e-6)
    p.add_argument("--output", type=Path)
    p.add_argument("bases", nargs="+", help="CUTOFF, or CUTOFF/E0[/WIDTH] in Ry for a constant cutoff")
    p = sub.add_parser("table")
    p.add_argument("runs", type=Path, nargs="+")
    p.add_argument("--discard-ps", type=float)
    p.add_argument("--min-discard-ps", type=float, default=0.5, help="of a run that starts from a box")
    p.add_argument("--branch-discard-ps", type=float, default=0.05, help="of a run that starts from a state")
    p.add_argument("--min-production-ps", type=float, default=0.1)
    p.add_argument("--window-ps", type=float, default=0.1)
    p.add_argument("--plot-window-ps", type=float, default=0.02)
    p.add_argument("--branch-step-ps", type=float, default=0.01)
    p.add_argument("--min-blocks", type=int, default=8)
    p.add_argument("--output", type=Path)
    p.add_argument("--json", type=Path)
    args = parser.parse_args()
    return {"start": start, "run": run, "pressure": pressure, "table": table}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
