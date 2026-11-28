#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Density of liquid water from NPT dynamics through CPMDCCalculator.

start writes the box a run begins from: WATERS molecules at DENSITY g/cm^3,
oxygens on the sites of a simple cubic lattice taken in checkerboard order,
seeded random orientations, then PS ps of Langevin dynamics (0.5 fs steps,
friction 0.01/fs) at TEMPERATURE in a classical flexible model: the TIP3P
charges and Lennard-Jones of ASE's TIP3P, with harmonic O-H bonds and H-O-H
angle that take the SPC/Fw force constants about the TIP3P geometry.

run integrates hydrostatic NPT at PRESSURE bar and TEMPERATURE K through
CPMDCCalculator: PBE, the Troullier-Martins O_MT_PBE and H_MT_PBE
pseudopotentials of the CPMD library or the files that --oxygen and
--hydrogen name (cpmdc reads a name that ends in .UPF as UPF), basis cutoff
CUTOFF, orbital convergence CONV, and orbitals predicted by ASPC of order
--aspc-order.
--rs-wfn-keep keeps the real-space orbitals between the SCF iterations
(REAL SPACE WFN KEEP); --odiis-vectors and --hamiltonian-cutoff, when not 0,
set the ODIIS vectors and HAMILTONIAN CUTOFF. The integrator is ASE's
LangevinBAOAB by default: the MTK barostat with a Langevin piston, and
Langevin atoms with the friction 1/T_TAU. --thermostat svr gives the atoms
stochastic velocity rescaling (Bussi, Donadio and Parrinello, J. Chem. Phys.
126, 014101 (2007)) with the time constant T_TAU in place of that friction
and noise. --barostat pr is ASE's MelchionnaNPT with the cell shape held: a
Parrinello-Rahman piston whose mass follows from P_TAU and --bulk-modulus,
and a Nose-Hoover thermostat with the time constant T_TAU.
--e0 adds the constant-cutoff functional at E0 (A = 50 x CUTOFF and sigma =
E0 / 10 unless --height and --width are given); without it each new cell is
a new setup at CUTOFF, whose stress is that of a fixed number of plane
waves. A --height so small that the penalty vanishes, such as 1e-6 with
--e0 50 --width 1 on an 85 Ry deck, and a CPMDC_CELL_IN_PLACE above any
strain of the run keep the plane waves of the first cell through the run
instead: a basis of a constant number of plane waves, whose effective
cutoff falls as the cell grows. Hydrogen carries --hydrogen-mass, deuterium
by default, so that 1 fs steps resolve the O-H stretch; configurational
averages such as the density do not depend on the masses. The density is that of H2O, the water count
times 18.015 u over the volume, and the pressure is the instantaneous one
that the barostat holds, with the kinetic term. A run stops when the density
leaves --density-range, as it does when the cell collapses. The MTK piston
has the friction 1/P_TAU and the mass of ASE's heuristic, which grows as
P_TAU squared. --start takes the box of start, whose momenta the run draws at
TEMPERATURE, or the state file of another run, whose positions, momenta, cell
and strain rate it continues from with its own settings; the record then
holds the state's step count times DT as start_ps, the time into the run that
wrote the state.

Each step appends the time, volume, density, pressure, temperature,
potential energy, wall time, SCF iterations and new setups to the record
OUT. The last two come from the CPMDC_SCF_STEPS lines of the error stream,
which rank 0 sends to OUT with the suffix .err. Every --checkpoint steps
rank 0 rewrites OUT, with the time to 1e-6 ps, the volume to 1e-6 nm^3, the
density to 1e-3 kg/m^3, the pressure to 0.01 bar, the temperature to
0.01 K, the energy to 1e-6 eV and the wall time to 1e-4 s, and the state
file (suffix .state.npz), and --resume continues the run from them with the
settings recorded in OUT. Every
--keep-every steps, when it is not 0, rank 0 also keeps the state as
OUT.STEP.state.npz. CPMDC_CELL_IN_PLACE applies as set.

pressure evaluates one configuration, a box or the state file of a run, at
each basis it is given, CUTOFF for a fixed cutoff or CUTOFF/E0[/WIDTH] for a
constant cutoff with the default A, and writes an Org table of the energy and
of two virial pressures without the kinetic term: that of the atoms, and that
of the molecular centres of mass (with the masses of the state file), which
leaves out the intramolecular forces. --density first scales the molecular
centres and the cell to that density with each molecule held rigid.
--oxygen and --hydrogen name other pseudopotential files, as for run.

pressures gathers tables of pressure: the molecular pressure of each
configuration, by its file and density, at each pair of pseudopotential
files and basis that the tables give, then the mean and standard error of
each over the configurations and of its difference from the --against
column, taken configuration by configuration. Every configuration must have
every column.

table names each run by the file of its record and writes Org tables of the
basis, P_tau and start of each run, with the pseudopotentials when the runs
differ in them and the thermostat and barostat of records that name them; of
the block-averaged volume, density,
pressure and temperature (bench/npt_block_average.py: the equilibration
detected in the first three discarded, at least --min-discard-ps of a run
that starts from a box and --branch-discard-ps of one that starts from a
state, Flyvbjerg-Petersen standard error); of the runs that stopped or hold
less than --min-production-ps after that, with the density at their start and
end; of the cost per step; of the standard error against the block length,
every --curve-step blocking levels; of the first blocking level whose blocks
outlast the correlation by the criterion of Lee et al. (Phys. Rev. E 83,
066706 (2011)), or the coarsest level tried; of the isothermal compressibility from
the volume fluctuations; and of the density every --branch-step-ps after the
start of the runs that start from the same state. The other settings must be
the same for every run. The run that wrote a state is the record with a
frame at start_ps of the volume the state starts at. --json also writes what
docs/render_constant_cutoff_npt.py draws: the means, the window means, the
start of each run from the start of the first of the runs that wrote the
states before it, and the density of every step of the runs that start from
the same state.

cost compares runs with the same start, physical settings and ranks that
differ in their library or SCF options: the options in which the runs
differ, the mean wall time and SCF iterations per step from
--from-step to --to-step, and the largest difference of the potential energy
over the steps both runs hold from the run listed before (--against above)
or from the first run (--against first). The header gives the options the
runs share.

capnp must be on PATH and CPMDC_PSEUDO_DIR must name the pseudopotential
directory. Under mpirun with CPMDC_RANK_CONSENSUS=1 every rank integrates the
same trajectory and rank 0 writes.

usage: constant_cutoff_npt.py start --out BOX.xyz [--waters 32] [--density 1.0]
                                    [--ps 10] [--temperature 350] [--seed 7]
       constant_cutoff_npt.py run --library LIB --out RUN.json
                                  (--start BOX.xyz|RUN.state.npz --cutoff RY [--e0 RY]
                                   [--height RY] [--width RY] [--conv 1e-6]
                                   [--dt 1] [--temperature 350] [--pressure 1]
                                   [--t-tau 100] [--p-tau 100] [--seed 11]
                                   [--rs-wfn-keep] [--aspc-order 4] [--odiis-vectors 0]
                                   [--hamiltonian-cutoff 0] [--thermostat langevin|svr]
                                   [--barostat mtk|pr] [--bulk-modulus 2.2]
                                   [--oxygen O_MT_PBE.psp] [--hydrogen H_MT_PBE.psp]
                                   | --resume) [--steps 3000] [--checkpoint 10]
                                  [--keep-every 0]
       constant_cutoff_npt.py pressure --library LIB --start BOX.xyz|RUN.state.npz
                                       [--conv 1e-6]
                                       [--density KG_M3 | --scale-volume FACTOR]
                                       [--oxygen O_MT_PBE.psp] [--hydrogen H_MT_PBE.psp]
                                       [--output FILE] BASIS [BASIS ...]
       constant_cutoff_npt.py pressures TABLE.org [TABLE.org ...] [--against 1]
                                        [--output FILE]
       constant_cutoff_npt.py pulay TABLE.org [TABLE.org ...] [--reference BASIS]
                                    [--output FILE]
       constant_cutoff_npt.py pv TABLE.org [TABLE.org ...] [--against BASIS]
                                 [--temperature 350] [--pressure 1] [--degree 2]
                                 [--no-kinetic] [--output FILE]
       constant_cutoff_npt.py table RUN.json [RUN.json ...] [--json OUT]
                                    [--discard-ps PS] [--min-discard-ps 0.5]
                                    [--branch-discard-ps 0.05] [--min-production-ps 0.1]
                                    [--window-ps 0.1] [--plot-window-ps 0.02]
                                    [--branch-step-ps 0.01] [--min-blocks 8]
                                    [--curve-step 2] [--prefix-check] [--output FILE]
       constant_cutoff_npt.py cost RUN.json [RUN.json ...] [--from-step 1]
                                   [--to-step N] [--against above|first]
                                   [--output FILE]
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
from typing import TYPE_CHECKING

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
from ase.md.velocitydistribution import Stationary  # noqa: E402

if TYPE_CHECKING:
    from ase.md.langevinbaoab import LangevinBAOAB

try:
    from ase.md.melchionna import MelchionnaNPT
except ImportError:
    class MelchionnaNPT:  # noqa: D101
        """Stand-in so an MTK run can ask isinstance without this class."""

from constant_cutoff_scan import cpu_model  # noqa: E402
from cpmdc import find_schema  # noqa: E402
from cpmdc.ase import CPMDCCalculator  # noqa: E402
from npt_block_average import TERMS, analyse, block_error, curve_rows, detail_rows, statistical_inefficiency  # noqa: E402


def thermalize_momenta(atoms: Atoms, temperature_K: float, rng) -> None:
    """Maxwell-Boltzmann momenta.

    ASE 3.29 names this ``thermalize_momenta``. Earlier releases expose the
    same draw as ``MaxwellBoltzmannDistribution``.
    """

    try:
        from ase.md.velocitydistribution import thermalize_momenta as draw
    except ImportError:
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as draw

    draw(atoms, temperature_K=temperature_K, rng=rng)


def langevin_baoab():
    """ASE's hydrostatic Langevin NPT integrator, added in ASE 3.27."""

    try:
        from ase.md.langevinbaoab import LangevinBAOAB
    except ImportError as exc:
        raise ImportError(
            "constant_cutoff_npt.py run needs ASE 3.27 or newer "
            "(ase.md.langevinbaoab.LangevinBAOAB)"
        ) from exc
    return LangevinBAOAB


WATER_MASS = float(atomic_masses[8] + 2.0 * atomic_masses[1])
BOND_K = 1059.162 * units.kcal / units.mol
ANGLE_K = 75.90 * units.kcal / units.mol
THETA0 = np.radians(angleHOH)
KB_J = 1.380649e-23
OXYGEN, HYDROGEN = "O_MT_PBE.psp", "H_MT_PBE.psp"
FIELDS = ("time_ps", "Volume", "Density", "Pressure", "Temperature", "Energy", "step_wall_s", "scf_steps",
          "setups")
# Decimals each record keeps: ps, nm^3, kg/m^3, bar, K, eV, s.
DIGITS = {"time_ps": 6, "Volume": 6, "Density": 3, "Pressure": 2, "Temperature": 2, "Energy": 6, "step_wall_s": 4}
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


def deck(cutoff: float, side: float, conv: float, keep: bool = False, order: int = 4, odiis: int = 0,
         hthrs: float = 0.0, oxygen: str = OXYGEN, hydrogen: str = HYDROGEN) -> bytes:
    text = f"""(
  functional = "PBE", cutOffRy = {cutoff:.2f}, charge = 0, multiplicity = 1, task = "gradient",
  inputSections = [
    ( cpmd = ( optimizeWavefunction = true, convergenceOrbitals = {conv:g}, maxStep = 200, maxIter = 200,
        extrapolateWfnOptions = "ASPC", extrapolateWfnOrder = {order}, realSpaceWfnKeep = {str(keep).lower()},
        odiisVectors = {odiis}, hamiltonianCutoff = {hthrs:g},
        directives = [ ( keyword = "CENTER MOLECULE OFF", args = [] ) ] ) ),
    ( dft = ( functional = "PBE", lsd = false, directives = [ ( keyword = "OLDCODE", args = [] ) ] ) ),
    ( system = ( symmetry = 1, angstrom = true, cutOffRy = {cutoff:.2f},
        cell = [{side:.6f}, 1.0, 1.0, 0.0, 0.0, 0.0] ) ),
    ( atoms = ( pseudopotentials = [
        ( element = "O", path = "{oxygen}", lmax = 1 ),
        ( element = "H", path = "{hydrogen}", lmax = 0 ) ] ) ) ] )
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


def call_and_finalize(call):
    """Run call and finalize MPI even when call raises."""
    try:
        return call()
    finally:
        finalize_mpi()


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


def _bind_mtk_classes() -> None:
    """Publish LangevinBAOAB and RescalingBAOAB once ASE provides them.

    Importing this module stays possible on ASE releases that predate
    LangevinBAOAB. The classes are bound on the first MTK step.
    """
    global LangevinBAOAB, RescalingBAOAB
    base = langevin_baoab()
    if globals().get("_mtk_base") is base:
        return

    class RescalingBAOAB(base):
        """LangevinBAOAB whose atoms take stochastic velocity rescaling (Bussi, Donadio and Parrinello, J. Chem. Phys.
        126, 014101 (2007)) with time constant T_tau over each O step instead of the Langevin friction and noise; the
        barostat drag on the velocities and the Langevin piston stay."""

        def set_temperature(self, temperature_K, from_init=False):
            super().set_temperature(temperature_K, from_init)
            self.gamma = 0.0

        def _BAOAB_OU(self, drag_gamma_mod=0.0):
            super()._BAOAB_OU(drag_gamma_mod)
            if self.temperature_K is None:
                return
            n = self.atoms.get_number_of_degrees_of_freedom()
            kinetic = self.atoms.get_kinetic_energy()
            if kinetic <= 0.0:
                return
            target = 0.5 * n * units.kB * self.temperature_K
            c = np.exp(-self.dt / self.T_tau)
            r = self.rng.normal()
            rest = self.rng.chisquare(n - 1)
            new = (c * kinetic + (1.0 - c) * target * (r * r + rest) / n
                   + 2.0 * r * np.sqrt(c * (1.0 - c) * kinetic * target / n))
            self.atoms.set_velocities(self.atoms.get_velocities() * np.sqrt(new / kinetic))

    RescalingBAOAB.__name__ = "RescalingBAOAB"
    RescalingBAOAB.__qualname__ = "RescalingBAOAB"
    LangevinBAOAB = base
    globals()["RescalingBAOAB"] = RescalingBAOAB
    globals()["_mtk_base"] = base


def _melchionna_npt():
    """ASE's Parrinello-Rahman integrator, or a clear failure when it is absent."""
    global MelchionnaNPT
    if getattr(MelchionnaNPT, "__module__", "").startswith("ase."):
        return MelchionnaNPT
    try:
        from ase.md.melchionna import MelchionnaNPT as kind
    except ImportError as exc:
        raise ImportError(
            "constant_cutoff_npt.py --barostat pr needs ase.md.melchionna.MelchionnaNPT"
        ) from exc
    MelchionnaNPT = kind
    return kind


def write_atomic(path: Path, write_to) -> None:
    tmp = path.with_name(path.name + ".tmp")
    write_to(tmp)
    os.replace(tmp, path)


def integrator(setup: argparse.Namespace, meta: dict, atoms: Atoms):
    """LangevinBAOAB, or its stochastic velocity rescaling, for the MTK barostat; ASE's MelchionnaNPT with the cell
    shape held for the Parrinello-Rahman barostat, whose Nose-Hoover thermostat has the time constant T_tau and
    whose piston mass follows from P_tau and --bulk-modulus."""
    dt = setup.dt * units.fs
    if meta.get("barostat", "mtk") == "pr":
        kind = _melchionna_npt()
        dyn = kind(atoms, dt, temperature_K=setup.temperature,
                   externalstress=setup.pressure * units.bar, ttime=setup.t_tau * units.fs,
                   pfactor=(setup.p_tau * units.fs) ** 2 * meta["bulk_modulus"] * units.GPa)
        dyn.set_fraction_traceless(0.0)
        return dyn
    _bind_mtk_classes()
    kind = RescalingBAOAB if meta.get("thermostat", "langevin") == "svr" else LangevinBAOAB
    return kind(atoms, dt, temperature_K=setup.temperature, externalstress=-setup.pressure * units.bar,
                hydrostatic=True, T_tau=setup.t_tau * units.fs, P_tau=setup.p_tau * units.fs,
                rng=np.random.default_rng(setup.seed + 1))


def strain_rate(saved) -> float:
    """d ln(L)/dt of a state file, per ASE time unit."""
    if "strain_rate" in saved:
        return float(saved["strain_rate"])
    return float(saved["p_eps"]) / float(saved["barostat_mass"])


def restore_dynamics(dyn, atoms: Atoms, saved) -> None:
    """Put a restarted integrator back on the saved state.

    MelchionnaNPT's constructor removes the centre-of-mass momentum, and its
    first step rebuilds the previous Verlet state. A resume restores the
    saved momenta and, once a step has been taken, that previous state.
    """
    atoms.set_momenta(np.asarray(saved["momenta"]))
    dyn.nsteps = int(saved["nsteps"])
    if isinstance(dyn, MelchionnaNPT):
        if "q" in saved.files:
            dyn.eta = np.asarray(saved["eta"])
            dyn.eta_past = np.asarray(saved["eta_past"])
            dyn.zeta = float(saved["zeta"])
            dyn.zeta_past = float(saved["zeta_past"])
            dyn.zeta_integrated = float(saved["zeta_integrated"])
            dyn.h = np.asarray(saved["h"])
            dyn.h_past = np.asarray(saved["h_past"])
            dyn.inv_h = np.linalg.inv(dyn.h)
            dyn.timeelapsed = float(saved["timeelapsed"])
            dyn.q = np.asarray(saved["q"])
            dyn.q_past = np.asarray(saved["q_past"])
            dyn.q_future = np.asarray(saved["q_future"])
            dyn.initialized = 1
        else:
            set_strain_rate(dyn, strain_rate(saved))
            dyn.zeta = float(saved["zeta"])
            dyn.zeta_integrated = float(saved["zeta_integrated"])
    else:
        dyn.p_eps = float(saved["p_eps"])
        dyn.rng.bit_generator.state = json.loads(str(saved["rng"]))


def set_strain_rate(dyn, rate: float) -> None:
    if isinstance(dyn, MelchionnaNPT):
        # ASE 3.29's set_strain_rate calls a method that no longer exists; before initialize() it only sets eta.
        assert not dyn.initialized
        dyn.eta = rate * np.identity(3)
    else:
        dyn.p_eps = rate * dyn.barostat_mass


def save(out: Path, record: dict, atoms: Atoms, dyn, keep: bool = False) -> None:
    rounded = {**record, **{name: [round(value, digits) for value in record[name]]
                            for name, digits in DIGITS.items() if name in record}}
    write_atomic(out, lambda p: p.write_text(json.dumps(rounded) + "\n"))
    if isinstance(dyn, MelchionnaNPT):
        extra = {"strain_rate": np.asarray(np.trace(dyn.get_strain_rate()) / 3.0), "zeta": np.asarray(dyn.zeta),
                 "zeta_integrated": np.asarray(dyn.zeta_integrated)}
        if dyn.initialized:
            extra.update(eta=np.asarray(dyn.eta), eta_past=np.asarray(dyn.eta_past),
                         zeta_past=np.asarray(dyn.zeta_past), h=np.asarray(dyn.h),
                         h_past=np.asarray(dyn.h_past), timeelapsed=np.asarray(dyn.timeelapsed),
                         q=np.asarray(dyn.q), q_past=np.asarray(dyn.q_past),
                         q_future=np.asarray(dyn.q_future))
    else:
        extra = {"p_eps": np.asarray(dyn.p_eps), "barostat_mass": np.asarray(dyn.barostat_mass),
                 "strain_rate": np.asarray(dyn.p_eps / dyn.barostat_mass),
                 "rng": json.dumps(dyn.rng.bit_generator.state)}

    def state(p: Path) -> None:
        with p.open("wb") as handle:
            np.savez(handle, numbers=atoms.numbers, positions=atoms.positions, momenta=atoms.get_momenta(),
                     masses=atoms.get_masses(), cell=atoms.cell.array, nsteps=dyn.nsteps, **extra)

    write_atomic(out.with_suffix(".state.npz"), state)
    if keep:
        write_atomic(out.with_suffix(f".{dyn.nsteps}.state.npz"), state)


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
                     "build": Path(args.library).resolve().parent.name,
                     "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
                     "cell_in_place": os.environ.get("CPMDC_CELL_IN_PLACE", ""), "ase": ase.__version__,
                     "rs_wfn_keep": args.rs_wfn_keep,
                     "aspc_order": args.aspc_order, "odiis_vectors": args.odiis_vectors,
                     "hamiltonian_cutoff": args.hamiltonian_cutoff,
                     "thermostat": "nose-hoover" if args.barostat == "pr" else args.thermostat,
                     "barostat": args.barostat, "bulk_modulus": args.bulk_modulus, "oxygen": args.oxygen,
                     "hydrogen": args.hydrogen},
            **{name: [] for name in FIELDS},
            "stopped": None,
        }
    waters = record["args"]["waters"]
    cc = None if setup.e0 is None else {"cutoff_ry": setup.e0, "height_ry": setup.height,
                                        "width_ry": setup.width}
    numerics = record["args"]
    calc = CPMDCCalculator(deck(setup.cutoff, float(atoms.cell[0, 0]), setup.conv,
                                bool(numerics.get("rs_wfn_keep", False)), int(numerics.get("aspc_order", 4)),
                                int(numerics.get("odiis_vectors", 0)), float(numerics.get("hamiltonian_cutoff", 0.0)),
                                numerics.get("oxygen", OXYGEN), numerics.get("hydrogen", HYDROGEN)),
                           library=str(args.library), constant_cutoff=cc)
    atoms.calc = calc

    def dynamics() -> int:
        clock = [time.perf_counter()]
        dyn = integrator(setup, numerics, atoms)
        if args.resume:
            restore_dynamics(dyn, atoms, saved)
        elif args.start.suffix == ".npz":
            set_strain_rate(dyn, strain_rate(saved))
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
            kept = args.keep_every > 0 and dyn.nsteps % args.keep_every == 0
            if rank == 0 and (dyn.nsteps % args.checkpoint == 0 or kept or not lo <= density <= hi):
                save(out, record, atoms, dyn, kept)
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
        return 0

    return call_and_finalize(dynamics)


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


def scale_molecules(atoms: Atoms, factor: float) -> None:
    """Scale the cell and the centres of mass of the O, H, H molecules by factor, keeping each molecule rigid."""
    cell = atoms.cell.array
    r = atoms.positions.reshape(-1, 3, 3)
    d = r - r[:, :1]
    d = d - np.round(d @ np.linalg.inv(cell)) @ cell
    m = atoms.get_masses().reshape(-1, 3, 1)
    offset = np.sum(m * d, axis=1, keepdims=True) / np.sum(m, axis=1, keepdims=True)
    atoms.set_cell(cell * factor, scale_atoms=False)
    atoms.positions = ((r[:, :1] + offset) * factor + (d - offset)).reshape(-1, 3)


def pressure(args: argparse.Namespace) -> int:
    rank = mpi_rank()
    atoms = load_configuration(args.start)
    if np.any(atoms.numbers.reshape(-1, 3) != [8, 1, 1]):
        raise SystemExit("pressure needs the molecules stored as O, H, H")
    waters = len(atoms) // 3
    stored = waters * WATER_MASS / units.mol / atoms.get_volume() * 1e27
    if args.density is not None and args.scale_volume is not None:
        raise SystemExit("pressure accepts --density or --scale-volume, not both")
    if args.density is not None:
        scale_molecules(atoms, (stored / args.density) ** (1.0 / 3.0))
    elif args.scale_volume is not None:
        scale_molecules(atoms, args.scale_volume ** (1.0 / 3.0))
    volume = atoms.get_volume()
    density = waters * WATER_MASS / units.mol / volume * 1e27
    scaled = ("" if args.density is None and args.scale_volume is None
              else f" with the molecules scaled rigidly from {stored:.1f} kg/m^3")
    lines = [
        "# Generated by bench/constant_cutoff_npt.py pressure; do not edit by hand.",
        f"# {waters} waters from {args.start.name}{scaled}, {volume:.2f} Angstrom^3, {density:.1f} kg/m^3; PBE, "
        f"{args.oxygen} and {args.hydrogen}, orbital convergence {args.conv:g}; {mpi_size()} rank(s) on "
        f"{cpu_model()}",
        "",
        "| Cutoff, Ry | E0, Ry | A, Ry | sigma, Ry | Energy, eV | Atomic pressure, kbar | "
        "Molecular pressure, kbar | Wall, s |",
        "|------------+--------+-------+-----------+------------+-----------------------+"
        "--------------------------+---------|",
    ]
    def pressures() -> int:
        for basis in args.bases:
            fields = [float(x) for x in basis.split("/")]
            cutoff, e0 = fields[0], fields[1] if len(fields) > 1 else None
            width = fields[2] if len(fields) > 2 else None
            cc = None if e0 is None else {"cutoff_ry": e0, "width_ry": width}
            probe = atoms.copy()
            clock = time.perf_counter()
            with CPMDCCalculator(deck(cutoff, float(atoms.cell[0, 0]), args.conv, oxygen=args.oxygen,
                                      hydrogen=args.hydrogen), library=str(args.library), constant_cutoff=cc) as calc:
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
        return 0

    return call_and_finalize(pressures)


PRESSURE_HEAD = re.compile(r"^# (\d+) waters from (\S+?)(?: with the molecules scaled rigidly from [\d.]+ kg/m\^3)?, "
                           r"([\d.]+) Angstrom\^3, ([\d.]+) kg/m\^3; PBE, (\S+) and (\S+), orbital convergence", re.M)


def parse_pressure_table(path: Path) -> tuple[tuple[int, str, float, str, str, str],
                                               list[tuple[str, float, float, float]]]:
    """Header metadata and (basis, energy_eV, atomic_kbar, molecular_kbar) rows of a pressure table."""
    text = path.read_text()
    head = PRESSURE_HEAD.search(text)
    if head is None:
        raise SystemExit(f"{path} is not a table of constant_cutoff_npt.py pressure")
    waters, start, volume, density, oxygen, hydrogen = head.groups()
    rows = []
    for line in text.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if line.startswith("|") and cells[0][:1].isdigit():
            if len(cells) < 7:
                raise SystemExit(f"{path} has a pressure row shorter than the energy and the two pressures")
            bases = "/".join(cell for cell in (cells[0], cells[1], cells[3]) if cell)
            rows.append((bases, float(cells[4]), float(cells[5]), float(cells[6])))
    return (int(waters), start, float(volume), density, oxygen, hydrogen), rows


def read_pressure_table(path: Path) -> tuple[tuple[str, str], list[tuple[str, float]]]:
    """The configuration and density of a pressure table, and the files, basis and molecular pressure of each row."""
    (_, start, _, density, oxygen, hydrogen), full_rows = parse_pressure_table(path)
    return (start, density), [(f"{oxygen}, {hydrogen}, {bases}", mol) for bases, _, _, mol in full_rows]


def pressures(args: argparse.Namespace) -> int:
    held: dict[tuple[str, str], dict[str, float]] = {}
    columns: list[str] = []
    for path in args.tables:
        configuration, rows = read_pressure_table(path)
        values = held.setdefault(configuration, {})
        for column, value in rows:
            if column in values:
                raise SystemExit(f"two tables give {configuration[0]} at {configuration[1]} kg/m^3 with {column}")
            values[column] = value
            if column not in columns:
                columns.append(column)
    for (start, density), values in held.items():
        for column in columns:
            if column not in values:
                raise SystemExit(f"no table gives {start} at {density} kg/m^3 with {column}")
    if not 1 <= args.against <= len(columns):
        raise SystemExit(f"--against takes a column from 1 to {len(columns)}")
    reference = columns[args.against - 1]
    grid = np.array([[values[column] for column in columns] for values in held.values()])
    differences = grid - grid[:, [args.against - 1]]
    count = len(grid)

    def error(samples: np.ndarray) -> str:
        return f"{np.std(samples, ddof=1) / np.sqrt(count):.3f}" if count > 1 else ""

    head = "| Configuration | Density, kg/m^3 | " + " | ".join(f"{column}, kbar" for column in columns) + " |"
    lines = ["# Generated by bench/constant_cutoff_npt.py pressures; do not edit by hand.",
             f"# Molecular pressures without the kinetic term at {count} configuration(s); the differences are "
             f"from {reference}, configuration by configuration",
             "", head, "|" + "+".join("-" * len(cell) for cell in head.split("|")[1:-1]) + "|"]
    for (start, density), row in zip(held, grid):
        lines.append(f"| {start} | {density} | " + " | ".join(f"{value:.2f}" for value in row) + " |")
    head = ("| Files and basis | Mean, kbar | Standard error, kbar | Difference, kbar | "
            "Standard error of the difference, kbar |")
    lines += ["", head, "|" + "+".join("-" * len(cell) for cell in head.split("|")[1:-1]) + "|"]
    for index, column in enumerate(columns):
        difference = ("|  |  |" if column == reference else
                      f"| {differences[:, index].mean():.3f} | {error(differences[:, index])} |")
        lines.append(f"| {column} | {grid[:, index].mean():.3f} | {error(grid[:, index])} {difference}")
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
    return 0


def pulay(args: argparse.Namespace) -> int:
    """Pulay correction P_ref - P_fixed(E_cut) and its step-to-step change across saved frames."""
    frames: dict[tuple[str, str], dict[str, tuple[float, float]]] = {}
    fixed_order: list[str] = []
    cc_bases: list[str] = []
    for path in args.tables:
        (_, start, _, density, _, _), full_rows = parse_pressure_table(path)
        entry = frames.setdefault((start, density), {})
        for basis_name, _, atomic, mol in full_rows:
            entry[basis_name] = (atomic, mol)
            if "/" in basis_name:
                if basis_name not in cc_bases:
                    cc_bases.append(basis_name)
            elif basis_name not in fixed_order:
                fixed_order.append(basis_name)
    fixed_order.sort(key=float)
    if not fixed_order:
        raise SystemExit("pulay needs at least one fixed-cutoff basis row")
    reference = args.reference or (cc_bases[-1] if cc_bases else fixed_order[-1])
    for (start, density), entry in frames.items():
        if reference not in entry:
            raise SystemExit(f"{start} at {density} kg/m^3 has no reference row {reference}")
        for cut in fixed_order:
            if cut not in entry:
                raise SystemExit(f"{start} at {density} kg/m^3 has no fixed-cutoff row {cut}")
    n_frames = len(frames)
    keys = list(frames)

    def sem(values: np.ndarray) -> str:
        return f"{np.std(values, ddof=1) / np.sqrt(n_frames):.3f}" if n_frames > 1 else ""

    lines = [
        "# Generated by bench/constant_cutoff_npt.py pulay; do not edit by hand.",
        f"# Pulay correction (reference {reference} minus fixed-cutoff pressure) over {n_frames} configuration(s)",
        "",
        "| Cutoff, Ry | Mol P_fixed mean, kbar | Mol Pulay, kbar | Mol Pulay SE, kbar | Mol dPulay, kbar | "
        "Max |Mol dPulay|, kbar | Atom Pulay, kbar | Atom Pulay SE, kbar | Max |Atom dPulay|, kbar |",
        "|------------+------------------------+-----------------+--------------------+------------------+"
        "------------------------+------------------+---------------------+-------------------------|",
    ]
    prev_mol: np.ndarray | None = None
    prev_atom: np.ndarray | None = None
    for cut in fixed_order:
        mol_fixed = np.array([frames[k][cut][1] for k in keys])
        mol_pulay = np.array([frames[k][reference][1] - frames[k][cut][1] for k in keys])
        atom_pulay = np.array([frames[k][reference][0] - frames[k][cut][0] for k in keys])
        if prev_mol is None or prev_atom is None:
            d_mol_str, max_mol_str, max_atom_str = "", "", ""
        else:
            dm = mol_pulay - prev_mol
            da = atom_pulay - prev_atom
            d_mol_str = f"{dm.mean():+.3f}"
            max_mol_str = f"{np.max(np.abs(dm)):.2f}"
            max_atom_str = f"{np.max(np.abs(da)):.2f}"
        lines.append(
            f"| {cut} | {mol_fixed.mean():.2f} | {mol_pulay.mean():+.3f} | {sem(mol_pulay)} | {d_mol_str} | "
            f"{max_mol_str} | {atom_pulay.mean():+.3f} | {sem(atom_pulay)} | {max_atom_str} |"
        )
        prev_mol, prev_atom = mol_pulay, atom_pulay
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
    return 0


def _solve_pv_root(
    coefs: np.ndarray, target_kbar: float, guess: float, lower: float, upper: float, origin: float = 0.0,
) -> float:
    """Real root of p(x) = target inside [lower, upper], nearest to guess.

    x is the density minus origin. lower and upper bound x. The real part of a
    complex pair is not a root. None inside the sampled densities means the fit
    has not located an equilibrium in the volumes it was given.
    """
    poly = coefs.copy()
    poly[-1] -= target_kbar
    roots = np.roots(poly)
    real = [float(np.real(z)) for z in roots if np.isfinite(z) and abs(float(np.imag(z))) < 1e-8]
    span = max(upper - lower, 1.0)
    tol = 1e-8 * span
    inside = [r for r in real if lower - tol <= r <= upper + tol]
    if not inside:
        raise SystemExit(
            f"no p(V) root lies in the sampled density range {lower + origin:.6g}..{upper + origin:.6g} kg/m^3"
        )
    return min(inside, key=lambda r: abs(r - guess))


def _pooled_pv_fit(
    frame_rows: list[list[tuple[float, float]]], degree: int, target_kbar: float
) -> tuple[float, float]:
    """Fixed-effects polynomial fit P_{k,j} = alpha_k + sum_{m=1..deg} c_m (rho - rho_ref)^m."""
    rho_ref = 850.0
    n_f = len(frame_rows)
    design, y = [], []
    for i, pts in enumerate(frame_rows):
        for rho, press in pts:
            dr = rho - rho_ref
            row = [0.0] * n_f + [dr**m for m in range(1, degree + 1)]
            row[i] = 1.0
            design.append(row)
            y.append(press)
    beta, *_ = np.linalg.lstsq(np.asarray(design), np.asarray(y), rcond=None)
    alpha_bar = float(np.mean(beta[:n_f]))
    poly = np.array([*beta[n_f:][::-1], alpha_bar])
    sampled = [rho for pts in frame_rows for rho, _ in pts]
    lo, hi = min(sampled), max(sampled)
    dr_eq = _solve_pv_root(
        poly, target_kbar, float(np.mean(sampled)) - rho_ref, lo - rho_ref, hi - rho_ref, origin=rho_ref,
    )
    rho_eq = rho_ref + dr_eq
    dp_drho = float(np.polyval(np.polyder(poly), dr_eq))
    kappa_bar = 1.0 / (rho_eq * dp_drho) * 1e-3
    return rho_eq, kappa_bar


def pv(args: argparse.Namespace) -> int:
    """Fit p(V) across scaled volumes per setup and report implied equilibrium density shifts."""
    by_setup: dict[str, dict[str, list[tuple[float, float]]]] = {}
    for path in args.tables:
        (waters, start, vol_a3, density, oxygen, hydrogen), full_rows = parse_pressure_table(path)
        rho_val = float(density)
        p_kin = 0.0 if args.no_kinetic else waters * KB_J * args.temperature / (vol_a3 * 1e-30) / 1e8
        for basis_name, _, _, mol in full_rows:
            setup = f"{oxygen}, {hydrogen}, {basis_name}"
            by_setup.setdefault(setup, {}).setdefault(start, []).append((rho_val, mol + p_kin))
    setups = list(by_setup)
    if not setups:
        raise SystemExit("pv needs at least one pressure table")
    frames = sorted({start for entries in by_setup.values() for start in entries})
    for setup in setups:
        for start in frames:
            pts = by_setup[setup].get(start, [])
            if len(pts) <= args.degree:
                raise SystemExit(f"{setup} on {start} has {len(pts)} volume(s), needs > {args.degree}")
    reference = args.against or setups[0]
    if reference not in by_setup:
        raise SystemExit(f"unknown --against setup {reference}")
    target_kbar = args.pressure * 1e-3
    n_f = len(frames)

    per_frame_rho: dict[str, np.ndarray] = {}
    per_frame_kap: dict[str, np.ndarray] = {}
    pooled_rho: dict[str, float] = {}
    pooled_kap: dict[str, float] = {}
    jk_rho: dict[str, np.ndarray] = {}

    for setup in setups:
        rhos, kaps = [], []
        ordered_pts = [sorted(by_setup[setup][start]) for start in frames]
        for pts in ordered_pts:
            r_arr = np.array([p[0] for p in pts])
            p_arr = np.array([p[1] for p in pts])
            coef = np.polyfit(r_arr, p_arr, args.degree)
            req = _solve_pv_root(
                coef, target_kbar, float(np.mean(r_arr)), float(r_arr.min()), float(r_arr.max()),
            )
            dp_drho = float(np.polyval(np.polyder(coef), req))
            rhos.append(req)
            kaps.append(1.0 / (req * dp_drho) * 1e-3)
        per_frame_rho[setup] = np.asarray(rhos)
        per_frame_kap[setup] = np.asarray(kaps)
        pooled_rho[setup], pooled_kap[setup] = _pooled_pv_fit(ordered_pts, args.degree, target_kbar)
        if n_f > 1:
            jk_rho[setup] = np.array([
                _pooled_pv_fit([ordered_pts[j] for j in range(n_f) if j != i], args.degree, target_kbar)[0]
                for i in range(n_f)
            ])
        else:
            jk_rho[setup] = np.array([pooled_rho[setup]])

    def paired_se(arr: np.ndarray) -> str:
        return f"{np.std(arr, ddof=1) / np.sqrt(n_f):.2f}" if n_f > 1 else ""

    def jk_se(arr: np.ndarray) -> str:
        if n_f <= 1:
            return ""
        return f"{np.sqrt((n_f - 1.0) / n_f * np.sum((arr - arr.mean()) ** 2)):.2f}"

    kin_note = "virial only" if args.no_kinetic else f"with N_mol k_B T / V at {args.temperature:g} K"
    lines = [
        "# Generated by bench/constant_cutoff_npt.py pv; do not edit by hand.",
        f"# Polynomial degree {args.degree} p(V) fit at {args.pressure:g} bar ({kin_note}) across {n_f} frame(s); "
        f"shifts from {reference}",
        "",
        "| Setup | Per-frame rho_eq, kg/m^3 | Paired SE, kg/m^3 | Paired shift, kg/m^3 | Shift SE, kg/m^3 | "
        "Pooled rho_eq, kg/m^3 | Jackknife SE, kg/m^3 | Pooled shift, kg/m^3 | Shift JK SE, kg/m^3 | "
        "kappa_rigid, 1e-5/bar |",
        "|-------+--------------------------+-------------------+----------------------+------------------+-"
        "---------------------+---------------------+----------------------+---------------------+-----------------------|",
    ]
    for setup in setups:
        pf = per_frame_rho[setup]
        pf_shift = pf - per_frame_rho[reference]
        pl_shift = pooled_rho[setup] - pooled_rho[reference]
        jk_shift = jk_rho[setup] - jk_rho[reference]
        shift_cols = ("|  |  |" if setup == reference else f"| {pf_shift.mean():+.2f} | {paired_se(pf_shift)} |")
        jk_shift_cols = ("|  |  |" if setup == reference else f"| {pl_shift:+.2f} | {jk_se(jk_shift)} |")
        lines.append(
            f"| {setup} | {pf.mean():.2f} | {paired_se(pf)} {shift_cols} {pooled_rho[setup]:.2f} | "
            f"{jk_se(jk_rho[setup])} {jk_shift_cols} {pooled_kap[setup] * 1e5:.2f} |"
        )
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
    return 0


def basis(meta: dict) -> str:
    if meta["e0"] is None:
        return f"{meta['cutoff']:g} Ry, no constant cutoff"
    height = "" if meta.get("height") is None else f", A {meta['height']:g} Ry"
    return f"{meta['cutoff']:g} Ry, E0 {meta['e0']:g} Ry{height}"


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
    lines = [f"| After the start, ps | {head} |", "|-" + "-+-" * len(runs) + "-|"]
    for c in centres:
        cells = [f"{column[c]:.1f}" if c in column else "" for column in columns]
        lines.append(f"| {c - 0.5 * width_ps:g}-{c + 0.5 * width_ps:g} | " + " | ".join(cells) + " |")
    return lines


def writers(runs: list[tuple[str, dict]]) -> dict[str, str]:
    """The run that wrote the state each run starts from: the one with a frame at start_ps that has the volume
    of the run's first frame."""
    found = {}
    for name, data in runs:
        start_ps = data["args"].get("start_ps")
        if not start_ps:
            continue
        for other, parent in runs:
            if other == name:
                continue
            at = np.flatnonzero(np.abs(np.asarray(parent["time_ps"]) - start_ps) < 1e-9)
            if at.size and parent["Volume"][at[0]] == data["Volume"][0]:
                found[name] = other
                break
    return found


def chain_starts(runs: list[tuple[str, dict]], writer: dict[str, str]) -> dict[str, float]:
    """When each run starts, from the start of the first run of the writers before it: start_ps holds only the
    time into the writer."""
    meta = {name: data["args"] for name, data in runs}
    starts = {}
    for name in meta:
        total, seen, here = 0.0, {name}, name
        while True:
            total += meta[here].get("start_ps") or 0.0
            if here not in writer or writer[here] in seen:
                break
            here = writer[here]
            seen.add(here)
        starts[name] = total
    return starts


def start_of(meta: dict, writer: str | None) -> str:
    if not meta.get("start_ps"):
        return meta["start"]
    if writer is None:
        return f"{meta['start']} at {meta['start_ps']:g} ps"
    return f"{meta['start']}, {meta['start_ps']:g} ps into {writer}"


def dynamics(meta: dict) -> str:
    if meta.get("barostat", "mtk") == "pr":
        return f"Nose-Hoover, Parrinello-Rahman at {meta['bulk_modulus']:g} GPa"
    thermostat = "stochastic velocity rescaling" if meta.get("thermostat") == "svr" else "Langevin"
    return f"{thermostat}, MTK"


def potentials(meta: dict) -> str:
    return f"{meta.get('oxygen', OXYGEN)} and {meta.get('hydrogen', HYDROGEN)}"


def in_place(meta: dict) -> str:
    return meta.get("cell_in_place") or "unset"


def run_rows(runs: list[tuple[str, dict]], writer: dict[str, str]) -> list[str]:
    metas = [data["args"] for _, data in runs]
    columns = [("Basis", basis)]
    if len({in_place(meta) for meta in metas}) > 1:
        columns.append(("CPMDC_CELL_IN_PLACE", in_place))
    if len({potentials(meta) for meta in metas}) > 1:
        columns.append(("Pseudopotentials", potentials))
    if any("barostat" in meta for meta in metas):
        columns.append(("Thermostat, barostat", dynamics))
    head = "| Run | " + " | ".join(title for title, _ in columns) + " | P_tau, fs | Start | Steps |"
    lines = [head, "|" + "+".join("-" * len(cell) for cell in head.split("|")[1:-1]) + "|"]
    for name, data in runs:
        meta = data["args"]
        cells = " | ".join(cell(meta) for _, cell in columns)
        lines.append(f"| {name} | {cells} | {meta['p_tau']:g} | {start_of(meta, writer.get(name))} | "
                     f"{len(data['time_ps']) - 1} |")
    return lines


def summary_rows(results: list[tuple[str, dict]]) -> list[str]:
    """npt_block_average's means of the volume, density and pressure, and the temperature."""
    terms = (*TERMS, ("Temperature", "K", 0))
    lines = ["| Run | Length, ps | Discarded, ps | Frames used | Volume, nm^3 | Density, kg/m^3 | Pressure, bar | "
             "Temperature, K | Wall, s |",
             "|-----+------------+---------------+-------------+--------------+-----------------+---------------+"
             "----------------+---------|"]
    for name, r in results:
        cells = [f"{r['terms'][n]['mean']:.{d}f} +- {r['terms'][n]['sem']:.{d + 1}f}" for n, _, d in terms]
        lines.append(f"| {name} | {r['length_ps']:g} | {r['discard_ps']:g} | {r['production_frames']} | "
                     + " | ".join(cells) + f" | {r['wall_s']:.0f} |")
    return lines


def cost_rows(runs: list[tuple[str, dict]]) -> list[str]:
    lines = ["| Run | Ranks | Wall per step, s | SCF iterations per step | New setups | Wall, h |",
             "|-----+-------+------------------+-------------------------+------------+---------|"]
    for name, data in runs:
        wall = np.asarray(data["step_wall_s"][1:])
        scf = np.asarray(data["scf_steps"][1:])
        lines.append(f"| {name} | {data['args']['ranks']} | {wall.mean():.1f} | {scf.mean():.1f} | "
                     f"{int(np.sum(data['setups'][1:]))} | {np.sum(data['step_wall_s']) / 3600.0:.2f} |")
    return lines


def plateau_rows(results: list[tuple[str, dict]]) -> list[str]:
    """The level npt_block_average.plateau finds for each observable, or the coarsest level it tried."""
    lines = ["| Run | Observable | Block, ps | Blocks | Error | Its uncertainty |",
             "|-----+------------+-----------+--------+-------+-----------------|"]
    for name, r in results:
        for term, unit, d in TERMS:
            s = r["terms"][term]
            if s["plateau"] is None:
                tried = s["block_frames"] * r["frame_ps"]
                lines.append(f"| {name} | {term}, {unit} | none up to {tried:g} |  |  |  |")
                continue
            blocks, sem, err = s["curve"][s["plateau"]]
            lines.append(f"| {name} | {term}, {unit} | {2 ** s['plateau'] * r['frame_ps']:g} | {blocks} | "
                         f"{sem:.{d + 1}f} | {err:.{d + 1}f} |")
    return lines


def branch_groups(runs: list[tuple[str, dict]]) -> list[list[tuple[str, dict]]]:
    """The runs that start from the same state, in groups of at least two."""
    groups: dict[tuple[str, float], list[tuple[str, dict]]] = {}
    for name, data in runs:
        meta = data["args"]
        if meta.get("start_ps"):
            groups.setdefault((meta["start"], meta["start_ps"]), []).append((name, data))
    return [group for group in groups.values() if len(group) >= 2]


def branch_rows(runs: list[tuple[str, dict]], writer: dict[str, str], step_ps: float) -> list[str]:
    """The density every step_ps after the start of the runs that start from the same state."""
    lines = []
    for branches in branch_groups(runs):
        first, meta = branches[0][0], branches[0][1]["args"]
        columns = [{round(float(t), 6): float(rho) for t, rho in zip(data["time_ps"], data["Density"])}
                   for _, data in branches]
        end = max(float(data["time_ps"][-1]) for _, data in branches)
        marks = [round(k * step_ps, 6) for k in range(int(np.floor(end / step_ps + 1e-9)) + 1)]
        head = " | ".join(f"{name}, kg/m^3" for name, _ in branches)
        lines += ["", f"Density after the start of the runs from {start_of(meta, writer.get(first))}:", "",
                  f"| After the start, ps | {head} |", "|-" + "-+-" * len(branches) + "-|"]
        for mark in marks:
            cells = [f"{column[mark]:.1f}" if mark in column else "" for column in columns]
            lines.append(f"| {mark:g} | " + " | ".join(cells) + " |")
    return lines


SHARED = ("waters", "conv", "dt", "temperature", "pressure", "t_tau", "hydrogen_mass")


def table(args: argparse.Namespace) -> int:
    names = [path.stem for path in args.runs]
    if len(set(names)) != len(names):
        raise SystemExit("two records share a name; the table names each run by its file")
    loaded = [json.loads(path.read_text()) for path in args.runs]
    meta = loaded[0]["args"]
    for key in SHARED:
        if len({json.dumps(data["args"][key]) for data in loaded}) > 1:
            raise SystemExit(f"the runs differ in {key}")
    named = list(zip(names, loaded))
    writer = writers(named)
    chain = chain_starts(named, writer)
    branched = {name for group in branch_groups(named) for name, _ in group}
    finished, unaveraged, plotted = [], [], []
    for name, data in named:
        data["wall_s"] = float(np.sum(data["step_wall_s"]))
        start_ps = data["args"].get("start_ps", 0.0)
        floor = args.branch_discard_ps if start_ps else args.min_discard_ps
        length = float(data["time_ps"][-1])
        entry = {"label": name, "basis": basis(data["args"]), "cutoff": data["args"]["cutoff"],
                 "e0": data["args"]["e0"], "p_tau": data["args"]["p_tau"], "start_ps": start_ps,
                 "chain_start_ps": chain[name], "length_ps": length, "windows": windows(data, args.plot_window_ps),
                 "stopped": data["stopped"]}
        if name in branched:
            entry["steps"] = [[float(t), float(rho)] for t, rho in zip(data["time_ps"], data["Density"])]
        if data["stopped"] is None and length - floor >= args.min_production_ps - 1e-9:
            result = analyse(data, args.discard_ps, floor, args.window_ps, args.min_blocks)
            temperature = np.asarray(data["Temperature"])[len(data["time_ps"]) - result["production_frames"]:]
            result["terms"]["Temperature"] = {"mean": float(temperature.mean()), "std": float(temperature.std()),
                                              **block_error(temperature, args.min_blocks)}
            volume = np.asarray(data["Volume"])[len(data["time_ps"]) - result["production_frames"]:]
            scale = 1e-27 * 1e5 / (KB_J * data["args"]["temperature"] * volume.mean())
            square = (volume - volume.mean()) ** 2
            result["kappa"] = (float(square.mean()) * scale, block_error(square, args.min_blocks)["sem"] * scale)
            finished.append((name, result))
            entry.update(mean=result["terms"]["Density"]["mean"], sem=result["terms"]["Density"]["sem"],
                         discard_ps=result["discard_ps"])
        else:
            why = data["stopped"] or (f"{length:g} ps, less than {floor:g} ps to discard and "
                                      f"{args.min_production_ps:g} ps to average")
            unaveraged.append((name, data, why))
        plotted.append(entry)
    shared = {potentials(data["args"]) for data in loaded}
    cells = {in_place(data["args"]) for data in loaded}
    lines = [
        "# Generated by bench/constant_cutoff_npt.py table; do not edit by hand.",
        f"# {meta['waters']} waters, PBE, "
        + (shared.pop() if len(shared) == 1 else "the pseudopotentials of the runs table")
        + f", orbital convergence {meta['conv']:g}, "
        f"ASPC order {meta.get('aspc_order', 4)}; "
        + ("LangevinBAOAB hydrostatic" if all("barostat" not in data["args"] for data in loaded)
           else "hydrostatic NPT as the runs table gives")
        + f", {meta['temperature']:g} K, {meta['pressure']:g} bar, "
        f"T_tau {meta['t_tau']:g} fs, {meta['dt']:g} fs steps, "
        f"hydrogen mass {meta['hydrogen_mass']:g} u; "
        + (f"CPMDC_CELL_IN_PLACE={cells.pop()}" if len(cells) == 1 else "CPMDC_CELL_IN_PLACE as the runs table gives"),
        f"# {meta['cpu']}, ASE {meta['ase']}; the cost table gives the ranks of each run",
    ]
    lines += ["", *run_rows(named, writer)]
    if finished:
        lines += ["", *summary_rows(finished)]
    if unaveraged:
        lines += ["", "| Run | No mean | Density at the start, kg/m^3 | Density at the end, kg/m^3 | "
                  "Mean pressure, bar |",
                  "|-----+---------+------------------------------+----------------------------+--------------------|"]
        for name, data, why in unaveraged:
            lines.append(f"| {name} | {why} | {data['Density'][0]:.1f} | {data['Density'][-1]:.1f} | "
                         f"{np.mean(data['Pressure'][1:]):.0f} |")
    lines += ["", *cost_rows(named)]
    if finished:
        lines += ["", *detail_rows(finished), "", "Standard error of the mean against block length:", "",
                  *curve_rows(finished, args.min_blocks, range(0, 16, args.curve_step))]
        lines += ["", "The first block of B frames, of n, with B^3 > 2 n (SE_B / SE_1)^4 (Lee et al.), among the "
                  "levels with at least --min-blocks blocks:", "", *plateau_rows(finished)]
        lines += ["", "Isothermal compressibility var(V)/(kT <V>) over the averaged frames, its error from blocking "
                  "(V - <V>)^2:", "", "| Run | kappa_T, 1e-5/bar | Error |", "|-----+-------------------+-------|"]
        lines += [f"| {name} | {r['kappa'][0] * 1e5:.2f} | {r['kappa'][1] * 1e5:.2f} |" for name, r in finished]
        if args.prefix_check:
            by_name = dict(named)
            lines += ["", "Equilibration start stability and production statistical inefficiency:", "",
                      "| Run | Discard, ps | Half-run discard, ps | Three-quarter discard, ps | "
                      "Density g_prod | Settled (< 1 ps) |",
                      "|-----+-------------+----------------------+---------------------------+"
                      "----------------+------------------|"]
            for name, r in finished:
                data = by_name[name]
                floor = args.branch_discard_ps if data["args"].get("start_ps", 0.0) else args.min_discard_ps
                n_total = len(data["time_ps"])
                half_n = max(3, n_total // 2)
                three_n = max(3, (3 * n_total) // 4)
                if min(half_n, three_n) < args.min_blocks:
                    raise SystemExit(
                        f"{name}: prefix check has no blocking level with at least "
                        f"{args.min_blocks} blocks"
                    )
                half = {k: v[:half_n] for k, v in data.items() if isinstance(v, list)}
                three_q = {k: v[:three_n] for k, v in data.items() if isinstance(v, list)}
                d_half = analyse(half, args.discard_ps, floor, args.window_ps, args.min_blocks)["discard_ps"]
                d_3q = analyse(three_q, args.discard_ps, floor, args.window_ps, args.min_blocks)["discard_ps"]
                prod_rho = np.asarray(data["Density"])[n_total - r["production_frames"] :]
                g_prod = statistical_inefficiency(prod_rho)
                settled = "yes" if max(abs(d_half - r["discard_ps"]), abs(d_3q - r["discard_ps"])) <= 1.0 else "no"
                lines.append(f"| {name} | {r['discard_ps']:g} | {d_half:g} | {d_3q:g} | {g_prod:.1f} | {settled} |")
    lines += ["", "Mean density over consecutive windows:", "", *density_window_rows(named, args.window_ps)]
    lines += branch_rows(named, writer, args.branch_step_ps)
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
    if args.json:
        args.json.write_text(json.dumps({"settings": lines[1:3], "runs": plotted}, indent=1) + "\n")
    return 0


TIMED = ("waters", "start", "start_ps", *PHYSICS, "ranks", "cpu")
LEVERS = (("Build", lambda meta: meta.get("build") or "unrecorded"),
          ("REAL SPACE WFN KEEP", lambda meta: "yes" if meta.get("rs_wfn_keep") else "no"),
          ("ASPC order", lambda meta: f"{meta.get('aspc_order', 4)}"),
          ("ODIIS vectors", lambda meta: f"{meta.get('odiis_vectors', 0) or 'default'}"),
          ("HAMILTONIAN CUTOFF", lambda meta: f"{meta['hamiltonian_cutoff']:g} Ha" if meta.get("hamiltonian_cutoff")
           else "default"))


def cost(args: argparse.Namespace) -> int:
    runs = [(path.stem, json.loads(path.read_text())) for path in args.runs]
    for key in TIMED:
        if len({json.dumps(data["args"].get(key)) for _, data in runs}) > 1:
            raise SystemExit(f"the runs differ in {key}")
    metas = [data["args"] for _, data in runs]
    varied = [(title, value) for title, value in LEVERS if len({value(meta) for meta in metas}) > 1]
    held = ", ".join(f"{title} {value(metas[0])}" for title, value in LEVERS if (title, value) not in varied)
    against = "the run above" if args.against == "above" else runs[0][0]
    head = ("| Run | " + "".join(f"{title} | " for title, _ in varied) + "Steps | Wall per step, s | "
            f"SCF iterations per step | Largest energy difference from {against}, eV |")
    meta = metas[0]
    lines = ["# Generated by bench/constant_cutoff_npt.py cost; do not edit by hand.",
             f"# {meta['waters']} waters from {meta['start']}, {basis(meta)}, orbital convergence {meta['conv']:g}, "
             f"P_tau {meta['p_tau']:g} fs, seed {meta['seed']}; {held}; {meta['ranks']} rank(s) on {meta['cpu']}; "
             f"means over steps {args.from_step} to {'the last' if args.to_step is None else args.to_step}",
             "", head, "|" + "+".join("-" * len(cell) for cell in head.split("|")[1:-1]) + "|"]
    stop = None if args.to_step is None else args.to_step + 1
    reference = None
    for name, data in runs:
        wall = np.asarray(data["step_wall_s"][args.from_step:stop])
        scf = np.asarray(data["scf_steps"][args.from_step:stop])
        energy = np.asarray(data["Energy"][:stop])
        diff = ""
        if reference is not None:
            n = min(len(reference), len(energy))
            diff = f"{np.max(np.abs(energy[:n] - reference[:n])):.1e}"
        if reference is None or args.against == "above":
            reference = energy
        cells = "".join(f"{value(data['args'])} | " for _, value in varied)
        lines.append(f"| {name} | {cells}{len(wall)} | {wall.mean():.3f} | {scf.mean():.2f} | {diff} |")
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    sys.stdout.write(text)
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
    p.add_argument("--p-tau", type=float, default=100.0, help="fs")
    p.add_argument("--hydrogen-mass", type=float, default=2.014)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--density-range", type=float, nargs=2, default=(500.0, 1500.0), help="kg/m^3")
    p.add_argument("--rs-wfn-keep", action="store_true")
    p.add_argument("--aspc-order", type=int, default=4)
    p.add_argument("--odiis-vectors", type=int, default=0, help="0 keeps the CPMD default")
    p.add_argument("--hamiltonian-cutoff", type=float, default=0.0, help="Hartree; 0 keeps the CPMD default")
    p.add_argument("--thermostat", choices=("langevin", "svr"), default="langevin",
                   help="of the atoms under the MTK barostat")
    p.add_argument("--barostat", choices=("mtk", "pr"), default="mtk",
                   help="pr uses ASE's MelchionnaNPT and its Nose-Hoover thermostat")
    p.add_argument("--bulk-modulus", type=float, default=2.2,
                   help="GPa; sets only the mass of the pr piston, with --p-tau")
    p.add_argument("--oxygen", default=OXYGEN)
    p.add_argument("--hydrogen", default=HYDROGEN)
    p.add_argument("--steps", type=int, default=3000)
    p.add_argument("--checkpoint", type=int, default=10)
    p.add_argument("--keep-every", type=int, default=0, help="also keep the state of every Nth step as "
                   "OUT.N.state.npz")
    p = sub.add_parser("pressure")
    p.add_argument("--library", required=True)
    p.add_argument("--start", type=Path, required=True)
    p.add_argument("--conv", type=float, default=1e-6)
    p.add_argument("--density", type=float, help="kg/m^3; scales the molecules rigidly to it")
    p.add_argument("--scale-volume", type=float, help="volume factor V/V0; scales the molecules rigidly")
    p.add_argument("--oxygen", default=OXYGEN)
    p.add_argument("--hydrogen", default=HYDROGEN)
    p.add_argument("--output", type=Path)
    p.add_argument("bases", nargs="+", help="CUTOFF, or CUTOFF/E0[/WIDTH] in Ry for a constant cutoff")
    p = sub.add_parser("pressures")
    p.add_argument("tables", type=Path, nargs="+")
    p.add_argument("--against", type=int, default=1, help="the column the differences are taken from, from 1")
    p.add_argument("--output", type=Path)
    p = sub.add_parser("pulay")
    p.add_argument("tables", type=Path, nargs="+")
    p.add_argument("--reference", help="reference basis label (default: last constant-cutoff basis)")
    p.add_argument("--output", type=Path)
    p = sub.add_parser("pv")
    p.add_argument("tables", type=Path, nargs="+")
    p.add_argument("--against", help="reference setup string 'O, H, basis' (default: first setup)")
    p.add_argument("--temperature", type=float, default=350.0, help="K, for N_mol k_B T / V")
    p.add_argument("--pressure", type=float, default=1.0, help="bar, target pressure for p(V) root")
    p.add_argument("--degree", type=int, default=2, choices=(1, 2), help="polynomial degree in rho")
    p.add_argument("--no-kinetic", action="store_true", help="fit virial pressure without N_mol k_B T / V")
    p.add_argument("--output", type=Path)
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
    p.add_argument("--curve-step", type=int, default=2, help="blocking levels between the columns of the curve")
    p.add_argument("--prefix-check", action="store_true")
    p.add_argument("--output", type=Path)
    p.add_argument("--json", type=Path)
    p = sub.add_parser("cost")
    p.add_argument("runs", type=Path, nargs="+")
    p.add_argument("--from-step", type=int, default=1)
    p.add_argument("--to-step", type=int, help="the last step of the means; the last of each run by default")
    p.add_argument("--against", choices=("above", "first"), default="above",
                   help="the run whose potential energy each run's is compared with")
    p.add_argument("--output", type=Path)
    args = parser.parse_args()
    return {"start": start, "run": run, "pressure": pressure, "pressures": pressures, "pulay": pulay,
            "pv": pv, "table": table, "cost": cost}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
