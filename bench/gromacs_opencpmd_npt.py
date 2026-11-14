#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""GROMACS QM/MM NPT on OpenCPMD: one BLYP water in TIP3P water.

By default the system of bench/gromacs_opencpmd.py: 27 flexible TIP3P
waters on a 0.6 nm lattice in a 1.8 nm box (0.14 g/cm^3), the one at the
centre of the deck's cluster cell QM (BLYP, Hockney, CELL NONE, mechanical
embedding). --lattice 4 --spacing 0.31 --minimize gives 64 waters at
1.00 g/cm^3, relaxed by a classical steepest descent first. The lattice
is shifted to put a site, the QM water, at the centre of the cluster cell.
The deck predicts the orbitals of each step with ASPC of order 4.

mdrun integrates STEPS steps of DT ps from 300 K velocities with the
v-rescale thermostat at 300 K and isotropic C-rescale or
Parrinello-Rahman pressure coupling at 1 bar, once with the QM water and
once with every water classical. The pressure, and so the box, then carry
the QM and QM-MM virial of libcpmdc_gromacs.

The report gives, for both runs, the mean and standard deviation of the
volume, density and pressure over the run and over its second half, the
volume at the first and the last step, and for C-rescale the least-squares
slope of GROMACS' conserved energy, which includes the work of the
barostat. GMX is a GROMACS built with scripts/ci/build_gromacs_cpmdc.sh
against an OpenCPMD build of cpmdc, so it carries the virial patch.
--classical-only runs the control alone, which needs no pseudopotentials
and costs little enough to run far longer than the QM run.
bench/npt_block_average.py turns the --json file into block-averaged means.

usage: gromacs_opencpmd_npt.py --gmx GMX (--pseudo-dir DIR | --classical-only)
                               [--mpiexec MPIEXEC] [--ranks N]
                               [--pcoupl C-rescale|Parrinello-Rahman]
                               [--lattice N] [--spacing NM] [--minimize]
                               [--steps N] [--dt PS] [--tau-p PS]
                               [--convergence C] [--workdir DIR] [--json FILE]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gromacs_opencpmd import QM_CHARGES, ROOT, gc, qmi_input  # noqa: E402

TERMS = ("Volume", "Density", "Pressure", "Temperature", "Conserved-En.")
CLUSTER_CENTRE_NM = np.array([0.5, 0.5, 0.5])
QM_WATER_NM = CLUSTER_CENTRE_NM + np.array([[0.0, 0.0, 0.01173], [0.0, 0.07572, -0.04692],
                                            [0.0, -0.07572, -0.04692]])


def system(args) -> tuple[np.ndarray, np.ndarray, int]:
    """An n^3 lattice shifted to put a site at the cluster centre, which takes the QM water.

    3 and 0.6 nm is the system of gromacs_opencpmd.py.
    """
    n, spacing = args.lattice, args.spacing
    site = int(round(CLUSTER_CENTRE_NM[0] / spacing - 0.5))
    x = gc.water_lattice(n, spacing) + CLUSTER_CENTRE_NM[0] - (site + 0.5) * spacing
    qm = site * n * n + site * n + site
    x[3 * qm : 3 * qm + 3] = QM_WATER_NM
    return x, np.eye(3) * n * spacing, qm


def deck(workdir: Path, convergence: float) -> Path:
    text = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text()
    text, n = re.subn(r"convergenceOrbitals = [0-9.eE+-]+", f"convergenceOrbitals = {convergence:g}", text)
    text, m = re.subn(r"maxStep = \d+,", 'maxStep = 200, extrapolateWfnOptions = "ASPC", '
                      "extrapolateWfnOrder = 4,", text)
    if n != 1 or m != 1:
        raise RuntimeError("tests/water_blyp_params.capnp.txt changed shape")
    target = workdir / "water_blyp_aspc.bin"
    done = subprocess.run(
        ["capnp", "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode(), stdout=subprocess.PIPE, check=True,
    )
    target.write_bytes(done.stdout)
    return target


def mdp(args, **extra) -> dict:
    return {
        "rcoulomb": 0.55, "rvdw": 0.55, "nsteps": args.steps, "dt": args.dt,
        "nstxout": 0, "nstvout": 0, "nstfout": 0, "nstenergy": 10, "nstcalcenergy": 10,
        "tcoupl": "v-rescale", "tc-grps": "System", "tau-t": 0.1, "ref-t": 300,
        "pcoupl": args.pcoupl, "pcoupltype": "isotropic", "tau-p": args.tau_p, "ref-p": 1.0,
        "compressibility": 4.5e-5, "nstpcouple": 10, "refcoord-scaling": "com",
        "comm-mode": "Linear", "gen-vel": "yes", "gen-temp": 300, "gen-seed": 4928, **extra,
    }


def summary(workdir: Path, gmx: str, conserved: bool) -> dict:
    columns = gc.energy_columns(gmx, workdir, *TERMS)
    t = columns["Time"]
    half = t >= t[-1] / 2
    out = {"time_ps": t.tolist()}
    for name in TERMS:
        column = columns[name]
        out[name] = column.tolist()
        out[f"{name} mean"] = float(column.mean())
        out[f"{name} std"] = float(column.std())
        out[f"{name} mean second half"] = float(column[half].mean())
        out[f"{name} std second half"] = float(column[half].std())
    if conserved:
        energy = columns["Conserved-En."]
        out["conserved drift kj_mol_ps"] = float(np.polyfit(t, energy, 1)[0])
        out["conserved max excursion kj_mol"] = float(np.abs(energy - energy[0]).max())
    out["volume first nm3"] = float(columns["Volume"][0])
    out["volume last nm3"] = float(columns["Volume"][-1])
    log = (workdir / "md.log").read_text()
    wall = re.search(r"Time:\s+\S+\s+(\S+)", log)
    out["wall_s"] = float(wall.group(1)) if wall else None
    return out


def mdrun(args, workdir: Path, ranks: int = 1) -> None:
    env = {"CPMDC_PSEUDO_DIR": str(args.pseudo_dir), "OMP_NUM_THREADS": "1"}
    if args.mpiexec:
        gc.run(args.gmx, workdir, launcher=args.mpiexec.split(), ranks=ranks,
               mdrun_args=("-npme", "0"), env_extra=env)
    else:
        gc.run(args.gmx, workdir, env_extra=env)


def minimized(args, x: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Classical steepest descent of every water, nm."""
    workdir = args.workdir / "minimize"
    gc.write_system(workdir, positions=x, box=box, qm_molecules=[], qmi="", qm_lj=True,
                    mdp={"integrator": "steep", "nsteps": 2000, "emtol": 100.0, "rcoulomb": 0.55,
                         "rvdw": 0.55, "nstxout": 0, "nstvout": 0, "nstfout": 0,
                         "qmmm-cp2k-active": "false"})
    mdrun(args, workdir)
    lines = (workdir / "md.gro").read_text().splitlines()[2 : 2 + len(x)]
    return np.array([[float(v) for v in line[20:].split()[:3]] for line in lines])


def run(args, x: np.ndarray, box: np.ndarray, qm: int, params: Path | None) -> dict:
    workdir = args.workdir / ("qmmm" if params else "classical")
    if params:
        qmi = qmi_input(params.resolve(), embedding="MECHANICAL", qm_charges=QM_CHARGES,
                        qm_ranks="ALL")
        setup = {"qm_molecules": [qm], "qmi": qmi, "mdp": mdp(args)}
        ranks = args.ranks
    else:
        setup = {"qm_molecules": [], "qmi": "", "mdp": mdp(args, **{"qmmm-cp2k-active": "false"})}
        ranks = 1
    gc.write_system(workdir, positions=x, box=box, qm_lj=True, **setup)
    mdrun(args, workdir, ranks)
    out = summary(workdir, args.gmx, args.pcoupl == "C-rescale")
    if params:
        text = (workdir / "md_cp2k.out").read_text()
        backend = re.search(r"first step: backend (\S+)", text)
        out["backend"] = backend.group(1) if backend else None
        out["stress"] = "first stress request: available 1" in text
        out["missing stress logged"] = "returned no stress tensor" in (workdir / "md.log").read_text()
    return out


def report(label: str, r: dict, conserved: bool) -> None:
    print(f"{label}: wall {r['wall_s']} s")
    for name, unit in (("Volume", "nm^3"), ("Density", "kg/m^3"), ("Pressure", "bar"),
                       ("Temperature", "K")):
        print(f"  {name:11s} mean {r[f'{name} mean']:.6g} std {r[f'{name} std']:.4g} {unit}; "
              f"second half mean {r[f'{name} mean second half']:.6g} "
              f"std {r[f'{name} std second half']:.4g}")
    print(f"  volume first {r['volume first nm3']:.6f} last {r['volume last nm3']:.6f} nm^3")
    if conserved:
        print(f"  conserved energy drift {r['conserved drift kj_mol_ps']:.4g} kJ/mol/ps, "
              f"max excursion {r['conserved max excursion kj_mol']:.4g} kJ/mol")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--pseudo-dir", type=Path)
    parser.add_argument("--classical-only", action="store_true")
    parser.add_argument("--mpiexec")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--pcoupl", choices=["C-rescale", "Parrinello-Rahman"], default="C-rescale")
    parser.add_argument("--lattice", type=int, default=3)
    parser.add_argument("--spacing", type=float, default=0.6)
    parser.add_argument("--minimize", action="store_true")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--dt", type=float, default=0.0005)
    parser.add_argument("--tau-p", type=float, default=0.5)
    parser.add_argument("--convergence", type=float, default=1e-6)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.pseudo_dir is None and not args.classical_only:
        parser.error("the QM run needs --pseudo-dir")
    args.pseudo_dir = args.pseudo_dir.resolve() if args.pseudo_dir else ""
    if args.mpiexec is None and args.ranks != 1:
        parser.error("--ranks other than 1 needs --mpiexec and an MPI GROMACS")
    args.workdir = (args.workdir or Path(tempfile.mkdtemp(prefix="gromacs-opencpmd-npt-"))).resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)

    x, box, qm = system(args)
    if args.minimize:
        x = minimized(args, x, box)
    control = run(args, x, box, qm, None)
    qmmm = None if args.classical_only else run(args, x, box, qm, deck(args.workdir, args.convergence))
    conserved = args.pcoupl == "C-rescale"
    print(f"workdir {args.workdir}")
    print(f"{len(x)} atoms, 3 QM, box {box[0, 0]:g} nm; {args.steps} steps of {args.dt * 1000:g} fs; "
          f"{args.pcoupl} tau-p {args.tau_p:g} ps, v-rescale 300 K; orbital convergence "
          f"{args.convergence:g}; QM run on {args.ranks} rank(s)")
    if qmmm:
        print(f"backend {qmmm['backend']}; stress available {qmmm['stress']}; "
              f"missing stress logged {qmmm['missing stress logged']}")
        report("QM/MM (one BLYP water)", qmmm, conserved)
    report("classical control", control, conserved)
    if args.json:
        runs = {"classical": control} if qmmm is None else {"qmmm": qmmm, "classical": control}
        args.json.write_text(json.dumps({"args": {k: str(v) for k, v in vars(args).items()}, **runs},
                                        indent=1))
    if qmmm is None:
        return 0
    ok = qmmm["backend"] == "opencpmd" and qmmm["stress"] and not qmmm["missing stress logged"]
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
