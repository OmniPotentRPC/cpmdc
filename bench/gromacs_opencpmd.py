#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""GROMACS QM/MM on OpenCPMD: one QM water in MM water, a step and NVE.

The system: 27 flexible TIP3P waters on a 0.6 nm lattice in a 1.8 nm box;
the one at (0.5, 0.5, 0.5) nm, the centre of the 1 nm cluster cell of
tests/water_blyp_params.capnp.txt, is QM (BLYP, 50 Ry, Hockney). CELL NONE,
mechanical embedding with the TIP3P charges on the QM atoms, TIP3P oxygen
and a repulsive hydrogen for the QM-MM Lennard-Jones, PME for MM-MM.

For each rank count mdrun integrates STEPS steps of DT ps NVE from 300 K
velocities. The report gives the backend libcpmdc_gromacs ran (it must be
opencpmd), the step-0 QM energy and forces against the first rank count,
and the drift of GROMACS' total energy, the conserved quantity of NVE: the
least-squares slope over the run, and its largest excursion from the first
step. The same numbers for the run with every water classical, on one rank,
give the integrator's share: the flexible O-H bonds make leap-frog's total
energy oscillate with the time step squared.

GMX is a GROMACS built against an OpenCPMD build of cpmdc, for example
scripts/ci/build_gromacs_cpmdc.sh 2025.5 PREFIX build-cpmd/tools/gromacs mpi;
libgromacs records the library's path, so no environment selects it.

usage: gromacs_opencpmd.py --gmx GMX --pseudo-dir DIR
                           [--mpiexec MPIEXEC] [--ranks 1,2] [--qm-ranks ALL]
                           [--steps N] [--dt PS] [--convergence C]
                           [--workdir DIR] [--json FILE]
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

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests" / "python"))

import gromacs_case as gc  # noqa: E402
from cpmdc.gromacs import qmi_input  # noqa: E402

BOX_NM = np.eye(3) * 1.8
QM = 0
QM_CHARGES = [gc.TIP3P["O"], gc.TIP3P["H"], gc.TIP3P["H"]]
NATOMS = 81
KJ_PER_HARTREE = 2625.4996394798254


def positions() -> np.ndarray:
    x = gc.water_lattice(3, 0.6) + 0.2
    x[3 * QM : 3 * QM + 3] = [
        [0.50, 0.50, 0.51173],
        [0.50, 0.57572, 0.45308],
        [0.50, 0.42428, 0.45308],
    ]
    return x


def deck(workdir: Path, convergence: float) -> Path:
    text = (ROOT / "tests" / "water_blyp_params.capnp.txt").read_text()
    text, n = re.subn(r"convergenceOrbitals = [0-9.eE+-]+", f"convergenceOrbitals = {convergence:g}", text)
    text, m = re.subn(r"maxStep = \d+", "maxStep = 200", text)
    if n != 1 or m != 1:
        raise RuntimeError("tests/water_blyp_params.capnp.txt changed shape")
    target = workdir / "water_blyp.bin"
    done = subprocess.run(
        ["capnp", "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=text.encode(), stdout=subprocess.PIPE, check=True,
    )
    target.write_bytes(done.stdout)
    return target


def mdp(args, **extra) -> dict:
    # No pair-list buffer and a cutoff off the lattice spacing, so every
    # decomposition sees the same pairs.
    return {"rcoulomb": 0.55, "rvdw": 0.55, "rlist": 0.55, "nstlist": 1,
            "verlet-buffer-tolerance": -1, "nsteps": args.steps, "dt": args.dt,
            "nstxout": 0, "nstvout": 0, "nstfout": args.steps, "nstenergy": 1,
            "gen-vel": "yes", "gen-temp": 300, "gen-seed": 4928, **extra}


def drift(table: np.ndarray) -> dict:
    t, kinetic, total = table[:, 0], table[:, -2], table[:, -1]
    slope = np.polyfit(t, total, 1)[0]
    return {
        "kinetic_en_kj_mol": kinetic.tolist(),
        "total_energy_kj_mol": total.tolist(),
        "time_ps": t.tolist(),
        "drift_kj_mol_ps": float(slope),
        "drift_kj_mol_ps_atom": float(slope / NATOMS),
        "max_excursion_kj_mol": float(np.max(np.abs(total - total[0]))),
        "kinetic_std_kj_mol": float(np.std(kinetic)),
    }


def classical(args) -> dict:
    """The same run with the QM water a TIP3P water: the integrator's share."""
    workdir = args.workdir / "classical"
    gc.write_system(workdir, positions=positions(), box=BOX_NM, qm_molecules=[], qmi="",
                    qm_lj=True, mdp=mdp(args, **{"qmmm-cp2k-active": "false"}))
    mdrun(args, workdir, 1)
    return drift(gc.energies(args.gmx, workdir, "Kinetic-En.", "Total-Energy"))


def mdrun(args, workdir: Path, ranks: int) -> None:
    env = {"CPMDC_PSEUDO_DIR": str(args.pseudo_dir), "OMP_NUM_THREADS": "1"}
    launcher = None if args.mpiexec is None else [args.mpiexec]
    extra = () if launcher is None else ("-npme", "0")
    gc.run(args.gmx, workdir, launcher=launcher, ranks=ranks, mdrun_args=extra, env_extra=env)


def run(args, params: Path, ranks: int) -> dict:
    workdir = args.workdir / f"ranks{ranks}"
    gc.write_system(
        workdir,
        positions=positions(),
        box=BOX_NM,
        qm_molecules=[QM],
        qmi=qmi_input(params.resolve(), embedding="MECHANICAL", qm_charges=QM_CHARGES,
                      qm_ranks=args.qm_ranks),
        qm_lj=True,
        mdp=mdp(args),
    )
    mdrun(args, workdir, ranks)
    out = (workdir / "md_cp2k.out").read_text()
    first = re.search(r"first step: backend (\S+) qm_energy_hartree (\S+)", out)
    table = gc.energies(args.gmx, workdir, "Quantum-En.", "Kinetic-En.", "Total-Energy")
    frames = gc.read_trr(workdir / "md.trr")
    wall = re.search(r"Time:\s+\S+\s+(\S+)", (workdir / "md.log").read_text())
    return {
        "ranks": ranks,
        "backend": first.group(1) if first else None,
        "qm_energy_hartree_step0": float(first.group(2)) if first else None,
        "quantum_en_kj_mol": table[:, 1].tolist(),
        **drift(table),
        "forces_step0": frames[0]["f"].tolist(),
        "forces_last": frames[-1]["f"].tolist(),
        "wall_s": float(wall.group(1)) if wall else None,
    }


def report(label: str, r: dict) -> None:
    print(
        f"  {label}total energy drift {r['drift_kj_mol_ps']:.4g} kJ/mol/ps "
        f"({r['drift_kj_mol_ps_atom']:.3g} per atom), max excursion "
        f"{r['max_excursion_kj_mol']:.4g} kJ/mol, kinetic energy std "
        f"{r['kinetic_std_kj_mol']:.4g} kJ/mol"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--pseudo-dir", type=Path, required=True)
    parser.add_argument("--mpiexec")
    parser.add_argument("--ranks", default="1")
    parser.add_argument("--qm-ranks", default="ALL")
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--dt", type=float, default=0.00025)
    parser.add_argument("--convergence", type=float, default=1e-7)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    args.pseudo_dir = args.pseudo_dir.resolve()
    rank_counts = [int(r) for r in args.ranks.split(",")]
    if args.mpiexec is None and rank_counts != [1]:
        parser.error("--ranks other than 1 needs --mpiexec and an MPI GROMACS")
    args.workdir = (args.workdir or Path(tempfile.mkdtemp(prefix="gromacs-opencpmd-"))).resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)
    params = deck(args.workdir, args.convergence)

    control = classical(args)
    results = [run(args, params, n) for n in rank_counts]
    ref = results[0]
    print(f"workdir {args.workdir}")
    print(f"{NATOMS} atoms, 3 QM; {args.steps} steps of {args.dt * 1000:g} fs; "
          f"orbital convergence {args.convergence:g}; qm_ranks {args.qm_ranks}")
    ok = True
    for r in results:
        if r["backend"] != "opencpmd":
            ok = False
        f, f_ref = np.array(r["forces_step0"]), np.array(ref["forces_step0"])
        g, g_ref = np.array(r["forces_last"]), np.array(ref["forces_last"])
        e, e_ref = np.array(r["total_energy_kj_mol"]), np.array(ref["total_energy_kj_mol"])
        print(
            f"ranks {r['ranks']}: backend {r['backend']}; step 0 E_QM "
            f"{r['qm_energy_hartree_step0']:.12f} Ha, |dE_QM| vs ranks {ref['ranks']} "
            f"{abs(r['qm_energy_hartree_step0'] - ref['qm_energy_hartree_step0']) * KJ_PER_HARTREE:.3e} kJ/mol, "
            f"max|dF| {np.max(np.abs(f - f_ref)):.3e} kJ/mol/nm (max|F| {np.max(np.abs(f_ref)):.1f}); "
            f"over the run max|dE_total| {np.max(np.abs(e - e_ref)):.3e} kJ/mol, "
            f"last step max|dF| {np.max(np.abs(g - g_ref)):.3e} kJ/mol/nm; wall {r['wall_s']} s"
        )
        report("", r)
    report("classical control: ", control)
    if args.json:
        args.json.write_text(json.dumps({"args": {k: str(v) for k, v in vars(args).items()},
                                         "classical": control, "results": results}, indent=1))
    if not ok:
        print("libcpmdc_gromacs did not run OpenCPMD: GMX is linked to a build without it")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
