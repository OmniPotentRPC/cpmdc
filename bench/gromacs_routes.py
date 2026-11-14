#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Per-step overhead of coupling GROMACS to cpmdc, beside the LAMMPS plugin.

Every route integrates the system of lammps_routes.py: two water molecules
in a 0.9 nm (9 Angstrom) cubic box, NVE at 0.25 fs from 300 K velocities,
on the reference evaluator with the periodic deck. The routes:

  direct         Session.calculate_au in a Python loop, no integrator
  lammps-zero    LAMMPS with pair_style zero (LAMMPS alone)
  plugin         LAMMPS with pair_style cpmdc
  gromacs-zero   GROMACS mdrun, QM/MM off, atoms without charges, LJ or
                 bonds (GROMACS alone)
  gromacs-cpmdc  GROMACS mdrun, both waters in the QM group, libcpmdc_gromacs
                 through the CP2K QM/MM interface, CELL PERIODIC

The LAMMPS routes are timed as in lammps_routes.py. A GROMACS route runs
mdrun twice, for SHORT and SHORT + STEPS steps, and divides the difference
of the wall times md.log reports by STEPS, so process startup, grompp,
input parsing and the first evaluation cancel. mdrun runs one thread-MPI
rank with one OpenMP thread and no output inside the loop. Each repeat
starts fresh; the summary is the median per step with every raw value.

usage: gromacs_routes.py [--build DIR] [--gmx GMX] [--steps N] [--short N]
                         [--repeats N] [--routes a,b,...] [--json FILE]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests" / "python"))
sys.path.insert(0, str(ROOT / "bench"))

import gromacs_case as gc  # noqa: E402
import lammps_routes  # noqa: E402
from cpmdc.gromacs import qmi_input  # noqa: E402

ROUTES = ("direct", "lammps-zero", "plugin", "gromacs-zero", "gromacs-cpmdc")
EDGE_NM = lammps_routes.EDGE / 10.0
WATERS_NM = (np.asarray(lammps_routes.WATERS) / 10.0).tolist()

BARE = """[ defaults ]
1 2 no 1.0 1.0

[ atomtypes ]
OZ 8 15.9994 0.0 A 0.0 0.0
HZ 1 1.008 0.0 A 0.0 0.0

[ moleculetype ]
BARE 0

[ atoms ]
1 OZ 1 BARE O 1 0.0 15.9994
2 HZ 1 BARE H1 1 0.0 1.008
3 HZ 1 BARE H2 1 0.0 1.008

[ system ]
bare waters

[ molecules ]
BARE {n}
"""

WALL = re.compile(r"^\s*Time:\s+(\S+)\s+(\S+)", re.M)


def gromacs_wall(args, route: str, steps: int, work: Path) -> float:
    """md.log's wall time of one mdrun of ``steps`` steps."""
    mdp = {
        "nsteps": steps,
        "coulombtype": "Cut-off",
        "rcoulomb": 0.4,
        "rvdw": 0.4,
        "nstxout": 0,
        "nstvout": 0,
        "nstfout": 0,
        "nstenergy": 0,
        "nstlog": 0,
        "nstcalcenergy": 100,
        "gen-vel": "yes",
        "gen-temp": 300,
        "gen-seed": 4928,
    }
    if route == "gromacs-zero":
        mdp["qmmm-cp2k-active"] = "false"
    qmi = qmi_input(args.params, cell="PERIODIC")
    gc.write_system(
        work, positions=WATERS_NM, box=np.eye(3) * EDGE_NM, qm_molecules=[0, 1],
        qmi=qmi, mdp=mdp,
    )
    (work / "topol.top").write_text(BARE.format(n=len(WATERS_NM) // 3))
    gc.run(args.gmx, work, with_qmi=route == "gromacs-cpmdc")
    match = WALL.search((work / "md.log").read_text())
    if match is None:
        raise RuntimeError(f"{work / 'md.log'} has no Time: line")
    return float(match.group(2))


def time_gromacs(args, route: str) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        short = gromacs_wall(args, route, args.short, Path(tmp) / "short")
        long = gromacs_wall(args, route, args.short + args.steps, Path(tmp) / "long")
    return (long - short) / args.steps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--gmx", default=os.environ.get("CPMDC_GMX") or shutil.which("gmx_d"))
    parser.add_argument("--steps", type=int, default=20000)
    parser.add_argument("--short", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--routes", default=",".join(ROUTES))
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    args.build = args.build.resolve()
    args.params = args.build / "water_blyp_periodic_params.bin"
    args.library = args.build / "libcpmdc.so"
    args.plugin = args.build / "cpmdcplugin.so"
    os.environ.setdefault("CPMDC_LIBRARY", str(args.library))

    results = {}
    for route in args.routes.split(","):
        if route not in ROUTES:
            parser.error(f"unknown route {route!r}; choose from {', '.join(ROUTES)}")
        if route.startswith("gromacs") and not args.gmx:
            parser.error("no gmx_d: pass --gmx or set CPMDC_GMX")
        raw = []
        for _ in range(args.repeats):
            if route == "direct":
                raw.append(lammps_routes.time_direct(args))
            elif route.startswith("gromacs"):
                raw.append(time_gromacs(args, route))
            else:
                raw.append(lammps_routes.time_lammps(args, route))
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
        meta = {
            "steps": args.steps,
            "short": args.short,
            "repeats": args.repeats,
            "atoms": len(WATERS_NM),
        }
        args.json.write_text(json.dumps({"setup": meta, "routes": results}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
