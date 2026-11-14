"""Build one LAMMPS system on cpmdc, in process or under mpiexec.

``test_lammps.py`` imports :func:`build`, :func:`attach`, and
:func:`measure`. Run as ``python lammps_case.py SPEC.json`` (also under
``mpiexec``), the script builds the system the JSON spec describes, runs
its commands, and prints one JSON line of results from world rank 0.

A spec holds ``units``, ``cell`` (three rows in LAMMPS units, or ``None``
for a box without periodic axes), ``positions`` in LAMMPS units inside the
box, ``types``, ``masses``, ``elements``, ``route`` (``"fix"`` for
``cpmdc.lammps.FixExternal`` or ``"plugin"`` for ``pair_style cpmdc``),
``params``, ``plugin``, and optional ``setup`` and ``run`` command lists.
"""

from __future__ import annotations

import json
import sys
from typing import Any

import numpy as np

CMDARGS = ["-log", "none", "-screen", "none", "-nocite"]


def build(spec: dict[str, Any]) -> Any:
    from lammps import lammps

    lmp = lammps(cmdargs=CMDARGS)
    lines = [f"units {spec['units']}", "atom_modify map array sort 0 0"]
    lines += spec.get("setup", [])
    cell = spec.get("cell")
    if cell is None:
        lo, hi = spec.get("block", (-20.0, 20.0))
        lines += ["boundary f f f", f"region box block {lo} {hi} {lo} {hi} {lo} {hi}"]
    else:
        (lx, _, _), (xy, ly, _), (xz, yz, lz) = cell
        lines += [
            "boundary p p p",
            f"region box prism 0 {lx!r} 0 {ly!r} 0 {lz!r} {xy!r} {xz!r} {yz!r}",
        ]
    lines.append(f"create_box {len(spec['masses'])} box")
    lines += [f"mass {index + 1} {mass!r}" for index, mass in enumerate(spec["masses"])]
    lmp.commands_string("\n".join(lines))
    positions = np.asarray(spec["positions"], dtype=float)
    count = len(positions)
    created = lmp.create_atoms(
        count,
        list(range(1, count + 1)),
        list(spec["types"]),
        positions.ravel().tolist(),
    )
    if created != count:
        raise RuntimeError(f"LAMMPS created {created} of {count} atoms")
    return lmp


def attach(lmp: Any, spec: dict[str, Any]) -> Any:
    """Couple cpmdc; return the FixExternal driver, or None for the plugin."""
    if spec["route"] == "fix":
        from cpmdc.lammps import FixExternal

        return FixExternal(lmp, spec["params"], spec["elements"])
    elements = spec["elements"]
    if not isinstance(elements, str):
        elements = " ".join(str(item) for item in elements)
    lmp.command(f"plugin load {spec['plugin']}")
    lmp.command(f"pair_style cpmdc {spec['params']}")
    lmp.command(f"pair_coeff * * {elements}")
    return None


def measure(lmp: Any) -> dict[str, Any]:
    """Potential energy, forces by atom ID, and the virial ``xx yy zz xy xz yz``.

    The virial is the ``compute pressure NULL virial`` tensor times the
    volume over ``nktv2p``, in LAMMPS energy units. Call after a run whose
    thermo output includes ``c_cpmdc_vir``.
    """
    pressure = np.array(lmp.numpy.extract_compute("cpmdc_vir", 0, 1), dtype=float)
    volume = lmp.get_thermo("vol")
    return {
        "pe": lmp.get_thermo("pe"),
        "forces": np.array(lmp.gather_atoms("f", 1, 3), dtype=float)
        .reshape(-1, 3)
        .tolist(),
        "virial": (pressure * volume / lmp.extract_global("nktv2p")).tolist(),
    }


THERMO = [
    "compute cpmdc_vir all pressure NULL virial",
    "compute cpmdc_proc all property/atom proc",
    "thermo_style custom step pe ke etotal econserve vol press c_cpmdc_vir[*]",
]


def main(path: str) -> int:
    with open(path, encoding="utf-8") as handle:
        spec = json.load(handle)
    lmp = build(spec)
    driver = attach(lmp, spec)
    lmp.commands_string("\n".join(THERMO))
    for command in spec.get("run", ["run 0"]):
        if driver is not None and command.startswith("run "):
            driver.run(*command.split()[1:])
        else:
            lmp.command(command)
    out = measure(lmp)
    procs = np.asarray(lmp.gather("c_cpmdc_proc", 1, 1), dtype=int)
    size = lmp.extract_setting("world_size")
    out["world_size"] = size
    out["nlocal"] = np.bincount(procs, minlength=size).tolist()
    if driver is not None:
        driver.close()
    if lmp.extract_setting("world_rank") == 0:
        print(json.dumps(out))
    lmp.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
