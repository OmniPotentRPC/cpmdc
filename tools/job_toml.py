#!/usr/bin/env python3
"""Compile one job TOML into the CPMDParams text both routes read.

The in-process route encodes this text (`capnp encode schema/Potentials.capnp
CPMDParams`) and passes the binary as `params_path`. The file route renders
that same message with `cpmdc_params_render_input_deck`. `file_route_text`
names the keywords that renderer writes for the fields this TOML sets.

A `[constant_cutoff]` table sets `CONSTANT CUTOFF` in `&SYSTEM`, the
constant-cutoff kinetic functional of Bernasconi et al. (1995) whose stress
carries no Pulay stress: `cutoff_ry` is E0, `height_ry` is A (50 times
`cutoff_ry` of the job when omitted), and `width_ry` is sigma (E0 / 10 when
omitted), all in Ry.
"""

from __future__ import annotations

import argparse
import math
import sys
import tomllib
from pathlib import Path
from typing import Any, Mapping


def load_job(path: Path) -> dict[str, Any]:
    """Read a job TOML. `functional` and `cutoff_ry` are required."""
    with path.open("rb") as handle:
        job = tomllib.load(handle)
    if not isinstance(job, dict):
        raise TypeError("job TOML must be a table")
    missing = [key for key in ("functional", "cutoff_ry") if key not in job]
    if missing:
        raise KeyError(",".join(missing))
    constant_cutoff(job)
    return job


_CONSTANT_CUTOFF_KEYS = ("cutoff_ry", "height_ry", "width_ry")


def constant_cutoff(job: Mapping[str, Any]) -> tuple[float, float, float] | None:
    """A, sigma, and E0 in Ry from `[constant_cutoff]`, or None without one.

    Raises ValueError for an unknown key, a value that is not a finite
    positive number, or E0 + 2 sigma above `cutoff_ry` of the job.
    """
    table = job.get("constant_cutoff")
    if table is None:
        return None
    if not isinstance(table, Mapping):
        raise ValueError("constant_cutoff must be a table")
    unknown = sorted(set(table) - set(_CONSTANT_CUTOFF_KEYS))
    if unknown:
        raise ValueError(f"unknown constant_cutoff keys: {unknown}")
    if "cutoff_ry" not in table:
        raise ValueError("constant_cutoff needs cutoff_ry")
    values = {}
    for key in _CONSTANT_CUTOFF_KEYS:
        if key not in table:
            continue
        value = table[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"constant_cutoff.{key} must be a number")
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"constant_cutoff.{key} must be finite and positive")
        values[key] = float(value)
    basis = float(job["cutoff_ry"])
    e0 = values["cutoff_ry"]
    width = values.get("width_ry", 0.1 * e0)
    height = values.get("height_ry", 50.0 * basis)
    if e0 + 2.0 * width > basis:
        raise ValueError(
            f"constant_cutoff.cutoff_ry {e0:g} Ry plus twice its width "
            f"{width:g} Ry exceeds cutoff_ry {basis:g} Ry"
        )
    return height, width, e0


def _fmt(value: float) -> str:
    text = f"{float(value):.10g}"
    if "e" in text:
        mant, exp = text.split("e")
        if "." not in mant:
            mant += ".0"
        return f"{mant}e{int(exp)}"
    if "." not in text:
        text += ".0"
    return text


def _lmax_letter(lmax: int) -> str:
    if lmax <= 0:
        return "S"
    if lmax == 1:
        return "P"
    if lmax == 2:
        return "D"
    return "F"


def _directives(items: list[Mapping[str, Any]] | None) -> str:
    if not items:
        return ""
    lines = ["directives = ["]
    for item in items:
        keyword = str(item["keyword"])
        args = item.get("args") or []
        arg_text = ", ".join(f'"{arg}"' for arg in args)
        lines.append(f'  ( keyword = "{keyword}", args = [{arg_text}] ),')
    lines.append("]")
    return "\n        ".join(lines)


def capnp_text(job: Mapping[str, Any]) -> str:
    """CPMDParams text for `capnp encode`. Field names are the schema names."""
    wave = job.get("wavefunction") or {}
    cell = job.get("cell") or {}
    dft = job.get("dft") or {}
    lines = [
        "(",
        f'  functional = "{job["functional"]}",',
        f'  cutOffRy = {_fmt(float(job["cutoff_ry"]))},',
        f'  charge = {int(job.get("charge", 0))},',
        f'  multiplicity = {int(job.get("multiplicity", 1))},',
        f'  task = "{job.get("task", "gradient")}",',
    ]
    if job.get("title"):
        lines.append(f'  title = "{job["title"]}",')
    blocks = job.get("blocks") or []
    if blocks:
        lines.append("  inputBlocks = [")
        for block in blocks:
            escaped = str(block).replace("\\", "\\\\").replace('"', '\\"')
            lines.append(f'    "{escaped}",')
        lines.append("  ],")
    lines.append("  inputSections = [")
    cpmd_bits = [
        f'optimizeWavefunction = {"true" if wave.get("optimize", True) else "false"}',
        f'convergenceOrbitals = {_fmt(float(wave.get("convergence_orbitals", 1.0e-5)))}',
    ]
    if int(wave.get("max_iter", 0) or 0) > 0:
        cpmd_bits.append(f'maxIter = {int(wave["max_iter"])}')
    directive = _directives(wave.get("directives"))
    if directive:
        cpmd_bits.append(directive)
    lines.append("    ( cpmd = ( " + ", ".join(cpmd_bits) + " ) ),")
    dft_bits = [
        f'functional = "{job["functional"]}"',
        f'newCode = {"true" if dft.get("newcode", False) else "false"}',
    ]
    gc = float(dft.get("gc_cutoff", 0.0) or 0.0)
    if gc > 0.0:
        dft_bits.append(f"gcCutoff = {_fmt(gc)}")
    lines.append("    ( dft = ( " + ", ".join(dft_bits) + " ) ),")
    system_bits = [
        f'symmetry = {int(cell.get("symmetry", 0))}',
        f'angstrom = {"true" if cell.get("angstrom", True) else "false"}',
        f'cutOffRy = {_fmt(float(job["cutoff_ry"]))}',
    ]
    if cell.get("absolute"):
        system_bits.append("cellAbsolute = true")
    if cell.get("degree"):
        system_bits.append("cellDegree = true")
    if cell.get("vectors"):
        system_bits.append(
            "cell = [" + ", ".join(_fmt(float(v)) for v in cell["vectors"]) + "]"
        )
    kinetic = constant_cutoff(job)
    if kinetic:
        system_bits.append("constantCutoff = [" + ", ".join(_fmt(v) for v in kinetic) + "]")
    sys_dir = _directives(cell.get("directives"))
    if sys_dir:
        system_bits.append(sys_dir)
    lines.append("    ( system = ( " + ", ".join(system_bits) + " ) ),")
    psps = job.get("pseudopotentials") or []
    if psps:
        lines.append("    ( atoms = ( pseudopotentials = [")
        for psp in psps:
            bits = [
                f'element = "{psp["element"]}"',
                f'path = "{psp["file"]}"',
                f'lmax = {int(psp["lmax"])}',
                f'loc = {int(psp["loc"])}',
                "kleinmanBylander = "
                + ("true" if psp.get("kleinman_bylander", False) else "false"),
            ]
            lines.append("      ( " + ", ".join(bits) + " ),")
        lines.append("    ] ) ),")
    disp = job.get("dispersion") or {}
    if disp:
        vdw = [
            "empiricalCorrection = true",
            f'grimme = "{disp.get("grimme", "D2")}"',
        ]
        lines.append("    ( vdwParams = ( " + ", ".join(vdw) + " ) ),")
    lines.append("  ]")
    lines.append(")")
    return "\n".join(lines) + "\n"


def file_route_text(job: Mapping[str, Any]) -> str:
    """Keywords `cpmdc_params_render_input_deck` writes for this job.

    Coordinate counts stay at zero until a geometry is merged. The letters
    follow `lmax_letter` in `src/cpmdc_params.c`.
    """
    dft = job.get("dft") or {}
    cell = job.get("cell") or {}
    lines = [
        f'FUNCTIONAL {job["functional"]}',
        "CUTOFF",
        _fmt(float(job["cutoff_ry"])),
    ]
    gc = float(dft.get("gc_cutoff", 0.0) or 0.0)
    if gc > 0.0:
        lines.extend(["GC-CUTOFF", _fmt(gc)])
    if dft.get("newcode"):
        lines.append("NEWCODE")
    lines.append(f'SYMMETRY\n{int(cell.get("symmetry", 0))}')
    kinetic = constant_cutoff(job)
    if kinetic:
        lines.extend(["CONSTANT CUTOFF", " ".join(f"{v:.10g}" for v in kinetic)])
    if int(job.get("charge", 0) or 0) != 0:
        lines.extend(["CHARGE", str(int(job["charge"]))])
    for psp in job.get("pseudopotentials") or []:
        mark = " KLEINMAN-BYLANDER" if psp.get("kleinman_bylander") else ""
        lines.append(f'*{psp["file"]}{mark}')
        lines.append(
            f'LMAX={_lmax_letter(int(psp["lmax"]))} LOC={_lmax_letter(int(psp["loc"]))}'
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    text = capnp_text(load_job(args.job))
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
