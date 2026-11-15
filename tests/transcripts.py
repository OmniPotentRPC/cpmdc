"""Match tutorial transcripts against a run of the built library.

A tutorial shows one transcript per backend, headed by the
``cpmdc/<version> available=<0|1>`` line the example prints first. The
reference well is deterministic, so its transcript must match exactly. An
OpenCPMD run interleaves SCF output and converges to a threshold, so only
its ``energy_h=`` and ``force[i] =`` lines are compared, within ``tol``.
"""

from __future__ import annotations

import re
import subprocess
import sys

NUMBER = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")


def _lines(text: str) -> list[str]:
    return [line.rstrip() for line in text.splitlines() if line.strip()]


def _results(lines: list[str]) -> list[tuple[str, list[float]]]:
    found = []
    for line in lines:
        key = line.split("=", 1)[0].strip()
        if key == "energy_h" or key.startswith("force["):
            found.append((key, [float(v) for v in NUMBER.findall(line.split("=", 1)[1])]))
    return found


def opencpmd_linked(library: str) -> bool:
    """True when ``library`` reports ``cpmdc_available() == 1``."""
    proc = subprocess.run(
        [sys.executable, "-c", "import ctypes, sys; print(ctypes.CDLL(sys.argv[1]).cpmdc_available())", library],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    return proc.returncode == 0 and proc.stdout.strip().endswith("1")


def ase_mismatch(stdout: str, transcript: str, opencpmd: bool) -> str | None:
    """Exact match on the reference well; the printed fields under OpenCPMD."""
    got = _lines(stdout)
    if not opencpmd:
        return None if got == _lines(transcript) else "ase-energy transcript differs"
    fields = {}
    for line in got:
        key, sep, value = line.partition("=")
        if sep and key in ("energy_ev", "natoms", "fmax"):
            fields[key] = value
    if fields.get("natoms") != "3" or "energy_ev" not in fields or "fmax" not in fields:
        return "ase-energy printed no energy_ev, natoms=3 and fmax"
    if float(fields["energy_ev"]) >= 0.0:
        return "OpenCPMD water energy is not negative"
    return None


def mismatch(stdout: str, transcripts: list[str], tol: float = 1.0e-5) -> str | None:
    """Return why ``stdout`` matches none of ``transcripts``, or None."""
    got = _lines(stdout)
    if not got:
        return "no output"
    for transcript in transcripts:
        expect = _lines(transcript)
        if not expect or expect[0] != got[0]:
            continue
        if got[0].endswith("available=0"):
            return None if got == expect else "reference transcript differs"
        want, have = _results(expect), _results(got)
        if [k for k, _ in want] != [k for k, _ in have]:
            return "OpenCPMD transcript lists different results"
        for (key, a), (_, b) in zip(want, have):
            if len(a) != len(b) or any(abs(x - y) > tol for x, y in zip(a, b)):
                return f"{key} differs beyond {tol:g}"
        return None
    return f"no transcript for {got[0]!r}"
