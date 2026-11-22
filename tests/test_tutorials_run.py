#!/usr/bin/env python3
"""Run the narrative tutorials against the built library."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from transcripts import ase_mismatch, mismatch, opencpmd_linked


def blocks(text: str, lang: str) -> list[str]:
    pattern = rf"#\+begin_src {lang}\n(.*?)\n#\+end_src"
    return re.findall(pattern, text, flags=re.S | re.I)


def run(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=merged,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def first_energy(root: Path, library_dir: Path) -> int:
    text = (root / "docs/orgmode/tutorials/first-energy.org").read_text(encoding="utf-8")
    capnp_blocks = blocks(text, "capnp")
    sources = blocks(text, "c")
    transcripts = [block for block in blocks(text, "text") if "force[0]" in block]
    if len(capnp_blocks) < 2 or not sources or not transcripts:
        return fail("first-energy tutorial is missing a runnable block")
    schema = root / "schema" / "Potentials.capnp"
    with tempfile.TemporaryDirectory(prefix="cpmdc-first-") as tmp:
        work = Path(tmp)
        (work / "water.params.txt").write_text(capnp_blocks[0] + "\n", encoding="utf-8")
        (work / "water.step.txt").write_text(capnp_blocks[1] + "\n", encoding="utf-8")
        (work / "first_forces.c").write_text(sources[0] + "\n", encoding="utf-8")
        for struct, src, dst in (
            ("CPMDParams", "water.params.txt", "water.params.bin"),
            ("ForceInput", "water.step.txt", "water.step.bin"),
        ):
            proc = subprocess.run(
                ["capnp", "encode", str(schema), struct],
                input=(work / src).read_bytes(),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            if proc.returncode != 0:
                return fail(proc.stderr.decode())
            (work / dst).write_bytes(proc.stdout)
        compiled = run(
            ["mpicc", "-c", "first_forces.c", "-I", str(root / "include")],
            work,
        )
        if compiled.returncode != 0:
            return fail(compiled.stderr)
        linked = run(
            [
                "mpif90",
                "first_forces.o",
                "-L",
                str(library_dir),
                "-lcpmdc",
                "-Wl,-rpath," + str(library_dir),
                "-o",
                "first_forces",
            ],
            work,
        )
        if linked.returncode != 0:
            return fail(linked.stderr)
        ran = run(
            [str(work / "first_forces"), "water.params.bin", "water.step.bin"],
            work,
        )
        if ran.returncode != 0:
            return fail(ran.stderr + ran.stdout)
        why = mismatch(ran.stdout, transcripts)
        if why:
            return fail(f"first-energy transcript mismatch ({why}):\n" + ran.stdout)
    return 0


def build_and_test(root: Path) -> int:
    text = (root / "docs/orgmode/howto/quickstart.org").read_text(encoding="utf-8")
    bash = [block for block in blocks(text, "bash") if "capnp encode" in block]
    if not bash:
        return fail("the build-and-test how-to has no capnp encode example")
    with tempfile.TemporaryDirectory(prefix="cpmdc-build-test-") as tmp:
        script = bash[0].replace("tutorial-work", tmp)
        proc = subprocess.run(
            ["bash", "-c", script],
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            return fail(proc.stderr)
        if not Path(tmp, "cpmd_params.bin").is_file() or not Path(tmp, "force_input.bin").is_file():
            return fail("the build-and-test encode did not write both messages")
    return 0


def python_examples(root: Path, library: Path, params: Path, schema: Path) -> int:
    env = {
        "PYTHONPATH": str(root / "python"),
        "CPMDC_LIBRARY": str(library),
        "CPMDC_PARAMS": str(params),
        "CPMDC_SCHEMA": str(schema),
    }
    energy = run(
        [sys.executable, str(root / "examples/python/energy.py"), str(params), str(library)],
        root,
        env,
    )
    if energy.returncode != 0 or "energy_h=" not in energy.stdout:
        return fail(energy.stderr + energy.stdout)
    ase = run(
        [sys.executable, str(root / "examples/python/ase_energy.py")],
        root,
        env,
    )
    if ase.returncode != 0 or "energy_ev=" not in ase.stdout:
        return fail(ase.stderr + ase.stdout)
    tutorial = (root / "docs/orgmode/tutorials/ase-energy.org").read_text(encoding="utf-8")
    transcripts = [block for block in blocks(tutorial, "text") if "energy_ev=" in block]
    if not transcripts:
        return fail("ase-energy tutorial has no transcript")
    why = ase_mismatch(ase.stdout, transcripts[0], opencpmd_linked(str(library)))
    if why:
        return fail(f"ase-energy transcript mismatch ({why}):\n" + ase.stdout)
    if "examples/python/ase_energy.py" not in tutorial:
        return fail("ase-energy tutorial does not run the example")
    eon = (root / "docs/orgmode/tutorials/eon-rgpot.org").read_text(encoding="utf-8")
    if "backend = cpmdc" not in eon or "5.000000 5.000000 5.117300" not in eon:
        return fail("eon tutorial lost the water geometry or the cpmdc backend")
    return 0


def section(text: str, title: str) -> str:
    """Return the body of the Org heading ``title``, up to the next heading."""
    pattern = rf"^\*+ {re.escape(title)}\n(.*?)(?=^\*+ |\Z)"
    found = re.search(pattern, text, flags=re.S | re.M)
    return found.group(1) if found else ""


def results_dat(text: str) -> dict[str, str]:
    """Parse eOn's results.dat, whose lines read ``value key``."""
    found = {}
    for line in text.splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2:
            found[parts[1].strip()] = parts[0]
    return found


def matter_rows(text: str) -> list[list[float]]:
    """Return the iteration, step, force norm, and energy of each ``[Matter]`` row."""
    rows = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 5 and fields[0] == "[Matter]" and fields[1].isdigit():
            rows.append([float(v) for v in fields[1:]])
    return rows


def close(got: float, want: float, rel: float) -> bool:
    return abs(got - want) <= rel * max(1.0, abs(want))


def eon_minimisation(root: Path, library: Path) -> int:
    """Run the eOn tutorial and compare its numbers; skip without eonclient."""
    eonclient = os.environ.get("EONCLIENT") or shutil.which("eonclient")
    if not eonclient:
        print("eon-rgpot: skipped, no eonclient on PATH and EONCLIENT is unset")
        return 0
    live = opencpmd_linked(str(library))
    if live and not os.environ.get("CPMDC_PSEUDO_DIR"):
        print("eon-rgpot: skipped, the library links OpenCPMD and CPMDC_PSEUDO_DIR is unset")
        return 0
    page = (root / "docs/orgmode/tutorials/eon-rgpot.org").read_text(encoding="utf-8")
    first = (root / "docs/orgmode/tutorials/first-energy.org").read_text(encoding="utf-8")
    build = "OpenCPMD build" if live else "Default build"
    expected = section(page, build)
    want_console = [block for block in blocks(expected, "text") if "[Matter] Iter" in block]
    want_results = [block for block in blocks(expected, "text") if "termination_reason" in block]
    con = [block for block in blocks(page, "text") if block.startswith("Generated for the cpmdc eOn tutorial")]
    ini = blocks(page, "ini")
    if not (want_console and want_results and con and ini):
        return fail(f"eon tutorial is missing a runnable block or the {build} transcript")
    want = results_dat(want_results[0])
    with tempfile.TemporaryDirectory(prefix="cpmdc-eon-") as tmp:
        work = Path(tmp)
        (work / "pos.con").write_text(con[0] + "\n", encoding="utf-8")
        (work / "config.ini").write_text(ini[0] + "\n", encoding="utf-8")
        encoded = subprocess.run(
            ["capnp", "encode", str(root / "schema" / "Potentials.capnp"), "CPMDParams"],
            input=(blocks(first, "capnp")[0] + "\n").encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if encoded.returncode != 0:
            return fail(encoded.stderr.decode())
        (work / "water.params.bin").write_bytes(encoded.stdout)
        ran = run(
            [eonclient],
            work,
            {"CPMDC_LIBRARY": str(library), "OMP_NUM_THREADS": "1"},
        )
        (work / "run.log").write_text(ran.stdout, encoding="utf-8")
        if ran.returncode != 0:
            return fail("eonclient failed:\n" + ran.stderr + ran.stdout[-2000:])
        if not (work / "results.dat").is_file():
            return fail("eonclient wrote no results.dat:\n" + ran.stdout[-2000:])
        got = results_dat((work / "results.dat").read_text(encoding="utf-8"))
        for key in ("termination_reason", "termination_reason_text", "job_type", "potential_type"):
            if got.get(key) != want[key]:
                return fail(f"eon results.dat {key} is {got.get(key)!r}, the page shows {want[key]!r}")
        energy_tol = 1.0e-4 if live else 1.0e-6
        if not close(float(got["potential_energy"]), float(want["potential_energy"]), energy_tol):
            return fail(
                f"eon potential_energy {got['potential_energy']} is not within {energy_tol:g} of "
                f"{want['potential_energy']}"
            )
        have_rows = matter_rows(ran.stdout)
        if not have_rows:
            return fail("eonclient printed no [Matter] rows:\n" + ran.stdout[-2000:])
        if have_rows[-1][2] >= 0.01:
            return fail(f"the final force norm {have_rows[-1][2]} is not below 0.01")
        if not live:
            if got["total_force_calls"] != want["total_force_calls"]:
                return fail(f"eon made {got['total_force_calls']} force calls, the page shows {want['total_force_calls']}")
            want_rows = matter_rows(want_console[0])
            if len(have_rows) != len(want_rows):
                return fail(f"eon printed {len(have_rows)} [Matter] rows, the page shows {len(want_rows)}")
            for have, wanted in zip(have_rows, want_rows):
                if any(not close(a, b, 1.0e-4) for a, b in zip(have, wanted)):
                    return fail(f"[Matter] row {int(wanted[0])} is {have}, the page shows {wanted}")
    return 0


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "usage: test_tutorials_run.py SOURCE_ROOT LIBRARY_DIR LIBRARY PARAMS",
            file=sys.stderr,
        )
        return 2
    root = Path(sys.argv[1])
    library_dir = Path(sys.argv[2])
    library = Path(sys.argv[3])
    params = Path(sys.argv[4])
    schema = root / "schema" / "Potentials.capnp"
    for step in (
        lambda: first_energy(root, library_dir),
        lambda: build_and_test(root),
        lambda: python_examples(root, library, params, schema),
        lambda: eon_minimisation(root, library),
    ):
        status = step()
        if status != 0:
            return status
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
