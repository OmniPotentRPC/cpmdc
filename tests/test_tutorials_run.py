#!/usr/bin/env python3
"""Run the narrative tutorials against the built library."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


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
        got = [line.rstrip() for line in ran.stdout.splitlines() if line.strip()]
        expect = [line.rstrip() for line in transcripts[0].splitlines() if line.strip()]
        if got != expect:
            return fail("first-energy transcript mismatch:\n" + ran.stdout + "\nexpected:\n" + transcripts[0])
    return 0


def quickstart(root: Path) -> int:
    text = (root / "docs/orgmode/tutorials/quickstart.org").read_text(encoding="utf-8")
    bash = [block for block in blocks(text, "bash") if "capnp encode" in block]
    if not bash:
        return fail("quickstart has no capnp encode example")
    with tempfile.TemporaryDirectory(prefix="cpmdc-quick-") as tmp:
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
            return fail("quickstart encode did not write both messages")
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
    got = [line.rstrip() for line in ase.stdout.splitlines() if line.strip()]
    expect = [line.rstrip() for line in transcripts[0].splitlines() if line.strip()]
    if got != expect:
        return fail("ase-energy transcript mismatch:\n" + ase.stdout)
    if "examples/python/ase_energy.py" not in tutorial:
        return fail("ase-energy tutorial does not run the example")
    eon = (root / "docs/orgmode/tutorials/eon-rgpot.org").read_text(encoding="utf-8")
    if "backend = cpmdc" not in eon or "5.000000 5.000000 5.117300" not in eon:
        return fail("eon tutorial lost the water geometry or the cpmdc backend")
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
        lambda: quickstart(root),
        lambda: python_examples(root, library, params, schema),
    ):
        status = step()
        if status != 0:
            return status
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
