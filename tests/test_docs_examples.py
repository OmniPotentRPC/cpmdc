#!/usr/bin/env python3
"""Run every documentation example the stub library can execute.

Each fenced block under docs/orgmode is either executed here or named in
SKIP with the tool it needs. A new block fails until it is one or the other.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

BLOCK = re.compile(r"#\+begin_src (\w+)\n(.*?)\n#\+end_src", re.S | re.I)

# file relative to docs/orgmode, block index -> why this process cannot run it
SKIP = {
    "contributing/index.org:0": "reconfigures the tree",
    "contributing/index.org:1": "expects a failing focused suite before the code exists",
    "contributing/index.org:2": "docs build needs the docs pixi environment",
    "contributing/index.org:3": "docs build needs the docs pixi environment",
    "contributing/index.org:4": "docs export needs emacs and ox-rst",
    "contributing/index.org:5": "link check needs lychee",
    "howto/debug-deck.org:0": "needs the caller's params.bin and step.bin",
    "howto/debug-deck.org:1": "diffs a caller-supplied reference deck",
    "howto/embedding.org:0": "fragment, not a complete program",
    "howto/embedding.org:1": "expected transcript, not a command",
    "howto/embedding.org:2": "fragment, not a complete program",
    "howto/embedding.org:3": "fragment, not a complete program",
    "howto/embedding.org:4": "fragment, not a complete program",
    "howto/embedding.org:5": "points at a caller-built library and pseudopotentials",
    "howto/extrapolate-orbitals.org:0": "section fragment, not a CPMDParams message",
    "howto/extrapolate-orbitals.org:1": "log transcript",
    "howto/install.org:0": "downloads subprojects",
    "howto/install.org:1": "clones the repository",
    "howto/install.org:2": "reconfigures the tree",
    "howto/install.org:3": "needs an OpenCPMD build",
    "howto/install.org:4": "needs an OpenCPMD build",
    "howto/install.org:5": "installs into the caller prefix",
    "howto/install.org:6": "installs from PyPI",
    "howto/install.org:7": "builds and installs a wheel into the caller environment",
    "howto/mpi.org:0": "needs a multi-rank host and pseudopotentials",
    "howto/mpi.org:1": "fragment, not a complete program",
    "howto/mpi.org:2": "fragment, not a complete program",
    "howto/mpi.org:3": "needs a multi-rank host",
    "howto/mpi.org:4": "fragment, not a complete program",
    "howto/opencpmd-archive.org:0": "clones OpenCPMD",
    "howto/opencpmd-archive.org:1": "configures an OpenCPMD source tree",
    "howto/opencpmd-archive.org:2": "patches an OpenCPMD source tree",
    "howto/opencpmd-archive.org:3": "needs an OpenCPMD source tree",
    "howto/opencpmd-archive.org:4": "needs an OpenCPMD build",
    "howto/opencpmd-archive.org:5": "needs an OpenCPMD build",
    "howto/opencpmd-archive.org:6": "configures an OpenCPMD source tree",
    "howto/pwscale-ranks.org:0": "benchmark, not a correctness example",
    "howto/pwscale-ranks.org:1": "fragment, not a complete program",
    "howto/socket-drivers.org:0": "needs a listening i-PI or ASE server",
    "howto/socket-drivers.org:2": "i-PI input fragment",
    "howto/socket-drivers.org:3": "needs i-PI",
    "howto/socket-drivers.org:4": "i-PI input fragment",
    "howto/socket-drivers.org:5": "LAMMPS input; untested",
    "howto/socket-drivers.org:6": "needs i-PI and LAMMPS",
    "howto/troubleshooting.org:0": "reads a caller build directory",
    "howto/wavefunction-optimiser.org:0": "section fragment, not a CPMDParams message",
    "howto/wavefunction-optimiser.org:1": "section fragment, not a CPMDParams message",
    "howto/wavefunction-optimiser.org:2": "patches an OpenCPMD source tree",
    "howto/write-cpmdparams.org:1": "reads cluster.params.txt from the caller's directory",
    "howto/write-cpmdparams.org:2": "needs the caller's step.bin",
    "index.org:0": "reconfigures the tree",
    "reference/cpmd-options.org:0": "re-runs the inventory suite",
    "reference/feature-inventory.org:0": "fragment, not a complete program",
    "reference/feature-inventory.org:1": "re-runs the inventory suite",
    "reference/feature-inventory.org:2": "needs an OpenCPMD source tree",
    "tutorials/ase-energy.org:1": "expected transcript",
    "tutorials/eon-rgpot.org:0": "prepares a directory for eonclient",
    "tutorials/eon-rgpot.org:1": "eOn structure file",
    "tutorials/eon-rgpot.org:2": "eOn config.ini",
    "tutorials/eon-rgpot.org:3": "needs libcpmdc from an OpenCPMD build and pseudopotentials",
    "tutorials/eon-rgpot.org:4": "reads the caller's water.params.txt",
    "tutorials/eon-rgpot.org:5": "needs eonclient",
    "tutorials/first-energy.org:0": "clones the repository",
    "tutorials/first-energy.org:1": "opens a pixi shell outside the checkout",
    "tutorials/first-energy.org:5": "runs example_host_step; the compile step below is the lesson",
    "tutorials/first-energy.org:6": "expected transcript",
    "tutorials/first-energy.org:9": "expected transcript",
    "tutorials/first-energy.org:10": "needs an OpenCPMD build",
    "tutorials/first-energy.org:11": "expected transcript from an OpenCPMD run",
    "tutorials/quickstart.org:0": "pixi run test-stub would re-enter this suite",
    "tutorials/quickstart.org:1": "reconfigures the tree",
    "tutorials/quickstart.org:2": "example-host-step is its own meson test",
    "tutorials/quickstart.org:3": "expected transcript",
    "tutorials/quickstart.org:4": "re-runs the e2e suite",
    "tutorials/quickstart.org:5": "needs an OpenCPMD build and pseudopotentials",
}

# The bash block encodes the messages that the Python blocks after it read.
PYTHON_HOWTO = [f"howto/python.org:{index}" for index in range(6)]


def blocks(root: Path) -> dict[str, tuple[str, str]]:
    found: dict[str, tuple[str, str]] = {}
    base = root / "docs" / "orgmode"
    for path in sorted(base.rglob("*.org")):
        rel = path.relative_to(base).as_posix()
        for index, (lang, body) in enumerate(BLOCK.findall(path.read_text(encoding="utf-8"))):
            found[f"{rel}:{index}"] = (lang, body)
    return found


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def encode(schema: Path, struct: str, text: str) -> None:
    proc = subprocess.run(
        ["capnp", "encode", str(schema), struct],
        input=text.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0 or not proc.stdout:
        raise RuntimeError(proc.stderr.decode() or "capnp encode failed")


def run_python(root: Path, library: Path, params: Path, schema: Path, body: str) -> None:
    work = Path(os.environ["CPMDC_EXAMPLE_WORK"])
    (work / "water.params.bin").write_bytes(params.read_bytes())
    build = work / "build"
    build.mkdir(exist_ok=True)
    target = build / "libcpmdc.so"
    if not target.exists():
        target.symlink_to(library)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "python")
    env["CPMDC_LIBRARY"] = str(library)
    env["CPMDC_SCHEMA"] = str(schema)
    proc = subprocess.run(
        [sys.executable, "-c", body],
        cwd=work,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr + proc.stdout)


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "usage: test_docs_examples.py SOURCE_ROOT BUILD LIBRARY PARAMS FD_CHECK",
            file=sys.stderr,
        )
        return 2
    root, build, library, params, fd_check = (Path(arg).resolve() for arg in sys.argv[1:])
    schema = root / "schema" / "Potentials.capnp"
    found = blocks(root)
    run_keys = {
        "howto/check-forces.org:0",
        *PYTHON_HOWTO,
        "howto/socket-drivers.org:1",
        "howto/write-cpmdparams.org:0",
        "reference/cpmd-options.org:1",
        "tutorials/ase-energy.org:0",
        "tutorials/first-energy.org:2",
        "tutorials/first-energy.org:3",
        "tutorials/first-energy.org:4",
        "tutorials/first-energy.org:7",
        "tutorials/first-energy.org:8",
        "tutorials/quickstart.org:6",
    }
    missing = sorted(set(found) - set(SKIP) - run_keys)
    extra = sorted((set(SKIP) | run_keys) - set(found))
    if missing or extra:
        return fail("unclassified docs blocks:\n" + "\n".join(missing + extra))

    import tempfile

    with tempfile.TemporaryDirectory(prefix="cpmdc-docs-") as tmp:
        os.environ["CPMDC_EXAMPLE_WORK"] = tmp
        work = Path(tmp)
        try:
            encoded = subprocess.run(
                ["bash", "-c", found[PYTHON_HOWTO[0]][1]],
                cwd=work,
                env={
                    **os.environ,
                    "CPMDC": str(root),
                    "PYTHONPATH": str(root / "python"),
                    "CPMDC_SCHEMA": str(schema),
                    "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"],
                },
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if encoded.returncode != 0:
                return fail(encoded.stderr)
            for key in PYTHON_HOWTO[1:]:
                run_python(root, library, params, schema, found[key][1])
            run_python(root, library, params, schema, found["howto/socket-drivers.org:1"][1])
            encode(schema, "CPMDParams", found["howto/write-cpmdparams.org:0"][1])
            encode(schema, "CPMDParams", found["reference/cpmd-options.org:1"][1])
            encode(schema, "CPMDParams", found["tutorials/first-energy.org:2"][1])
            encode(schema, "ForceInput", found["tutorials/first-energy.org:3"][1])
            (work / "water.params.txt").write_text(
                found["tutorials/first-energy.org:2"][1] + "\n", encoding="utf-8"
            )
            (work / "water.step.txt").write_text(
                found["tutorials/first-energy.org:3"][1] + "\n", encoding="utf-8"
            )
            env = os.environ.copy()
            env["CPMDC"] = str(root)
            env["PYTHONPATH"] = str(root / "python")
            env["CPMDC_LIBRARY"] = str(library)
            env["CPMDC_PARAMS"] = str(params)
            env["CPMDC_SCHEMA"] = str(schema)
            encoded = subprocess.run(
                ["bash", "-c", found["tutorials/first-energy.org:4"][1]],
                cwd=work,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if encoded.returncode != 0:
                return fail(encoded.stderr)
            quick = found["tutorials/quickstart.org:6"][1].replace("tutorial-work", tmp)
            proc = subprocess.run(
                ["bash", "-c", quick],
                cwd=root,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                return fail(proc.stderr)
            proc = subprocess.run(
                [str(fd_check), str(params), str(build / "force_input_water.bin")],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                return fail(proc.stdout + proc.stderr)
            ase = subprocess.run(
                ["bash", "-c", found["tutorials/ase-energy.org:0"][1]],
                cwd=root,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if ase.returncode != 0 or "energy_ev=" not in ase.stdout:
                return fail(ase.stderr + ase.stdout)
            (work / "first_forces.c").write_text(
                found["tutorials/first-energy.org:7"][1] + "\n", encoding="utf-8"
            )
            ran = subprocess.run(
                ["bash", "-c", found["tutorials/first-energy.org:8"][1]],
                cwd=work,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            if ran.returncode != 0:
                return fail(ran.stderr + ran.stdout)
            expect = [
                line.rstrip()
                for line in found["tutorials/first-energy.org:9"][1].splitlines()
                if line.strip()
            ]
            got = [line.rstrip() for line in ran.stdout.splitlines() if line.strip()]
            if got != expect:
                return fail("first-energy example mismatch:\n" + ran.stdout)
        except RuntimeError as exc:
            return fail(str(exc))
    print(f"ran {len(run_keys)} docs examples, skipped {len(SKIP)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
