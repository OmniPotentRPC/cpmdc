#!/usr/bin/env python3
"""Run every documentation example the stub library can execute.

Each fenced block in README.md and under docs/orgmode is classified once:

- run here against the stub build,
- compiled here against the public headers (C fragments),
- named in OPENCPMD, because it needs an OpenCPMD build, an OpenCPMD source
  tree, or its pseudopotentials; the text around it must say so,
- or named in SKIP with the other tool or input it needs.

A new block fails until it is one of these. Every shell block, run or not, is
also checked against the tree: meson -D options (except in the blocks of
OTHER_MESON_PROJECT, which configure another project), pixi tasks and environments,
meson test names and suites, repository paths, build outputs, and i-PI
species strings must exist.
"""

from __future__ import annotations

import configparser
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import tomllib
from pathlib import Path

from transcripts import ase_mismatch, mismatch, opencpmd_linked

ORG_BLOCK = re.compile(
    r"^[ \t]*#\+begin_src (\w+)[^\n]*\n(.*?)\n[ \t]*#\+end_src", re.S | re.M | re.I
)
ORG_HEADING = re.compile(r"^\*+ ", re.M)
MD_FENCE = re.compile(r"^```(\w*)[^\n]*\n(.*?)\n```[ \t]*$", re.S | re.M)
MD_HEADING = re.compile(r"^#+ ", re.M)
SHELL = {"bash", "sh"}
OPENCPMD_WORDS = re.compile(
    r"OpenCPMD|\bCPMD\b|available=1|build-cpmd|with_cpmd|CPMD_ROOT|PP_LIBRARY|pseudopotential"
)
ORG_TITLE = re.compile(r"^#\+TITLE:.*$", re.M)

# key -> what of OpenCPMD the block needs
OPENCPMD = {
    "README.md:9": "an OpenCPMD build",
    "contributing/index.org:1": "an OpenCPMD source tree",
    "howto/embedding.org:4": "libcpmdc from an OpenCPMD build and pseudopotentials",
    "explanation/md-engines.org:0": "an OpenCPMD build and an MPI GROMACS built against it",
    "explanation/md-engines.org:1": "an OpenCPMD build and an MPI GROMACS built against it",
    "howto/install.org:3": "an OpenCPMD build",
    "howto/install.org:4": "an OpenCPMD build",
    "howto/mpi.org:0": "a multi-rank host on an OpenCPMD build and pseudopotentials",
    "howto/opencpmd-archive.org:0": "clones OpenCPMD",
    "howto/opencpmd-archive.org:1": "an OpenCPMD source tree",
    "howto/opencpmd-archive.org:2": "an OpenCPMD source tree",
    "howto/opencpmd-archive.org:3": "an OpenCPMD source tree",
    "howto/opencpmd-archive.org:4": "an OpenCPMD build",
    "howto/opencpmd-archive.org:5": "an OpenCPMD build",
    "howto/opencpmd-archive.org:6": "an OpenCPMD source tree",
    "howto/wavefunction-optimiser.org:2": "an OpenCPMD source tree",
    "reference/feature-inventory.org:2": "an OpenCPMD source tree",
    "tutorials/eon-rgpot.org:7": "libcpmdc from an OpenCPMD build, pseudopotentials, and eonclient",
    "tutorials/eon-rgpot.org:8": "transcript of an OpenCPMD run with eonclient",
    "tutorials/eon-rgpot.org:10": "libcpmdc from an OpenCPMD build, pseudopotentials, and eonclient",
    "tutorials/first-energy.org:10": "an OpenCPMD build",
    "tutorials/first-energy.org:11": "the OpenCPMD pseudopotentials",
    "tutorials/first-energy.org:12": "an OpenCPMD build and pseudopotentials",
    "tutorials/first-energy.org:13": "transcript of an OpenCPMD run",
}

# key -> why this process cannot run it
SKIP = {
    "README.md:0": "installs from PyPI",
    "README.md:1": "reconfigures the tree",
    "README.md:2": "pixi run test-stub would re-enter this suite",
    "README.md:4": "expected output",
    "README.md:5": "example-host-step is its own meson test",
    "README.md:7": "needs -Dwith_readcon=true; pixi run test-readcon runs restart-con-golden",
    "README.md:8": "re-runs the inventory suite",
    "README.md:10": "docs build needs the docs pixi environment",
    "README.md:11": "docs export needs emacs and ox-rst",
    "contributing/index.org:0": "reconfigures the tree",
    "contributing/index.org:2": "docs build needs the docs pixi environment",
    "contributing/index.org:3": "docs build needs the docs pixi environment",
    "contributing/index.org:4": "docs export needs emacs and ox-rst",
    "contributing/index.org:5": "link check needs lychee",
    "tutorials/quickstart.org:0": "pixi run test-stub would re-enter this suite",
    "tutorials/quickstart.org:1": "reconfigures the tree",
    "tutorials/quickstart.org:2": "example-host-step is its own meson test",
    "tutorials/quickstart.org:3": "expected transcript",
    "explanation/md-engines.org:2": "benchmark, not a correctness example",
    "explanation/md-engines.org:3": "benchmark, not a correctness example",
    "howto/debug-deck.org:1": "diffs a caller-supplied reference deck",
    "howto/embedding.org:1": "expected transcript, not a command",
    "howto/extrapolate-orbitals.org:0": "section fragment, not a CPMDParams message",
    "howto/extrapolate-orbitals.org:1": "log transcript",
    "howto/install.org:0": "downloads subprojects",
    "howto/install.org:1": "clones the repository",
    "howto/install.org:2": "reconfigures the tree",
    "howto/install.org:5": "installs into the caller prefix",
    "howto/install.org:6": "installs from PyPI",
    "howto/install.org:7": "builds and installs a wheel into the caller environment",
    "howto/install.org:8": "docs build needs the docs pixi environment",
    "howto/gromacs.org:0": "reconfigures the tree and builds GROMACS",
    "howto/gromacs.org:1": "installs into build/prefix and builds GROMACS; CI's gromacs-qmmm (installed) runs it",
    "howto/gromacs.org:2": "GROMACS mdp fragment",
    "howto/gromacs.org:4": "needs gmx_d built against libcpmdc_gromacs; tests/python/test_gromacs.py runs these steps",
    "howto/gromacs.org:5": "GROMACS mdp fragment",
    "howto/lammps.org:0": "needs the lammps Python module; tests/python/test_lammps.py runs these steps",
    "howto/lammps.org:1": "reconfigures the tree",
    "howto/lammps.org:2": "LAMMPS input script for lmp and the plugin",
    "howto/mpi.org:3": "needs the caller's MPI host on 16 ranks",
    "howto/pwscale-ranks.org:0": "benchmark, not a correctness example",
    "howto/pwscale-ranks.org:1": "needs -Dwith_pwscale=true; pwscale-agreement makes the same calls",
    "howto/socket-drivers.org:0": "needs a listening i-PI or ASE server",
    "howto/socket-drivers.org:2": "i-PI input fragment",
    "howto/socket-drivers.org:3": "needs i-PI",
    "howto/socket-drivers.org:4": "i-PI input fragment",
    "howto/socket-drivers.org:5": "LAMMPS input; untested",
    "howto/socket-drivers.org:6": "needs i-PI and LAMMPS",
    "howto/socket-drivers.org:7": "needs a listening i-PI or ASE server",
    "howto/wavefunction-optimiser.org:0": "section fragment, not a CPMDParams message",
    "howto/wavefunction-optimiser.org:1": "section fragment, not a CPMDParams message",
    "index.org:0": "reconfigures the tree",
    "reference/cpmd-options.org:0": "re-runs the inventory suite",
    "reference/feature-inventory.org:1": "re-runs the inventory suite",
    "tutorials/ase-energy.org:1": "expected transcript",
    "tutorials/eon-rgpot.org:0": "installs Eigen with pixi global",
    "tutorials/eon-rgpot.org:1": "sets the environment of the eonclient build",
    "tutorials/eon-rgpot.org:2": "clones and builds rgpot",
    "tutorials/eon-rgpot.org:3": "clones and builds eOn",
    "tutorials/eon-rgpot.org:11": "error transcript, not a command",
    "tutorials/first-energy.org:0": "clones the repository",
    "tutorials/first-energy.org:1": "opens a pixi shell outside the checkout",
    "tutorials/first-energy.org:6": "expected transcript",
    "tutorials/first-energy.org:9": "expected transcript",
}

# key -> the other project whose meson setup the block runs. Its -D options
# are that project's, so they are not checked against meson_options.txt.
OTHER_MESON_PROJECT = {
    "tutorials/eon-rgpot.org:2": "rgpot",
    "tutorials/eon-rgpot.org:3": "eOn",
}

# C fragments: hoisted includes, then the body inside a function after PRELUDE.
COMPILE = {
    "README.md:3",
    "howto/embedding.org:0",
    "howto/embedding.org:2",
    "howto/embedding.org:3",
    "howto/mpi.org:1",
    "howto/mpi.org:2",
    "howto/mpi.org:4",
    "reference/feature-inventory.org:0",
}

# The names the fragments take from the host program around them.
PRELUDE = """\
#include <cpmdc.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N_ATOMS 3
extern const void *params_bytes;
extern size_t params_size;
extern const void *force_input_bytes;
extern size_t force_input_size;
extern int n_atoms;
extern const double *positions_ang;
extern const int *atomic_numbers;
extern double *gradient_h_bohr;
extern double *forces_h_bohr;
extern CPMDCSession *session;
extern int world_rank;
extern CPMDCResult result;
extern double *forces;
"""

# The bash block of each page encodes the messages that its Python blocks read.
PYTHON_HOWTO = [f"howto/python.org:{index}" for index in range(3)]
ASE_HOWTO = [f"howto/ase.org:{index}" for index in range(7)]

RUN = {
    "README.md:6",
    "tutorials/quickstart.org:4",
    "howto/check-forces.org:0",
    "howto/debug-deck.org:0",
    "howto/gromacs.org:3",
    *PYTHON_HOWTO,
    *ASE_HOWTO,
    "howto/socket-drivers.org:1",
    "howto/troubleshooting.org:0",
    "howto/write-cpmdparams.org:0",
    "howto/write-cpmdparams.org:1",
    "howto/write-cpmdparams.org:2",
    "reference/cpmd-options.org:1",
    "tutorials/ase-energy.org:0",
    "tutorials/eon-rgpot.org:4",
    "tutorials/eon-rgpot.org:5",
    "tutorials/eon-rgpot.org:6",
    "tutorials/eon-rgpot.org:9",
    "tutorials/first-energy.org:2",
    "tutorials/first-energy.org:3",
    "tutorials/first-energy.org:4",
    "tutorials/first-energy.org:5",
    "tutorials/first-energy.org:7",
    "tutorials/first-energy.org:8",
}

# meson test options that take a value
TEST_VALUE_OPTIONS = {
    "-C",
    "--suite",
    "--no-suite",
    "--setup",
    "--num-processes",
    "--repeat",
    "-t",
    "--timeout-multiplier",
    "--wrapper",
    "--wrap",
    "--test-args",
}
REPO_DIRS = ("tools", "scripts", "bench", "examples", "tests", "schema", "include", "python", "subprojects", "docs")


class Block:
    def __init__(self, lang: str, body: str, context: str) -> None:
        self.lang = lang
        self.body = body
        self.context = context


def _scan(text: str, pattern: re.Pattern[str], heading: re.Pattern[str]) -> list[Block]:
    found = []
    previous = 0
    title = ORG_TITLE.search(text)
    for match in pattern.finditer(text):
        before = text[previous : match.start()]
        headings = list(heading.finditer(before))
        if headings:
            before = before[headings[-1].start() :]
        if title:
            before = title.group(0) + "\n" + before
        found.append(Block(match.group(1), textwrap.dedent(match.group(2)), before))
        previous = match.end()
    return found


def blocks(root: Path) -> dict[str, Block]:
    found: dict[str, Block] = {}
    base = root / "docs" / "orgmode"
    for path in sorted(base.rglob("*.org")):
        rel = path.relative_to(base).as_posix()
        for index, block in enumerate(_scan(path.read_text(encoding="utf-8"), ORG_BLOCK, ORG_HEADING)):
            found[f"{rel}:{index}"] = block
    readme = (root / "README.md").read_text(encoding="utf-8")
    for index, block in enumerate(_scan(readme, MD_FENCE, MD_HEADING)):
        found[f"README.md:{index}"] = block
    return found


def fail(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


class Tree:
    """What the shell blocks may name: options, tasks, tests, paths."""

    def __init__(self, root: Path, build: Path) -> None:
        self.root = root
        self.build = build
        info = build / "meson-info"
        options = json.loads((info / "intro-buildoptions.json").read_text(encoding="utf-8"))
        self.options = {option["name"] for option in options}
        tests = json.loads((info / "intro-tests.json").read_text(encoding="utf-8"))
        self.tests = {test["name"] for test in tests}
        self.suites = {suite.split(":", 1)[-1] for test in tests for suite in test["suite"]}
        pixi = tomllib.loads((root / "pixi.toml").read_text(encoding="utf-8"))
        self.default_tasks = set(pixi.get("tasks", {}))
        self.env_tasks: dict[str, set[str]] = {}
        features = pixi.get("feature", {})
        for env, spec in pixi.get("environments", {}).items():
            names = spec.get("features", []) if isinstance(spec, dict) else spec
            tasks = set(self.default_tasks)
            for feature in names:
                tasks |= set(features.get(feature, {}).get("tasks", {}))
            self.env_tasks[env] = tasks
        sys.path.insert(0, str(root / "python"))
        from cpmdc.ipi import species_numbers

        self.species_numbers = species_numbers

    def path_problem(self, token: str, loops: dict[str, list[str]], made: list[str]) -> str | None:
        for prefix in ('"$CPMDC"/', "$CPMDC/", "/path/to/cpmdc/"):
            if token.startswith(prefix):
                token = token[len(prefix) :]
                break
        else:
            if token.startswith("/") or token.startswith("$"):
                return None
        candidates = [token]
        for var, words in loops.items():
            if f"${var}" in token:
                candidates = [c.replace(f"${var}", word) for c in candidates for word in words]
        for candidate in candidates:
            if "$" in candidate or "*" in candidate or (candidate + "/").startswith(tuple(made)):
                continue
            head = candidate.split("/", 1)[0]
            if head == "build" and "/" in candidate:
                if not (self.build / candidate.split("/", 1)[1]).exists():
                    return f"build output {candidate} does not exist"
            elif head in REPO_DIRS and "/" in candidate:
                if not (self.root / candidate).exists():
                    return f"path {candidate} does not exist"
        return None

    def problems(self, body: str, own_meson: bool = True) -> list[str]:
        found = []
        text = body.replace("\\\n", " ")
        loops = {
            var: words.split()
            for var, words in re.findall(r"for (\w+) in ([^;]*?); do", text)
        }
        made = [path.rstrip("/") + "/" for path in re.findall(r"--prefix[= ]\"?(?:\$PWD/)?([^\s\"]+)", text)]
        for line in text.splitlines():
            try:
                tokens = shlex.split(line, comments=True)
            except ValueError:
                continue
            while tokens and re.match(r"^[A-Za-z_]\w*=", tokens[0]):
                tokens = tokens[1:]
            if tokens[:1] == ["export"]:
                tokens = [t.split("=", 1)[1] for t in tokens[1:] if "=" in t]
            for index, token in enumerate(tokens):
                problem = self.path_problem(token, loops, made)
                if problem:
                    found.append(problem)
                if token == "--species" and index + 1 < len(tokens):
                    try:
                        self.species_numbers(tokens[index + 1])
                    except ValueError as exc:
                        found.append(f"--species {tokens[index + 1]!r}: {exc}")
            for start in [i for i, t in enumerate(tokens) if t == "meson"]:
                if own_meson:
                    found += self.meson_problems(tokens[start + 1 :])
            for start in [i for i, t in enumerate(tokens) if t == "pixi"]:
                found += self.pixi_problems(tokens[start + 1 :])
        return found

    def meson_problems(self, args: list[str]) -> list[str]:
        if not args:
            return []
        found = []
        args = args[: args.index("&&")] if "&&" in args else args
        if args[0] == "setup":
            for arg in args[1:]:
                if arg.startswith("-D"):
                    name = arg[2:].split("=", 1)[0]
                    if name not in self.options:
                        found.append(f"meson option {name} does not exist")
        if args[0] == "test":
            rest = iter(args[1:])
            for arg in rest:
                if arg in TEST_VALUE_OPTIONS:
                    value = next(rest, "")
                    if arg == "--suite" and value not in self.suites:
                        found.append(f"meson test suite {value} does not exist")
                elif arg.startswith("-"):
                    continue
                elif arg not in self.tests:
                    found.append(f"meson test {arg} does not exist")
        return found

    def pixi_problems(self, args: list[str]) -> list[str]:
        if not args or args[0] != "run":
            return []
        env = None
        rest = iter(args[1:])
        for arg in rest:
            if arg in ("-e", "--environment"):
                env = next(rest, "")
            elif arg == "--manifest-path":
                next(rest, "")
            else:
                if env is not None and env not in self.env_tasks:
                    return [f"pixi environment {env} does not exist"]
                tasks = self.env_tasks[env] if env else self.default_tasks
                if arg not in tasks:
                    return [f"pixi task {arg} does not exist" + (f" in -e {env}" if env else "")]
                return []
        return []


def encode(schema: Path, struct: str, text: str) -> bytes:
    proc = subprocess.run(
        ["capnp", "encode", str(schema), struct],
        input=text.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0 or not proc.stdout:
        raise RuntimeError(proc.stderr.decode() or "capnp encode failed")
    return proc.stdout


def shell(key: str, body: str, cwd: Path, env: dict[str, str]) -> str:
    proc = subprocess.run(
        ["bash", "-e", "-c", body],
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{key} failed:\n{proc.stderr}{proc.stdout}")
    return proc.stdout


def at_build(body: str, build: Path) -> str:
    """Point a block's ``build/`` and ``$CPMDC/build`` paths at the meson build."""
    body = body.replace('"$CPMDC/build"', f'"{build}"').replace("$CPMDC/build/", f"{build}/")
    return re.sub(r"(?<![\w/.$-])build/", f"{build}/", body)


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


def compile_fragment(key: str, body: str, root: Path, work: Path) -> None:
    if re.search(r"^\s*int main\(", body, re.M):
        source = body + "\n"
    else:
        includes = [line for line in body.splitlines() if line.startswith("#include")]
        rest = [line for line in body.splitlines() if not line.startswith("#include")]
        source = (
            PRELUDE
            + "\n".join(includes)
            + "\nvoid docs_fragment(void);\nvoid docs_fragment(void) {\n"
            + "\n".join(rest)
            + "\n}\n"
        )
    path = work / (re.sub(r"\W", "_", key) + ".c")
    path.write_text(source, encoding="utf-8")
    compiler = shutil.which("mpicc") or os.environ.get("CC", "cc")
    proc = subprocess.run(
        [
            compiler,
            "-std=c11",
            "-fsyntax-only",
            "-Werror=implicit-function-declaration",
            "-Werror=incompatible-pointer-types",
            "-Werror=int-conversion",
            "-I",
            str(root / "include"),
            str(path),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{key} does not compile against include/:\n{proc.stderr}")


def meson_on_path(work: Path, env: dict[str, str]) -> None:
    if shutil.which("meson", path=env["PATH"]):
        return
    shim = work / "bin" / "meson"
    shim.parent.mkdir(exist_ok=True)
    shim.write_text(f'#!/bin/sh\nexec "{sys.executable}" -m mesonbuild.mesonmain "$@"\n', encoding="utf-8")
    shim.chmod(0o755)
    env["PATH"] = str(shim.parent) + os.pathsep + env["PATH"]


def classify(found: dict[str, Block]) -> list[str]:
    errors = []
    buckets = [set(OPENCPMD), set(SKIP), COMPILE, RUN]
    known = set().union(*buckets)
    for key in sorted(set(found) - known):
        errors.append(f"unclassified docs block: {key}")
    for key in sorted((known | set(OTHER_MESON_PROJECT)) - set(found)):
        errors.append(f"classified block does not exist: {key}")
    for index, bucket in enumerate(buckets):
        for other in buckets[index + 1 :]:
            for key in sorted(bucket & other):
                errors.append(f"block classified twice: {key}")
    for key in sorted(set(OPENCPMD) & set(found)):
        block = found[key]
        if not OPENCPMD_WORDS.search(block.context + block.body):
            errors.append(f"{key} needs {OPENCPMD[key]}, but its text does not say OpenCPMD")
    return errors


def main() -> int:
    if len(sys.argv) != 7:
        print(
            "usage: test_docs_examples.py SOURCE_ROOT BUILD LIBRARY PARAMS FD_CHECK RESTART_SAMPLE",
            file=sys.stderr,
        )
        return 2
    root, build, library, params, fd_check, restart_sample = (
        Path(arg).resolve() for arg in sys.argv[1:]
    )
    schema = root / "schema" / "Potentials.capnp"
    found = blocks(root)
    errors = classify(found)
    tree = Tree(root, build)
    for key, block in sorted(found.items()):
        if block.lang in SHELL:
            own = key not in OTHER_MESON_PROJECT
            errors += [f"{key}: {problem}" for problem in tree.problems(block.body, own)]
    if errors:
        return fail("\n".join(errors))

    stub = not opencpmd_linked(str(library))
    with tempfile.TemporaryDirectory(prefix="cpmdc-docs-") as tmp:
        os.environ["CPMDC_EXAMPLE_WORK"] = tmp
        work = Path(tmp)
        env = os.environ.copy()
        env["CPMDC"] = str(root)
        env["PYTHONPATH"] = str(root / "python")
        env["CPMDC_LIBRARY"] = str(library)
        env["CPMDC_PARAMS"] = str(params)
        env["CPMDC_SCHEMA"] = str(schema)
        env["PATH"] = os.pathsep.join([str(build), str(Path(sys.executable).parent), env["PATH"]])
        try:
            for key in (PYTHON_HOWTO[0], ASE_HOWTO[0]):
                shell(key, found[key].body, work, env)
            for key in PYTHON_HOWTO[1:] + ASE_HOWTO[1:]:
                run_python(root, library, params, schema, found[key].body)
            run_python(root, library, params, schema, found["howto/socket-drivers.org:1"].body)
            run_python(root, library, params, schema, found["howto/gromacs.org:3"].body)
            if "&CPMDC" not in (work / "cpmdc.inp").read_text(encoding="utf-8"):
                return fail("howto/gromacs.org:3 wrote no &CPMDC section")
            encode(schema, "CPMDParams", found["reference/cpmd-options.org:1"].body)

            cluster = work / "write-cpmdparams"
            cluster.mkdir()
            (cluster / "cluster.params.txt").write_text(
                found["howto/write-cpmdparams.org:0"].body + "\n", encoding="utf-8"
            )
            (cluster / "schema").symlink_to(root / "schema")
            (cluster / "step.bin").write_bytes((build / "force_input_water.bin").read_bytes())
            shell("howto/write-cpmdparams.org:1", found["howto/write-cpmdparams.org:1"].body, cluster, env)
            shell(
                "howto/write-cpmdparams.org:2",
                at_build(found["howto/write-cpmdparams.org:2"].body, build),
                cluster,
                env,
            )
            if "&CPMD" not in (cluster / "cluster.inp").read_text(encoding="utf-8"):
                return fail("howto/write-cpmdparams.org:2 rendered no &CPMD section")

            deck = work / "debug-deck"
            deck.mkdir()
            (deck / "params.bin").write_bytes(params.read_bytes())
            (deck / "step.bin").write_bytes((build / "force_input_water.bin").read_bytes())
            out = shell("howto/debug-deck.org:0", at_build(found["howto/debug-deck.org:0"].body, build), deck, env)
            if "&CPMD" not in out:
                return fail("howto/debug-deck.org:0 printed no &CPMD section:\n" + out)

            checkout = work / "checkout"
            checkout.mkdir()
            (checkout / "build").symlink_to(build)
            meson_on_path(work, env)
            out = shell("howto/troubleshooting.org:0", found["howto/troubleshooting.org:0"].body, checkout, env)
            version = (root / "version.txt").read_text(encoding="utf-8").strip()
            if f"cpmdc {version}" not in out or "with_cpmd" not in out:
                return fail("howto/troubleshooting.org:0 printed no version and options:\n" + out)

            tutorial = work / "cpmdc-tutorial"
            tutorial.mkdir()
            encode(schema, "CPMDParams", found["tutorials/first-energy.org:2"].body)
            encode(schema, "ForceInput", found["tutorials/first-energy.org:3"].body)
            (tutorial / "water.params.txt").write_text(
                found["tutorials/first-energy.org:2"].body + "\n", encoding="utf-8"
            )
            (tutorial / "water.step.txt").write_text(
                found["tutorials/first-energy.org:3"].body + "\n", encoding="utf-8"
            )
            shell("tutorials/first-energy.org:4", found["tutorials/first-energy.org:4"].body, tutorial, env)
            out = shell(
                "tutorials/first-energy.org:5",
                at_build(found["tutorials/first-energy.org:5"].body, build),
                tutorial,
                env,
            )
            expected = found["tutorials/first-energy.org:6"].body
            if stub and out.split() != expected.split():
                return fail(f"tutorials/first-energy.org:5 printed\n{out}\nnot\n{expected}")
            (tutorial / "first_forces.c").write_text(
                found["tutorials/first-energy.org:7"].body + "\n", encoding="utf-8"
            )
            out = shell(
                "tutorials/first-energy.org:8",
                at_build(found["tutorials/first-energy.org:8"].body, build),
                tutorial,
                env,
            )
            why = mismatch(
                out,
                [found["tutorials/first-energy.org:9"].body, found["tutorials/first-energy.org:13"].body],
            )
            if why:
                return fail(f"first-energy example mismatch ({why}):\n" + out)

            shell(
                "tutorials/eon-rgpot.org:4",
                found["tutorials/eon-rgpot.org:4"].body,
                work,
                env,
            )
            eon = work / "eon-water"
            (eon / "pos.con").write_text(found["tutorials/eon-rgpot.org:5"].body + "\n", encoding="utf-8")
            (eon / "config.ini").write_text(found["tutorials/eon-rgpot.org:6"].body + "\n", encoding="utf-8")
            from ase.io import read

            symbols = read(eon / "pos.con", format="eon").get_chemical_symbols()
            if symbols != ["O", "H", "H"]:
                return fail(f"tutorials/eon-rgpot.org:5 reads as {symbols}")
            config = configparser.ConfigParser()
            config.read(eon / "config.ini")
            if not (eon / config["RgpotPot"]["params_path"]).is_file():
                return fail("tutorials/eon-rgpot.org:6 params_path names no file in eon-water")
            shell("tutorials/eon-rgpot.org:9", found["tutorials/eon-rgpot.org:9"].body, eon, env)

            quick = found["tutorials/quickstart.org:4"].body.replace("tutorial-work", tmp)
            shell("tutorials/quickstart.org:4", quick, root, env)
            shell("howto/check-forces.org:0", at_build(found["howto/check-forces.org:0"].body, build), root, env)
            ase_out = shell(
                "tutorials/ase-energy.org:0",
                at_build(found["tutorials/ase-energy.org:0"].body, build),
                root,
                env,
            )
            why = ase_mismatch(ase_out, found["tutorials/ase-energy.org:1"].body, not stub)
            if why:
                return fail(f"tutorials/ase-energy.org:0 ({why}):\n" + ase_out)

            restart = work / "restart"
            restart.mkdir()
            sample = subprocess.run(
                [str(restart_sample), "--write-sample", str(restart / "RESTART.1")],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
            if sample.returncode != 0:
                return fail(sample.stdout)
            shell("README.md:6", found["README.md:6"].body, restart, env)
            if (restart / "RESTART.1.next").read_bytes() != (restart / "RESTART.1").read_bytes():
                return fail("README.md:6 changed the positions it read back")

            for key in sorted(COMPILE):
                compile_fragment(key, found[key].body, root, work)
        except RuntimeError as exc:
            return fail(str(exc))
    print(
        f"ran {len(RUN)} docs examples, compiled {len(COMPILE)} C fragments, "
        f"{len(OPENCPMD)} need OpenCPMD, skipped {len(SKIP)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
