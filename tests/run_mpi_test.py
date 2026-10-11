#!/usr/bin/env python3
"""Run a program under an MPI launcher.

Exit 77, Meson's skip code, when no mpirun or mpiexec is on PATH.
Open MPI needs --oversubscribe before it will place more ranks than
the host has cores. Launchers whose help text lacks that flag are
started without it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys


def supports_oversubscribe(launcher: str) -> bool:
    """Open MPI 5's ``mpirun --help`` does not list options.

    A one-rank probe is enough: the flag is accepted, or the launcher
    rejects it. ``--help`` text is checked first for launchers that
    print the option and refuse an unknown probe.
    """
    try:
        listed = subprocess.run(
            [launcher, "--help"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        listed = None
    if listed is not None:
        blob = (listed.stdout or "") + (listed.stderr or "")
        if "--oversubscribe" in blob:
            return True
    try:
        probed = subprocess.run(
            [launcher, "--oversubscribe", "-n", "1", "true"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return probed.returncode == 0


def launch_env() -> dict[str, str]:
    """Environment for the launcher and its ranks.

    Open MPI 5's memory patcher rewrites instructions. A binary built
    with AddressSanitizer dies with SIGILL unless that patcher is off.
    UCX opens the same patcher during MPI_Init on some CPUs, so those
    components are excluded while a sanitizer is active.
    Open MPI selects AVX-512 reductions from CPUID and does not check
    that the OS enables them. A host that advertises AVX-512F uses the
    scalar op component so a faulting instruction cannot kill the rank.
    LeakSanitizer also counts allocations that belong to Open MPI, so
    leak detection is disabled for this process only.
    """
    env = os.environ.copy()
    env["OMPI_MCA_memory"] = "^patcher"
    if env.get("ASAN_OPTIONS"):
        env["OMPI_MCA_pml"] = "^ucx"
        env["OMPI_MCA_osc"] = "^ucx"
        env["OMPI_MCA_btl"] = "^uct"
    if host_advertises_avx512():
        env["OMPI_MCA_op"] = "^avx"
    parts = [
        part
        for part in env.get("ASAN_OPTIONS", "").split(":")
        if part and not part.startswith("detect_leaks=")
    ]
    parts.append("detect_leaks=0")
    env["ASAN_OPTIONS"] = ":".join(parts)
    return env


def host_advertises_avx512() -> bool:
    """True when ``/proc/cpuinfo`` lists the AVX-512F flag."""
    try:
        text = open("/proc/cpuinfo", encoding="utf-8", errors="replace").read()
    except OSError:
        return False
    for line in text.splitlines():
        if line.startswith("flags") and "avx512f" in line.split():
            return True
    return False


def mca_accepted(launcher: str, flags: list[str]) -> bool:
    try:
        probed = subprocess.run(
            [launcher, *flags, "-n", "1", "true"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return probed.returncode == 0


def openmpi_mca_args(launcher: str, asan: bool) -> list[str]:
    """Command-line MCA flags, or an empty list when the launcher rejects them.

    ``OMPI_MCA_*`` in the environment is not enough on every runner.
    ``mpirun --mca`` reaches the ranks directly. MPICH has no ``--mca``.
    """
    flags = ["--mca", "memory", "^patcher"]
    if asan:
        flags.extend(
            [
                "--mca",
                "pml",
                "^ucx",
                "--mca",
                "osc",
                "^ucx",
                "--mca",
                "btl",
                "^uct",
            ]
        )
    if not mca_accepted(launcher, flags):
        return []
    if host_advertises_avx512():
        with_scalar_ops = [*flags, "--mca", "op", "^avx"]
        if mca_accepted(launcher, with_scalar_ops):
            return with_scalar_ops
    return flags


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(
            "usage: run_mpi_test.py N executable [args...]",
            file=sys.stderr,
        )
        return 2
    nprocs = argv[1]
    executable = argv[2]
    rest = argv[3:]
    launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if not launcher:
        print("skipped: no MPI launcher (mpirun or mpiexec) on PATH")
        return 77
    command = [launcher]
    command.extend(openmpi_mca_args(launcher, bool(os.environ.get("ASAN_OPTIONS"))))
    command.extend(["-n", nprocs])
    if supports_oversubscribe(launcher):
        command.append("--oversubscribe")
    command.append(executable)
    command.extend(rest)
    print("running:", " ".join(command), flush=True)
    completed = subprocess.run(command, check=False, env=launch_env())
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
