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
import signal
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


def btl_exclusion(asan: bool) -> str:
    """Open MPI 4's openib BTL probes RDMA devices it has no parameters
    for, such as the Azure MANA NIC of GitHub runners, and mpirun can then
    abort with ORTE_ERROR_LOG: Out of resource in show_help.c. Open MPI 5
    has no openib and ignores the name.
    """
    return "^uct,openib" if asan else "^openib"


def launch_env() -> dict[str, str]:
    """Environment for the launcher and its ranks.

    A traced CI failure reached libfabric's memhooks code while it was
    replacing AddressSanitizer's mprotect interceptor during MPI_Init.
    A separate Open MPI 5 reproduction still failed with that monitor off
    and passed with Open MPI's memory patcher off. The sanitize job and this
    launcher disable both observed patching paths. UCX is also excluded under
    a sanitizer because selecting another transport avoided the reproduction.
    LeakSanitizer also counts allocations that belong to Open MPI, so
    leak detection is disabled for this process only.
    """
    env = os.environ.copy()
    env["OMPI_MCA_memory"] = "^patcher"
    env["OMPI_MCA_btl"] = btl_exclusion(bool(env.get("ASAN_OPTIONS")))
    if env.get("ASAN_OPTIONS"):
        env["OMPI_MCA_pml"] = "^ucx"
        env["OMPI_MCA_osc"] = "^ucx"
    parts = [
        part
        for part in env.get("ASAN_OPTIONS", "").split(":")
        if part and not part.startswith("detect_leaks=")
    ]
    parts.append("detect_leaks=0")
    env["ASAN_OPTIONS"] = ":".join(parts)
    return env


def openmpi_mca_args(launcher: str, asan: bool) -> list[str]:
    """Command-line MCA flags, or an empty list when the launcher rejects them.

    ``OMPI_MCA_*`` in the environment is not enough on every runner.
    ``mpirun --mca`` reaches the ranks directly. MPICH has no ``--mca``.
    """
    flags = ["--mca", "memory", "^patcher", "--mca", "btl", btl_exclusion(asan)]
    if asan:
        flags.extend(["--mca", "pml", "^ucx", "--mca", "osc", "^ucx"])
    try:
        probed = subprocess.run(
            [launcher, *flags, "-n", "1", "true"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    if probed.returncode != 0:
        return []
    return flags


def run_with_optional_timeout(
    command: list[str], env: dict[str, str], timeout_s: float | None
) -> int:
    if timeout_s is None:
        completed = subprocess.run(command, check=False, env=env)
        return completed.returncode
    proc = subprocess.Popen(command, env=env, start_new_session=True)
    try:
        return proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        print(
            f"FAIL: timed out after {timeout_s:g}s: {' '.join(command)}",
            file=sys.stderr,
            flush=True,
        )
        return 124


def main(argv: list[str]) -> int:
    args = argv[1:]
    timeout_s: float | None = None
    if len(args) >= 2 and args[0] == "--timeout":
        try:
            timeout_s = float(args[1])
        except ValueError:
            print("invalid --timeout value", file=sys.stderr)
            return 2
        if timeout_s <= 0:
            print("--timeout must be positive", file=sys.stderr)
            return 2
        args = args[2:]
    if len(args) < 2:
        print(
            "usage: run_mpi_test.py [--timeout SECONDS] N|--direct executable [args...]",
            file=sys.stderr,
        )
        return 2
    nprocs = args[0]
    executable = args[1]
    rest = args[2:]
    if nprocs == "--direct":
        return run_with_optional_timeout([executable, *rest], launch_env(), timeout_s)
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
    return run_with_optional_timeout(command, launch_env(), timeout_s)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
