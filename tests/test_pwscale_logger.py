#!/usr/bin/env python3
"""Regression for the scaffold logger.

A timed row must carry that sample's gaps and visit flag, and a poisoned
sample must not exit as a pass.
"""

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, env):
    return subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)


def rows_of(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_pwscale_logger.py HARNESS", file=sys.stderr)
        return 2
    harness = sys.argv[1]
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / "raw.csv"
        env = os.environ.copy()
        env["OMP_DYNAMIC"] = "false"
        proc = run(
            [
                harness,
                "--mode",
                "force-strong",
                "--threads",
                "2,3",
                "--runs",
                "2",
                "--natoms",
                "4",
                "--npw",
                "27",
                "--out",
                str(raw),
            ],
            env,
        )
        if proc.returncode != 0:
            failures.append(f"clean harness exit {proc.returncode}\n{proc.stderr}")
        else:
            rows = rows_of(raw)
            times = [row for row in rows if row["kind"] == "time"]
            checks = [row for row in rows if row["kind"] == "check"]
            if len(checks) != 2 or len(times) != 4:
                failures.append(f"expected 2 checks and 4 times, got {len(checks)} {len(times)}")
            for row in times:
                if row["max_abs_dE"] == "" or row["max_abs_dF"] == "" or row["counts_match"] == "":
                    failures.append(f"timed row left a gap blank: {row}")
                if row["counts_match"] != "1":
                    failures.append(f"clean timed visit flag {row['counts_match']}")
                if row["team_observed"] != "1":
                    failures.append(f"team was not observed: {row}")
                if row["omp_threads"] != row["workers"]:
                    failures.append(f"workers {row['workers']} omp {row['omp_threads']}")
                if row["omp_threads"] not in {"2", "3"}:
                    failures.append(f"unexpected team {row['omp_threads']}")
                if row["mpi_ranks"] != "1":
                    failures.append(f"mpi ranks {row['mpi_ranks']}")
        poisoned = Path(tmp) / "poison.csv"
        env["CPMDC_PWSCALE_POISON_RUN"] = "2"
        proc = run(
            [
                harness,
                "--mode",
                "force-strong",
                "--threads",
                "2",
                "--runs",
                "2",
                "--natoms",
                "4",
                "--npw",
                "27",
                "--out",
                str(poisoned),
            ],
            env,
        )
        if proc.returncode == 0:
            failures.append("poisoned timed sample exited 0")
        if poisoned.is_file():
            bad = [row for row in rows_of(poisoned) if row["kind"] == "time" and row["run"] == "2"]
            if len(bad) != 1:
                failures.append(f"poisoned log missing run 2: {poisoned.read_text()}")
            elif float(bad[0]["max_abs_dE"]) < 0.5:
                failures.append(f"poisoned dE stayed at the preflight gap: {bad[0]['max_abs_dE']}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("pwscale logger ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
