#!/usr/bin/env python3
"""Median speedups from the scaffold CSV. Five-run medians, hardware copied through."""

import csv
import statistics
import sys
from collections import defaultdict


def median(values):
    return statistics.median(values)


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: pwscale_summarize.py raw.csv [summary.csv]", file=sys.stderr)
        sys.exit(2)
    groups = defaultdict(list)
    checks = {}
    hardware = None
    with open(sys.argv[1], newline="") as handle:
        for row in csv.DictReader(handle):
            if hardware is None and row.get("cpu_model"):
                hardware = row
            key = (
                row["mode"],
                row["dist"],
                row["scaling"],
                row["backend"],
                row.get("parallel") or "threads",
                row["workers"],
                row["natoms"],
                row["npw"],
            )
            if row["kind"] == "time" and row["seconds"]:
                groups[key].append(float(row["seconds"]))
            if row["kind"] == "check":
                checks[key] = row

    baselines = {}
    for key, samples in groups.items():
        mode, dist, scaling, backend, parallel, workers, _natoms, _npw = key
        if workers == "1":
            baselines[(mode, dist, scaling, backend, parallel)] = median(samples)

    summary = []
    for key in sorted(groups, key=lambda item: (item[0], item[4], int(item[5]))):
        samples = groups[key]
        mode, dist, scaling, backend, parallel, workers, natoms, npw = key
        med = median(samples)
        base = baselines.get((mode, dist, scaling, backend, parallel))
        speedup = (base / med) if base and med > 0 else None
        n = int(workers)
        if speedup is None:
            efficiency = None
        elif scaling == "strong":
            efficiency = speedup / n if n else None
        else:
            efficiency = speedup
        check = checks.get(key, {})
        summary.append(
            {
                "mode": mode,
                "dist": dist,
                "scaling": scaling,
                "backend": backend,
                "parallel": parallel,
                "workers": workers,
                "natoms": natoms,
                "npw": npw,
                "runs": len(samples),
                "median_s": f"{med:.6g}",
                "min_s": f"{min(samples):.6g}",
                "max_s": f"{max(samples):.6g}",
                "speedup": "" if speedup is None else f"{speedup:.4f}",
                "efficiency": "" if efficiency is None else f"{efficiency:.4f}",
                "max_abs_dE": check.get("max_abs_dE", ""),
                "max_abs_dF": check.get("max_abs_dF", ""),
                "counts_match": check.get("counts_match", ""),
                "alloc_bytes": check.get("alloc_bytes", ""),
                "plan_s": check.get("plan_s", ""),
            }
        )

    if len(sys.argv) >= 3:
        with open(sys.argv[2], "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()) if summary else ["mode"])
            writer.writeheader()
            writer.writerows(summary)

    if hardware:
        print(
            f"cpu_model={hardware['cpu_model']} nproc={hardware['nproc']} "
            f"uname_m={hardware['uname_m']} gpu={hardware['gpu']}"
        )
    for row in summary:
        print(
            f"{row['mode']} {row['dist']} {row['parallel']} workers={row['workers']} "
            f"natoms={row['natoms']} npw={row['npw']} median_s={row['median_s']} "
            f"speedup={row['speedup']} efficiency={row['efficiency']} "
            f"dE={row['max_abs_dE']} dF={row['max_abs_dF']} counts={row['counts_match']}"
        )


if __name__ == "__main__":
    main()
