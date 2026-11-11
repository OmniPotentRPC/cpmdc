#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Medians, interquartile ranges, paired wins, and rank efficiency for ABBA runs.

Reads b<block>-s<slot>-<side>.csv from run_pwscale_abba.sh.
Each block pairs slot 1 with slot 2 and slot 4 with slot 3, so every pair
is one A pass next to one B pass. A pair is a B win when the B pass median
is below the A pass median.
Strong efficiency is T(1) / (n T(n)) and weak efficiency is T(1) / T(n),
with T(1) the one-rank median of the same side.
"""

import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

NAME = re.compile(r"b(\d+)-s(\d+)-([AB])\.csv$")
PAIRS = ((1, 2), (4, 3))


def quartiles(values):
    if len(values) < 2:
        return values[0], values[0]
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[0], q[2]


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: pwscale_abba.py BENCH_DIR [summary.csv] [--metric call_s|seconds]", file=sys.stderr)
        return 2
    metric = "call_s"
    args = sys.argv[1:]
    if "--metric" in args:
        i = args.index("--metric")
        metric = args[i + 1]
        del args[i : i + 2]
    bench = Path(args[0])
    pooled = defaultdict(list)
    passes = defaultdict(list)
    for path in sorted(bench.glob("b*-s*-*.csv")):
        m = NAME.search(path.name)
        if not m:
            continue
        block, slot, side = int(m.group(1)), int(m.group(2)), m.group(3)
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["kind"] != "time" or not row.get(metric):
                    continue
                key = (row["mode"], row["scaling"], row["parallel"], int(row["workers"]))
                value = float(row[metric])
                pooled[(side,) + key].append(value)
                passes[(block, slot, side) + key].append(value)

    keys = sorted({k[1:] for k in pooled}, key=lambda k: (k[0], k[2], k[3]))
    blocks = sorted({k[0] for k in passes})
    out = []
    for key in keys:
        mode, scaling, parallel, workers = key
        row = {"mode": mode, "scaling": scaling, "parallel": parallel, "workers": workers}
        for side in "AB":
            samples = pooled.get((side,) + key, [])
            if not samples:
                continue
            q1, q3 = quartiles(samples)
            med = statistics.median(samples)
            row[f"{side}_n"] = len(samples)
            row[f"{side}_median_s"] = med
            row[f"{side}_iqr_s"] = q3 - q1
            base = pooled.get((side, mode, scaling, parallel, 1))
            if parallel == "ranks" and base:
                t1 = statistics.median(base)
                eff = t1 / (workers * med) if scaling == "strong" else t1 / med
                row[f"{side}_efficiency"] = eff
            serial = pooled.get((side, mode, scaling, "threads", 1))
            if parallel == "ranks" and serial and scaling == "strong":
                row[f"{side}_vs_serial"] = statistics.median(serial) / med
        wins = total = 0
        for block in blocks:
            for a_slot, b_slot in PAIRS:
                a = passes.get((block, a_slot, "A") + key)
                b = passes.get((block, b_slot, "B") + key)
                if not a or not b:
                    continue
                total += 1
                if statistics.median(b) < statistics.median(a):
                    wins += 1
        row["b_wins"] = wins
        row["pairs"] = total
        if row.get("A_median_s") and row.get("B_median_s"):
            row["b_over_a"] = row["B_median_s"] / row["A_median_s"]
        out.append(row)

    fields = [
        "mode", "scaling", "parallel", "workers",
        "A_n", "A_median_s", "A_iqr_s", "A_efficiency", "A_vs_serial",
        "B_n", "B_median_s", "B_iqr_s", "B_efficiency", "B_vs_serial",
        "b_over_a", "b_wins", "pairs",
    ]
    if len(args) >= 2:
        with open(args[1], "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in out:
                writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in row.items()})

    def fmt(row, name):
        value = row.get(name)
        return "-" if value is None else (f"{value:.6g}" if isinstance(value, float) else str(value))

    print(f"metric {metric}")
    for row in out:
        print(
            f"{row['mode']} {row['parallel']} n={row['workers']}"
            f" A median={fmt(row, 'A_median_s')} iqr={fmt(row, 'A_iqr_s')} eff={fmt(row, 'A_efficiency')}"
            f" B median={fmt(row, 'B_median_s')} iqr={fmt(row, 'B_iqr_s')} eff={fmt(row, 'B_efficiency')}"
            f" B/A={fmt(row, 'b_over_a')} wins={row['b_wins']}/{row['pairs']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
