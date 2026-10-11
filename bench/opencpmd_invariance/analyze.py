#!/usr/bin/env python3
"""Build the bitwise, attribution, cost, and scaling tables from the CSV shards."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SKIP = {"ngw", "nrxpl_lo", "nrxpl_hi"}


def load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def key_of(row: dict) -> tuple:
    return (row["quantity"], int(row["index"]))


def group_runs(rows: list[dict]) -> dict:
    runs = defaultdict(dict)
    for row in rows:
        if row["quantity"] in SKIP:
            continue
        run = (int(row["ranks"]), int(row["threads"]), int(row["rep"]))
        runs[run][key_of(row)] = (row["value_hex"], float(row["value"]), int(row["iters"]), float(row["wall_step_s"]))
    return runs


def bitwise_table(system: str, rows: list[dict]) -> list[str]:
    runs = group_runs(rows)
    ref_run = (1, 1, 0)
    if ref_run not in runs:
        return [f"missing {system} reference 1x1 rep 0"]
    ref = runs[ref_run]
    cells = sorted({(r, t) for r, t, _rep in runs})
    quantities = sorted({q for q, _i in ref})
    lines = [
        "| ranks | threads | quantity | n | equal to 1x1 rep0 | equal to each other | max abs | unit |",
        "|-------+---------+----------+---+-------------------+---------------------+---------+------|",
    ]
    for ranks, threads in cells:
        reps = sorted(rep for r, t, rep in runs if r == ranks and t == threads)
        for quantity in quantities:
            indexes = sorted(i for q, i in ref if q == quantity)
            n = len(reps)
            equal_ref = 0
            max_abs = 0.0
            hexes = []
            for rep in reps:
                run = runs.get((ranks, threads, rep), {})
                same = True
                for index in indexes:
                    got = run.get((quantity, index))
                    base = ref.get((quantity, index))
                    if got is None or base is None or got[0] != base[0]:
                        same = False
                    if got is not None and base is not None:
                        max_abs = max(max_abs, abs(got[1] - base[1]))
                if same:
                    equal_ref += 1
                hexes.append(tuple(run.get((quantity, index), ("",))[0] for index in indexes))
            counts = defaultdict(int)
            for item in hexes:
                counts[item] += 1
            equal_each = max(counts.values()) if counts else 0
            unit = "Ha/Bohr" if quantity == "force" else "GPa" if quantity == "stress_gpa" else "Ha"
            if quantity == "stress_paiu":
                unit = "Ha"
            lines.append(
                f"| {ranks} | {threads} | {quantity} | {n} | {equal_ref} | {equal_each} | {max_abs:.3e} | {unit} |"
            )
    return lines


def cost_table(by_variant: dict[str, list[dict]]) -> list[str]:
    lines = [
        "| variant | system | ranks | threads | median s | min s | max s | overhead % | oversub |",
        "|---------+--------+-------+---------+----------+-------+-------+------------+---------|",
    ]
    stock = defaultdict(list)
    if "stock" in by_variant:
        for row in by_variant["stock"]:
            if row["quantity"] != "etot":
                continue
            stock[(row["system"], int(row["ranks"]), int(row["threads"]))].append(float(row["wall_step_s"]))
    for variant, rows in by_variant.items():
        buckets = defaultdict(list)
        for row in rows:
            if row["quantity"] != "etot":
                continue
            buckets[(row["system"], int(row["ranks"]), int(row["threads"]))].append(float(row["wall_step_s"]))
        for key in sorted(buckets):
            vals = buckets[key]
            med = statistics.median(vals)
            base_vals = stock.get(key, [])
            base = statistics.median(base_vals) if base_vals else float("nan")
            overhead = 0.0 if not base_vals or base == 0 else 100.0 * (med - base) / base
            system, ranks, threads = key
            over = "yes" if ranks * threads > 8 else "no"
            lines.append(
                f"| {variant} | {system} | {ranks} | {threads} | {med:.4f} | {min(vals):.4f} | {max(vals):.4f} | {overhead:.1f} | {over} |"
            )
    return lines


def scaling_table(by_variant: dict[str, list[dict]], variants: list[str]) -> list[str]:
    lines = [
        "| variant | system | ranks | threads | median s | min s | max s | efficiency |",
        "|---------+--------+-------+---------+----------+-------+-------+------------|",
    ]
    for variant in variants:
        rows = by_variant.get(variant, [])
        buckets = defaultdict(list)
        for row in rows:
            if row["quantity"] != "etot" or int(row["threads"]) != 1:
                continue
            buckets[(row["system"], int(row["ranks"]))].append(float(row["wall_step_s"]))
        for system in sorted({k[0] for k in buckets}):
            base = statistics.median(buckets[(system, 1)]) if (system, 1) in buckets else float("nan")
            for ranks in sorted(r for s, r in buckets if s == system):
                vals = buckets[(system, ranks)]
                med = statistics.median(vals)
                eff = base / (ranks * med) if med and base == base else float("nan")
                lines.append(
                    f"| {variant} | {system} | {ranks} | 1 | {med:.4f} | {min(vals):.4f} | {max(vals):.4f} | {eff:.3f} |"
                )
    return lines


def _hex_groups(runs: dict, quantity: str, indexes: list[int], selected: list[tuple]) -> tuple[int, int]:
    groups = defaultdict(int)
    for key in selected:
        run = runs.get(key, {})
        hexes = tuple(run.get((quantity, index), ("",))[0] for index in indexes)
        if any(item == "" for item in hexes):
            continue
        groups[hexes] += 1
    return len(groups), sum(groups.values())


def attribution_table(by_system_variant: dict[tuple[str, str], list[dict]]) -> list[str]:
    """Split discrepancies into run-to-run, thread, and rank axes.

    equal_grid is yes only when every rank, thread, and repeat shares one hex.
    thread_groups is the number of distinct hexes across threads at rank 1.
    rank_groups is the number of distinct hexes across ranks at 1 thread.
    A fix that removes one axis shows a single group on that axis.
    """
    lines = [
        "| system | variant | quantity | runs | cells r2r | thread groups at r1 | rank groups at t1 | r8 with r1 | grid groups | equal grid | max abs vs 1x1 |",
        "|--------+---------+----------+------+-----------+---------------------+-------------------+------------+-------------+------------+----------------|",
    ]
    systems = sorted({system for system, _variant in by_system_variant})
    for system in systems:
        variants = [variant for s, variant in sorted(by_system_variant) if s == system]
        for variant in variants:
            runs = group_runs(by_system_variant[(system, variant)])
            if not runs:
                continue
            ref = runs.get((1, 1, 0), next(iter(runs.values())))
            quantities = sorted({q for q, _i in ref})
            nruns = len(runs)
            cells = sorted({(ranks, threads) for ranks, threads, _rep in runs})
            for quantity in quantities:
                indexes = sorted(i for q, i in ref if q == quantity)
                max_abs = 0.0
                base_run = runs.get((1, 1, 0))
                r2r_cells = 0
                r2r_seen = 0
                for ranks, threads in cells:
                    reps = [rep for r, t, rep in runs if r == ranks and t == threads]
                    if len(reps) < 2:
                        continue
                    r2r_seen += 1
                    groups_cell, _n = _hex_groups(
                        runs, quantity, indexes, [(ranks, threads, rep) for rep in reps]
                    )
                    if groups_cell == 1:
                        r2r_cells += 1
                if base_run is not None:
                    for run in runs.values():
                        for index in indexes:
                            got = run.get((quantity, index))
                            base = base_run.get((quantity, index))
                            if got is not None and base is not None:
                                max_abs = max(max_abs, abs(got[1] - base[1]))
                thread_keys = [
                    key for key in runs if key[0] == 1 and (quantity, indexes[0]) in runs[key]
                ]
                rank_keys = [
                    key
                    for key in runs
                    if key[1] == 1 and key[0] in (1, 2, 3, 4) and (quantity, indexes[0]) in runs[key]
                ]
                r8_keys = [
                    key
                    for key in runs
                    if key[1] == 1 and key[0] in (1, 8) and (quantity, indexes[0]) in runs[key]
                ]
                thread_groups, _tn = _hex_groups(runs, quantity, indexes, thread_keys)
                rank_groups, _rn = _hex_groups(runs, quantity, indexes, rank_keys)
                r8_groups, _r8n = _hex_groups(runs, quantity, indexes, r8_keys)
                has_r8 = any(key[0] == 8 for key in r8_keys)
                grid_groups, _gn = _hex_groups(runs, quantity, indexes, list(runs))
                equal = "yes" if grid_groups == 1 and nruns > 0 else "no"
                lines.append(
                    "| "
                    + " | ".join(
                        [
                            system,
                            variant,
                            quantity,
                            str(nruns),
                            f"{r2r_cells}/{r2r_seen}",
                            str(thread_groups),
                            str(rank_groups),
                            str(r8_groups) if has_r8 else "na",
                            str(grid_groups),
                            equal,
                            f"{max_abs:.3e}",
                        ]
                    )
                    + " |"
                )
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "tables.org")
    args = parser.parse_args()
    shards = sorted(args.results.glob("*_*.csv"))
    by_variant = defaultdict(list)
    by_system_variant = defaultdict(list)
    for path in shards:
        rows = load(path)
        if not rows:
            continue
        variant = rows[0]["variant"]
        system = rows[0]["system"]
        by_variant[variant].extend(rows)
        by_system_variant[(system, variant)].extend(rows)
    parts = ["* Bitwise"]
    for system, variant in sorted(by_system_variant):
        parts.append(f"** {system} {variant}")
        parts.extend(bitwise_table(system, by_system_variant[(system, variant)]))
    parts.append("* Attribution")
    parts.extend(attribution_table(by_system_variant))
    parts.append("* Cost")
    parts.extend(cost_table(by_variant))
    parts.append("* Strong scaling")
    parts.extend(scaling_table(by_variant, [v for v in ("stock", "full") if v in by_variant]))
    text = "\n".join(parts) + "\n"
    args.out.write_text(text, encoding="utf-8")
    print(args.out)


if __name__ == "__main__":
    main()
