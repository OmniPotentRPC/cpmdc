#!/usr/bin/env python3
# Copyright (c) cpmdc developers
"""Equilibration and block-averaged means of gromacs_opencpmd_npt.py runs.

Reads the --json file of bench/gromacs_opencpmd_npt.py, which holds every
energy frame of the QM/MM run and of the classical control, and optionally
the --json file of a longer --classical-only run as a reference, and writes
Org tables of the volume, density and pressure from those frames alone.

Equilibration: for each observable the start t0 in the first half of the
run that maximises the number of uncorrelated frames after it, (N - t0) / g,
with g the statistical inefficiency of the frames from t0 on (Chodera,
J. Chem. Theory Comput. 12, 1799 (2016)). A run discards the latest t0 of
its three observables, and at least --min-discard-ps (default four times
tau-p); --discard-ps replaces both.

Error bars: the standard error of the mean of the production frames from
Flyvbjerg-Petersen blocking (J. Chem. Phys. 91, 461 (1989)): the frames
are averaged in pairs level by level, and the error is that of the coarsest
level with at least --min-blocks blocks. The tables give that error's own
uncertainty, the error at half the block length, and the blocking curve:
while the error still grows with the block length, the blocks are not yet
longer than the correlation time and the error is a lower bound. plateau()
finds the first level whose blocks outlast the correlation by the criterion
of Lee et al. (Phys. Rev. E 83, 066706 (2011)).

The window tables give the means over consecutive windows, which show
whether the volume and density have levelled off.

usage: npt_block_average.py RUN.json [--reference LONG.json]
                            [--discard-ps PS | --min-discard-ps PS]
                            [--window-ps PS] [--reference-window-ps PS]
                            [--min-blocks N] [--output FILE]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

TERMS = (("Volume", "nm^3", 4), ("Density", "kg/m^3", 1), ("Pressure", "bar", 0))
CURVE_LEVELS = range(0, 16, 2)


def statistical_inefficiency(x: np.ndarray) -> float:
    """1 + 2 sum (1 - t/N) C(t), summed until the autocorrelation C(t) first drops to zero or below."""
    n = len(x)
    dx = x - x.mean()
    var = float(dx @ dx) / n
    if n < 3 or var == 0.0:
        return 1.0
    spectrum = np.fft.rfft(dx, 2 * n)
    acov = np.fft.irfft(spectrum * np.conj(spectrum))[:n] / np.arange(n, 0, -1)
    c = acov[1:] / var
    stop = np.flatnonzero(c <= 0.0)
    c = c[: stop[0]] if stop.size else c
    t = np.arange(1, len(c) + 1)
    return max(1.0, 1.0 + 2.0 * float(np.sum((1.0 - t / n) * c)))


def equilibration_start(x: np.ndarray, stride: int = 1) -> tuple[int, float]:
    """The frame t0 in the first half that maximises (N - t0) / g(x[t0:]), and that g.

    Later starts leave so few frames that g is underestimated and a start
    near the end can win.
    """
    n = len(x)
    best = (0, 1.0, -1.0)
    for t0 in range(0, n // 2 + 1, stride):
        g = statistical_inefficiency(x[t0:])
        if (n - t0) / g > best[2]:
            best = (t0, g, (n - t0) / g)
    return best[0], best[1]


def blocking(x: np.ndarray) -> list[tuple[int, float, float]]:
    """Flyvbjerg-Petersen levels: (blocks, standard error of the mean, its uncertainty)."""
    levels = []
    b = np.asarray(x, dtype=float)
    while len(b) >= 2:
        n = len(b)
        sem = float(np.sqrt(np.mean((b - b.mean()) ** 2) / (n - 1)))
        levels.append((n, sem, sem / np.sqrt(2.0 * (n - 1))))
        b = 0.5 * (b[: n - n % 2 : 2] + b[1 : n - n % 2 : 2])
    return levels


def plateau(levels: list[tuple[int, float, float]], min_blocks: int) -> int | None:
    """The first blocking level, of those with at least min_blocks blocks, whose blocks outlast the correlation.

    Blocks of B frames, of n, outlast it when B^3 > 2 n (SE_B / SE_1)^4, where
    (SE_B / SE_1)^2 is the correlation length that they show (Lee et al.,
    Phys. Rev. E 83, 066706 (2011), Eq. 14 and the text after it). Lee et al.
    scan down from the longest block, take the last B that meets it, and ask
    for B < n/50, which gives the error to 10%. Scanning up gives the same
    level while every longer block meets it too, and a coarse level of few
    blocks, whose error is uncertain, cannot hide it. None means that no level
    with enough blocks meets it: the run is too short for the blocking curve
    to level off.
    """
    n, first = levels[0][0], levels[0][1]
    for k, (blocks, sem, _) in enumerate(levels):
        if blocks < min_blocks:
            break
        if first == 0.0 or 2.0 ** (3 * k) > 2.0 * n * (sem / first) ** 4:
            return k
    return None


def block_error(x: np.ndarray, min_blocks: int) -> dict:
    """The coarsest blocking level with at least min_blocks blocks, the level before it, and the plateau level."""
    levels = blocking(x)
    k = max(i for i, lv in enumerate(levels) if lv[0] >= min_blocks)
    return {"blocks": levels[k][0], "block_frames": 2**k, "sem": levels[k][1], "sem_err": levels[k][2],
            "sem_half": levels[k - 1][1] if k else None, "curve": levels, "plateau": plateau(levels, min_blocks)}


def analyse(run: dict, discard_ps: float | None, min_discard_ps: float, window_ps: float,
            min_blocks: int) -> dict:
    t = np.asarray(run["time_ps"])
    dt = float(t[1] - t[0])
    stride = max(1, len(t) // 2000)
    starts = {name: equilibration_start(np.asarray(run[name]), stride)[0] for name, _, _ in TERMS}
    if discard_ps is None:
        discard_ps = max(min_discard_ps, max(float(t[v]) for v in starts.values()))
    t0 = int(np.searchsorted(t, discard_ps - 1e-6 * dt))
    out = {
        "frame_ps": dt, "length_ps": float(t[-1]), "discard_ps": float(t[t0]),
        "production_frames": len(t) - t0, "detected_start_ps": {k: float(t[v]) for k, v in starts.items()},
        "wall_s": run.get("wall_s"), "terms": {}, "windows": [],
    }
    for name, _, _ in TERMS:
        x = np.asarray(run[name])[t0:]
        out["terms"][name] = {"mean": float(x.mean()), "std": float(x.std()), **block_error(x, min_blocks)}
    edges = np.arange(0.0, t[-1] + 0.5 * dt, window_ps)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (t > lo - 0.5 * dt) & ((t < hi - 0.5 * dt) | (hi == edges[-1]))
        out["windows"].append((float(lo), float(hi), *(float(np.mean(np.asarray(run[n])[sel])) for n, _, _ in TERMS)))
    return out


def summary_rows(results: list[tuple[str, dict]]) -> list[str]:
    lines = [
        "| Run | Length, ps | Discarded, ps | Frames used | Volume, nm^3 | Density, kg/m^3 | Pressure, bar | Wall, s |",
        "|-----+------------+---------------+-------------+--------------+-----------------+---------------+---------|",
    ]
    for label, r in results:
        cells = [f"{r['terms'][n]['mean']:.{d}f} +- {r['terms'][n]['sem']:.{d + 1}f}" for n, _, d in TERMS]
        wall = f"{r['wall_s']:.0f}" if r["wall_s"] is not None else ""
        lines.append(f"| {label} | {r['length_ps']:g} | {r['discard_ps']:g} | {r['production_frames']} | "
                     + " | ".join(cells) + f" | {wall} |")
    return lines


def detail_rows(results: list[tuple[str, dict]]) -> list[str]:
    lines = [
        "| Run | Observable | Detected start, ps | Std | Blocks | Block, ps | Error | Its uncertainty | Error at half the block |",
        "|-----+------------+--------------------+-----+--------+-----------+-------+-----------------+-------------------------|",
    ]
    for label, r in results:
        for name, unit, d in TERMS:
            s = r["terms"][name]
            half = f"{s['sem_half']:.{d + 1}f}" if s["sem_half"] is not None else ""
            lines.append(f"| {label} | {name}, {unit} | {r['detected_start_ps'][name]:g} | {s['std']:.{d + 1}f} | "
                         f"{s['blocks']} | {s['block_frames'] * r['frame_ps']:g} | {s['sem']:.{d + 1}f} | "
                         f"{s['sem_err']:.{d + 1}f} | {half} |")
    return lines


def curve_rows(results: list[tuple[str, dict]], min_blocks: int, levels=CURVE_LEVELS) -> list[str]:
    dt = results[0][1]["frame_ps"]
    head = " | ".join(f"{2**k * dt:g} ps" for k in levels)
    lines = [f"| Run | Observable | {head} |",
             "|-----+------------+" + "+".join("-" * (len(f"{2**k * dt:g} ps") + 2) for k in levels) + "|"]
    for label, r in results:
        for name, unit, d in TERMS:
            curve = r["terms"][name]["curve"]
            cells = [f"{curve[k][1]:.{d + 1}f}" if k < len(curve) and curve[k][0] >= min_blocks else ""
                     for k in levels]
            lines.append(f"| {label} | {name}, {unit} | " + " | ".join(cells) + " |")
    return lines


def window_rows(results: list[tuple[str, dict]]) -> list[str]:
    head = " | ".join(f"{label} {name.lower()}, {unit}" for label, _ in results for name, unit, _ in TERMS)
    lines = [f"| Window, ps | {head} |",
             "|------------+" + "+".join("-" * (len(f"{label} {name.lower()}, {unit}") + 2)
                                        for label, _ in results for name, unit, _ in TERMS) + "|"]
    for i, w in enumerate(results[0][1]["windows"]):
        cells = []
        for _, r in results:
            win = r["windows"][i]
            cells += [f"{win[2]:.4f}", f"{win[3]:.1f}", f"{win[4]:.0f}"]
        lines.append(f"| {w[0]:g}-{w[1]:g} | " + " | ".join(cells) + " |")
    return lines


def describe(args: dict) -> str:
    return (f"{args['pcoupl']} tau-p {args['tau_p']} ps, {args['steps']} steps of {float(args['dt']) * 1000:g} fs, "
            f"lattice {args['lattice']} spacing {args['spacing']} nm")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("json", type=Path)
    parser.add_argument("--reference", type=Path)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--discard-ps", type=float)
    group.add_argument("--min-discard-ps", type=float)
    parser.add_argument("--window-ps", type=float, default=1.0)
    parser.add_argument("--reference-window-ps", type=float, default=100.0)
    parser.add_argument("--min-blocks", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    def load(path: Path, window_ps: float) -> tuple[dict, dict]:
        data = json.loads(path.read_text())
        floor = args.min_discard_ps if args.min_discard_ps is not None else 4.0 * float(data["args"]["tau_p"])
        runs = {key: analyse(data[key], args.discard_ps, floor, window_ps, args.min_blocks)
                for key in ("qmmm", "classical") if key in data}
        return data["args"], runs

    meta, runs = load(args.json, args.window_ps)
    main_runs = [("QM water", runs["qmmm"]), ("classical", runs["classical"])]
    results = list(main_runs)
    lines = [
        "# Generated by bench/npt_block_average.py from the energy frames of",
        "# bench/gromacs_opencpmd_npt.py; do not edit by hand.",
        f"# QM water and classical: {describe(meta)}, QM run on {meta['ranks']} rank(s)",
    ]
    reference = None
    if args.reference:
        ref_meta, ref_runs = load(args.reference, args.reference_window_ps)
        reference = [("classical, long", ref_runs["classical"])]
        results += reference
        lines.append(f"# classical, long: {describe(ref_meta)}")
    lines += ["", *summary_rows(results), "", *detail_rows(results), "",
              "Standard error of the mean against block length:", "", *curve_rows(results, args.min_blocks),
              "", *window_rows(main_runs)]
    if reference:
        lines += ["", *window_rows(reference)]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
