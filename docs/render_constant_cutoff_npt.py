#!/usr/bin/env python3
"""Draw the density of the constant-cutoff NPT runs of water.

Reads the --json file of bench/constant_cutoff_npt.py table, by default
bench/results/constant_cutoff_npt.json, and writes an SVG. The left panel
draws the density of each run over windows of its record against the time
from the start of the first run before it, with the run's block-averaged mean
and one standard error over the frames it averages; the right panel draws
every step of the runs that start from the same state against the time after
that state. Both draw the PBE density of Miceli, de Gironcoli and
Pasquarello, J. Chem. Phys. 142, 034501 (2015), doi:10.1063/1.4905333, about
0.87 g/cm^3 at 349 K. Colour follows the cutoff and a dashed line is a run
without the constant cutoff; the legend gives each run's name, basis, and
mean. The script refuses a file without runs, a run with fewer than two
windows or steps, a value that is not finite, and a file where no run has a
mean.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = ROOT / "bench" / "results" / "constant_cutoff_npt.json"
REFERENCE_G_CM3 = 0.87
REFERENCE_LABEL = ("PBE, Miceli et al. (2015)", "about 0.87 g/cm^3 at 349 K")
COLOURS = ("#5b21b6", "#0f766e", "#b45309", "#be123c")
WIDTH, HEIGHT = 900, 440
LEFT, RIGHT, TOP, BOTTOM, GAP = 64, 230, 24, 52, 36
BRANCH_SHARE = 0.3


def load(path: Path) -> list[dict]:
    runs = json.loads(path.read_text(encoding="utf-8"))["runs"]
    for run in runs:
        start = float(run.get("chain_start_ps") or 0.0)
        if "steps" in run:
            run["points"] = [(t, rho / 1000.0) for t, rho in run["steps"]]
        else:
            run["points"] = [(start + t, rho / 1000.0) for t, rho in run["windows"]]
        if run.get("mean") is not None:
            offset = 0.0 if "steps" in run else start
            run["band"] = (offset + run["discard_ps"], offset + run["length_ps"],
                           run["mean"] / 1000.0, run["sem"] / 1000.0)
    return runs


def check(runs: list[dict]) -> None:
    if not runs:
        raise SystemExit("the file holds no runs")
    for run in runs:
        if len(run["points"]) < 2:
            raise SystemExit(f"{run['label']} has fewer than two {'steps' if 'steps' in run else 'windows'}")
        values = [v for point in run["points"] for v in point] + list(run.get("band", ()))
        if not all(math.isfinite(v) for v in values):
            raise SystemExit(f"{run['label']} has a value that is not finite")
    if not any("band" in run for run in runs):
        raise SystemExit("no run has a mean")


def ticks(lo: float, hi: float) -> list[float]:
    """At most eight ticks at a step of 1, 2, or 5 times a power of ten."""
    scale = 10.0 ** math.floor(math.log10(hi - lo))
    step = next(m * scale / 10.0 for m in (1, 2, 5, 10, 20, 50) if (hi - lo) / (m * scale / 10.0) <= 8)
    first = math.ceil(lo / step - 1e-9)
    return [round(k * step, 10) for k in range(first, int(math.floor(hi / step + 1e-9)) + 1)]


def svg(runs: list[dict]) -> str:
    ys = [rho for run in runs for _, rho in run["points"]] + [REFERENCE_G_CM3]
    ys += [b[2] + s * b[3] for run in runs if "band" in run for b in [run["band"]] for s in (-1, 1)]
    pad = 0.05 * (max(ys) - min(ys))
    y0, y1 = min(ys) - pad, max(ys) + pad
    plot_w, plot_h = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM
    panels = [[run for run in runs if "steps" not in run], [run for run in runs if "steps" in run]]
    if not panels[1]:
        boxes = [(LEFT, plot_w)]
    else:
        branch_w = BRANCH_SHARE * (plot_w - GAP)
        boxes = [(LEFT, plot_w - GAP - branch_w), (LEFT + plot_w - branch_w, branch_w)]
    labels = ("time, ps", "time after the state, ps")

    def py(y: float) -> float:
        return TOP + (y1 - y) / (y1 - y0) * plot_h

    cutoffs = sorted({run["cutoff"] for run in runs})
    colour = {cutoff: COLOURS[i % len(COLOURS)] for i, cutoff in enumerate(cutoffs)}
    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="sans-serif" font-size="12">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    for (x_left, width), drawn, label in zip(boxes, panels, labels):
        x1 = max(t for run in drawn for t, _ in run["points"]) * 1.02

        def px(x: float, x_left: float = x_left, width: float = width, x1: float = x1) -> float:
            return x_left + x / x1 * width

        out.append(f'<rect x="{x_left:.2f}" y="{TOP}" width="{width:.2f}" height="{plot_h}" fill="none" '
                   'stroke="#444444"/>')
        for x in ticks(0.0, x1):
            out.append(f'<line x1="{px(x):.2f}" y1="{TOP + plot_h}" x2="{px(x):.2f}" y2="{TOP + plot_h + 5}" '
                       'stroke="#444444"/>')
            out.append(f'<text x="{px(x):.2f}" y="{TOP + plot_h + 18}" text-anchor="middle">{x:g}</text>')
        for y in ticks(y0, y1):
            out.append(f'<line x1="{x_left - 5:.2f}" y1="{py(y):.2f}" x2="{x_left + width:.2f}" y2="{py(y):.2f}" '
                       'stroke="#e5e5e5"/>')
            if x_left == LEFT:
                out.append(f'<text x="{LEFT - 8}" y="{py(y) + 4:.2f}" text-anchor="end">{y:.2f}</text>')
        out.append(f'<text x="{x_left + width / 2:.2f}" y="{HEIGHT - 14}" text-anchor="middle">{label}</text>')
        out.append(f'<line x1="{x_left:.2f}" y1="{py(REFERENCE_G_CM3):.2f}" x2="{x_left + width:.2f}" '
                   f'y2="{py(REFERENCE_G_CM3):.2f}" stroke="#111111" stroke-width="1.5" stroke-dasharray="2 4"/>')
        for run in drawn:
            stroke = colour[run["cutoff"]]
            dash_attr = "" if run.get("e0") is not None else ' stroke-dasharray="7 4"'
            if "band" in run:
                lo, hi, mean, sem = run["band"]
                out.append(f'<rect x="{px(lo):.2f}" y="{py(mean + sem):.2f}" width="{px(hi) - px(lo):.2f}" '
                           f'height="{py(mean - sem) - py(mean + sem):.2f}" fill="{stroke}" fill-opacity="0.18"/>')
                out.append(f'<line x1="{px(lo):.2f}" y1="{py(mean):.2f}" x2="{px(hi):.2f}" y2="{py(mean):.2f}" '
                           f'stroke="{stroke}" stroke-width="1"{dash_attr}/>')
            poly = " ".join(f"{px(t):.2f},{py(rho):.2f}" for t, rho in run["points"])
            out.append(f'<polyline fill="none" stroke="{stroke}" stroke-width="2"{dash_attr} points="{poly}"/>')
    out.append(f'<text x="16" y="{TOP + plot_h / 2:.2f}" text-anchor="middle" '
               f'transform="rotate(-90 16 {TOP + plot_h / 2:.2f})">density, g/cm^3</text>')
    legend = [(REFERENCE_LABEL, "#111111", "2 4")]
    for run in runs:
        value = (f"{run['band'][2]:.3f} +- {run['band'][3]:.3f} g/cm^3" if "band" in run
                 else str(run.get("stopped") or "no mean"))
        legend.append(((run["label"], run["basis"], value), colour[run["cutoff"]],
                       "" if run.get("e0") is not None else "7 4"))
    y = TOP + 14
    for lines, stroke, dash in legend:
        x = LEFT + plot_w + 14
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<line x1="{x}" y1="{y - 4}" x2="{x + 24}" y2="{y - 4}" stroke="{stroke}" '
                   f'stroke-width="2"{dash_attr}/>')
        for k, text in enumerate(lines):
            fill = "" if k == 0 else ' fill="#444444"'
            out.append(f'<text x="{x + 30}" y="{y + 14 * k}" font-size="11"{fill}>{text}</text>')
        y += 14 * len(lines) + 8
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    runs = load(args.json)
    check(runs)
    text = svg(runs)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
