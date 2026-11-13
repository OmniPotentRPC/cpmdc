#!/usr/bin/env python3
"""Draw the reference harmonic well from the evaluator's formula.

The curve is one hydrogen atom, cutoff 70 Ry, and no deck offset:
k = 1e-3 * max(0.1, cutoff/70), energy = 0.5 * k * (x/bohr)^2.
The script writes an SVG. It refuses a file that does not rise on both
sides of the origin.
"""

from __future__ import annotations

import argparse
from pathlib import Path

BOHR = 0.529177210903


def samples(count: int = 41) -> list[tuple[float, float]]:
    cutoff = 70.0
    k = 1.0e-3 * max(0.1, cutoff / 70.0)
    points = []
    for step in range(count):
        x = -1.0 + 2.0 * step / (count - 1)
        r = x / BOHR
        energy = 0.5 * k * r * r
        points.append((x, energy))
    return points


def svg(points: list[tuple[float, float]]) -> str:
    width, height = 640, 360
    pad = 36
    xs = [item[0] for item in points]
    ys = [item[1] for item in points]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)

    def px(x: float) -> float:
        return pad + (x - x0) / (x1 - x0) * (width - 2 * pad)

    def py(y: float) -> float:
        return height - pad - (y - y0) / (y1 - y0) * (height - 2 * pad)

    poly = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in points)
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        f'<polyline fill="none" stroke="#5b21b6" stroke-width="2" points="{poly}"/>\n'
        f'<text x="{pad}" y="24" font-size="14">reference well, Hartree</text>\n'
        "</svg>\n"
    )


def check(points: list[tuple[float, float]]) -> None:
    if len(points) < 8:
        raise SystemExit("the well needs at least 8 samples")
    mid = min(points, key=lambda item: abs(item[0]))
    if mid[1] > points[0][1] or mid[1] > points[-1][1]:
        raise SystemExit("the well does not have its minimum at the origin")
    text = svg(points)
    if "placeholder" in text.lower():
        raise SystemExit("the figure is a placeholder")
    if text.count(",") < 16:
        raise SystemExit("the figure has too few points")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    points = samples()
    check(points)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(svg(points), encoding="utf-8")
    if args.check and args.out is None:
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
