#!/usr/bin/env python3
"""The constant-cutoff NPT figure draws every run of bench/results and refuses bad input."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path


def draw(script: Path, data: Path, out: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(script), "--json", str(data), "--out", str(out)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_constant_cutoff_npt_plot.py SOURCE_ROOT", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    script = root / "docs" / "render_constant_cutoff_npt.py"
    data = root / "bench" / "results" / "constant_cutoff_npt.json"
    runs = json.loads(data.read_text(encoding="utf-8"))["runs"]
    with tempfile.TemporaryDirectory(prefix="cpmdc-npt-plot-") as tmp:
        out = Path(tmp) / "constant-cutoff-npt.svg"
        proc = draw(script, data, out)
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr)
            return proc.returncode
        text = out.read_text(encoding="utf-8")
        if text.count("<polyline") != len(runs):
            print(f"{text.count('<polyline')} curves for {len(runs)} runs", file=sys.stderr)
            return 1
        means = sum(1 for run in runs if run.get("mean") is not None)
        if text.count('fill-opacity="0.18"') != means:
            print(f"{text.count('fill-opacity')} bands for {means} runs with a mean", file=sys.stderr)
            return 1
        for run in runs:
            if run["label"] not in text:
                print(f"{run['label']} is not in the legend", file=sys.stderr)
                return 1
        if "Miceli" not in text:
            print("the reference density is not drawn", file=sys.stderr)
            return 1
        branches = [run for run in runs if "steps" in run]
        if bool(branches) != ("time after the state" in text):
            print(f"{len(branches)} runs from a shared state, but the panel for them is "
                  f"{'missing' if branches else 'drawn'}", file=sys.stderr)
            return 1

        first = next(dict(run) for run in runs if run.get("mean") is not None)
        bad = {
            "no runs": [],
            "one window": [{**first, "windows": first["windows"][:1]}],
            "a value that is not finite": [{**first, "windows": [[0.0, float("nan")], *first["windows"][1:]]}],
            "no mean": [{**first, "mean": None}],
        }
        if branches:
            bad["one step"] = [first, {**branches[0], "steps": branches[0]["steps"][:1]}]
            bad["a step that is not finite"] = [first, {**branches[0], "steps": [[0.0, float("inf")],
                                                                                   *branches[0]["steps"][1:]]}]
        for case, broken in bad.items():
            path = Path(tmp) / "broken.json"
            path.write_text(json.dumps({"runs": broken}), encoding="utf-8")
            if draw(script, path, Path(tmp) / "broken.svg").returncode == 0:
                print(f"the figure accepted a file with {case}", file=sys.stderr)
                return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
