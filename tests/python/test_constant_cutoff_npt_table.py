"""bench/constant_cutoff_npt.py table: the runs that wrote each state, and the tables the explanation quotes."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("constant_cutoff_npt", ROOT / "bench" / "constant_cutoff_npt.py")
npt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(npt)

RESULTS = ROOT / "bench" / "results"
EXPLANATION = ROOT / "docs" / "orgmode" / "explanation" / "constant-cutoff-npt.org"
QUOTED = ("| Run | Basis |", "| Run | Length, ps |", "| Run | No mean |", "| Run | Ranks |")


def record(start: str, start_ps: float, volumes: list[float]) -> dict:
    return {"args": {"start": start, "start_ps": start_ps},
            "time_ps": [k * 1e-3 for k in range(len(volumes))], "Volume": volumes}


def chain() -> list[tuple[str, dict]]:
    box = record("box.xyz", 0.0, [1.0, 1.1, 1.2, 1.3])
    cont = record("box.state.npz", 0.003, [1.3, 1.25, 1.2])
    return [("box", box), ("cont", cont),
            ("left", record("cont.final.state.npz", 0.001, [1.25, 1.4])),
            ("right", record("cont.final.state.npz", 0.001, [1.25, 1.1])),
            ("other", record("elsewhere.state.npz", 0.002, [2.0, 2.1]))]


def test_each_state_belongs_to_the_run_with_its_volume_at_start_ps():
    assert npt.writers(chain()) == {"cont": "box", "left": "cont", "right": "cont"}


def test_a_run_starts_after_the_runs_that_wrote_its_state():
    runs = chain()
    starts = npt.chain_starts(runs, npt.writers(runs))
    assert starts == pytest.approx({"box": 0.0, "cont": 0.003, "left": 0.004, "right": 0.004, "other": 0.002})


def test_the_runs_from_one_state_form_a_branch_group():
    groups = npt.branch_groups(chain())
    assert [[name for name, _ in group] for group in groups] == [["left", "right"]]


def tables(path: Path) -> list[tuple[str, list[str]]]:
    """Each Org table of path, with the last line of text before it."""
    found, caption, rows = [], "", []
    for line in path.read_text().splitlines() + [""]:
        if line.startswith("|"):
            rows.append(line)
            continue
        if rows:
            found.append((caption, rows))
            rows = []
        if line.strip():
            caption = line
    return found


def test_explanation_carries_the_generated_tables():
    quoted = [rows for caption, rows in tables(RESULTS / "constant_cutoff_npt.org")
              if rows[0].startswith(QUOTED) or caption.startswith("Density after the start of the runs from")]
    quoted += [rows for _, rows in tables(RESULTS / "constant_cutoff_pressure.org")]
    assert len(quoted) == 6
    page = EXPLANATION.read_text()
    missing = [row for rows in quoted for row in rows if row not in page]
    assert not missing, missing
