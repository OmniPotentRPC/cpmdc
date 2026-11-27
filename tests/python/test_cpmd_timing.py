"""bench/cpmd_timing.py: the cost per evaluation between two CPMD timing tables."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("cpmd_timing", ROOT / "bench" / "cpmd_timing.py")
timing = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(timing)


def table(rows: list[tuple[str, int, float]]) -> str:
    lines = [" " + "*" * 80, " *" + " " * 35 + " TIMING " + " " * 35 + "*", " " + "*" * 80,
             " SUBROUTINE                       CALLS             SELF TIME          TOTAL TIME",
             "                                            AVERAGE   MAXIMUM   AVERAGE   MAXIMUM"]
    lines += [f" {name:<30}{calls:>8}{0.0:>12.2f}{0.0:>10.2f}{total:>10.2f}{total + 0.01:>10.2f}"
              for name, calls, total in rows]
    return "\n".join(lines + [" " + "*" * 80, " CPMD output after the table"]) + "\n"


def output(tmp_path: Path, name: str, *tables: str) -> Path:
    path = tmp_path / name
    path.write_text("".join(tables), encoding="utf-8")
    return path


def test_the_cost_is_the_change_between_the_last_two_tables_over_the_evaluations(tmp_path):
    path = output(tmp_path, "run.out",
                  table([("embed_scf", 1, 9.0), ("invfftn", 100, 4.0)]),
                  table([("embed_scf", 3, 29.0), ("invfftn", 300, 10.0), ("vpsi", 4, 3.0)]),
                  table([("embed_scf", 7, 49.0), ("invfftn", 1300, 26.0), ("vpsi", 12, 7.0)]))
    assert timing.per_evaluation(path) == (4, {"embed_scf": (1.0, 5.0), "invfftn": (250.0, 4.0),
                                               "vpsi": (2.0, 1.0)})


def test_the_rows_follow_the_first_run_and_each_run_has_its_columns(tmp_path):
    first = output(tmp_path, "base.out", table([("embed_scf", 1, 9.0), ("vpsi", 10, 2.0), ("invfftn", 10, 1.0)]),
                   table([("embed_scf", 3, 27.0), ("vpsi", 30, 8.0), ("invfftn", 30, 5.0)]))
    other = output(tmp_path, "keep.out", table([("embed_scf", 1, 9.0), ("invfftn", 10, 1.0)]),
                   table([("embed_scf", 2, 15.0), ("invfftn", 15, 2.0)]))
    lines = timing.report([("base", first), ("keep", other)], None, 1).splitlines()
    assert lines[1] == ("# Per evaluation between the last two timing tables of base (base.out, 2 evaluations), "
                        "keep (keep.out, 1 evaluation)")
    assert lines[3:] == ["| Routine | base calls | base, s | keep calls | keep, s |",
                         "|---------+------------+---------+------------+---------|",
                         "| embed_scf | 1.0 | 9.000 | 1.0 | 6.000 |",
                         "| vpsi | 10.0 | 3.000 | 0.0 | 0.000 |"]
    named = timing.report([("keep", other)], ["invfftn"], 12).splitlines()
    assert named[3:] == ["| Routine | keep calls | keep, s |", "|---------+------------+---------|",
                         "| invfftn | 5.0 | 1.000 |"]


def test_an_output_without_two_tables_or_without_evaluations_stops(tmp_path):
    with pytest.raises(SystemExit, match="1 timing table"):
        timing.per_evaluation(output(tmp_path, "one.out", table([("embed_scf", 1, 9.0)])))
    with pytest.raises(SystemExit, match="no embed_scf calls"):
        timing.per_evaluation(output(tmp_path, "none.out", table([("vpsi", 1, 9.0)]), table([("vpsi", 2, 9.5)])))
