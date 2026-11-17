"""bench/constant_cutoff_npt.py table: the runs that wrote each state, and the tables the explanations quote."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("constant_cutoff_npt", ROOT / "bench" / "constant_cutoff_npt.py")
npt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(npt)

RESULTS = ROOT / "bench" / "results"
EXPLANATION = ROOT / "docs" / "orgmode" / "explanation" / "constant-cutoff-npt.org"
PBE_WATER = ROOT / "docs" / "orgmode" / "explanation" / "pbe-water-density.org"
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


def test_the_runs_table_names_the_pseudopotentials_only_when_the_runs_differ_in_them():
    meta = {"cutoff": 100.0, "e0": 85.0, "p_tau": 100.0, "start": "box.xyz", "start_ps": 0.0}
    ours = {"args": {**meta, "thermostat": "langevin", "barostat": "mtk"}, "time_ps": [0.0, 0.001]}
    fhi = {"args": {**ours["args"], "oxygen": "O.pbe-mt_fhi.UPF", "hydrogen": "H.pbe-mt_fhi.UPF"},
           "time_ps": [0.0, 0.001, 0.002]}
    rows = npt.run_rows([("ours", ours), ("fhi", fhi)], {})
    assert rows[0] == "| Run | Basis | Pseudopotentials | Thermostat, barostat | P_tau, fs | Start | Steps |"
    assert rows[1] == "|-----+-------+------------------+----------------------+-----------+-------+-------|"
    assert rows[2] == "| ours | 100 Ry, E0 85 Ry | O_MT_PBE.psp and H_MT_PBE.psp | Langevin, MTK | 100 | box.xyz | 1 |"
    assert rows[3].startswith("| fhi | 100 Ry, E0 85 Ry | O.pbe-mt_fhi.UPF and H.pbe-mt_fhi.UPF | Langevin, MTK |")
    plain = npt.run_rows([("ours", {"args": meta, "time_ps": [0.0]})], {})
    assert plain[:2] == ["| Run | Basis | P_tau, fs | Start | Steps |", "|-----+-------+-----------+-------+-------|"]


def test_the_runs_table_gives_the_in_place_limit_and_the_height_of_runs_that_differ_in_them():
    meta = {"cutoff": 100.0, "e0": 85.0, "height": None, "p_tau": 100.0, "start": "box.xyz", "start_ps": 0.0,
            "cell_in_place": "0.01"}
    fixed = {**meta, "cutoff": 85.0, "e0": 50.0, "height": 1e-6, "cell_in_place": "0.25"}
    rows = npt.run_rows([("cc", {"args": meta, "time_ps": [0.0]}), ("fixed", {"args": fixed, "time_ps": [0.0]})], {})
    assert rows[0] == "| Run | Basis | CPMDC_CELL_IN_PLACE | P_tau, fs | Start | Steps |"
    assert rows[2] == "| cc | 100 Ry, E0 85 Ry | 0.01 | 100 | box.xyz | 0 |"
    assert rows[3] == "| fixed | 85 Ry, E0 50 Ry, A 1e-06 Ry | 0.25 | 100 | box.xyz | 0 |"


def test_the_plateau_table_gives_the_level_or_the_coarsest_level_tried():
    curve = [(64, 1.0, 0.1), (32, 2.0, 0.2), (16, 2.1, 0.3), (8, 2.2, 0.4)]
    terms = {"Volume": {"plateau": None, "block_frames": 8, "curve": curve},
             "Density": {"plateau": 2, "block_frames": 8, "curve": curve},
             "Pressure": {"plateau": 0, "block_frames": 8, "curve": curve}}
    rows = npt.plateau_rows([("run", {"frame_ps": 0.001, "terms": terms})])
    assert rows[:2] == ["| Run | Observable | Block, ps | Blocks | Error | Its uncertainty |",
                        "|-----+------------+-----------+--------+-------+-----------------|"]
    assert rows[2:] == ["| run | Volume, nm^3 | none up to 0.008 |  |  |  |",
                        "| run | Density, kg/m^3 | 0.004 | 16 | 2.10 | 0.30 |",
                        "| run | Pressure, bar | 0.001 | 64 | 1.0 | 0.1 |"]


def test_cost_compares_each_run_with_the_run_above_it(tmp_path, capsys):
    meta = {"waters": 32, "start": "box.xyz", "start_ps": 0.0, "cutoff": 100.0, "e0": 85.0, "height": None,
            "width": 4.25, "conv": 1e-5, "dt": 1.0, "temperature": 350.0, "pressure": 1.0, "t_tau": 100.0,
            "p_tau": 500.0, "hydrogen_mass": 2.014, "seed": 11, "density_range": [500.0, 1500.0], "ranks": 8,
            "cpu": "cpu", "build": "build-cpmd", "rs_wfn_keep": False}
    runs = (("base", {}, [9.0, 3.0, 5.0], [-1.0, -2.0, -3.0]),
            ("fast", {"build": "build-fftw", "rs_wfn_keep": True}, [8.0, 2.0, 3.0],
             [-1.0, -2.0, -3.5]),
            ("order6", {"aspc_order": 6}, [9.0, 3.0, 7.0], [-1.0, -2.25, -3.0]),
            ("slow", {"p_tau": 100.0}, [9.0, 3.0, 5.0], [-1.0, -2.0, -3.0]))
    paths = []
    for name, extra, wall, energy in runs:
        paths.append(tmp_path / f"{name}.json")
        paths[-1].write_text(json.dumps({"args": {**meta, **extra}, "step_wall_s": wall, "scf_steps": [20, 9, 7],
                                         "Energy": energy}))

    def cost(paths, **options):
        assert npt.cost(argparse.Namespace(**{"runs": paths, "from_step": 1, "to_step": None, "against": "above",
                                              "output": None, **options})) == 0
        return capsys.readouterr().out.splitlines()

    lines = cost(paths[:2])
    assert "; ASPC order 4, ODIIS vectors default, HAMILTONIAN CUTOFF default; 8 rank(s) on cpu;" in lines[1]
    assert lines[3] == ("| Run | Build | REAL SPACE WFN KEEP | Steps | Wall per step, s | "
                        "SCF iterations per step | Largest energy difference from the run above, eV |")
    assert lines[5:] == ["| base | build-cpmd | no | 2 | 4.000 | 8.00 |  |",
                         "| fast | build-fftw | yes | 2 | 2.500 | 8.00 | 5.0e-01 |"]
    lines = cost([paths[0], paths[2], paths[0]], to_step=1, against="first")
    assert lines[3].startswith("| Run | ASPC order | Steps |") and lines[3].endswith(" from base, eV |")
    assert lines[5:] == ["| base | 4 | 1 | 3.000 | 9.00 |  |", "| order6 | 6 | 1 | 3.000 | 9.00 | 2.5e-01 |",
                         "| base | 4 | 1 | 3.000 | 9.00 | 0.0e+00 |"]
    with pytest.raises(SystemExit, match="the runs differ in p_tau"):
        npt.cost(argparse.Namespace(runs=paths[::3], from_step=1, to_step=None, against="above", output=None))


def pressure_table(start: str, oxygen: str, rows: list[str]) -> str:
    return "\n".join(["# Generated by bench/constant_cutoff_npt.py pressure; do not edit by hand.",
                      f"# 32 waters from {start}, 1100.31 Angstrom^3, 870.0 kg/m^3; PBE, {oxygen} and H.upf, "
                      "orbital convergence 1e-06; 8 rank(s) on cpu", "",
                      "| Cutoff, Ry | E0, Ry | A, Ry | sigma, Ry | Energy, eV | Atomic pressure, kbar | "
                      "Molecular pressure, kbar | Wall, s |",
                      "|------------+--------+-------+-----------+------------+-----------------------+"
                      "--------------------------+---------|", *rows]) + "\n"


def test_pressures_takes_the_differences_configuration_by_configuration(tmp_path, capsys):
    paths = []
    for start, fixed, constant, other in (("a.state.npz", -1.0, 1.5, 0.5), ("b.state.npz", -3.0, -0.3, -1.5)):
        cc = "| 100 | 85 | 5000 | 4.25 | -1.0 | 2.00 |"
        one = tmp_path / f"{start}-one.org"
        one.write_text(pressure_table(start, "O.upf", [f"| 85 |  |  |  | -1.0 | 2.00 | {fixed:.2f} | 1.0 |",
                                                       f"{cc} {constant:.2f} | 1.0 |"]))
        two = tmp_path / f"{start}-two.org"
        two.write_text(pressure_table(start, "O.psp", [f"{cc} {other:.2f} | 1.0 |"]))
        paths += [one, two]
    npt.pressures(argparse.Namespace(tables=paths, against=2, output=None))
    lines = capsys.readouterr().out.splitlines()
    assert lines[1].endswith("the differences are from O.upf, H.upf, 100/85/4.25, configuration by configuration")
    assert lines[3] == ("| Configuration | Density, kg/m^3 | O.upf, H.upf, 85, kbar | "
                        "O.upf, H.upf, 100/85/4.25, kbar | O.psp, H.upf, 100/85/4.25, kbar |")
    assert lines[5:7] == ["| a.state.npz | 870.0 | -1.00 | 1.50 | 0.50 |",
                          "| b.state.npz | 870.0 | -3.00 | -0.30 | -1.50 |"]
    assert lines[10:] == ["| O.upf, H.upf, 85 | -2.000 | 1.000 | -2.600 | 0.100 |",
                          "| O.upf, H.upf, 100/85/4.25 | 0.600 | 0.900 |  |  |",
                          "| O.psp, H.upf, 100/85/4.25 | -0.500 | 1.000 | -1.100 | 0.100 |"]
    with pytest.raises(SystemExit, match="no table gives b.state.npz at 870.0 kg/m\\^3 with O.psp"):
        npt.pressures(argparse.Namespace(tables=paths[:3], against=1, output=None))
    with pytest.raises(SystemExit, match="not a table of constant_cutoff_npt.py pressure"):
        npt.pressures(argparse.Namespace(tables=[RESULTS / "constant_cutoff_npt.org"], against=1, output=None))


def tables(source: Path | str) -> list[tuple[str, list[str]]]:
    """Each Org table of a file or of text, with the last line of text before it."""
    text = source.read_text() if isinstance(source, Path) else source
    found, caption, rows = [], "", []
    for line in text.splitlines() + [""]:
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


def sections(path: Path) -> dict[str, str]:
    """The text under each top-level heading of an Org file, and under "" the text before the first."""
    parts = re.split(r"^\* (.+)$", path.read_text(), flags=re.M)
    return {"": parts[0], **dict(zip(parts[1::2], parts[2::2]))}


def test_the_pbe_water_page_carries_the_generated_tables():
    long_run = sections(RESULTS / "pbe_water_npt.org")
    eight = ("| Run | Length, ps |", "| Run | Observable | Detected start", "| Run | Observable | Block, ps |",
             "| Run | kappa_T", "| After the start, ps |")
    four = ("| Run | Observable | 0.001 ps |", "| Run | Observable | Block, ps |")
    quoted = [rows for _, rows in tables(long_run["At least 8 blocks"]) if rows[0].startswith(eight)]
    quoted += [rows for _, rows in tables(long_run["At least 4 blocks"]) if rows[0].startswith(four)]
    quoted += [rows for _, rows in tables(RESULTS / "pbe_water_cost.org")]
    quoted += [rows for name, text in sections(RESULTS / "pbe_water_pressure.org").items()
               if name and not name.startswith("mesh-") for _, rows in tables(text)]
    quoted += [rows for _, rows in tables(sections(RESULTS / "pbe_water_states.org")[""])]
    assert len(quoted) == 20
    page = PBE_WATER.read_text()
    missing = [row for rows in quoted for row in rows if row not in page]
    assert not missing, missing
