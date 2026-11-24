"""pair_style cpmdc electrostatic embedding: NULL types, empty MM set, width.

The plugin is ``CPMDC_LAMMPS_PLUGIN`` or ``build/cpmdcplugin.so``; every
case skips with the reason when LAMMPS or the plugin is missing.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest
from conftest import ROOT
from cpmdc import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Session

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lammps_case import CMDARGS

WATER = [(1.0, 1.0, 1.1173), (1.0, 1.7572, 0.5308), (1.0, 0.2428, 0.5308)]
MM_SITE = (3.0, 1.0, 1.1)


def _plugin() -> Path | None:
    named = os.environ.get("CPMDC_LAMMPS_PLUGIN")
    candidate = Path(named) if named else ROOT / "build" / "cpmdcplugin.so"
    return candidate if candidate.is_file() else None


def _lammps_problem() -> str:
    try:
        import lammps

        lammps.lammps(cmdargs=list(CMDARGS)).close()
    except Exception as exc:
        return f"LAMMPS does not start: {exc}"
    return ""


PLUGIN = _plugin()
PROBLEM = _lammps_problem()
pytestmark = [
    pytest.mark.skipif(PLUGIN is None, reason="cpmdcplugin.so is not built"),
    pytest.mark.skipif(bool(PROBLEM), reason=PROBLEM or "LAMMPS"),
]


def _system(params: Path, *, mm: bool, width: str = ""):
    from lammps import lammps

    lmp = lammps(cmdargs=list(CMDARGS))
    lines = [
        "units metal",
        "atom_style charge",
        "atom_modify map array sort 0 0",
        "boundary f f f",
        "region box block -20 20 -20 20 -20 20",
        "create_box 3 box",
        "mass 1 15.999",
        "mass 2 1.008",
        "mass 3 1.0",
        f"plugin load {PLUGIN}",
    ]
    for kind, (x, y, z) in zip((1, 2, 2), WATER):
        lines.append(f"create_atoms {kind} single {x} {y} {z}")
    if mm:
        x, y, z = MM_SITE
        lines += [f"create_atoms 3 single {x} {y} {z}", "set type 3 charge 1.0"]
    lines += [
        f"pair_style cpmdc {params}{width}",
        "pair_coeff * * O H NULL",
        "thermo_style custom step pe",
        "thermo_modify format float %.15g",
    ]
    lmp.commands_string("\n".join(lines))
    return lmp


def _energy(lmp) -> float:
    lmp.command("run 0")
    return float(lmp.get_thermo("pe"))


def test_mm_charge_energy_and_forces_match_the_session(embedding_params):
    lmp = _system(embedding_params, mm=False)
    bare = _energy(lmp)
    lmp.close()

    lmp = _system(embedding_params, mm=True)
    embedded = _energy(lmp)
    forces = np.array(lmp.numpy.extract_atom("f")[:4], dtype=float)
    lmp.close()

    with Session(embedding_params) as session:
        session.set_external_charges(np.array(MM_SITE), [1.0], [0.30])
        result = session.energy_forces(np.array(WATER), [8, 1, 1])
        mm_force_au = np.array(session.external_charge_forces_au(1))

    assert embedded != pytest.approx(bare, abs=1e-6)
    assert embedded == pytest.approx(result.energy_h * HARTREE_TO_EV, rel=1e-9)
    np.testing.assert_allclose(
        forces[3], mm_force_au * HARTREE_TO_EV / BOHR_TO_ANGSTROM, rtol=1e-9, atol=1e-12
    )
    assert np.abs(forces[3]).max() > 1e-6


def test_deleting_every_mm_atom_restores_the_bare_energy(embedding_params):
    lmp = _system(embedding_params, mm=False)
    bare = _energy(lmp)
    lmp.close()

    lmp = _system(embedding_params, mm=True)
    embedded = _energy(lmp)
    lmp.commands_string("group mm type 3\ndelete_atoms group mm")
    after = _energy(lmp)
    lmp.close()

    assert embedded != pytest.approx(bare, abs=1e-6)
    assert after == pytest.approx(bare, rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("width", ["0", "-0.1", "nan", "inf", "1e-13"])
def test_unusable_width_is_refused(embedding_params, width):
    with pytest.raises(Exception, match="width must be finite"):
        _system(embedding_params, mm=True, width=f" width {width}")
