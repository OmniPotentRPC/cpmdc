"""Cold and warm steps of a live OpenCPMD build, one at a time.

The library keeps the orbitals of the last step after a ``Session``
closes, so a step is set up from the deck or starts from those orbitals
depending on what ran before it in the process. These cases fix each
path on purpose instead of leaving it to the order of the suite.
``CPMDC_SCF_STEPS=1`` makes the library report where each step started:
``-2`` is the deck, nothing was stored, and ``0`` is the stored orbitals.
The reference well keeps nothing, so the file is skipped there.
"""

from __future__ import annotations

import re

import pytest
from conftest import WATER, park_orbitals

from cpmdc import CPMDCError, Session

NOTHING_STORED = -2
STORED = 0
WARM_START = re.compile(r"^cpmdc_warm_start (-?\d+)", re.MULTILINE)
OUTSIDE_THE_DECK = [[c * 1e200 for c in row] for row in WATER]


@pytest.fixture(autouse=True)
def live(opencpmd, monkeypatch):
    if not opencpmd:
        pytest.skip("the reference well keeps no orbitals")
    monkeypatch.setenv("CPMDC_SCF_STEPS", "1")


def water_step(params, capfd) -> list[int]:
    capfd.readouterr()
    with Session(params) as session:
        session.energy_forces(WATER, [8, 1, 1])
    return [int(source) for source in WARM_START.findall(capfd.readouterr().err)]


def test_parked_orbitals_make_the_next_session_start_from_the_deck(params, capfd):
    park_orbitals(params)
    assert water_step(params, capfd) == [NOTHING_STORED]
    assert water_step(params, capfd) == [STORED]
    park_orbitals(params)
    assert water_step(params, capfd) == [NOTHING_STORED]


def test_the_cold_start_fixture_starts_from_the_deck(cold_start, params, capfd):
    assert water_step(params, capfd) == [NOTHING_STORED]


@pytest.mark.parametrize(
    ("path", "reason"),
    [
        ("cold", r"^CPMD stopgm during embed SCF$"),
        ("warm", r"^an atom is \*+ Angstrom from a face of the isolated cell"),
    ],
)
def test_a_geometry_beyond_the_deck_format_is_refused_on_each_path(
    params, capfd, path, reason
):
    """Written into the deck with ``3F14.6`` the coordinates read ``****``.

    A cold step writes them and CPMD stops reading them. A warm step does
    not rewrite the deck, and the recentring refuses the step first.
    """
    park_orbitals(params)
    if path == "warm":
        assert water_step(params, capfd) == [NOTHING_STORED]
    with Session(params) as session:
        with pytest.raises(CPMDCError, match=reason):
            session.energy_forces(OUTSIDE_THE_DECK, [8, 1, 1])
