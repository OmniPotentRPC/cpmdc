"""Fixtures for the pytest suite of python/cpmdc.

``CPMDC_TEST_PARAMS`` and ``CPMDC_TEST_PERIODIC_PARAMS`` name the encoded
cluster and periodic water messages. Without them the text messages in
``tests/`` are encoded with ``capnp``. ``CPMDC_TEST_INSTALLED=1`` tests the
importable package as installed; otherwise ``python/`` of this checkout is
imported and ``build/libcpmdc.so`` is used when ``CPMDC_LIBRARY`` is unset.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
INSTALLED = os.environ.get("CPMDC_TEST_INSTALLED") == "1"

if not INSTALLED:
    sys.path.insert(0, str(ROOT / "python"))
    if not os.environ.get("CPMDC_LIBRARY"):
        for name in ("libcpmdc.so", "libcpmdc.dylib"):
            candidate = ROOT / "build" / name
            if candidate.is_file():
                os.environ["CPMDC_LIBRARY"] = str(candidate)
                break

WATER = [
    [0.0, 0.0, 0.1173],
    [0.0, 0.7572, -0.4692],
    [0.0, -0.7572, -0.4692],
]

# No test evaluates this molecule, so the orbitals parked on it match no test.
PARKING = ([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]], [1, 1])


def _params(env: str, source: str, tmp: pytest.TempPathFactory) -> Path:
    named = os.environ.get(env)
    if named:
        return Path(named)
    capnp = shutil.which("capnp")
    if capnp is None:
        pytest.skip(f"set {env} or put capnp on PATH")
    target = tmp.mktemp("params") / (Path(source).stem + ".bin")
    done = subprocess.run(
        [capnp, "encode", str(ROOT / "schema" / "Potentials.capnp"), "CPMDParams"],
        input=(ROOT / "tests" / source).read_bytes(),
        stdout=subprocess.PIPE,
        check=True,
    )
    target.write_bytes(done.stdout)
    return target


@pytest.fixture(scope="session")
def params(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Cluster water: SYMMETRY 0 and the Hockney Poisson solver."""
    return _params("CPMDC_TEST_PARAMS", "water_blyp_params.capnp.txt", tmp_path_factory)


@pytest.fixture(scope="session")
def periodic_params(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Periodic water: SYMMETRY 1, so a cell gives a stress tensor."""
    return _params(
        "CPMDC_TEST_PERIODIC_PARAMS", "water_blyp_periodic_params.capnp.txt", tmp_path_factory
    )


@pytest.fixture(scope="session")
def opencpmd(params: Path) -> bool:
    """True when the library is linked to OpenCPMD, not the reference well."""
    from cpmdc import Session

    with Session(params) as session:
        return session.available()


def park_orbitals(params: Path) -> None:
    """Drop the orbitals OpenCPMD keeps for the process, on a live build.

    The library keeps them after a ``Session`` closes and reuses them for
    any later session of the same elements, cell, and deck, so the first
    step of a test starts warm or cold depending on what ran before it.
    Orbitals of another composition replace them, and the next step of
    water then sets CPMD up again from its deck. The reference well keeps
    nothing, so this does nothing there.
    """
    from cpmdc import Session

    with Session(params) as session:
        if session.available():
            session.energy_forces(*PARKING)


@pytest.fixture
def cold_start(params: Path) -> None:
    """The first step of the test sets CPMD up from its deck."""
    park_orbitals(params)
