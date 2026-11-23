"""Fixtures for the pytest suite of python/cpmdc.

``CPMDC_TEST_PARAMS`` and ``CPMDC_TEST_PERIODIC_PARAMS`` name the encoded
cluster and periodic water messages. Without them the text messages in
``tests/`` are encoded with ``capnp``. ``CPMDC_TEST_INSTALLED=1`` tests the
importable package as installed; otherwise ``python/`` of this checkout is
imported and ``build/libcpmdc.so`` is used when ``CPMDC_LIBRARY`` is unset.
"""

from __future__ import annotations

import inspect
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import blas_probe

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


def library_is_opencpmd() -> bool:
    """True when this process's libcpmdc is the OpenCPMD archive.

    ``cpmdc_available()`` is 1 only for that archive. The reference well
    returns 0. :func:`opencpmd` reports the same bit.
    ``pytest_pycollect_makeitem`` omits a ``cpmd_backend`` test whose name
    is the other one, and ``pytest_generate_tests`` parametrizes the test
    that remains.
    """
    from cpmdc import _load, find_library

    lib = _load(find_library())
    return bool(lib.cpmdc_available())


_LINKED_BACKEND: str | None = None


def linked_backend() -> str:
    """``opencpmd`` or ``reference``, from :func:`library_is_opencpmd`."""
    global _LINKED_BACKEND
    if _LINKED_BACKEND is None:
        _LINKED_BACKEND = "opencpmd" if library_is_opencpmd() else "reference"
    return _LINKED_BACKEND


def _cpmd_backend_mark(obj: object):
    func = getattr(obj, "__func__", obj)
    if not inspect.isfunction(func):
        return None
    raw = getattr(func, "pytestmark", [])
    marks = raw if isinstance(raw, list) else [raw]
    found = None
    for decorator in marks:
        mark = getattr(decorator, "mark", decorator)
        if getattr(mark, "name", None) != "cpmd_backend":
            continue
        if found is not None:
            raise pytest.UsageError(f"{func.__name__}: cpmd_backend is set twice")
        found = mark
    return found


def _checked_backend(name: str, mark: pytest.Mark) -> str:
    if len(mark.args) != 1 or mark.args[0] not in ("reference", "opencpmd"):
        raise pytest.UsageError(
            f"{name}: cpmd_backend takes reference or opencpmd"
        )
    reason = mark.kwargs.get("reason")
    if not isinstance(reason, str) or not reason.strip():
        raise pytest.UsageError(f"{name}: cpmd_backend needs a reason")
    if set(mark.kwargs) - {"reason"}:
        raise pytest.UsageError(f"{name}: cpmd_backend takes only reason=")
    return str(mark.args[0])


def pytest_sessionstart(session: pytest.Session) -> None:
    """Stop before collection when ``@`` returns wrong numbers (blas_probe.py)."""
    problem = blas_probe.environment_problem()
    if problem is not None:
        pytest.exit(problem, returncode=3)


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "cpmd_backend(name, *, reason): collect this test only when "
        "linked_backend() is name (reference or opencpmd); reason says why "
        "the other library does not run it",
    )


def pytest_pycollect_makeitem(
    collector: pytest.Collector, name: str, obj: object
):
    """Omit a cpmd_backend test the linked library does not run.

    Returning an empty list means the test is not collected. That is not a
    skip and not a deselection. The matching tests fall through to pytest's
    own collector, and ``pytest_generate_tests`` gives them the fixture.
    """
    del collector
    mark = _cpmd_backend_mark(obj)
    if mark is None:
        return None
    if _checked_backend(name, mark) != linked_backend():
        return []
    return None


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    mark = metafunc.definition.get_closest_marker("cpmd_backend")
    if mark is None:
        return
    wanted = _checked_backend(metafunc.definition.nodeid, mark)
    metafunc.parametrize("cpmd_backend", [wanted], indirect=True, ids=[wanted])


@pytest.fixture
def cpmd_backend(request: pytest.FixtureRequest, opencpmd: bool) -> str:
    """The backend this test was collected for.

    Collection keeps the test only when ``linked_backend()`` equals the
    marker. This fixture repeats that check with the ``opencpmd`` fixture,
    which is the session's reading of the same ``cpmdc_available()`` bit.
    """
    wanted = request.param
    active = "opencpmd" if opencpmd else "reference"
    if wanted != active:
        raise AssertionError(
            f"collected for the {wanted} backend, but the opencpmd fixture "
            f"says {active}"
        )
    return active


@pytest.fixture(scope="session")
def opencpmd(params: Path) -> bool:
    """True when the library is linked to OpenCPMD, not the reference well."""
    from cpmdc import Session

    with Session(params) as session:
        available = bool(session.available())
    if available != library_is_opencpmd():
        raise RuntimeError(
            "cpmdc_available() changed between collection and the opencpmd fixture"
        )
    return available


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
