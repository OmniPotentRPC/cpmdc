"""An installed wheel uses its own package, library, and schema."""

from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms
from conftest import INSTALLED, ROOT, WATER

import cpmdc

pytestmark = pytest.mark.skipif(not INSTALLED, reason="CPMDC_TEST_INSTALLED=1 is not set")


def test_package_comes_from_site_packages():
    package = Path(cpmdc.__file__).resolve().parent
    assert not package.is_relative_to(ROOT)
    assert cpmdc.find_library().resolve().parent == package.parent / ".cpmdc.mesonpy.libs"
    assert cpmdc.find_schema().resolve() == package / "schema" / "Potentials.capnp"


def test_wheel_library_is_the_reference_build(params):
    with cpmdc.Session(params) as session:
        assert session.version().startswith("cpmdc/")
        assert not session.available()


def test_calculator_runs_from_the_wheel(params, periodic_params):
    from cpmdc.ase import CPMDCCalculator

    cluster = Atoms("OHH", positions=WATER)
    with CPMDCCalculator(params) as calc:
        cluster.calc = calc
        assert cluster.get_forces().shape == (3, 3)
    periodic = Atoms("OHH", positions=WATER, cell=[10.0, 10.0, 10.0], pbc=True)
    with CPMDCCalculator(periodic_params) as calc:
        periodic.calc = calc
        assert periodic.get_stress().shape == (6,)
