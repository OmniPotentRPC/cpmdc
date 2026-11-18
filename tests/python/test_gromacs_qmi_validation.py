"""qmi_input refuses values that would write a file the library reads differently."""

from __future__ import annotations

import pytest

from cpmdc.gromacs import qmi_input, write_qmi


@pytest.mark.parametrize(
    ("options", "error", "fragment"),
    [
        ({"ewald_rcut": float("nan")}, ValueError, "ewald_rcut must be finite"),
        ({"ewald_rcut": float("inf")}, ValueError, "ewald_rcut must be finite"),
        ({"ewald_rcut": -0.5}, ValueError, "ewald_rcut must be positive"),
        ({"ewald_rcut": 0.0}, ValueError, "ewald_rcut must be positive"),
        ({"ewald_precision": -1.0e-8}, ValueError, "ewald_precision must be positive"),
        ({"ewald_rcut": "0.7"}, TypeError, "ewald_rcut must be a real number, got str"),
        ({"ewald_rcut": True}, TypeError, "ewald_rcut must be a real number, got bool"),
        ({"charge": 1.5}, ValueError, "charge must be an integer"),
        ({"charge": "1"}, TypeError, "charge must be an integer, got str"),
        ({"multiplicity": 0}, ValueError, "multiplicity must be 1 or more"),
        ({"multiplicity": 2.5}, ValueError, "multiplicity must be an integer"),
        (
            {"embedding": "MECHANICAL", "qm_charges": [0.0, float("nan")]},
            ValueError,
            "qm_charges must be finite",
        ),
    ],
)
def test_bad_options_are_refused(options, error, fragment):
    with pytest.raises(error, match=fragment):
        qmi_input("water.bin", **options)


def test_whole_floats_and_numpy_integers_are_charges():
    numpy = pytest.importorskip("numpy")
    text = qmi_input("water.bin", charge=numpy.int64(-1), multiplicity=2.0)
    assert "CHARGE -1" in text
    assert "MULTIPLICITY 2\n" in text


def test_write_qmi_is_qmi_input_in_a_file(tmp_path):
    target = tmp_path / "water.qmi"
    assert write_qmi(target, "water.bin", ewald_rcut=0.7, charge=1) == target
    assert target.read_text(encoding="utf-8") == qmi_input("water.bin", ewald_rcut=0.7, charge=1)
    with pytest.raises(ValueError, match="multiplicity must be 1 or more"):
        write_qmi(tmp_path / "bad.qmi", "water.bin", multiplicity=0)
    assert not (tmp_path / "bad.qmi").exists()
