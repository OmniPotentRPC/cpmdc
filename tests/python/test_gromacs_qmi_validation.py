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
        ({"epsilon_r": float("nan")}, ValueError, "epsilon_r must be finite"),
        ({"epsilon_r": float("inf")}, ValueError, "epsilon_r must be finite"),
        ({"epsilon_r": -2.0}, ValueError, "epsilon_r must be positive"),
        ({"epsilon_r": 0.0}, ValueError, "epsilon_r must be positive"),
        ({"epsilon_r": "2.0"}, TypeError, "epsilon_r must be a real number, got str"),
        ({"epsilon_r": True}, TypeError, "epsilon_r must be a real number, got bool"),
        ({"epsilon_r": 4.0}, ValueError, "epsilon_r needs embedding='MECHANICAL'"),
        (
            {"embedding": "NONE", "epsilon_r": 4.0},
            ValueError,
            "epsilon_r needs embedding='MECHANICAL'",
        ),
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
        ({"links": [(0, 4)]}, ValueError, "links are 1-based atom numbers"),
        ({"links": [(1, -4)]}, ValueError, "links are 1-based atom numbers"),
        ({"links": [(2, 2)]}, ValueError, "a link joins two atoms"),
        ({"links": [(1, 2, 3)]}, ValueError, "each link is a"),
        ({"links": [(1.5, 4)]}, ValueError, "links qm_index must be an integer"),
        ({"links": [("1", 4)]}, TypeError, "links qm_index must be an integer, got str"),
        ({"links": [(1, 4)], "link_alpha": 1.0}, ValueError, "link_alpha must exceed 1"),
        ({"links": [(1, 4)], "link_alpha": float("nan")}, ValueError, "link_alpha must be finite"),
        ({"link_alpha": 1.38}, ValueError, "link_alpha needs links"),
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


def test_links_are_link_sections_of_qmmm():
    text = qmi_input("water.bin", links=[(1, 4), (numpy_index(2), 7)], link_alpha=1.5)
    assert text.count("&LINK") == 2
    assert text.count("&END LINK") == 2
    assert text.count("ALPHA_IMOMM 1.5") == 2
    lines = [line.strip() for line in text.splitlines()]
    assert lines.index("QM_INDEX 1") < lines.index("MM_INDEX 4") < lines.index("&END LINK")
    assert "QM_INDEX 2" in lines and "MM_INDEX 7" in lines
    assert lines.index("&QMMM") < lines.index("&LINK") < lines.index("&END QMMM")
    assert lines.index("&END QMMM") < lines.index("&SUBSYS")
    assert "&QMMM" not in qmi_input("water.bin")
    assert "ALPHA_IMOMM" not in qmi_input("water.bin", links=[(1, 4)])


def numpy_index(value):
    numpy = pytest.importorskip("numpy")
    return numpy.int64(value)


def test_write_qmi_is_qmi_input_in_a_file(tmp_path):
    target = tmp_path / "water.qmi"
    assert (
        write_qmi(
            target,
            "water.bin",
            embedding="MECHANICAL",
            qm_charges=[-0.834, 0.417, 0.417],
            ewald_rcut=0.7,
            epsilon_r=78.4,
            charge=1,
        )
        == target
    )
    assert target.read_text(encoding="utf-8") == qmi_input(
        "water.bin",
        embedding="MECHANICAL",
        qm_charges=[-0.834, 0.417, 0.417],
        ewald_rcut=0.7,
        epsilon_r=78.4,
        charge=1,
    )
    assert "EPSILON_R 78.4\n" in target.read_text(encoding="utf-8")
    electrostatic = tmp_path / "water_es.qmi"
    assert (
        write_qmi(
            electrostatic,
            "water.bin",
            embedding="ELECTROSTATIC",
            cell="NONE",
            mm_width=0.025,
            epsilon_r=2.5,
        )
        == electrostatic
    )
    assert electrostatic.read_text(encoding="utf-8") == qmi_input(
        "water.bin",
        embedding="ELECTROSTATIC",
        cell="NONE",
        mm_width=0.025,
        epsilon_r=2.5,
    )
    assert "MM_WIDTH 0.025\n" in electrostatic.read_text(encoding="utf-8")
    assert "EPSILON_R 2.5\n" in electrostatic.read_text(encoding="utf-8")
    linked = write_qmi(tmp_path / "linked.qmi", "water.bin", links=[(1, 4)], link_alpha=1.4)
    assert linked.read_text(encoding="utf-8") == qmi_input("water.bin", links=[(1, 4)], link_alpha=1.4)
    with pytest.raises(ValueError, match="multiplicity must be 1 or more"):
        write_qmi(tmp_path / "bad.qmi", "water.bin", multiplicity=0)
    assert not (tmp_path / "bad.qmi").exists()


@pytest.mark.parametrize(
    ("params", "error", "fragment"),
    [
        (b"water.bin", TypeError, "params must be a path string or path-like object"),
        (bytearray(b"water.bin"), TypeError, "params must be a path string or path-like object"),
        ("", ValueError, "params must not be empty"),
    ],
)
def test_params_path_is_validated(params, error, fragment):
    with pytest.raises(error, match=fragment):
        qmi_input(params)


@pytest.mark.parametrize(
    ("path", "error", "fragment"),
    [
        (b"water.qmi", TypeError, "path must be a path string or path-like object"),
        ("", ValueError, "path must not be empty"),
    ],
)
def test_write_qmi_target_path_is_validated(tmp_path, path, error, fragment):
    with pytest.raises(error, match=fragment):
        write_qmi(path, tmp_path / "water.bin")
