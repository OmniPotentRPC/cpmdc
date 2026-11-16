"""One job TOML feeds the params message and the file-route keywords."""

import importlib.util
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _tool():
    path = ROOT / "tools" / "job_toml.py"
    assert path.is_file()
    spec = importlib.util.spec_from_file_location("job_toml", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_si3n4_toml_feeds_both_routes():
    tool = _tool()
    job = tool.load_job(ROOT / "examples" / "si3n4-isomer1.toml")
    message = tool.capnp_text(job)
    deck = tool.file_route_text(job)
    assert 'functional = "BLYP"' in message
    assert "cutOffRy = 70.0" in message
    assert "gcCutoff = 1.0e-7" in message
    assert "newCode = true" in message
    assert 'path = "Si_MT_BLYP.psp"' in message
    assert "lmax = 2" in message
    assert "lmax = 1" in message
    assert "FUNCTIONAL BLYP" in deck
    assert "70.0" in deck
    assert "GC-CUTOFF" in deck
    assert "1.0e-7" in deck
    assert "NEWCODE" in deck
    assert "*Si_MT_BLYP.psp KLEINMAN-BYLANDER" in deck
    assert "LMAX=D LOC=D" in deck
    assert "LMAX=P LOC=P" in deck
    assert job["functional"] in message and job["functional"] in deck


def test_asin_toml_keeps_charge_and_cell_on_both_routes():
    tool = _tool()
    job = tool.load_job(ROOT / "examples" / "asin-tls.toml")
    message = tool.capnp_text(job)
    deck = tool.file_route_text(job)
    assert "charge = 1" in message
    assert "symmetry = 8" in message
    assert "9.984" in message
    assert 'grimme = "D2"' in message
    assert "CHARGE" in deck
    assert "\n1\n" in deck
    assert "O_SPRIK_BLYP.psp" in message
    assert "O_SPRIK_BLYP.psp" in deck
    assert "LMAX=S LOC=S" in deck


def test_constant_cutoff_table_feeds_both_routes():
    tool = _tool()
    job = tool.load_job(ROOT / "examples" / "si3n4-isomer1.toml")
    job["constant_cutoff"] = {"cutoff_ry": 50}
    assert tool.constant_cutoff(job) == (3500.0, 5.0, 50.0)
    assert "constantCutoff = [3500.0, 5.0, 50.0]" in tool.capnp_text(job)
    assert "CONSTANT CUTOFF\n3500 5 50\n" in tool.file_route_text(job)
    job["constant_cutoff"] = {"cutoff_ry": 50, "height_ry": 400, "width_ry": 2.5}
    assert "constantCutoff = [400.0, 2.5, 50.0]" in tool.capnp_text(job)
    del job["constant_cutoff"]
    assert "constantCutoff" not in tool.capnp_text(job)
    assert "CONSTANT CUTOFF" not in tool.file_route_text(job)


@pytest.mark.parametrize(
    "table, message",
    [
        ({"cutoff_ry": 65}, "exceeds cutoff_ry 70 Ry"),
        ({"height_ry": 100}, "needs cutoff_ry"),
        ({"cutoff_ry": 50, "sigma": 2}, "unknown constant_cutoff keys"),
        ({"cutoff_ry": -1}, "finite and positive"),
        ({"cutoff_ry": True}, "must be a number"),
    ],
)
def test_constant_cutoff_table_is_checked(tmp_path, table, message):
    tool = _tool()
    job = tomllib.loads((ROOT / "examples" / "si3n4-isomer1.toml").read_text(encoding="utf-8"))
    job["constant_cutoff"] = table
    with pytest.raises(ValueError, match=message):
        tool.constant_cutoff(job)
    lines = [f"{key} = {str(value).lower()}" for key, value in table.items()]
    path = tmp_path / "job.toml"
    text = (ROOT / "examples" / "si3n4-isomer1.toml").read_text(encoding="utf-8")
    path.write_text(text + "\n[constant_cutoff]\n" + "\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        tool.load_job(path)
