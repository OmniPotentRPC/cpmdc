"""The PBE page's core radii and projector cutoffs are the recorded headers."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "bench" / "results" / "pbe_pseudopotential_headers.org"


def test_the_header_record_names_the_sources_and_the_checksums():
    text = RECORD.read_text(encoding="utf-8")
    assert "https://github.com/OpenCPMD/pseudopotentials-extended-library/raw/97f98cb25cbee1cce67302addd2bf06e21a17978/pseudo.tar.gz" in text
    assert "https://pseudopotentials.quantum-espresso.org/upf_files/O.pbe-mt_fhi.UPF" in text
    assert "d26bb0bfa19ddaabe95d21169f50fa25915631e84df0d157e93bb32a7a9298ad" in text
    assert "8ff874f26491f2b6e482540342ae0ad85c9c7023936f3b050baa4ddced459988" in text
    assert "747670be5d6f008d921dc74551dc0705f2c4b34132e79590c482b532171f8cd1" in text
    assert "b6ed642482602f477eb8c167dce4a2c6d254b89522685b5b9ad8f1e904c50db4" in text
    assert "2    S    1.1200" in text
    assert "1    S     .5000" in text
    assert 'cutoff_radius="1.434117653854000E+000"' in text
    assert 'cutoff_radius="1.276346260702400E+000"' in text
