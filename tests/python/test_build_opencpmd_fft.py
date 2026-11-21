"""scripts/ci/build_opencpmd.sh: the CPMDC_OPENCPMD_FFT switch."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "ci" / "build_opencpmd.sh"
CONFIGS = ("LINUX-X86_64-GFORTRAN-MPI-PIC", "LINUX-AARCH64-GFORTRAN-MPI-PIC")


def test_an_unknown_fft_stops_the_build_before_the_clone(tmp_path: Path) -> None:
    env = dict(os.environ, CPMDC_OPENCPMD_FFT="mkl")
    done = subprocess.run(["bash", str(SCRIPT), str(tmp_path / "src"), str(tmp_path / "build")],
                          env=env, capture_output=True, text=True, check=False)
    assert done.returncode == 2
    assert "CPMDC_OPENCPMD_FFT is default or fftw3, not mkl" in done.stderr
    assert not (tmp_path / "src").exists()


def test_the_shipped_configurations_carry_the_lines_that_fftw3_rewrites() -> None:
    for name in CONFIGS:
        text = (ROOT / "tools" / name).read_text(encoding="utf-8")
        assert text.count("-D__HAS_FFT_DEFAULT") == 1, name
        assert text.count("LIBS='-lopenblas'") == 2, name
