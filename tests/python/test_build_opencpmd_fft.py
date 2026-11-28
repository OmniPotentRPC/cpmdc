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
    stderr = done.stderr.strip()
    earlier = f"no shipped OpenCPMD configuration for {os.uname().machine}"
    unknown = "CPMDC_OPENCPMD_FFT is default or fftw3, not mkl"
    assert done.returncode == 2
    assert stderr == earlier or stderr == unknown
    assert not (tmp_path / "src").exists()


def test_the_shipped_configurations_carry_the_lines_that_fftw3_rewrites() -> None:
    script = SCRIPT.read_text(encoding="utf-8")
    fft = "s/-D__HAS_FFT_DEFAULT/-D__HAS_FFT_FFTW3/"
    libs = "s/LIBS='-lopenblas'/LIBS='-lopenblas -lfftw3'/"
    assert fft in script
    assert libs in script
    for name in CONFIGS:
        source = ROOT / "tools" / name
        text = source.read_text(encoding="utf-8")
        assert text.count("-D__HAS_FFT_DEFAULT") == 1, name
        assert text.count("LIBS='-lopenblas'") == 2, name
        rewritten = subprocess.run(
            ["sed", "-e", fft, "-e", libs, str(source)],
            check=True, capture_output=True, text=True,
        ).stdout
        assert "__HAS_FFT_DEFAULT" not in rewritten, name
        assert rewritten.count("-D__HAS_FFT_FFTW3") == 1, name
        assert rewritten.count("-lfftw3'") == 2, name
