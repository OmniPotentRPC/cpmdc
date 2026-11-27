`scripts/ci/build_opencpmd.sh` builds OpenCPMD on FFTW 3 when
`CPMDC_OPENCPMD_FFT=fftw3`: it configures from a copy of the shipped
configuration that defines `__HAS_FFT_FFTW3` in place of `__HAS_FFT_DEFAULT`
and adds `-lfftw3` to `LIBS`, stops when the shipped file loses either line,
and rejects any other value before it clones.
`howto/opencpmd-archive.org` (section Build on FFTW) describes that build.
`tests/python/test_build_opencpmd_fft.py` checks the switch and the lines it
rewrites.
