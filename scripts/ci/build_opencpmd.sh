#!/usr/bin/env bash
# Clone OpenCPMD at the commit the patches in tools/ apply to, apply them in
# the order of tools/README.md, and build lib/libcpmd.a and obj/timetag.o
# under BUILD with the shipped -fPIC configuration of this machine. BUILD
# is the cpmd_root of meson setup -Dwith_cpmd=true.
#
# Usage: build_opencpmd.sh SRC BUILD
#   SRC is the clone (created when missing); BUILD is the configure.sh DEST.
#   mpif90 and gcc come from PATH; LIBS is -lopenblas, found on
#   LIBRARY_PATH, which gains $CONDA_PREFIX/lib inside a pixi environment.
#   CPMDC_OPENCPMD_FFT=fftw3 builds on FFTW 3 instead of the FFT shipped
#   with CPMD, and adds -lfftw3 to LIBS.
set -euo pipefail

commit=062582b7cfd832d36f88f504cd08e4ead42eb404
src=$1
build=$2
tools=$(cd "$(dirname "$0")/../../tools" && pwd)
fft=${CPMDC_OPENCPMD_FFT:-default}

case $(uname -m) in
  x86_64) config=LINUX-X86_64-GFORTRAN-MPI-PIC ;;
  aarch64) config=LINUX-AARCH64-GFORTRAN-MPI-PIC ;;
  *) echo "no shipped OpenCPMD configuration for $(uname -m)" >&2; exit 2 ;;
esac
case $fft in
  default|fftw3) ;;
  *) echo "CPMDC_OPENCPMD_FFT is default or fftw3, not $fft" >&2; exit 2 ;;
esac

if [ ! -d "$src/.git" ]; then
  git clone -q https://github.com/OpenCPMD/CPMD.git "$src"
fi
git -C "$src" checkout -q "$commit"
for p in embed_rinitwf embed_geometry embed_rwfopt converged_state \
         kpoints_inputfile stopgm_return tistopgm c_mem_addrs \
         embed_teardown embed_extrapolate omp_info pcg_step \
         constant_cutoff_precon; do
  patch -p1 -s -N -d "$src" < "$tools/opencpmd_$p.patch"
done
cp "$tools/$config" "$src/configure/"
if [ "$fft" = fftw3 ]; then
  sed -e 's/-D__HAS_FFT_DEFAULT/-D__HAS_FFT_FFTW3/' -e "s/LIBS='-lopenblas'/LIBS='-lopenblas -lfftw3'/" \
    "$tools/$config" > "$src/configure/$config-FFTW3"
  if grep -q __HAS_FFT_DEFAULT "$src/configure/$config-FFTW3" ||
     [ "$(grep -c -- "-lfftw3'" "$src/configure/$config-FFTW3")" -ne 2 ]; then
    echo "$config no longer has the FFT lines that fftw3 replaces" >&2
    exit 2
  fi
  config=$config-FFTW3
fi

if [ -n "${CONDA_PREFIX:-}" ]; then
  export LIBRARY_PATH="$CONDA_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
fi
mkdir -p "$build"
build=$(cd "$build" && pwd)
(cd "$src" && ./configure.sh -DEST="$build" "$config")
make -C "$build" -j "$(nproc 2>/dev/null || echo 2)"
make -C "$build/obj" -f "$build/Makefile" timetag.o
ls -l "$build/lib/libcpmd.a" "$build/obj/timetag.o"
