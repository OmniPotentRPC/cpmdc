#!/usr/bin/env bash
# Build a double-precision GROMACS whose CP2K QM/MM interface links
# libcpmdc_gromacs instead of libcp2k, and install it under PREFIX.
# CP2K_DIR holds pkgconfig/libcp2k.pc: BUILD/tools/gromacs for a meson
# build tree, or LIBDIR/cpmdc/gromacs under a meson install prefix. The
# directory goes first on PKG_CONFIG_PATH, so it wins over a real CP2K.
#
# Usage: build_gromacs_cpmdc.sh VERSION PREFIX CP2K_DIR [thread-mpi|mpi]
#   VERSION is 2024.6, 2025.5 or 2026.4 (the pinned tarballs).
#   thread-mpi (default) installs gmx_d; mpi builds against the MPI of
#   mpicc, mpicxx and mpifort on PATH and installs gmx_mpi_d. That MPI must
#   be the one libcpmdc_gromacs links.
# GROMACS links the BLAS libcpmdc links, when it links one (an OpenCPMD
# build): GROMACS' bundled BLAS and LAPACK export dgemm_ and the rest, and
# would take the place of OpenCPMD's in the process.
set -euo pipefail

version=$1
prefix=$2
cp2k_dir=$(cd "$3" && pwd)
mpi=${4:-thread-mpi}

case $version in
  2024.6) sha256=7cbad81f51c71a144d646515a7249aa74940b3f68071f51410e3a9473f05b339 ;;
  2025.5) sha256=7d0c7791675d0171af71dbfb35bfbb90bc4deb08266e62fbdcdac996f1f60174 ;;
  2026.4) sha256=58eda60979b124fcf9dff8b1b63cedfbc72fdb36d99e921dc1bdc310571098b6 ;;
  *) echo "no pinned tarball for GROMACS $version" >&2; exit 2 ;;
esac

if [ ! -f "$cp2k_dir/pkgconfig/libcp2k.pc" ]; then
  echo "$cp2k_dir has no pkgconfig/libcp2k.pc" >&2
  exit 2
fi
export PKG_CONFIG_PATH="$cp2k_dir/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
libdir=$(pkg-config --variable=libdir libcp2k)
blas=$(ldd "$libdir/libcpmdc.so" | awk '$1 ~ /^lib(openblas|blas|flexiblas|mkl_rt)[.-]/ {print $3; exit}')
if [ -n "$blas" ]; then
  blas_args=(-DGMX_EXTERNAL_BLAS=ON -DGMX_EXTERNAL_LAPACK=ON
             "-DGMX_BLAS_USER=$blas" "-DGMX_LAPACK_USER=$blas")
else
  blas_args=(-DGMX_EXTERNAL_BLAS=OFF -DGMX_EXTERNAL_LAPACK=OFF)
fi

case $mpi in
  thread-mpi) mpi_args=(-DGMX_MPI=OFF); gmx=gmx_d ;;
  mpi)
    mpi_args=(-DGMX_MPI=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
              -DMPI_Fortran_COMPILER=mpifort)
    gmx=gmx_mpi_d
    ;;
  *) echo "the fourth argument is thread-mpi or mpi, not $mpi" >&2; exit 2 ;;
esac

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
tarball="$work/gromacs-$version.tar.gz"
curl -fsSL --retry 3 -o "$tarball" "https://ftp.gromacs.org/gromacs/gromacs-$version.tar.gz"
echo "$sha256  $tarball" | sha256sum -c -
tar xzf "$tarball" -C "$work"

cmake -S "$work/gromacs-$version" -B "$work/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DGMX_DOUBLE=ON \
  "${mpi_args[@]}" \
  -DGMX_GPU=OFF \
  -DGMX_FFT_LIBRARY=fftpack \
  "${blas_args[@]}" \
  -DGMX_CP2K=ON \
  -DCP2K_DIR="$cp2k_dir" \
  -DGMX_USE_PLUMED=OFF \
  -DGMX_NNPOT=OFF \
  -DGMX_HWLOC=OFF \
  -DGMX_BUILD_HELP=OFF \
  -DGMX_DEFAULT_SUFFIX=ON \
  -DBUILD_TESTING=OFF
cmake --build "$work/build" --target install -j "$(nproc 2>/dev/null || echo 2)"
"$prefix/bin/$gmx" --version | grep -E 'GROMACS version|Precision|MPI library'
resolved=$(ldd "$prefix"/lib*/libgromacs*_d.so | awk '$1 == "libcpmdc_gromacs.so" {print $3}')
if [ "$resolved" != "$libdir/libcpmdc_gromacs.so" ]; then
  echo "libgromacs resolves libcpmdc_gromacs.so to '$resolved', not $libdir" >&2
  exit 1
fi
echo "libcpmdc_gromacs.so: $resolved"
