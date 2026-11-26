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
set -euo pipefail

commit=062582b7cfd832d36f88f504cd08e4ead42eb404
src=$1
build=$2
tools=$(cd "$(dirname "$0")/../../tools" && pwd)

case $(uname -m) in
  x86_64) config=LINUX-X86_64-GFORTRAN-MPI-PIC ;;
  aarch64) config=LINUX-AARCH64-GFORTRAN-MPI-PIC ;;
  *) echo "no shipped OpenCPMD configuration for $(uname -m)" >&2; exit 2 ;;
esac

if [ ! -d "$src/.git" ]; then
  git clone -q https://github.com/OpenCPMD/CPMD.git "$src"
fi
git -C "$src" checkout -q "$commit"
for p in embed_rinitwf embed_geometry embed_rwfopt converged_state \
         kpoints_inputfile stopgm_return tistopgm c_mem_addrs \
         embed_teardown embed_extrapolate embed_dipole embed_charges omp_info pcg_step \
         constant_cutoff_precon atomwf_kpts_ortho external_charges; do
  patch -p1 -s -N -d "$src" < "$tools/opencpmd_$p.patch"
done
cp "$tools/$config" "$src/configure/"

if [ -n "${CONDA_PREFIX:-}" ]; then
  export LIBRARY_PATH="$CONDA_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
fi
mkdir -p "$build"
build=$(cd "$build" && pwd)
(cd "$src" && ./configure.sh -DEST="$build" "$config")
make -C "$build" -j "$(nproc 2>/dev/null || echo 2)"
make -C "$build/obj" -f "$build/Makefile" timetag.o
ls -l "$build/lib/libcpmd.a" "$build/obj/timetag.o"
