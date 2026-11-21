#!/usr/bin/env bash
# Install rgpot.lammps into the selected Python.
#
# The frontend lives on OmniPotentRPC/rgpot at the pinned commit below
# (branch feat/lammps-fix-external-v3). No released wheel ships it, and a
# pip install of that commit builds the full project, including torch.
# This script builds the Python extension with that frontend and without
# the torch, Fortran, RPC, and Highway pieces.
#
# Usage: install_rgpot_lammps.sh
#        PYTHON=/path/to/python install_rgpot_lammps.sh
set -euo pipefail

pin=d2418b01502cf8e904d414d854e912b23c7aa6d5
python=${PYTHON:-python3}

if "$python" -c 'import rgpot.lammps' >/dev/null 2>&1; then
  "$python" -c 'import rgpot.lammps as m; print("rgpot.lammps", m.__file__)'
  exit 0
fi

if ! "$python" -m pip --version >/dev/null 2>&1; then
  "$python" -m ensurepip --upgrade \
    || curl -fsSL https://bootstrap.pypa.io/get-pip.py | "$python"
fi

"$python" -m pip install 'nanobind>=2.2' 'numpy>=1.24' 'meson>=1.5.2,<2' 'ninja>=1.11'

if [[ -z ${CXX:-} ]] && command -v g++ >/dev/null 2>&1; then
  export CXX=g++
fi
if [[ -z ${CC:-} ]] && command -v gcc >/dev/null 2>&1; then
  export CC=gcc
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --filter=blob:none --no-checkout https://github.com/OmniPotentRPC/rgpot.git "$work/src"
git -C "$work/src" fetch --depth 1 origin "$pin"
git -C "$work/src" checkout --detach FETCH_HEAD
test "$(git -C "$work/src" rev-parse HEAD)" = "$pin"

bindir=$(dirname "$("$python" -c 'import sys; print(sys.executable)')")
prefix=$("$python" -c 'import sys; print(sys.prefix)')
meson=${bindir}/meson
export PATH="${bindir}:${PATH}"
"$meson" setup "$work/build" "$work/src" \
  -Dwith_python=true \
  -Dwith_tests=false \
  -Dwith_examples=false \
  -Dwith_fortran_pots=disabled \
  -Dwith_highway=disabled \
  -Dwith_ipi=disabled \
  -Dwith_metatomic=false \
  -Dwith_expr=false \
  -Dwith_rpc=false \
  -Dwith_fast_nlist=disabled \
  -Dwith_dftd3=false \
  -Dwith_dftd4=false \
  -Dwith_xckernel=false \
  -Dpython.bytecompile=-1 \
  --prefix "$prefix"
"$meson" compile -C "$work/build"
"$meson" install -C "$work/build"

# _core NEEDs librgpot.so.<soversion>. Its RUNPATH searches $ORIGIN, the
# directory of the extension, which is not the prefix libdir. Match the
# rgpot package only: other wheels ship a _core extension too.
core=$(find "$prefix" -path '*/rgpot/_core*.so' -print -quit)
real=$(find "$prefix" -name 'librgpot.so.*.*' -print -quit)
test -n "$core"
test -n "$real"
dest=$(dirname "$core")
cp -a "$real" "$dest/"
soname=$(readelf -d "$real" | awk '/SONAME/ {gsub(/[][]/, "", $NF); print $NF; exit}')
test -n "$soname"
ln -sfn "$(basename "$real")" "$dest/$soname"
ln -sfn "$soname" "$dest/librgpot.so"

"$python" -c 'import rgpot.lammps as m; print("rgpot.lammps", m.__file__)'
