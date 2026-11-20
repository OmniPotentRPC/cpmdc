#!/usr/bin/env bash
# Install the rgpot 3.7.0 LAMMPS frontend into the selected Python.
set -euo pipefail

pin=d2418b01502cf8e904d414d854e912b23c7aa6d5
python=${PYTHON:-python3}

if "$python" -c 'import rgpot.lammps' >/dev/null 2>&1; then
  "$python" -c 'import rgpot.lammps as m; print("rgpot 3.7.0", m.__file__)'
  exit 0
fi

"$python" -m pip install 'nanobind>=2.2' 'numpy>=1.24' 'meson>=1.5.2,<2' 'ninja>=1.11'

if [[ -z ${CXX:-} ]] && command -v g++ >/dev/null 2>&1; then
  export CXX=g++
fi
if [[ -z ${CC:-} ]] && command -v gcc >/dev/null 2>&1; then
  export CC=gcc
fi

work=$(mktemp -d)
trap '"$python" -c '"'"'import shutil, sys; shutil.rmtree(sys.argv[1])'"'"' "$work"' EXIT
git clone --filter=blob:none --no-checkout https://github.com/OmniPotentRPC/rgpot.git "$work/src"
git -C "$work/src" fetch --depth 1 origin "$pin"
git -C "$work/src" checkout --detach FETCH_HEAD
test "$(git -C "$work/src" rev-parse HEAD)" = "$pin"

bindir=$(dirname "$("$python" -c 'import sys; print(sys.executable)')")
userbin=$("$python" -c 'import site; print(site.USER_BASE + "/bin")')
prefix=$(
  "$python" -c 'import site, sys; print(site.USER_BASE if sys.prefix == sys.base_prefix else sys.prefix)'
)
export PATH="${bindir}:${userbin}:${PATH}"
meson=$(command -v meson)
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

"$python" -c 'import rgpot.lammps as m; print("rgpot 3.7.0", m.__file__)'
