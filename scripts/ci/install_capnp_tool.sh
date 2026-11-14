#!/usr/bin/env bash
# Build the Cap'n Proto compiler (capnp, capnpc-c++) from a pinned release
# and install it under PREFIX. The manylinux images carry no capnproto
# package; the cpmdc build and the wheel tests only need the command.
#
# Usage: install_capnp_tool.sh [PREFIX]   (default /usr/local)
set -euo pipefail

version=1.0.2
sha256=9057dbc0223366b74bbeca33a05de164a229b0377927f1b7ef3828cdd8cb1d7e
prefix=${1:-/usr/local}

if command -v capnp >/dev/null 2>&1; then
  capnp --version
  exit 0
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
tarball="$work/capnproto-c++-$version.tar.gz"
curl -fsSL --retry 3 -o "$tarball" "https://capnproto.org/capnproto-c++-$version.tar.gz"
echo "$sha256  $tarball" | sha256sum -c -
tar xzf "$tarball" -C "$work"

cmake -S "$work/capnproto-c++-$version" -B "$work/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF \
  -DWITH_OPENSSL=OFF \
  -DWITH_ZLIB=OFF
cmake --build "$work/build" --target capnp_tool capnpc_cpp -j "$(nproc 2>/dev/null || echo 2)"
install -d "$prefix/bin"
install -m 0755 "$work/build/src/capnp/capnp" "$prefix/bin/capnp"
install -m 0755 "$work/build/src/capnp/capnpc-c++" "$prefix/bin/capnpc-c++"
ln -sf capnp "$prefix/bin/capnpc"
"$prefix/bin/capnp" --version
