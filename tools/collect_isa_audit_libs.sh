#!/bin/bash
# Copy unstripped libcpmdc and the top-level ELF test binaries.
# Nested Meson and CMake probes under the build tree are not cpmdc.
set -euo pipefail

build=${1:-build}
dest=${2:-isa-audit-libs}

if [ ! -d "$build" ]; then
  echo "build directory not found: $build" >&2
  exit 1
fi

rm -rf "$dest"
mkdir -p "$dest"

shopt -s nullglob
copied=0
for f in "$build"/libcpmdc.so* "$build"/test_* "$build"/example_host_step "$build"/cpmdc-restart; do
  if [ -L "$f" ] || [ ! -f "$f" ]; then
    continue
  fi
  if file "$f" | grep -q ELF; then
    cp -a "$f" "$dest/$(basename "$f")"
    copied=1
  fi
done

if [ "$copied" -ne 1 ]; then
  echo "no ELF libraries or test binaries copied from $build" >&2
  exit 1
fi

for f in "$dest"/*; do
  if ! readelf -S "$f" | grep -q '\.symtab'; then
    echo "NO SYMTAB: $f" >&2
    exit 1
  fi
  echo "symtab: $f"
done
