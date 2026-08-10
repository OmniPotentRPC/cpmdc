#!/usr/bin/env bash
set -euo pipefail

cpmd_root=${1:?usage: rebuild_opencpmd_embed.sh CPMD_ROOT [JOBS]}
jobs=${2:-${CPMDC_BUILD_JOBS:-1}}
make_program=${MAKE:-make}

if [[ ! -f "$cpmd_root/Makefile" ]]; then
  printf 'missing OpenCPMD Makefile: %s\n' "$cpmd_root/Makefile" >&2
  exit 2
fi
if [[ ! -d "$cpmd_root/obj" || ! -d "$cpmd_root/lib" ]]; then
  printf 'incomplete OpenCPMD build tree: %s\n' "$cpmd_root" >&2
  exit 2
fi

"$make_program" \
  -C "$cpmd_root/obj" \
  -f "$cpmd_root/Makefile" \
  -j "$jobs" \
  rwfopt_utils.mod.o \
  updwf_utils.mod.o \
  "$cpmd_root/lib/libcpmd.a"
