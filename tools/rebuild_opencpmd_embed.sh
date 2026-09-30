#!/usr/bin/env bash
set -euo pipefail

cpmd_root=${1:?usage: rebuild_opencpmd_embed.sh CPMD_ROOT}
make_program=${MAKE:-make}

if [[ ! -f "$cpmd_root/Makefile" ]]; then
  printf 'missing OpenCPMD Makefile: %s\n' "$cpmd_root/Makefile" >&2
  exit 2
fi
if [[ ! -d "$cpmd_root/obj" || ! -d "$cpmd_root/lib" ]]; then
  printf 'incomplete OpenCPMD build tree: %s\n' "$cpmd_root" >&2
  exit 2
fi

# gfortran stores a copy of Scex_t in rwfopt_utils.mod. One job writes
# scex_utils.mod before rwfopt_utils reads it. A parallel rebuild leaves
# the two module files naming different vtable components, and the next
# USE rwfopt_utils stops.
objects=(
  error_handling.mod.o
  scex_utils.mod.o
  rwfopt_utils.mod.o
  updwf_utils.mod.o
  rkpnt_utils.mod.o
)
"$make_program" \
  -C "$cpmd_root/obj" \
  -f "$cpmd_root/Makefile" \
  -j1 \
  "${objects[@]}"
# The generated Makefile sets AR to the archiver with no operation
# letter, so its libcpmd.a rule does not replace members.
ar_program=${AR:-ar}
"$ar_program" r "$cpmd_root/lib/libcpmd.a" \
  "${objects[@]/#/$cpmd_root/obj/}"
ranlib "$cpmd_root/lib/libcpmd.a"
