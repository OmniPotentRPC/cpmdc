#!/usr/bin/env bash
set -euo pipefail
cpmd_root=${1:?usage: rebuild_opencpmd_constant_cutoff.sh CPMD_ROOT}
test -f "$cpmd_root/Makefile"
test -f "$cpmd_root/lib/libcpmd.a"
"${MAKE:-make}" -C "$cpmd_root/obj" -f "$cpmd_root/Makefile" -j1 hesele_utils.mod.o
"${AR:-ar}" r "$cpmd_root/lib/libcpmd.a" "$cpmd_root/obj/hesele_utils.mod.o"
ranlib "$cpmd_root/lib/libcpmd.a"
