#!/bin/bash
# The sanitize MPI tests SIGILL only on runners that advertise avx_vnni and
# avx512_vp2intersect. AVX-512F alone does not. This runs those two
# instructions under a SIGILL handler so the log names the one that faults.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
pixi_bin="${root}/.pixi/envs/default/bin"
export PATH="${pixi_bin}:${PATH}"
cc=${CC:-${pixi_bin}/gcc}
"${cc}" -O0 -g -o /tmp/probe_extra_isa "${root}/tools/probe_extra_isa.c"
/tmp/probe_extra_isa
