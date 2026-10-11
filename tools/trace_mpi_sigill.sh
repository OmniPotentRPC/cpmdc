#!/bin/bash
# Rerun the two ASan MPI tests with a SIGILL handler in the ranks.
# prterun only reports the signal. This prints the bytes and the symbol.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
pixi_bin="${root}/.pixi/envs/default/bin"
export PATH="${pixi_bin}:${PATH}"
export ASAN_OPTIONS="${ASAN_OPTIONS:-halt_on_error=1:abort_on_error=1:print_summary=1:print_stacktrace=1:detect_leaks=0}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:abort_on_error=1:print_summary=1:print_stacktrace=1}"
cc=${CC:-${pixi_bin}/gcc}
so=/tmp/libtrace_sigill.so
asan=$("${cc}" -print-file-name=libasan.so)
"${cc}" -shared -fPIC -O0 -g -o "${so}" "${root}/tools/trace_sigill.c" -ldl
# ASan refuses to start unless its runtime is the first preloaded library.
export LD_PRELOAD="${asan}:${so}${LD_PRELOAD:+:${LD_PRELOAD}}"
cd "${root}"
for spec in "4 ${root}/build-sanitize/test_adopt_comm_stub" "4 ${root}/build-sanitize/test_adopt_comm_shell ${root}/build-sanitize/cpmd_params_parser.bin"; do
  echo "===== trace ${spec}"
  set +e
  timeout 90 python3 "${root}/tests/run_mpi_test.py" ${spec}
  echo "trace exit $?"
  set -e
done
