#!/usr/bin/env bash
set -euo pipefail
source=${1:?usage: check_opencpmd_omp_info.sh CPMD_SOURCE BUILD_DIRECTORY}
build=${2:?Expected fixture build directory}
tools=$(cd -- "$(dirname -- "$0")" && pwd)
fixture="$tools/../tests/opencpmd_omp_info_fixture.F90"
mkdir -p "$build"
build=$(cd "$build" && pwd)
sed -n '/^  SUBROUTINE print_omp_info()/,/^  END SUBROUTINE print_omp_info/p' \
  "$source/src/rggen_utils.mod.F90" >"$build/native_print_omp_info.inc"
test -s "$build/native_print_omp_info.inc"
cd "$build"
"${FC:-mpifort}" -cpp -O0 -g -I"$build" "$fixture" -o check-native-omp-info
"${MPIEXEC:-mpirun}" -np "${MPI_NP:-3}" --bind-to none ./check-native-omp-info
