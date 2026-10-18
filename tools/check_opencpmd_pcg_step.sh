#!/usr/bin/env bash
set -euo pipefail
source=${1:?usage: check_opencpmd_pcg_step.sh CPMD_SOURCE BUILD_DIRECTORY}
build=${2:?Expected isolated fixture build directory}
tools=$(cd -- "$(dirname -- "$0")" && pwd)
mkdir -p "$build"
build=$(cd "$build" && pwd)
sed -n '/^  PURE FUNCTION scaled_pcg_step(/,/^  END FUNCTION scaled_pcg_step/p' \
  "$source/src/pcgrad_driver.mod.F90" >"$build/native_scaled_pcg_step.inc"
test -s "$build/native_scaled_pcg_step.inc"
"${FC:-gfortran}" -cpp -O0 -g -I"$build" -J"$build" -c \
  "$tools/../tests/opencpmd_pcg_step_fixture.F90" -o "$build/step.o"
"${CXX:-c++}" -std=c++17 -O0 -g "$tools/../tests/test_opencpmd_pcg_step.cpp" \
  "$build/step.o" -o "$build/check-native-pcg-step"
"$build/check-native-pcg-step"
