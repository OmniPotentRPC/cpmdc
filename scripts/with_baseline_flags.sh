#!/bin/bash
# Drop ISA selectors from the compiler environment.
# conda-forge exports CFLAGS=-march=nocona, and Meson appends that
# environment after its own arguments, so the later -march would win.
# Meson then passes -march=x86-64 or -march=armv8-a itself.
strip_march() {
  local out="" tok
  for tok in $1; do
    case "$tok" in
      -march=*|-mtune=*|-mavx*|-mfma|-msse3|-mssse3|-msse4*|-mpopcnt|-mbmi*|-mf16c|-mlzcnt)
        ;;
      *)
        out="${out} ${tok}"
        ;;
    esac
  done
  printf '%s' "${out# }"
}

export CFLAGS="$(strip_march "${CFLAGS-}")"
export CXXFLAGS="$(strip_march "${CXXFLAGS-}")"
export FFLAGS="$(strip_march "${FFLAGS-}")"
export FCFLAGS="$(strip_march "${FCFLAGS-}")"
