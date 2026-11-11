#!/bin/sh
# Copyright (c) cpmdc developers
# Interleave two harness builds in ABBA order over the rank sweeps.
# usage: run_pwscale_abba.sh A_HARNESS B_HARNESS OUT
# Each block runs A, B, B, A one after another, never at the same time.
# A_ARGS and B_ARGS add harness flags to one side only.
# Ranks are bound with mpirun --bind-to "$BIND" and never oversubscribed.
set -eu
A=${1:?A harness}
B=${2:?B harness}
OUT=${3:?output directory}
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
RUNS=${RUNS:-9}
BLOCKS=${BLOCKS:-1}
RANKS=${RANKS:-1 2 4 8}
NATOMS=${NATOMS:-96}
NPW=${NPW:-32768}
ATOMS_PER=${ATOMS_PER:-96}
PW_PER=${PW_PER:-16384}
BIND=${BIND:-core}
A_ARGS=${A_ARGS:-}
B_ARGS=${B_ARGS:-}
MODES="force-strong force-weak-atoms force-weak-pw"
mkdir -p "$OUT/logs" "$OUT/benches"
unset OMP_NUM_THREADS || true
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMPI_MCA_btl_base_warn_component_unused=0
COMMON="--runs $RUNS --natoms $NATOMS --npw $NPW --atoms-per $ATOMS_PER --pw-per $PW_PER --cell 12 --alpha 0.4 --sigma 0.9"
{
  echo "nproc $(nproc)"
  grep -m1 'model name' /proc/cpuinfo || true
  echo "uname $(uname -srm)"
  echo "mpirun $(mpirun --version 2>&1 | head -n 1)"
  echo "bind $BIND runs $RUNS blocks $BLOCKS ranks $RANKS"
  echo "natoms $NATOMS npw $NPW atoms_per $ATOMS_PER pw_per $PW_PER"
  echo "a_args $A_ARGS"
  echo "b_args $B_ARGS"
} >"$OUT/logs/host.txt"

pass() {
  side=$1
  bin=$2
  args=$3
  csv=$4
  for mode in $MODES; do
    echo "serial mode $mode"
    "$bin" --mode "$mode" --threads 1 $COMMON --out "$csv"
  done
  for np in $RANKS; do
    for mode in $MODES; do
      echo "side $side mode $mode np $np"
      mpirun --bind-to "$BIND" -np "$np" "$bin" --mode "$mode" --ranks $args $COMMON --out "$csv"
    done
  done
}

block=1
while [ "$block" -le "$BLOCKS" ]; do
  slot=1
  for side in A B B A; do
    if [ "$side" = A ]; then
      bin=$A
      args=$A_ARGS
    else
      bin=$B
      args=$B_ARGS
    fi
    csv="$OUT/benches/b$block-s$slot-$side.csv"
    rm -f "$csv"
    echo "block $block slot $slot side $side start $(date -u +%Y-%m-%dT%H:%M:%SZ) load $(cat /proc/loadavg)" >>"$OUT/logs/passes.txt"
    pass "$side" "$bin" "$args" "$csv" >>"$OUT/logs/sweep.txt" 2>&1
    echo "block $block slot $slot side $side end $(date -u +%Y-%m-%dT%H:%M:%SZ) load $(cat /proc/loadavg)" >>"$OUT/logs/passes.txt"
    slot=$((slot + 1))
  done
  block=$((block + 1))
done
python3 "$ROOT/bench/pwscale_abba.py" "$OUT/benches" "$OUT/benches/abba.csv" >"$OUT/logs/abba.txt"
echo "wrote $OUT/benches/abba.csv"
