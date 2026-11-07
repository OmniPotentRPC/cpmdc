#!/bin/sh
# Run the CPU strong and weak sweeps and the CUDA request.
# The CUDA request records an empty timing when this build has no device.
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
OUT=${1:?output directory}
RUNS=${RUNS:-5}
THREADS=${THREADS:-1,2,4,8}
RANKS=${RANKS:-1 2 4 8}
NATOMS=${NATOMS:-96}
NPW=${NPW:-32768}
ATOMS_PER=${ATOMS_PER:-96}
PW_PER=${PW_PER:-16384}
FFT_N=${FFT_N:-128}
FFT_N_WEAK=${FFT_N_WEAK:-96}
mkdir -p "$OUT/logs" "$OUT/benches"
make -C "$ROOT/bench"
unset OMP_NUM_THREADS || true
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMPI_MCA_btl_base_warn_component_unused=0
export OMPI_MCA_btl=^openib
BIN="$ROOT/bench/build/cpmdc-pwscale"
RAW="$OUT/benches/raw.csv"
rm -f "$RAW"
{
  echo "agreement"
  "$ROOT/bench/build/test_pwscale_agreement"
  mpirun --bind-to none --oversubscribe -np 2 "$ROOT/bench/build/test_pwscale_agreement"
  mpirun --bind-to none --oversubscribe -np 4 "$ROOT/bench/build/test_pwscale_agreement"
  mpirun --bind-to none --oversubscribe -np 8 "$ROOT/bench/build/test_pwscale_agreement"
} >"$OUT/logs/agreement.txt" 2>&1

COMMON="--runs $RUNS --natoms $NATOMS --npw $NPW --atoms-per $ATOMS_PER --pw-per $PW_PER --fft-n $FFT_N --cell 12 --alpha 0.4 --sigma 0.9 --out $RAW"
{
  echo "threads $THREADS"
  for mode in force-strong force-weak-atoms force-weak-pw fft-strong fft-weak; do
    echo "mode $mode"
    if [ "$mode" = fft-weak ]; then
      "$BIN" --mode "$mode" --threads "$THREADS" $COMMON --fft-n "$FFT_N_WEAK"
    else
      "$BIN" --mode "$mode" --threads "$THREADS" $COMMON
    fi
  done
  echo "ranks $RANKS"
  for np in $RANKS; do
    for mode in force-strong force-weak-atoms force-weak-pw; do
      echo "mode $mode np $np"
      mpirun --bind-to none --oversubscribe -np "$np" "$BIN" --mode "$mode" --ranks $COMMON
    done
  done
} >"$OUT/logs/sweep.txt" 2>&1

set +e
"$BIN" --mode cuda --runs "$RUNS" --natoms "$NATOMS" --npw "$NPW" --out "$RAW" >"$OUT/logs/cuda.txt" 2>&1
CUDA_RC=$?
set -e
echo "cuda_exit=$CUDA_RC" >>"$OUT/logs/cuda.txt"
python3 "$ROOT/bench/pwscale_summarize.py" "$RAW" "$OUT/benches/summary.csv" >"$OUT/logs/summary.txt"
echo "wrote $RAW"
