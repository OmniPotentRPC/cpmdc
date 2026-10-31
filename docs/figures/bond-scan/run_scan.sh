#!/bin/bash
# Rerun both scans of docs/orgmode/tutorials/bond-scan.org and write scan.csv.
# Run inside `pixi shell` with the OpenCPMD build in $CPMDC/build-cpmd and
# CPMDC_PSEUDO_DIR set. The program is taken from the tutorial page itself.
# Usage: run_scan.sh water.params.bin [log-dir]
# With no log-dir, scratch files go to a new temporary directory.
# The script prints that path.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
CPMDC=${CPMDC:-$(cd "$here/../../.." && pwd)}
params=$(realpath "${1:?usage: run_scan.sh water.params.bin [log-dir]}")
if [ -n "${2:-}" ]; then
  mkdir -p "$2"
  logs=$(cd "$2" && pwd)
else
  logs=$(mktemp -d)
  printf '%s\n' "$logs"
fi

sed -n '/^#+begin_src c$/,/^#+end_src$/p' \
  "$CPMDC/docs/orgmode/tutorials/bond-scan.org" | sed '1d;$d' \
  > "$logs/scan_bond.c"
mpicc -c "$logs/scan_bond.c" -I "$CPMDC/include" -o "$logs/scan_bond.o"
mpif90 "$logs/scan_bond.o" -L "$CPMDC/build-cpmd" -lcpmdc \
  -Wl,-rpath,"$CPMDC/build-cpmd" -o "$logs/scan_bond"

export OMP_NUM_THREADS=1 CPMDC_SCF_STEPS=1
cd "$logs"
./scan_bond "$params" > cpmd-warm.out 2> scan-warm.log
CPMDC_WARM_RMSD=0 ./scan_bond "$params" > cpmd-cold.out 2> scan-cold.log

# Host lines and cpmdc lines go through separate buffers, so pair them by order.
columns() {
  paste -d, \
    <(sed -n 's/^scan r_oh=\([^ ]*\) energy_h=\([^ ]*\) f_bond=\(.*\)/\1,\2,\3/p' "$1") \
    <(awk '/^cpmdc_warm_start/ {print $2}' "$1") \
    <(awk '/^cpmdc_scf_steps/ {print $2}' "$1")
}
{
  echo "r_oh_ang,energy_h,f_bond_h_bohr,warm_start,scf_steps,energy_h_cold,f_bond_cold,warm_start_cold,scf_steps_cold"
  paste -d, <(columns scan-warm.log) <(columns scan-cold.log | cut -d, -f2-)
} > "$here/scan.csv"
