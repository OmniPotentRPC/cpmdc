#!/bin/sh
# Compare exported cpmdc_* symbols with the committed baseline.
# A removed or renamed symbol is an unannounced break and fails the job.
# Pass "write" as the third argument to refresh the baseline.
#
# When abidiff is on PATH and scripts/abi_baseline.abi exists, the library
# is also compared type by type. Added functions pass. Any other change
# fails. Without the tool or the corpus, that step is skipped and says so.
# The corpus is written without build or source locations:
#   abidw --headers-dir include --drop-private-types --drop-undefined-syms \
#     --no-corpus-path --no-comp-dir-path --no-show-locs \
#     --no-load-undefined-interfaces -o scripts/abi_baseline.abi libcpmdc.so
set -eu
lib=${1:?library}
baseline=${2:?baseline}
mode=${3:-check}
corpus=$(dirname "$baseline")/abi_baseline.abi
dump=$(mktemp)
trap 'rm -f "$dump"' EXIT

{
  echo "# cpmdc exported C ABI"
  nm -D --defined-only "$lib" | awk '{print $3}' | awk '/^cpmdc_/ {print}' | sort -u
} > "$dump"

if [ "$mode" = "write" ]; then
  cp "$dump" "$baseline"
  exit 0
fi
if [ ! -f "$baseline" ]; then
  echo "missing ABI baseline $baseline" >&2
  exit 1
fi
if ! diff -u "$baseline" "$dump"; then
  echo "C ABI differs from the committed baseline" >&2
  exit 1
fi
echo "exported symbols match $baseline"

abidiff=$(command -v abidiff || true)
if [ -z "$abidiff" ] || [ ! -f "$corpus" ]; then
  echo "abidiff corpus not used; type-level check skipped"
  exit 0
fi
include=$(CDPATH= cd "$(dirname "$baseline")/../include" && pwd)
status=0
"$abidiff" --no-added-syms --drop-private-types --headers-dir2 "$include" \
  "$corpus" "$lib" > "${dump}.abi" 2>&1 || status=$?
if [ "$status" -ne 0 ]; then
  cat "${dump}.abi"
  echo "libcpmdc differs from the ABI corpus beyond added functions (abidiff status $status)" >&2
  exit 1
fi
echo "libcpmdc keeps the ABI corpus; additions only"
