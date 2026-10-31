#!/bin/bash
# run_scan.sh with no log directory must not write scratch files here.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
script=$here/run_scan.sh
root=$(cd "$here/../../.." && pwd)

stub=$(mktemp -d)
printf '#!/bin/sh\nexit 127\n' > "$stub/mpicc"
cp "$stub/mpicc" "$stub/mpif90"
chmod +x "$stub/mpicc" "$stub/mpif90"

fail() {
  printf '%s\n' "$1" >&2
  exit 1
}

work=$(mktemp -d)
cd "$work"
touch params.bin
set +e
printed=$(PATH="$stub:/usr/bin:/bin" CPMDC="$root" "$script" "$work/params.bin" 2>"$work/err")
status=$?
set -e
[ "$status" -ne 0 ] || fail "expected the compile to fail"
[ -n "$printed" ] || fail "no log directory was printed"
[ -f "$printed/scan_bond.c" ] || fail "the program was not written to the printed directory"
[ ! -e "$work/scan_bond.c" ] || fail "scan_bond.c landed in the working directory"
[ ! -e "$work/scan_bond.o" ] || fail "an object file landed in the working directory"
grep -q '#include <cpmdc.h>' "$printed/scan_bond.c" || fail "the log directory does not hold the page program"
case "$printed" in
  "$work" | "$work"/*) fail "the default log directory is inside the working directory" ;;
esac

explicit=$(mktemp -d)
set +e
PATH="$stub:/usr/bin:/bin" CPMDC="$root" "$script" "$work/params.bin" "$explicit" >"$work/explicit.out" 2>"$work/explicit.err"
status=$?
set -e
[ "$status" -ne 0 ] || fail "expected the explicit run to stop at the compile"
[ ! -s "$work/explicit.out" ] || fail "an explicit log directory should not be printed"
[ -f "$explicit/scan_bond.c" ] || fail "the explicit log directory was not used"
[ ! -e "$work/scan_bond.c" ] || fail "the explicit run wrote scan_bond.c into the working directory"
printf 'ok\n'
