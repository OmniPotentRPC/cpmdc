#!/usr/bin/env bash
# Build OpenCPMD at the pinned SHA, apply tools/*.patch in tools/README.md
# order, and run the linked Meson suites. A cache hit restores
# lib/libcpmd.a, obj/timetag.o, and obj/*.mod. The Fortran bridge compiles
# against those module files, so they travel with the archive.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
opencpmd_sha=062582b7cfd832d36f88f504cd08e4ead42eb404
regtests_sha=b2b2c64d39da2ecb96857ab9a007f34332db9d1b
# Bust the cache when the build sequence changes.
recipe=v2-patch-configure-make-rebuild-timetag

patches=(
  opencpmd_embed_rinitwf.patch
  opencpmd_embed_geometry.patch
  opencpmd_embed_rwfopt.patch
  opencpmd_converged_state.patch
  opencpmd_kpoints_inputfile.patch
  opencpmd_stopgm_return.patch
  opencpmd_tistopgm.patch
  opencpmd_c_mem_addrs.patch
  opencpmd_embed_teardown.patch
)
config_name=LINUX-X86_64-GFORTRAN-MPI-PIC

# OpenCPMD/Regtests tests/PP_LIBRARY at regtests_sha. Si lives under pseudo/.
pp_files=(
  "H_CVB_BLYP.psp 86995e5ec76d0fbeca6f1e543eef2697e38fb3934e0befcdd331fefb79861f42 H_CVB_BLYP.psp"
  "N_MT_BLYP.psp fdf235554b60e48e2373487aa1e32c5ea889b604c3bb7718ee7d9ca5df25b158 N_MT_BLYP.psp"
  "O_MT_BLYP.psp 9d59ac632cd95b51a1d6ad6de8fe0ee4849d7973419492363b2adc1e2dad562c O_MT_BLYP.psp"
  "Si_MT_BLYP.psp 41ffd6dcf119648267245c7b61cc2cd57fe72f4829e58beb7b843fa9cd54c3a4 pseudo/Si_MT_BLYP.psp"
)

cache_key() {
  {
    printf '%s\n%s\n' "$opencpmd_sha" "$recipe"
    printf '%s %s\n' "$config_name" "$(sha256sum "$root/tools/$config_name" | awk '{print $1}')"
    local patch
    for patch in "${patches[@]}"; do
      printf '%s %s\n' "$patch" "$(sha256sum "$root/tools/$patch" | awk '{print $1}')"
    done
  } | sha256sum | awk '{print $1}'
}

if [[ "${1:-}" == "key" ]]; then
  cache_key
  exit 0
fi

key=$(cache_key)
work=${CPMDC_OPENCPMD_WORK:-/tmp/cpmdc-opencpmd-e2e}
src=$work/CPMD
pp=$work/pp
if [[ -n "${CPMDC_OPENCPMD_CACHE:-}" ]]; then
  cache=$CPMDC_OPENCPMD_CACHE
else
  cache=${HOME}/.cache/cpmdc-opencpmd/$key
fi
build=$root/build-cpmd
jobs=${CPMDC_OPENCPMD_JOBS:-$(nproc)}
pristine=${CPMDC_OPENCPMD_PRISTINE:-/tmp/opencpmd/CPMD}
started=$(date +%s)
archive_s=0
mode=cold

die() {
  printf 'ci_opencpmd_e2e: %s\n' "$*" >&2
  exit 1
}

need() {
  command -v "$1" >/dev/null 2>&1 || die "missing $1"
}

apply_patches() {
  local patch
  for patch in "${patches[@]}"; do
    patch -p1 --forward --batch -d "$src" -i "$root/tools/$patch"
  done
}

save_archive() {
  mkdir -p "$cache/lib" "$cache/obj"
  cp -a "$src/lib/libcpmd.a" "$cache/lib/libcpmd.a"
  find "$cache/obj" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  find "$src/obj" -maxdepth 1 \( -name '*.mod' -o -name 'timetag.o' \) -exec cp -a {} "$cache/obj/" \;
  printf '%s\n' "$key" >"$cache/CACHE_KEY"
  test -f "$cache/lib/libcpmd.a"
  test -f "$cache/obj/timetag.o"
  test -n "$(find "$cache/obj" -maxdepth 1 -name '*.mod' -print -quit)"
}

restore_archive() {
  test -f "$cache/CACHE_KEY" || return 1
  [[ $(cat "$cache/CACHE_KEY") == "$key" ]] || return 1
  test -f "$cache/lib/libcpmd.a" || return 1
  test -f "$cache/obj/timetag.o" || return 1
  mkdir -p "$src/lib" "$src/obj"
  cp -a "$cache/lib/libcpmd.a" "$src/lib/libcpmd.a"
  find "$src/obj" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  find "$cache/obj" -maxdepth 1 \( -name '*.mod' -o -name 'timetag.o' \) -exec cp -a {} "$src/obj/" \;
}

fetch_pseudopotentials() {
  local row name sum rel dest got
  mkdir -p "$pp"
  for row in "${pp_files[@]}"; do
    name=${row%% *}
    sum=${row#* }
    sum=${sum%% *}
    rel=${row##* }
    dest=$pp/$name
    if [[ -f $dest ]] && [[ $(sha256sum "$dest" | awk '{print $1}') == "$sum" ]]; then
      continue
    fi
    curl -fsSL \
      "https://raw.githubusercontent.com/OpenCPMD/Regtests/${regtests_sha}/tests/PP_LIBRARY/${rel}" \
      -o "$dest"
    got=$(sha256sum "$dest" | awk '{print $1}')
    [[ $got == "$sum" ]] || die "$name sha256 $got, expected $sum"
  done
}

reject_skips() {
  python3 - "$1" <<'PY'
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace")
bad = []
for chunk in text.split("===================================="):
    name = ""
    result = ""
    for line in chunk.splitlines():
        if line.startswith("test:"):
            name = line.split(":", 1)[1].strip()
        elif line.startswith("result:"):
            result = line.split(":", 1)[1].strip()
    if not name or not result:
        continue
    low = result.lower()
    if "skip" in low or low.endswith(" 77") or "exit status 77" in low:
        bad.append(f"{name}: {result}")
if bad:
    print("OpenCPMD-linked suite skipped:", file=sys.stderr)
    print("\n".join(bad), file=sys.stderr)
    sys.exit(1)
PY
}

run_rank() {
  local ranks=$1
  local threads=$2
  local log=$3
  local launcher
  printf 'cpmdc-rank-run mpi=%s omp=%s\n' "$ranks" "$threads"
  launcher=(python3 "$root/tests/run_mpi_test.py" "$ranks" "$build/test_e2e_single_point"
    "$build/water_blyp_params.bin" "$build/force_input_step_a.bin")
  OMP_NUM_THREADS=$threads timeout 900 "${launcher[@]}" | tee "$log"
  grep -E '^cpmdc-rank-report ' "$log" >/dev/null || die "missing cpmdc-rank-report for mpi=$ranks omp=$threads"
}

need python3
need curl
need patch
need sha256sum
need meson
need mpif90
need gcc

mkdir -p "$work"
if [[ -d $pristine/.git ]]; then
  got=$(git -C "$pristine" rev-parse HEAD)
  [[ $got == "$opencpmd_sha" ]] || die "pristine $pristine is $got, expected $opencpmd_sha"
  rm -rf "$src"
  cp -a "$pristine" "$src"
else
  rm -rf "$src"
  mkdir -p "$src"
  git -C "$src" init
  git -C "$src" remote add origin https://github.com/OpenCPMD/CPMD.git
  git -C "$src" fetch --depth 1 origin "$opencpmd_sha"
  git -C "$src" checkout --detach FETCH_HEAD
fi

archive_begin=$(date +%s)
if restore_archive; then
  mode=warm
  apply_patches
  printf 'opencpmd archive cache hit %s\n' "$key"
else
  mode=cold
  cp "$root/tools/$config_name" "$src/configure/$config_name"
  # GCC 14 rejects the unpatched cGetMemAddrs return, and the new embed
  # modules have to be in the Makefile. Apply the patches, then configure.
  apply_patches
  (cd "$src" && ./configure.sh -DEST="$src" "$config_name")
  make -C "$src" -j"$jobs"
  "$root/tools/rebuild_opencpmd_embed.sh" "$src"
  make -C "$src/obj" -f "$src/Makefile" timetag.o
  save_archive
  printf 'opencpmd archive cache saved %s\n' "$key"
fi
archive_s=$(($(date +%s) - archive_begin))

fetch_pseudopotentials
export CPMDC_PSEUDO_DIR=$pp
export CPMD_PP_LIBRARY_PATH=${pp}/
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMPI_MCA_btl_base_warn_component_unused=0
export OMPI_MCA_rmaps_base_oversubscribe=1

meson_begin=$(date +%s)
if [[ -f $build/build.ninja ]]; then
  meson setup --reconfigure "$build" -Dwith_cpmd=true -Dwith_tests=true -Dcpmd_root="$src"
else
  meson setup "$build" -Dwith_cpmd=true -Dwith_tests=true -Dcpmd_root="$src"
fi
meson compile -C "$build"
meson_s=$(($(date +%s) - meson_begin))

log=$build/meson-logs/testlog.txt
rm -f "$log"
suite_begin=$(date +%s)
meson test -C "$build" \
  --suite e2e --suite embed --suite water --suite opencpmd \
  --num-processes 1 --print-errorlogs
reject_skips "$log"
suite_s=$(($(date +%s) - suite_begin))

rank_begin=$(date +%s)
mkdir -p "$work/rank"
run_rank 1 1 "$work/rank/mpi1-omp1.log"
run_rank 2 1 "$work/rank/mpi2-omp1.log"
run_rank 1 2 "$work/rank/mpi1-omp2.log"
run_rank 1 "$jobs" "$work/rank/mpi1-ompn.log"
python3 - "$work/rank/mpi1-omp1.log" "$work/rank/mpi2-omp1.log" "$work/rank/mpi1-omp2.log" "$work/rank/mpi1-ompn.log" <<'PY'
import sys
from pathlib import Path


def report(path):
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("cpmdc-rank-report "):
            fields = {}
            for token in line.split()[1:]:
                key, value = token.split("=", 1)
                fields[key] = value
            return fields
    raise SystemExit(f"no report in {path}")


one, two, omp2, ompn = (report(path) for path in sys.argv[1:])
energy = abs(float(one["energy_h"]) - float(two["energy_h"]))
force = abs(float(one["max_abs_force"]) - float(two["max_abs_force"]))
print(
    "cpmdc-rank-deviation mpi=1-vs-2 omp=1 "
    f"energy_h={energy:.17g} max_abs_force={force:.17g} "
    f"fnv_equal={one['fnv'] == two['fnv']}"
)
print(
    "cpmdc-thread-report omp=1 "
    f"energy_h={one['energy_h']} max_abs_force={one['max_abs_force']} fnv={one['fnv']}"
)
print(
    "cpmdc-thread-report omp=2 "
    f"energy_h={omp2['energy_h']} max_abs_force={omp2['max_abs_force']} fnv={omp2['fnv']}"
)
print(
    "cpmdc-thread-report omp=n "
    f"energy_h={ompn['energy_h']} max_abs_force={ompn['max_abs_force']} fnv={ompn['fnv']}"
)
print(
    "cpmdc-thread-deviation report-only "
    f"fnv_1_eq_2={one['fnv'] == omp2['fnv']} fnv_1_eq_n={one['fnv'] == ompn['fnv']}"
)
PY
rank_s=$(($(date +%s) - rank_begin))
total_s=$(($(date +%s) - started))
printf 'opencpmd-ci mode=%s archive_s=%s meson_s=%s suite_s=%s rank_s=%s total_s=%s key=%s\n' \
  "$mode" "$archive_s" "$meson_s" "$suite_s" "$rank_s" "$total_s" "$key"
