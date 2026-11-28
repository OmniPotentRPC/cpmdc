#!/usr/bin/env bash
# Assemble the pseudopotential directory of the OpenCPMD test suite from
# OpenCPMD's Regtests at a fixed commit. Regtests keeps the files in two
# places: tests/PP_LIBRARY holds the water decks' O_MT_BLYP.psp and
# H_CVB_BLYP.psp, and tests/PP_LIBRARY/pseudo holds Si_MT_BLYP.psp and
# Si_MT_PBE.psp, which the Si3N4 and silylamine decks name. OpenCPMD reads
# every file of a run from one directory, CPMDC_PSEUDO_DIR, so DEST gets
# the files of PP_LIBRARY and, for names PP_LIBRARY lacks, those of
# PP_LIBRARY/pseudo.
#
# Usage: fetch_pseudopotentials.sh [--fhi] [--pbe] DEST [CLONE]
#   DEST is the directory to export as CPMDC_PSEUDO_DIR; CLONE (default
#   DEST.regtests) is the shallow Regtests checkout, reused when present.
#   --fhi also downloads O.pbe-mt_fhi.UPF and H.pbe-mt_fhi.UPF from the
#   Quantum ESPRESSO table and verifies their SHA-256 checksums.
#   --pbe also downloads O_MT_PBE.psp and H_MT_PBE.psp from the OpenCPMD
#   extended pseudopotential library and verifies their SHA-256 checksums.
#   A checksum that does not match prints the file name and sha256sum's
#   own report.
set -euo pipefail

fetch_fhi=0
fetch_pbe=0
while [ $# -gt 0 ]; do
  case "$1" in
    --fhi) fetch_fhi=1; shift ;;
    --pbe) fetch_pbe=1; shift ;;
    --) shift; break ;;
    -*) echo "fetch_pseudopotentials.sh: unknown option $1" >&2; exit 2 ;;
    *) break ;;
  esac
done

commit=b2b2c64d39da2ecb96857ab9a007f34332db9d1b
dest=${1:?usage: fetch_pseudopotentials.sh [--fhi] [--pbe] DEST [CLONE]}
clone=${2:-$1.regtests}

if [ ! -d "$clone/.git" ]; then
  git init -q "$clone"
fi
if [ "$(git -C "$clone" rev-parse -q --verify HEAD 2>/dev/null)" != "$commit" ]; then
  git -C "$clone" fetch -q --depth 1 https://github.com/OpenCPMD/Regtests.git "$commit"
  git -C "$clone" checkout -q FETCH_HEAD
fi
lib="$clone/tests/PP_LIBRARY"
mkdir -p "$dest"
find "$lib" -maxdepth 1 -type f -exec cp -f {} "$dest/" \;
for file in "$lib"/pseudo/*; do
  if [ -f "$file" ] && [ ! -e "$dest/$(basename "$file")" ]; then
    cp "$file" "$dest/"
  fi
done
for name in O_MT_BLYP.psp H_CVB_BLYP.psp Si_MT_BLYP.psp N_MT_BLYP.psp; do
  test -f "$dest/$name" || { echo "fetch_pseudopotentials.sh: $name missing" >&2; exit 1; }
done

check_sha256() {
  local path=$1 expected=$2
  if ! echo "$expected  $path" | sha256sum -c -; then
    echo "fetch_pseudopotentials.sh: checksum mismatch for $(basename "$path")" >&2
    exit 1
  fi
}

if [ "$fetch_fhi" -eq 1 ]; then
  fhi_base="https://pseudopotentials.quantum-espresso.org/upf_files"
  check_fhi() {
    local fname=$1 expected=$2
    local path="$dest/$fname"
    if [ ! -f "$path" ] || [ "$(sha256sum "$path" | awk '{print $1}')" != "$expected" ]; then
      curl -fsSL --retry 3 -o "$path" "$fhi_base/$fname"
    fi
    check_sha256 "$path" "$expected"
  }
  check_fhi O.pbe-mt_fhi.UPF 747670be5d6f008d921dc74551dc0705f2c4b34132e79590c482b532171f8cd1
  check_fhi H.pbe-mt_fhi.UPF b6ed642482602f477eb8c167dce4a2c6d254b89522685b5b9ad8f1e904c50db4
fi

if [ "$fetch_pbe" -eq 1 ]; then
  pbe_commit=97f98cb25cbee1cce67302addd2bf06e21a17978
  pbe_url="https://github.com/OpenCPMD/pseudopotentials-extended-library/raw/${pbe_commit}/pseudo.tar.gz"
  o_sum=d26bb0bfa19ddaabe95d21169f50fa25915631e84df0d157e93bb32a7a9298ad
  h_sum=8ff874f26491f2b6e482540342ae0ad85c9c7023936f3b050baa4ddced459988
  o_path="$dest/O_MT_PBE.psp"
  h_path="$dest/H_MT_PBE.psp"
  need=0
  if [ ! -f "$o_path" ] || [ "$(sha256sum "$o_path" | awk '{print $1}')" != "$o_sum" ]; then
    need=1
  fi
  if [ ! -f "$h_path" ] || [ "$(sha256sum "$h_path" | awk '{print $1}')" != "$h_sum" ]; then
    need=1
  fi
  if [ "$need" -eq 1 ]; then
    tmp=$(mktemp -d)
    curl -fsSL --retry 3 -o "$tmp/pseudo.tar.gz" "$pbe_url"
    tar -xf "$tmp/pseudo.tar.gz" -C "$tmp" ./pseudo/O_MT_PBE.psp ./pseudo/H_MT_PBE.psp
    cp "$tmp/pseudo/O_MT_PBE.psp" "$o_path"
    cp "$tmp/pseudo/H_MT_PBE.psp" "$h_path"
    rm -rf "$tmp"
  fi
  check_sha256 "$o_path" "$o_sum"
  check_sha256 "$h_path" "$h_sum"
fi
