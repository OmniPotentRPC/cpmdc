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
# Usage: fetch_pseudopotentials.sh DEST [CLONE]
#   DEST is the directory to export as CPMDC_PSEUDO_DIR; CLONE (default
#   DEST.regtests) is the shallow Regtests checkout, reused when present.
set -euo pipefail

commit=b2b2c64d39da2ecb96857ab9a007f34332db9d1b
dest=$1
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
