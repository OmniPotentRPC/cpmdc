#!/usr/bin/env bash
# Check or rewrap tracked Org and Markdown prose with snapper at 80 columns.
# --check goes through tools/snapper_gate.py, which fails on files that would
# reformat, fused lines and render-check refusals, and reports the other
# diagnostic classes; SNAPPER_FAIL_ON=long,wrap enforces those too.
# --write goes through tools/snapper_write.py, which refuses symlinks and
# paths outside the repository and replaces each file with a rename in its
# own directory.
set -euo pipefail

ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
CONFIG="$ROOT/.snapperrc.toml"
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if command -v snapper-fmt >/dev/null 2>&1; then
  SNAPPER=snapper-fmt
elif command -v snapper >/dev/null 2>&1; then
  SNAPPER=snapper
else
  echo "snapper-fmt (or snapper) not on PATH; install snapper-fmt 0.11.10" >&2
  exit 1
fi

mode="--check"
if [ "${1:-}" = "--write" ] || [ "${1:-}" = "--in-place" ]; then
  mode="--write"
  shift
elif [ "${1:-}" = "--check" ]; then
  shift
fi

if [ "$#" -gt 0 ]; then
  files=("$@")
else
  mapfile -d '' -t files < <(git -C "$ROOT" ls-files -z '*.org' '*.md')
  if [ "${#files[@]}" -eq 0 ]; then
    exit 0
  fi
  cd "$ROOT"
fi

if [ "$mode" = "--write" ]; then
  exec "${PYTHON:-python3}" "$HERE/snapper_write.py" --snapper "$SNAPPER" \
    --config "$CONFIG" --root "$ROOT" --max-width 80 -- "${files[@]}"
fi

exec "${PYTHON:-python3}" "$HERE/snapper_gate.py" --snapper "$SNAPPER" \
  --config "$CONFIG" --max-width 80 --fail-on "${SNAPPER_FAIL_ON:-}" \
  -- "${files[@]}"
