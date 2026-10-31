#!/bin/bash
# Redraw the bond-scan figures from the committed sources.
# Needs Graphviz `dot` and `uv`; the data come from run_scan.sh.
set -euo pipefail
cd "$(dirname "$0")"
dot -Tpng -Gdpi=150 session-loop.dot -o session-loop.png
uv run --quiet plot_scan.py
