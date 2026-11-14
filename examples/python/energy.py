#!/usr/bin/env python3
"""Print the energy in Hartree and the atom count of one isolated step.

Usage: energy.py PARAMS_BIN [LIBRARY]
"""

from __future__ import annotations

import sys
from pathlib import Path

from cpmdc import Session


def main(argv: list[str]) -> int:
    if len(argv) not in (2, 3):
        print(f"usage: {argv[0]} PARAMS_BIN [LIBRARY]", file=sys.stderr)
        return 2
    library = argv[2] if len(argv) == 3 else None
    with Session(Path(argv[1]).read_bytes(), library=library) as session:
        done = session.energy_forces(
            [[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]],
            [8, 1, 1],
        )
    print(f"energy_h={done.energy_h:.12f}")
    print(f"natoms={len(done.forces_h_bohr)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
