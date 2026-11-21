#!/usr/bin/env python3
"""Reject test-only fault controls in the installable shared library."""

from pathlib import Path
import sys


def main() -> int:
    library = Path(sys.argv[1])
    data = library.read_bytes()
    markers = (
        b"CPMDC_TEST_FAULT",
        b"CPMDC_TEST_ORBITAL_FAULT",
        b"energy-allocation",
        b"topology-allocation",
        b"drop-stale",
        b"session-energy-allocation",
        b"calculate-forces-reserve",
        b"calculate-result-reserve",
        b"calculate-result-allocation",
        b"embed-energy-gradient",
        b"potential-result-write",
        b"calculate-au-reserve",
        b"selector-allocation",
        b"embed-available",
    )
    found = [marker.decode() for marker in markers if marker in data]
    if found:
        print(f"{library}: test fault controls present: {', '.join(found)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
