#!/usr/bin/env python3
"""The optional nvcc line must see the scaffold header. It does not run CUDA."""

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_pwscale_cuda_include.py SOURCE_ROOT", file=sys.stderr)
        return 2
    text = (Path(sys.argv[1]) / "meson.build").read_text(encoding="utf-8")
    start = text.find("pwscale_cuda = custom_target(")
    end = text.find("pwscale_objects += [pwscale_cuda]", start)
    if start < 0 or end < 0:
        print("cuda custom_target is missing", file=sys.stderr)
        return 1
    block = text[start:end]
    if "'-I' + meson.current_source_dir() / 'include'" not in block:
        print("nvcc command has no scaffold include directory", file=sys.stderr)
        return 1
    if "src/cpmdc_pwscale_cuda.cu" not in block:
        print("nvcc command does not name the cuda source", file=sys.stderr)
        return 1
    print("pwscale cuda include ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
