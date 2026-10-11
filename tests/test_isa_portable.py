#!/usr/bin/env python3
"""Portable OpenCPMD configs, the empty ISA allowlist, and the golden digest."""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

FORBIDDEN = (
    "-xHost",
    "-march=native",
    "-axCORE-AVX512",
    "-mavx",
    "-march=x86-64-v3",
    "-march=x86-64-v2",
    "-march=armv8",
)

CONFIGS = (
    "tools/LINUX-X86_64-GFORTRAN-MPI-PIC",
    "tools/LINUX-AARCH64-GFORTRAN-MPI-PIC",
)


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def strip_comments(text: str) -> str:
    kept = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        kept.append(line.split("#", 1)[0])
    return "\n".join(kept)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_isa_portable.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    failures: list[str] = []

    for rel in CONFIGS:
        path = root / rel
        require(path.is_file(), f"missing {rel}", failures)
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        code = strip_comments(text)
        require(
            code.count("-ffp-contract=off") >= 4,
            f"{rel} does not put -ffp-contract=off on debug and optimised FFLAGS and CFLAGS",
            failures,
        )
        require(code.count("-fPIC") >= 4, f"{rel} does not pass -fPIC on FFLAGS and CFLAGS", failures)
        for flag in FORBIDDEN:
            require(flag not in code, f"{rel} contains {flag}", failures)
        require("-march" not in code, f"{rel} passes -march", failures)

    aarch64 = (root / CONFIGS[1]).read_text(encoding="utf-8") if (root / CONFIGS[1]).is_file() else ""
    require("armv8-a" in aarch64, "aarch64 config does not name the armv8-a default", failures)
    x86 = (root / CONFIGS[0]).read_text(encoding="utf-8") if (root / CONFIGS[0]).is_file() else ""
    require("x86-64" in x86, "x86_64 config does not name the x86-64 default", failures)

    meson = (root / "meson.build").read_text(encoding="utf-8")
    require("-ffp-contract=off" in meson, "meson.build does not pass -ffp-contract=off", failures)
    require("/fp:precise" in meson, "meson.build does not pass /fp:precise for MSVC", failures)
    require("-march=x86-64" in meson, "meson.build does not pin the x86-64 baseline", failures)
    require("-march=armv8-a" in meson, "meson.build does not pin the armv8-a baseline", failures)
    require("-march=native" not in meson, "meson.build passes -march=native", failures)
    require("-march=x86-64-v3" not in meson, "meson.build passes -march=x86-64-v3", failures)
    require(
        meson.find("add_global_arguments") < meson.find("subproject("),
        "fp-contract arguments are not registered before subproject()",
        failures,
    )

    allow = (root / "isa-allow.toml").read_text(encoding="utf-8")
    allow_code = strip_comments(allow)
    require('baseline = "x86-64"' in allow_code, "isa-allow.toml baseline is not x86-64", failures)
    require("[[allow]]" not in allow_code, "isa-allow.toml allows a cpmdc symbol", failures)
    require("[[vendored]]" not in allow_code, "isa-allow.toml vendors a library that CI links dynamically", failures)

    how = (root / "docs/orgmode/howto/opencpmd-archive.org").read_text(encoding="utf-8")
    how_rst = (root / "docs/source/howto/opencpmd-archive.rst").read_text(encoding="utf-8")
    for label, text in (("org", how), ("rst", how_rst)):
        require(
            "Building OpenCPMD for mixed clusters" in text,
            f"opencpmd-archive {label} lacks the mixed-cluster section",
            failures,
        )
        require("openblas_get_config" in text, f"opencpmd-archive {label} lacks openblas_get_config", failures)
        require("-ffp-contract=off" in text, f"opencpmd-archive {label} lacks -ffp-contract=off", failures)
        require(
            "LINUX-AARCH64-GFORTRAN-MPI-PIC" in text,
            f"opencpmd-archive {label} lacks the aarch64 config",
            failures,
        )
        require("10.1002/spe.3075" in text, f"opencpmd-archive {label} lacks the EESSI DOI", failures)

    simd = root / "docs/orgmode/explanation/simd-paths.org"
    simd_rst = root / "docs/source/explanation/simd-paths.rst"
    require(simd.is_file(), "missing docs/orgmode/explanation/simd-paths.org", failures)
    require(simd_rst.is_file(), "missing docs/source/explanation/simd-paths.rst", failures)
    if simd.is_file() and simd_rst.is_file():
        for label, text in (
            ("org", simd.read_text(encoding="utf-8")),
            ("rst", simd_rst.read_text(encoding="utf-8")),
        ):
            require("CPMDC_SIMD" in text, f"simd-paths {label} does not name CPMDC_SIMD", failures)
            require("-ffp-contract=off" in text, f"simd-paths {label} lacks -ffp-contract=off", failures)
            require("10.1016/j.cpc.2013.06.003" in text, f"simd-paths {label} lacks the Páll and Hess DOI", failures)

    index = (root / "docs/orgmode/index.org").read_text(encoding="utf-8")
    index_rst = (root / "docs/source/index.rst").read_text(encoding="utf-8")
    require("explanation/simd-paths" in index, "org index does not list simd-paths", failures)
    require("explanation/simd-paths" in index_rst, "rst index does not list simd-paths", failures)

    workflow = (root / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    require(
        "8957850fd15005ea460ba336781ed152f7ee3d92" in workflow,
        "ci.yml does not pin the I0 audit workflow by SHA",
        failures,
    )
    require("isa-allow.toml" in workflow, "ci.yml does not pass isa-allow.toml", failures)
    require("workflow_dispatch" in workflow, "ci.yml has no workflow_dispatch trigger", failures)
    require("ubuntu-24.04-arm" in workflow, "ci.yml has no ubuntu-24.04-arm job", failures)

    golden = root / "tests/golden/isa-hashes.toml"
    require(golden.is_file(), "missing tests/golden/isa-hashes.toml", failures)
    if golden.is_file():
        golden_text = golden.read_text(encoding="utf-8")
        digest = hashlib.sha256(b"").hexdigest()
        require(digest in golden_text, "golden file digest is not SHA-256 of the empty f64 buffer", failures)
        require("no hot numeric kernel" in golden_text, "golden file recipe does not name the empty kernel", failures)
        for thread in ("1", "2", "3", "8"):
            require(re.search(rf"\b{thread}\b", golden_text) is not None, f"golden file omits {thread} threads", failures)

    flags_script = root / "scripts/with_baseline_flags.sh"
    require(flags_script.is_file(), "missing scripts/with_baseline_flags.sh", failures)
    if flags_script.is_file():
        flags_text = flags_script.read_text(encoding="utf-8")
        require("-march=*" in flags_text, "baseline flag script does not drop -march", failures)
        check = subprocess.run(["bash", "-n", str(flags_script)], check=False, capture_output=True, text=True)
        require(check.returncode == 0, f"with_baseline_flags.sh: {check.stderr}", failures)
    pixi = (root / "pixi.toml").read_text(encoding="utf-8")
    require(pixi.count("with_baseline_flags.sh") >= 4, "pixi test tasks do not source the baseline flag script", failures)

    script = root / "tools/collect_isa_audit_libs.sh"
    require(script.is_file(), "missing tools/collect_isa_audit_libs.sh", failures)
    if script.is_file():
        check = subprocess.run(["bash", "-n", str(script)], check=False, capture_output=True, text=True)
        require(check.returncode == 0, f"collect_isa_audit_libs.sh: {check.stderr}", failures)

    flags = ""
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("flags"):
                flags = line
                break
    print(flags or "cpu flags: unavailable")
    print("path scalar: exercised (compiler default; no dispatched kernel to force)")
    print("path avx2: not exercised")
    print("path avx512: not exercised")
    print("path neon: not exercised")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
