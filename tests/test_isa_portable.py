#!/usr/bin/env python3
"""Portable OpenCPMD configs and the empty ISA allowlist."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

# Whole -march tokens. The substring "-march=armv8" also matches "-march=armv8-a".
FORBIDDEN = (
    "-xHost",
    "-march=native",
    "-axCORE-AVX512",
    "-mavx",
    "-march=x86-64-v3",
    "-march=x86-64-v2",
)

CONFIGS = {
    "tools/LINUX-X86_64-GFORTRAN-MPI-PIC": "-march=x86-64",
    "tools/LINUX-AARCH64-GFORTRAN-MPI-PIC": "-march=armv8-a",
}

AUDIT_SHA = "06638f2ece83b12731b387ba819b67dc2151556b"
OLD_AUDIT_SHA = "8957850fd15005ea460ba336781ed152f7ee3d92"


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


def march_tokens(code: str) -> list[str]:
    return [token.rstrip("'\"") for token in re.findall(r"-march=\S+", code)]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_isa_portable.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    failures: list[str] = []

    for rel, baseline in CONFIGS.items():
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
        tokens = march_tokens(code)
        require(
            len(tokens) >= 4 and all(token == baseline for token in tokens),
            f"{rel} march tokens {tokens} are not exactly {baseline}",
            failures,
        )
        for flag in FORBIDDEN:
            require(flag not in code, f"{rel} contains {flag}", failures)
        require("compiler default" not in text.lower(), f"{rel} claims a compiler default ISA", failures)
        require("does not pass -march" not in text, f"{rel} claims it leaves -march unset", failures)

    x86_path = root / "tools/LINUX-X86_64-GFORTRAN-MPI-PIC"
    if x86_path.is_file():
        x86 = x86_path.read_text(encoding="utf-8")
        require("x86-64-v2" in x86 and "x86-64-v3" in x86, "x86_64 config does not name the RHEL GCC defaults", failures)

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
    # mkrst writes the rst during the docs build. It is not committed.
    pages = [("org", how)]
    how_rst = root / "docs/source/howto/opencpmd-archive.rst"
    if how_rst.is_file():
        pages.append(("rst", how_rst.read_text(encoding="utf-8")))
    for label, text in pages:
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
        require("-march=x86-64" in text, f"opencpmd-archive {label} lacks -march=x86-64", failures)
        require("-march=armv8-a" in text, f"opencpmd-archive {label} lacks -march=armv8-a", failures)
        require("x86-64-v2" in text and "x86-64-v3" in text, f"opencpmd-archive {label} lacks the RHEL GCC defaults", failures)
        require("compiler default" not in text.lower(), f"opencpmd-archive {label} claims a compiler default ISA", failures)
        require("10.1002/spe.3075" in text, f"opencpmd-archive {label} lacks the EESSI DOI", failures)

    readme = (root / "tools/README.md").read_text(encoding="utf-8")
    require("-march=x86-64" in readme and "-march=armv8-a" in readme, "tools/README.md lacks the baseline march flags", failures)
    require("x86-64-v2" in readme and "x86-64-v3" in readme, "tools/README.md lacks the RHEL GCC defaults", failures)

    news = (root / "docs/newsfragments/+portable-fp-contract.changed.md").read_text(encoding="utf-8")
    require("-march=x86-64" in news and "-march=armv8-a" in news, "newsfragment lacks the baseline march flags", failures)
    require("x86-64-v2" in news and "x86-64-v3" in news, "newsfragment lacks the RHEL GCC defaults", failures)

    simd = root / "docs/orgmode/explanation/simd-paths.org"
    require(simd.is_file(), "missing docs/orgmode/explanation/simd-paths.org", failures)
    simd_pages = []
    if simd.is_file():
        simd_pages.append(("org", simd.read_text(encoding="utf-8")))
    simd_rst = root / "docs/source/explanation/simd-paths.rst"
    if simd_rst.is_file():
        simd_pages.append(("rst", simd_rst.read_text(encoding="utf-8")))
    for label, text in simd_pages:
        require("CPMDC_SIMD" in text, f"simd-paths {label} does not name CPMDC_SIMD", failures)
        require("-ffp-contract=off" in text, f"simd-paths {label} lacks -ffp-contract=off", failures)
        require("10.1016/j.cpc.2013.06.003" in text, f"simd-paths {label} lacks the Páll and Hess DOI", failures)
        require("compiler default" not in text.lower(), f"simd-paths {label} claims a compiler default ISA", failures)

    index = (root / "docs/orgmode/index.org").read_text(encoding="utf-8")
    require("explanation/simd-paths" in index, "org index does not list simd-paths", failures)
    require("-march=x86-64" in index and "-march=armv8-a" in index, "org index lacks the baseline march flags", failures)
    index_rst = root / "docs/source/index.rst"
    if index_rst.is_file():
        index_text = index_rst.read_text(encoding="utf-8")
        require("explanation/simd-paths" in index_text, "rst index does not list simd-paths", failures)
        require(
            "-march=x86-64" in index_text and "-march=armv8-a" in index_text,
            "rst index lacks the baseline march flags",
            failures,
        )

    workflow = (root / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    require(workflow.count(AUDIT_SHA) == 2, "ci.yml does not pin both ISA audit calls to the shipped SHA", failures)
    require(OLD_AUDIT_SHA not in workflow, "ci.yml still pins the I0 audit workflow", failures)
    require("actions: write" not in workflow, "ci.yml grants actions: write", failures)
    require("TODO(I0)" not in workflow, "ci.yml still has a TODO(I0) comment", failures)
    require("isa-allow.toml" in workflow, "ci.yml does not pass isa-allow.toml", failures)
    require("workflow_dispatch" in workflow, "ci.yml has no workflow_dispatch trigger", failures)
    require("ubuntu-24.04-arm" in workflow, "ci.yml has no ubuntu-24.04-arm job", failures)
    require(
        "tests/distro_cmocka_prefix.h" in workflow,
        "ci.yml does not include the distro cmocka prefix header",
        failures,
    )
    require("-include setjmp.h" not in workflow, "ci.yml includes setjmp.h before _GNU_SOURCE", failures)
    prefix = root / "tests/distro_cmocka_prefix.h"
    require(prefix.is_file(), "missing tests/distro_cmocka_prefix.h", failures)
    if prefix.is_file():
        prefix_text = prefix.read_text(encoding="utf-8")
        require(
            prefix_text.find("_GNU_SOURCE") < prefix_text.find("setjmp.h"),
            "distro cmocka prefix does not define _GNU_SOURCE before setjmp.h",
            failures,
        )
        require("stdarg.h" in prefix_text, "distro cmocka prefix does not include stdarg.h", failures)
    mpi_launch = (root / "tests/run_mpi_test.py").read_text(encoding="utf-8")
    require("--mca" in mpi_launch and "^patcher" in mpi_launch, "MPI launcher does not disable the memory patcher", failures)
    require("^ucx" in mpi_launch, "MPI launcher does not exclude UCX under ASan", failures)

    golden = root / "tests/golden/isa-hashes.toml"
    require(not golden.exists(), "tests/golden/isa-hashes.toml is the SHA-256 of zero bytes and must stay deleted", failures)

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
    print("path scalar: exercised (baseline march; no dispatched kernel to force)")
    print("path avx2: not exercised")
    print("path avx512: not exercised")
    print("path neon: not exercised")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
