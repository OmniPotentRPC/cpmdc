#!/usr/bin/env python3
"""Portable OpenCPMD configs and the ISA allowlist rules."""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
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

# The auditor applies each [[allow]] symbol with re.search, so an entry
# must start with an anchored Rust path. cpmdc's C and Fortran symbols
# contain no "::".
ALLOW_KEYS = {"symbol", "features", "reason"}
RUST_PATH = re.compile(r"\^<?([a-z_][a-z0-9_]*)::")
WILDCARDS = set(".*+[]{}\\")
DISPATCH_CHECKS = ("is_x86_feature_detected!", "is_available()", "CPUID")

# Entries the rules must reject: (entry, the reason it is rejected).
BAD_ALLOWS = (
    ({"symbol": "^minimage::.*", "features": ["avx"], "reason": "is_x86_feature_detected!(avx)"}, "wildcard"),
    ({"symbol": "minimage::fused::avx::", "features": ["avx"], "reason": "is_x86_feature_detected!(avx)"}, "unanchored"),
    ({"symbol": "^cpmdc_geometry_rmsd$", "features": ["avx2"], "reason": "is_x86_feature_detected!(avx2)"}, "cpmdc symbol"),
    ({"symbol": "^minimage::a::|^b", "features": ["avx"], "reason": "is_x86_feature_detected!(avx)"}, "top-level |"),
    ({"symbol": "^minimage::fused::avx::", "features": ["avx"], "reason": "minimage avx kernels"}, "no run-time check"),
    ({"symbol": "^minimage::fused::avx::", "features": ["avx", "fma"], "reason": "is_x86_feature_detected!(avx)"}, "feature not in reason"),
    ({"symbol": "^minimage::fused::avx::", "features": ["avx"]}, "no reason"),
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


def march_tokens(code: str) -> list[str]:
    return [token.rstrip("'\"") for token in re.findall(r"-march=\S+", code)]


def cpmdc_identifiers(root: Path) -> set[str]:
    names: set[str] = set()
    for pattern in ("src/*.c", "src/*.h", "src/*.f90", "src/*.F90", "include/*.h"):
        for path in root.glob(pattern):
            names.update(re.findall(r"\b[A-Za-z_]\w*\b", path.read_text(encoding="utf-8", errors="replace")))
    return names


def wildcard(symbol: str) -> str | None:
    """Return the first construct that widens a pattern past the names it spells out."""
    depth = 0
    previous = ""
    for char in symbol:
        if char in WILDCARDS or (char == "?" and previous != ")"):
            return char
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "|" and depth == 0:
            return char
        previous = char
    return None


def allow_errors(entry: dict, project: str, identifiers: set[str]) -> list[str]:
    symbol = entry.get("symbol")
    label = f"isa-allow.toml [[allow]] {symbol!r}"
    errors = []
    if set(entry) != ALLOW_KEYS:
        errors.append(f"{label} has keys {sorted(entry)}, not {sorted(ALLOW_KEYS)}")
    if not isinstance(symbol, str) or not symbol:
        return errors + [f"{label} names no symbol"]
    path = RUST_PATH.match(symbol)
    if path is None or symbol.count("^") != 1:
        errors.append(f"{label} does not start with one anchored Rust path")
    elif path.group(1).startswith(project):
        errors.append(f"{label} allows a {project} symbol")
    if not symbol.endswith(("$", "::")) or "$" in symbol[:-1]:
        errors.append(f"{label} does not end at a function name or a module path")
    construct = wildcard(symbol)
    if construct is not None:
        errors.append(f"{label} uses the wildcard {construct!r}")
    try:
        compiled = re.compile(symbol)
    except re.error as exc:
        errors.append(f"{label} is not a regex: {exc}")
    else:
        matched = sorted(name for name in identifiers if compiled.search(name))
        if matched:
            errors.append(f"{label} allows a {project} symbol: {matched[0]}")
    features = entry.get("features")
    if (
        not isinstance(features, list)
        or not features
        or not all(isinstance(item, str) and item for item in features)
        or len(set(features)) != len(features)
    ):
        errors.append(f"{label} features {features!r} are not distinct feature names")
        features = []
    reason = entry.get("reason")
    if not isinstance(reason, str) or not reason.strip():
        return errors + [f"{label} states no reason"]
    if not any(check in reason for check in DISPATCH_CHECKS):
        errors.append(f"{label} reason names none of the run-time checks {DISPATCH_CHECKS}")
    unstated = [item for item in features if not re.search(rf"\b{re.escape(item)}\b", reason)]
    if unstated:
        errors.append(f"{label} reason does not account for {unstated}")
    return errors


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

    allow = tomllib.loads((root / "isa-allow.toml").read_text(encoding="utf-8"))
    require(allow.get("baseline") == "x86-64", "isa-allow.toml baseline is not x86-64", failures)
    require("vendored" not in allow, "isa-allow.toml vendors a library that CI links dynamically", failures)
    project = re.search(r"project\(\s*'([^']+)'", meson)
    require(project is not None, "meson.build names no project", failures)
    project_name = project.group(1) if project else "cpmdc"
    identifiers = cpmdc_identifiers(root)
    require(f"{project_name}_geometry_rmsd" in identifiers, "no cpmdc identifiers found under src/ and include/", failures)
    allows = allow.get("allow", [])
    require(isinstance(allows, list) and len(allows) > 0, "isa-allow.toml has no [[allow]] entry for the minimage kernels", failures)
    for entry in allows if isinstance(allows, list) else []:
        for error in allow_errors(entry, project_name, identifiers):
            failures.append(error)
    for entry, why in BAD_ALLOWS:
        require(bool(allow_errors(entry, project_name, identifiers)), f"allow rules accept a bad entry ({why})", failures)

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
    simd_lines = re.findall(r'echo "cpmdc SIMD paths: ([^"]*)"', workflow)
    require(
        len(simd_lines) == workflow.count("- name: CPU flags"),
        "ci.yml CPU flags steps do not each print the cpmdc SIMD paths",
        failures,
    )
    for line in simd_lines:
        require("minimage" in line, f"ci.yml reports cpmdc SIMD paths {line!r} without the minimage kernels", failures)
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
    print("path scalar: exercised (baseline march; the suite calls the scalar minimage one-pair displacement)")
    print("minimage AVX, AVX2 and AVX-512 kernels: linked, picked at run time, not called by this suite")
    print("path avx2: not exercised")
    print("path avx512: not exercised")
    print("path neon: not exercised")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
