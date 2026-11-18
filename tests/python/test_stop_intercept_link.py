"""The serial stop intercept and MPI sentinel export list link on ELF and ld64.

ELF leaves an undefined weak cpmd_stopgm_hook unresolved and null. ld64, which
links the macOS wheel, reports that reference as an undefined symbol. The
serial object therefore does not mention the hook. An OpenCPMD build still
takes the weak reference, and a strong definition from the archive satisfies
it. ld64 also requires every literal name in an exported-symbols list to be
defined, while wildcard entries match only defined symbols: the serial list
names the literal public C ABI alone, and an OpenCPMD build adds wildcard
patterns for the Open MPI and MPICH Fortran sentinel common blocks.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "src" / "cpmd_stop_intercept.c"
ABI_INVENTORY = ROOT / "scripts" / "abi_inventory.py"
HOOK = "cpmd_stopgm_hook"
ELF_MAGIC = b"\x7fELF"
MPI_SENTINEL_SETS: dict[str, tuple[str, ...]] = {
    "openmpi": (
        "mpi_fortran_argv_null_",
        "mpi_fortran_argvs_null_",
        "mpi_fortran_bottom_",
        "mpi_fortran_errcodes_ignore_",
        "mpi_fortran_in_place_",
        "mpi_fortran_status_ignore_",
        "mpi_fortran_statuses_ignore_",
        "mpi_fortran_unweighted_",
        "mpi_fortran_weights_empty_",
    ),
    "mpich": ("mpipriv1_", "mpipriv2_", "mpiprivc_", "mpifcmb5_", "mpifcmb9_"),
}


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _nm_type(obj: Path, name: str) -> str | None:
    done = _run(["nm", str(obj)])
    assert done.returncode == 0, done.stderr
    for line in done.stdout.splitlines():
        fields = line.split()
        if fields and fields[-1] == name:
            return fields[-2]
    return None


def _require_elf_host(tmp_path: Path) -> None:
    if not sys.platform.startswith(("linux", "freebsd", "openbsd", "netbsd")):
        pytest.skip(f"ELF link check requires an ELF host, not {sys.platform}")
    if shutil.which("gcc") is None or shutil.which("nm") is None:
        pytest.skip("gcc and nm are required")
    probe_c = tmp_path / "elf_probe.c"
    probe_o = tmp_path / "elf_probe.o"
    probe_c.write_text("int cpmdc_elf_probe(void) { return 0; }\n", encoding="utf-8")
    done = _run(["gcc", "-c", "-fPIC", "-std=c11", str(probe_c), "-o", str(probe_o)])
    if done.returncode != 0 or probe_o.read_bytes()[:4] != ELF_MAGIC:
        pytest.skip("host gcc does not emit ELF objects")


def _compile_elf(obj: Path, *defines: str) -> None:
    done = _run(["gcc", "-c", "-fPIC", "-std=c11", str(SOURCE), "-o", str(obj), *defines])
    assert done.returncode == 0, done.stderr


def _ld64() -> Path | None:
    found = shutil.which("ld64.lld")
    if found:
        return Path(found)
    clang = shutil.which("clang")
    if clang is not None:
        probed = _run([clang, "-print-prog-name=ld64.lld"])
        if probed.returncode == 0:
            candidate = Path(probed.stdout.strip())
            if candidate.is_absolute() and os.access(candidate, os.X_OK):
                return candidate
    rustc = shutil.which("rustc")
    if rustc is None:
        return None
    sysroot = _run([rustc, "--print", "sysroot"])
    if sysroot.returncode != 0:
        return None
    matches = sorted(Path(sysroot.stdout.strip()).glob("lib/rustlib/*/bin/gcc-ld/ld64.lld"))
    return matches[0] if matches else None


@dataclass(frozen=True)
class MachOToolchain:
    clang: Path
    ld64: Path | None
    target: str | None

    def compile(self, source: Path, obj: Path, *defines: str) -> subprocess.CompletedProcess[str]:
        cmd = [str(self.clang), "-c", "-fPIC", "-std=c11", str(source), "-o", str(obj), *defines]
        if self.target is not None:
            cmd.insert(1, f"--target={self.target}")
        return _run(cmd)

    def link_dylib(
        self, out: Path, *objects: Path, export_list: Path | None = None
    ) -> subprocess.CompletedProcess[str]:
        if self.ld64 is not None:
            cmd = [
                str(self.ld64), "-dylib", "-arch", "arm64",
                "-platform_version", "macos", "15.0.0", "15.0.0",
                "-o", str(out), *(str(obj) for obj in objects),
            ]
            if export_list is not None:
                cmd[1:1] = ["-exported_symbols_list", str(export_list)]
            return _run(cmd)
        cmd = [str(self.clang), "-dynamiclib", "-o", str(out), *(str(obj) for obj in objects)]
        if export_list is not None:
            cmd[2:2] = [f"-Wl,-exported_symbols_list,{export_list}"]
        return _run(cmd)


def _probe_macho(tc: MachOToolchain, work: Path) -> str | None:
    work.mkdir(parents=True, exist_ok=True)
    src, obj, dylib = work / "empty.c", work / "probe.o", work / "probe.dylib"
    src.write_text("int cpmdc_macho_probe(void) { return 0; }\n", encoding="utf-8")
    compiled = tc.compile(src, obj)
    if compiled.returncode != 0:
        return compiled.stderr
    linked = tc.link_dylib(dylib, obj)
    return None if linked.returncode == 0 else linked.stderr


@pytest.fixture
def macho(tmp_path: Path) -> MachOToolchain:
    clang_found = shutil.which("clang")
    if clang_found is None:
        pytest.skip("clang is not on PATH")
    clang = Path(clang_found)
    linker = _ld64()
    if linker is not None:
        cross = MachOToolchain(clang=clang, ld64=linker, target="arm64-apple-macos15")
        err = _probe_macho(cross, tmp_path / "probe_lld")
        if err is None:
            return cross
    if sys.platform == "darwin":
        native = MachOToolchain(clang=clang, ld64=None, target=None)
        err = _probe_macho(native, tmp_path / "probe_native")
        if err is None:
            return native
        pytest.skip(f"clang -dynamiclib probe failed: {err}")
    pytest.skip("ld64.lld is not on PATH and host is not darwin")


def test_serial_elf_object_does_not_reference_the_stop_hook(tmp_path: Path) -> None:
    _require_elf_host(tmp_path)
    obj = tmp_path / "serial.o"
    _compile_elf(obj)
    assert _nm_type(obj, HOOK) is None
    shared = tmp_path / "libserial.so"
    linked = _run(["gcc", "-shared", "-Wl,-z,defs", "-o", str(shared), str(obj)])
    assert linked.returncode == 0, linked.stderr
    undefined = _run(["nm", "-D", "-u", str(shared)])
    assert undefined.returncode == 0, undefined.stderr
    names = {line.split()[-1] for line in undefined.stdout.splitlines() if line.split()}
    assert HOOK not in names
    assert not any(name.startswith(("mpifcmb", "mpipriv", "mpi_fortran_")) for name in names)


def test_opencpmd_elf_object_keeps_a_weak_hook(tmp_path: Path) -> None:
    _require_elf_host(tmp_path)
    obj = tmp_path / "opencpmd.o"
    _compile_elf(obj, "-DCPMDC_HAS_CPMD=1")
    assert _nm_type(obj, HOOK) == "w"
    provider = tmp_path / "provider.c"
    provider.write_text(
        "typedef int (*cpmd_stopgm_hook_fn)(int);\n"
        "cpmd_stopgm_hook_fn cpmd_stopgm_hook = 0;\n",
        encoding="utf-8",
    )
    provider_o = tmp_path / "provider.o"
    compiled = _run(["gcc", "-c", "-fPIC", "-std=c11", str(provider), "-o", str(provider_o)])
    assert compiled.returncode == 0, compiled.stderr
    shared = tmp_path / "libopencpmd.so"
    linked = _run(["gcc", "-shared", "-Wl,-z,defs", "-o", str(shared), str(obj), str(provider_o)])
    assert linked.returncode == 0, linked.stderr
    assert _nm_type(shared, HOOK) in {"B", "D", "T"}


def test_serial_macho_link_has_no_undefined_hook(tmp_path: Path, macho: MachOToolchain) -> None:
    obj = tmp_path / "serial.o"
    done = macho.compile(SOURCE, obj)
    assert done.returncode == 0, done.stderr
    linked = macho.link_dylib(tmp_path / "libcpmdc.dylib", obj)
    assert linked.returncode == 0, linked.stderr
    listed = tmp_path / "public.exp"
    listed.write_text("_cpmdc_stop_arm\n_cpmdc_stop_disarm\n", encoding="utf-8")
    listed_link = macho.link_dylib(tmp_path / "listed.dylib", obj, export_list=listed)
    assert listed_link.returncode == 0, listed_link.stderr
    sentinels = tmp_path / "sentinels.exp"
    sentinels.write_text(listed.read_text(encoding="utf-8") + "_mpifcmb5_\n_mpipriv1_\n", encoding="utf-8")
    rejected = macho.link_dylib(tmp_path / "rejected.dylib", obj, export_list=sentinels)
    assert rejected.returncode != 0
    assert "_mpifcmb5_" in rejected.stderr
    assert "_mpipriv1_" in rejected.stderr


def test_opencpmd_macho_link_resolves_a_strong_hook(tmp_path: Path, macho: MachOToolchain) -> None:
    obj = tmp_path / "opencpmd.o"
    done = macho.compile(SOURCE, obj, "-DCPMDC_HAS_CPMD=1")
    assert done.returncode == 0, done.stderr
    missing = macho.link_dylib(tmp_path / "missing.dylib", obj)
    assert missing.returncode != 0
    assert f"_{HOOK}" in missing.stderr
    provider = tmp_path / "provider.c"
    provider.write_text(
        "typedef int (*cpmd_stopgm_hook_fn)(int);\n"
        "cpmd_stopgm_hook_fn cpmd_stopgm_hook = 0;\n",
        encoding="utf-8",
    )
    provider_o = tmp_path / "provider.o"
    compiled = macho.compile(provider, provider_o)
    assert compiled.returncode == 0, compiled.stderr
    linked = macho.link_dylib(tmp_path / "libcpmdc.dylib", obj, provider_o)
    assert linked.returncode == 0, linked.stderr


@pytest.mark.parametrize("mpi_impl", sorted(MPI_SENTINEL_SETS))
def test_opencpmd_macho_export_list_accepts_each_mpi_symbol_set(
    tmp_path: Path, macho: MachOToolchain, mpi_impl: str
) -> None:
    serial_exp = tmp_path / "serial.exp"
    opencpmd_exp = tmp_path / "opencpmd.exp"
    for out, extra in ((serial_exp, ()), (opencpmd_exp, ("--mpi-fortran-sentinels",))):
        done = _run([sys.executable, str(ABI_INVENTORY), "--exported-symbols-list", str(out), *extra, "--root", str(ROOT)])
        assert done.returncode == 0, done.stderr
    public = [s.removeprefix("_") for s in serial_exp.read_text(encoding="utf-8").split()]
    sentinels = MPI_SENTINEL_SETS[mpi_impl]
    provider_c, provider_o = tmp_path / f"p_{mpi_impl}.c", tmp_path / f"p_{mpi_impl}.o"
    provider_c.write_text(
        "".join(f"void {s}(void) {{}}\n" for s in public)
        + "".join(f"int {s} = 0;\n" for s in sentinels)
        + "int cpmdc_internal_unlisted = 0;\n",
        encoding="utf-8",
    )
    assert macho.compile(provider_c, provider_o).returncode == 0
    dylib = tmp_path / f"libcpmdc_{mpi_impl}.dylib"
    linked = macho.link_dylib(dylib, provider_o, export_list=opencpmd_exp)
    assert linked.returncode == 0, linked.stderr

    refs = [public[0], *sentinels]
    use_c, use_o = tmp_path / f"u_{mpi_impl}.c", tmp_path / f"u_{mpi_impl}.o"
    use_c.write_text(
        "".join(f"extern int {s};\n" for s in refs)
        + "void *refs[] = {" + ", ".join(f"&{s}" for s in refs) + "};\n",
        encoding="utf-8",
    )
    assert macho.compile(use_c, use_o).returncode == 0
    assert macho.link_dylib(tmp_path / f"u_{mpi_impl}.dylib", use_o, dylib).returncode == 0

    hid_c, hid_o = tmp_path / f"h_{mpi_impl}.c", tmp_path / f"h_{mpi_impl}.o"
    hid_c.write_text("extern int cpmdc_internal_unlisted;\nvoid *r = &cpmdc_internal_unlisted;\n", encoding="utf-8")
    assert macho.compile(hid_c, hid_o).returncode == 0
    assert macho.link_dylib(tmp_path / f"h_{mpi_impl}.dylib", hid_o, dylib).returncode != 0


def test_darwin_meson_lists_sentinels_only_for_opencpmd() -> None:
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    start = text.index(".startswith('ld64')")
    window = text[start : text.index("cpmdc_shared = shared_library", start)]
    before, after = window.split("if has_cpmd", 1)
    assert "--exported-symbols-list" in before
    assert "--mpi-fortran-sentinels" not in before
    assert "--mpi-fortran-sentinels" in after
