"""scripts/abi_inventory.py fails on every kind of C ABI mismatch.

Each test runs the generator on a copy of the headers and the inventory,
so an injected mismatch never touches the checkout. A fake nm stands in for
the library; the last tests read the built libcpmdc named by CPMDC_LIBRARY.
"""

from __future__ import annotations

import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = "scripts/abi_inventory.py"
INVENTORY = "schema/inventory/surface_map_abi.txt"
COPIED = (
    SCRIPT,
    "include/cpmdc.h",
    "include/cpmdc_features.h",
    INVENTORY,
)
FAKE_NM = """\
import sys
from pathlib import Path

under = "_" if sys.platform == "darwin" else ""
names = Path(__file__).with_name("symbols.txt").read_text().split()
for i, name in enumerate(names):
    version = "@@CPMDC_0" if i == 0 and sys.platform != "darwin" else ""
    print(f"{0x1000 + i:016x} T {under}{name}{version}")
print(f"{0x9000:016x} t {under}cpmdc_local_helper")
print(f"{0x9100:016x} t {under}unrelated_local")
print("                 U cpmdc_undefined_import")
print("                 U unrelated_import")
"""


def load_generator(root: Path):
    spec = importlib.util.spec_from_file_location("abi_inventory", root / SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def tree(tmp_path: Path) -> Path:
    root = tmp_path / "src"
    for relative in COPIED:
        (root / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, root / relative)
    return root


@pytest.fixture
def nm(tmp_path: Path, tree: Path) -> Path:
    """Fake nm exporting exactly what the headers declare."""
    generator = load_generator(tree)
    symbols = generator.public_symbols(tree)
    fake = tmp_path / "nm" / "fake_nm.py"
    fake.parent.mkdir()
    fake.write_text(FAKE_NM)
    (fake.parent / "symbols.txt").write_text("\n".join(sorted(symbols)) + "\n")
    return fake


def run(tree: Path, mode: str, library: Path, nm: Path | None = None):
    env = dict(os.environ)
    if nm is not None:
        env["NM"] = f"{shlex.quote(sys.executable)} {shlex.quote(str(nm))}"
    return subprocess.run(
        [sys.executable, str(tree / SCRIPT), mode, str(library), "--root", str(tree)],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


def run_mode(tree: Path, mode: str, out: Path):
    return subprocess.run(
        [sys.executable, str(tree / SCRIPT), mode, str(out), "--root", str(tree)],
        capture_output=True,
        text=True,
        check=False,
    )


def inject_declaration(tree: Path, name: str) -> None:
    header = tree / "include/cpmdc_features.h"
    text = header.read_text()
    tail = "#ifdef __cplusplus\n}\n#endif\n"
    assert text.endswith(tail)
    header.write_text(text[: -len(tail)] + f"CPMDC_API int {name}(void);\n\n" + tail)


def add_symbol(nm: Path, name: str) -> None:
    symbols = nm.with_name("symbols.txt")
    symbols.write_text(symbols.read_text() + name + "\n")


def test_committed_inventory_matches_headers(tree: Path, nm: Path) -> None:
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 0, done.stdout + done.stderr


def test_write_is_sorted_and_deterministic(tree: Path, nm: Path) -> None:
    (tree / INVENTORY).write_text("")
    assert run(tree, "--write", tree / "libcpmdc.so", nm).returncode == 0
    first = (tree / INVENTORY).read_bytes()
    assert run(tree, "--write", tree / "libcpmdc.so", nm).returncode == 0
    assert (tree / INVENTORY).read_bytes() == first
    assert first == (ROOT / INVENTORY).read_bytes()
    text = first.decode()
    rows = [line.split("|")[1].strip() for line in text.splitlines()
            if line.startswith("| =")]
    assert rows and rows == sorted(rows)
    assert "Internal exports" not in text
    assert "cpmdc_local_helper" not in text
    assert "unrelated_local" not in text


def test_edited_inventory_fails_and_names_the_command(tree: Path, nm: Path) -> None:
    path = tree / INVENTORY
    path.write_text(path.read_text().replace("Parity: OK", "Parity: edited"))
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert "-Parity: edited" in done.stdout and "+Parity: OK" in done.stdout
    assert "is stale; run meson compile -C" in done.stderr
    assert "abi-inventory" in done.stderr


def test_declaration_without_export_fails(tree: Path, nm: Path) -> None:
    inject_declaration(tree, "cpmdc_injected_declaration")
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert (
        "declared in a public header but not exported: cpmdc_injected_declaration"
        in done.stderr
    )


def test_export_without_declaration_fails(tree: Path, nm: Path) -> None:
    add_symbol(nm, "cpmdc_injected_export")
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert "exported but declared in no public header: cpmdc_injected_export" in done.stderr


def test_write_still_fails_on_a_parity_gap(tree: Path, nm: Path) -> None:
    add_symbol(nm, "cpmdc_injected_export")
    done = run(tree, "--write", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    text = (tree / INVENTORY).read_text()
    assert "Parity: FAIL" in text and "cpmdc_injected_export" in text


def test_former_internal_export_fails(tree: Path, nm: Path) -> None:
    add_symbol(nm, "cpmdc_embed_init")
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert "exported but declared in no public header: cpmdc_embed_init" in done.stderr


def test_foreign_export_fails(tree: Path, nm: Path) -> None:
    add_symbol(nm, "capn_init_mem")
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert "exported but not a cpmdc_ symbol: capn_init_mem" in done.stderr


def test_declaration_without_export_macro_fails(tree: Path, nm: Path) -> None:
    header = tree / "include/cpmdc_features.h"
    header.write_text(
        header.read_text().replace(
            "CPMDC_API size_t cpmdc_feature_count", "size_t cpmdc_feature_count"
        )
    )
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert (
        "declared without CPMDC_API in cpmdc_features.h: cpmdc_feature_count"
        in done.stderr
    )


def test_export_lists_name_exactly_the_public_symbols(tree: Path, tmp_path: Path) -> None:
    public = load_generator(tree).public_symbols(tree)
    assert len(public) == len(set(public)) > 0
    script = tmp_path / "libcpmdc.map"
    done = run_mode(tree, "--version-script", script)
    assert done.returncode == 0, done.stderr
    text = script.read_text()
    globals_, locals_ = text.split("global:", 1)[1].split("local:", 1)
    assert [s.strip(" ;") for s in globals_.split()] == public
    assert locals_.split() == ["*;", "};"]
    exp = tmp_path / "libcpmdc.exp"
    done = run_mode(tree, "--exported-symbols-list", exp)
    assert done.returncode == 0, done.stderr
    assert exp.read_text().split() == [f"_{s}" for s in public]


def test_since_tag_is_recorded(tree: Path, nm: Path) -> None:
    header = tree / "include/cpmdc_features.h"
    header.write_text(
        header.read_text().replace(
            "/** @brief Number of entries", "/** @since ABI 7\n * @brief Number of entries"
        )
    )
    assert run(tree, "--write", tree / "libcpmdc.so", nm).returncode == 0
    rows = [line for line in (tree / INVENTORY).read_text().splitlines()
            if line.startswith("| =cpmdc_feature_")]
    assert [row.split("|")[3].strip() for row in rows] == ["7", "-", "-"]


def built_library() -> Path:
    named = os.environ.get("CPMDC_LIBRARY")
    if not named or not Path(named).is_file():
        pytest.skip("set CPMDC_LIBRARY to the built libcpmdc")
    if shutil.which(shlex.split(os.environ.get("NM", "nm"))[0]) is None:
        pytest.skip("no nm on PATH")
    return Path(named)


def test_built_library_matches_inventory(tree: Path) -> None:
    done = run(tree, "--check", built_library())
    assert done.returncode == 0, done.stdout + done.stderr


def test_built_library_exports_only_the_public_c_abi(tree: Path) -> None:
    library = built_library()
    generator = load_generator(tree)
    exported = generator.exported_symbols(library)
    public = set(generator.public_symbols(tree))
    assert sorted(exported - public) == []
    assert exported == public


def test_built_library_rejects_injected_declaration(tree: Path) -> None:
    library = built_library()
    inject_declaration(tree, "cpmdc_injected_declaration")
    done = run(tree, "--check", library)
    assert done.returncode == 1
    assert "not exported: cpmdc_injected_declaration" in done.stderr
    assert "abi-inventory" in done.stderr
