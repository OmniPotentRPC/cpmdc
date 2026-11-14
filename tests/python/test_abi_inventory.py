"""scripts/abi_inventory.py fails on every kind of C ABI mismatch.

Each test runs the generator on a copy of the headers, the internal list,
and the inventory, so an injected mismatch never touches the checkout. A
fake nm stands in for the library; the last tests read the built libcpmdc
named by CPMDC_LIBRARY.
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
    "scripts/abi_internal_exports.txt",
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
print(f"{0x9100:016x} T {under}unrelated_export")
print("                 U cpmdc_undefined_import")
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
    """Fake nm exporting exactly what the headers and the internal list name."""
    generator = load_generator(tree)
    declarations, _ = generator.header_declarations(tree)
    symbols = {d.symbol for d in declarations} | set(generator.internal_exports(tree))
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


def inject_declaration(tree: Path, name: str) -> None:
    header = tree / "include/cpmdc_features.h"
    text = header.read_text()
    tail = "#ifdef __cplusplus\n}\n#endif\n"
    assert text.endswith(tail)
    header.write_text(text[: -len(tail)] + f"int {name}(void);\n\n" + tail)


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
    for section in text.split("* Internal exports"):
        rows = [line.split("|")[1].strip() for line in section.splitlines()
                if line.startswith("| =")]
        assert rows and rows == sorted(rows)
    assert "cpmdc_local_helper" not in text
    assert "unrelated_export" not in text


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


def test_stale_internal_entry_fails(tree: Path, nm: Path) -> None:
    internal = tree / "scripts/abi_internal_exports.txt"
    internal.write_text(internal.read_text() + "cpmdc_injected_internal\n")
    done = run(tree, "--check", tree / "libcpmdc.so", nm)
    assert done.returncode == 1
    assert "but not exported: cpmdc_injected_internal" in done.stderr


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


def test_built_library_rejects_injected_declaration(tree: Path) -> None:
    library = built_library()
    inject_declaration(tree, "cpmdc_injected_declaration")
    done = run(tree, "--check", library)
    assert done.returncode == 1
    assert "not exported: cpmdc_injected_declaration" in done.stderr
    assert "abi-inventory" in done.stderr
