"""scripts/ci/fetch_pseudopotentials.sh --fhi and --pbe."""

from __future__ import annotations

import os
import stat
import subprocess
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "ci" / "fetch_pseudopotentials.sh"
OXYGEN = "747670be5d6f008d921dc74551dc0705f2c4b34132e79590c482b532171f8cd1"
HYDROGEN = "b6ed642482602f477eb8c167dce4a2c6d254b89522685b5b9ad8f1e904c50db4"
O_PBE = "d26bb0bfa19ddaabe95d21169f50fa25915631e84df0d157e93bb32a7a9298ad"
H_PBE = "8ff874f26491f2b6e482540342ae0ad85c9c7023936f3b050baa4ddced459988"
REGTESTS = (
    "O_MT_BLYP.psp",
    "H_CVB_BLYP.psp",
    "Si_MT_BLYP.psp",
    "N_MT_BLYP.psp",
)


def test_fhi_names_the_checksums_and_rejects_an_unknown_option(tmp_path: Path) -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "--fhi" in text
    assert "--pbe" in text
    assert OXYGEN in text
    assert HYDROGEN in text
    assert O_PBE in text
    assert H_PBE in text
    assert "sha256sum -c - >/dev/null" not in text
    assert "checksum mismatch" in text
    done = subprocess.run(
        ["bash", str(SCRIPT), "--not-fhi", str(tmp_path / "pp")],
        capture_output=True, text=True, check=False,
    )
    assert done.returncode == 2
    assert "unknown option" in done.stderr
    assert not (tmp_path / "pp").exists()


def _shim_dir(tmp_path: Path) -> Path:
    text = SCRIPT.read_text(encoding="utf-8")
    commit = text.split("commit=", 1)[1].splitlines()[0]
    bindir = tmp_path / "bin"
    bindir.mkdir()
    git = bindir / "git"
    git.write_text(
        "#!/bin/sh\n"
        "case \" $* \" in\n"
        f"  *rev-parse*) echo {commit} ;;\n"
        "  *) echo \"git shim refused: $*\" >&2; exit 1 ;;\n"
        "esac\n",
        encoding="utf-8",
    )
    curl = bindir / "curl"
    curl.write_text(
        textwrap.dedent(
            """\
            #!/bin/sh
            out=
            url=
            prev=
            for arg in "$@"; do
              if [ "$prev" = "-o" ]; then
                out=$arg
                prev=
                continue
              fi
              case "$arg" in
                -o) prev=-o ;;
                http*) url=$arg ;;
              esac
            done
            printf '%s\\n' "$url" >> "$CURL_LOG"
            case "$url" in
              *pseudo.tar.gz)
                python3 - "$out" << 'PY'
            import io, sys, tarfile
            data = b"not the pseudopotential\\n"
            with tarfile.open(sys.argv[1], "w") as tar:
                for name in ("./pseudo/O_MT_PBE.psp", "./pseudo/H_MT_PBE.psp"):
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))
            PY
                ;;
              *)
                printf 'not the file\\n' > "$out"
                ;;
            esac
            """
        ),
        encoding="utf-8",
    )
    git.chmod(git.stat().st_mode | stat.S_IEXEC)
    curl.chmod(curl.stat().st_mode | stat.S_IEXEC)
    return bindir


def _seed(tmp_path: Path) -> tuple[Path, Path]:
    clone = tmp_path / "regtests"
    library = clone / "tests" / "PP_LIBRARY"
    library.mkdir(parents=True)
    (clone / ".git").mkdir()
    for name in REGTESTS:
        (library / name).write_bytes(b"regtest\n")
    return tmp_path / "pp", clone


def _run(tmp_path: Path, flag: str) -> subprocess.CompletedProcess[str]:
    dest, clone = _seed(tmp_path)
    env = os.environ.copy()
    env["PATH"] = f"{_shim_dir(tmp_path)}:{env['PATH']}"
    env["CURL_LOG"] = str(tmp_path / "curl.log")
    return subprocess.run(
        ["bash", str(SCRIPT), flag, str(dest), str(clone)],
        capture_output=True, text=True, check=False, env=env,
    )


def test_fhi_checksum_mismatch_runs_the_downloader(tmp_path: Path) -> None:
    done = _run(tmp_path, "--fhi")
    log = (tmp_path / "curl.log").read_text(encoding="utf-8")
    assert "O.pbe-mt_fhi.UPF" in log
    assert done.returncode != 0
    assert "checksum mismatch for O.pbe-mt_fhi.UPF" in done.stderr
    assert "FAILED" in done.stdout + done.stderr


def test_pbe_checksum_mismatch_runs_the_downloader(tmp_path: Path) -> None:
    done = _run(tmp_path, "--pbe")
    log = (tmp_path / "curl.log").read_text(encoding="utf-8")
    assert "pseudo.tar.gz" in log
    assert done.returncode != 0
    assert "checksum mismatch for O_MT_PBE.psp" in done.stderr
    assert "FAILED" in done.stdout + done.stderr
