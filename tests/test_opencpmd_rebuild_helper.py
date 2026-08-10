#!/usr/bin/env python3
"""Verify that the OpenCPMD helper rebuilds real dependency targets."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    repo = Path(sys.argv[1]).resolve()
    helper = repo / "tools" / "rebuild_opencpmd_embed.sh"

    with tempfile.TemporaryDirectory() as scratch:
        root = Path(scratch) / "cpmd-root"
        (root / "obj").mkdir(parents=True)
        (root / "lib").mkdir()
        (root / "Makefile").write_text("# test fixture\n")

        log = Path(scratch) / "make-args.json"
        fake_make = Path(scratch) / "make"
        fake_make.write_text(
            "#!/usr/bin/env python3\n"
            "import json, os, sys\n"
            "open(os.environ['MAKE_LOG'], 'w').write(json.dumps(sys.argv[1:]))\n"
        )
        fake_make.chmod(0o755)

        env = os.environ.copy()
        env["MAKE"] = str(fake_make)
        env["MAKE_LOG"] = str(log)
        subprocess.run(
            [str(helper), str(root), "3"],
            cwd=repo,
            env=env,
            check=True,
        )

        actual = json.loads(log.read_text())
        expected = [
            "-C",
            str(root / "obj"),
            "-f",
            str(root / "Makefile"),
            "-j",
            "3",
            "rwfopt_utils.mod.o",
            "updwf_utils.mod.o",
            str(root / "lib" / "libcpmd.a"),
        ]
        if actual != expected:
            raise AssertionError(f"make arguments {actual!r}, expected {expected!r}")

    print("OK: OpenCPMD helper rebuilds patched objects and archive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
