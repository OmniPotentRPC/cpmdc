#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "usage: encode_capnp.py CAPNP SCHEMA TYPE INPUT OUTPUT",
            file=sys.stderr,
        )
        return 2
    capnp, schema, type_name, input_path, output_path = sys.argv[1:]
    text = Path(input_path).read_bytes()
    encoded = subprocess.run(
        [capnp, "encode", schema, type_name],
        input=text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if encoded.returncode != 0:
        sys.stderr.buffer.write(encoded.stderr)
        return encoded.returncode
    Path(output_path).write_bytes(encoded.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
