`pip install cpmdc` installs the Python package with its own `libcpmdc` and schema; wheels cover Linux x86_64 and aarch64 and macOS arm64, without OpenCPMD or MPI.
A periodic cell needs no other tool, because the package encodes each step as `ForceInput` itself; the `capnp` command is needed only to encode a `CPMDParams` text file.
