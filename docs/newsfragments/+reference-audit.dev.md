The `reference-audit` test, also run by `pixi run test-stub` with the rest of the Meson suite, fails when a public name has no contract entry in `reference/contracts.org`, or when an entry names a symbol that no longer exists.
It reads the C names from the headers that `meson install` places (the install plan of the build directory) and the Python names from the imported package.
A C function needs entries for its units, array layout, ownership, threads and MPI ranks, return and error codes, Python mapping, and defaults.
A Python name needs its signature, units, shape and dtype, ownership and threads, the exceptions it raises, and its defaults; the signature must equal what `inspect` reports, so a changed default or a renamed parameter fails until the page follows.
