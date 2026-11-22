AddressSanitizer leak detection stays on for serial C and Fortran tests
that do not preload the sanitizer runtime or run under
`tests/run_mpi_test.py`. `params-cp-dft-render` keeps its fixture paths
in `main` instead of leaving a heap copy for the process exit, and
`dlopen-leak` verifies that `dlopen` constructor leaks are reported.
`meson.build` sets `detect_leaks=0` for the four ASan preload tests
(`python-api`, `python-pytest`, `docs-examples`, and `tutorials-run`)
because CPython and compilers started from those tests do not free
their pools, and `tests/run_mpi_test.py` sets `detect_leaks=0` for the
MPI tests (`adopt-comm-stub`, `adopt-comm-shell`,
`orbital-limit-agree-mpi`, and `rank-consensus`, which is registered
only for the reference build) and the `pwscale` tests
(`pwscale-agreement`, `pwscale-logger`, `pwscale-fft-fail`,
`pwscale-agreement-mpi-*`, `pwscale-mpi-fault`,
`pwscale-collectives-*`). The pixi `test-stub-sanitize` task and the
GCC 13 CI job use that same split.

`tests/lsan.supp` is kept and narrowed to a single `leak:libibverbs.so`
rule, which `meson.build` passes through `LSAN_OPTIONS` only to the
`shared-dlopen-abi` test under AddressSanitizer on Linux. It covers the
driver that rdma-core's `libibverbs` registers while libfabric is loaded.
`shared-dlopen-abi` no longer suppresses the `ld-linux` allocation in
`_dl_catch_exception`.
