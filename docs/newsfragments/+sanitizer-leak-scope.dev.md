AddressSanitizer leak detection stays on for the C and Fortran tests.
`params-cp-dft-render` keeps its fixture paths in `main` instead of
leaving a heap copy for the process exit. Tests that preload `libasan`
into Python set `detect_leaks=0`, because CPython and a compiler started
from those tests do not free their pools. The pixi sanitizer task and
the GCC 13 CI job use that same split. `shared-dlopen-abi` suppresses
the 24 byte `ld-linux` allocation in `_dl_catch_exception` that GCC 14
LeakSanitizer reports during `dlopen`.
