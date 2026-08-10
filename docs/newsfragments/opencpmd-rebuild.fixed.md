Rebuild patched OpenCPMD module objects and `libcpmd.a` through an explicit
dependency-target helper, avoiding a `make lib` no-op against the existing
`lib/` directory.
