The first energy tutorial now compiles its host program with `mpicc` and links it with `mpif90`.
The wrappers are needed because the default `libcpmdc.so` links the Message Passing Interface (MPI) library whenever Meson finds one.
The tutorial transcripts now come from a run of version 0.3.0 against both the default build and a patched OpenCPMD archive.
