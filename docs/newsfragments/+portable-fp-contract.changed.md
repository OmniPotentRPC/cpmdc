The x86_64 OpenCPMD PIC file passes `-march=x86-64`, and the aarch64 file passes `-march=armv8-a`, on debug and release `FFLAGS` and `CFLAGS`.
GCC 11 on RHEL 9 defaults to x86-64-v2, and GCC 14 on RHEL 10 defaults to x86-64-v3 with AVX2 and FMA.
Both files also pass `-fPIC` and `-ffp-contract=off`.
cpmdc C and Fortran compile with `-ffp-contract=off` so a multiply and an add stay separate on every build target.
