The x86_64 and aarch64 OpenCPMD PIC configure files pass `-fPIC` and `-ffp-contract=off` and do not set `-march`.
cpmdc C and Fortran compile with `-ffp-contract=off` so a multiply and an add stay separate on every build target.
