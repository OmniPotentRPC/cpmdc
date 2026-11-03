`tests/test_extrapolate_wfn.c` declares `setenv`, `dup`, `dup2` and `fileno` through `_POSIX_C_SOURCE`, so it compiles under the `c11` standard of the build with GCC 14.
That compiler rejects the implicit declaration of `setenv` that the test used to rely on.
