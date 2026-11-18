`test_cpmd_restart` makes its scratch directory under `TMPDIR`, falling back to the C library's `P_tmpdir`, and removes it at the end.
