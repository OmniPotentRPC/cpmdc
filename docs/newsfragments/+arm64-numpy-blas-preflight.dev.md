The pytest session now checks NumPy's matrix product against a plain loop
before collecting any test, and stops with a message naming the remedies when
the two disagree. On aarch64 under `qemu-user` 8.2 with `-cpu max`, OpenBLAS
selects its `armv9sme` kernel and `@` returns wrong values, which made the
NumPy Ewald reference of `test_mechanical_embedding_matches_numpy_ewald` sum
to zero and fail on its energy assertion while the library was correct.
