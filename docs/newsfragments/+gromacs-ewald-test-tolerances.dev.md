The `rcut` and `alpha` assertions of `test_mechanical_embedding_matches_numpy_ewald`
now check `1e-14` relative with no absolute slack, where `rel=1e-15` had been an
absolute `1e-12` in effect, and the other tolerances of that test carry their
error bounds.
