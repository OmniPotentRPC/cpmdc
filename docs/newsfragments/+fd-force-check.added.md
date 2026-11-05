`cpmdc-fd-check` compares the forces of a session with central differences of its energy, after the CP2K `DEBUG_FORCES` check with `DX = 1e-3` Bohr.
`--stress` adds the symmetric-strain check of `PotentialResult.stress` after `DEBUG_STRESS_TENSOR`.
The tool is built with the tests and is not installed.
On the reference evaluator the forces agree to rounding and the toy stress is reported as a mismatch.
