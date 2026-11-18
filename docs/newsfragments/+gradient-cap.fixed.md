An evaluation no longer fails above 1365 atoms.
`cpmdc_session_calculate_au()` (LAMMPS, GROMACS, i-PI), the plain and session gradient and force calls, and the Cap'n Proto result calls all return the whole gradient above that size, where an earlier build of this branch rejected the call with `property gradient snapshot supports at most 1365 atoms`.
The deprecated `CPMDCPropertySnapshot.hessian[4096]` array still holds only 1365 atoms.
For a larger system it stays empty: `hessian_count` is 0 and `hessian` is zeroed, and `PotentialResult.gradient` is the negated `forces` in the units of the `ForceInput`.
`cpmdc_last_gradient()` and `cpmdc_session_last_gradient()` return -1 with `*count` set to the required length `n_atoms * 3`, also when `out` is null, so a caller can size a buffer.
The gradient itself is the one the evaluation call wrote to its own buffer. The library keeps no second, system-sized copy of every gradient per session.
