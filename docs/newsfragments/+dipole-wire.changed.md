For rgpot and other consumers of `PotentialResult` and `CPMDCPropertySnapshot`, the dipole fields now say whether a dipole exists.
`CPMDCPropertySnapshot.dipole_count` is 3 only when OpenCPMD recomputed the dipole for the current geometry and 0 otherwise; it was always 3, with zeros from the reference evaluator and the stale `ddip%pdipole` from OpenCPMD.
`PotentialResult.dipole` is an empty list in the same cases and was three zeros.
When present, it is the molecular dipole `dipomod%moment%dmom` in electron Bohr, with the ionic term, in place of `ddip%pdipole`.
A consumer must read an empty `dipole` list or `dipole_count == 0` as "no dipole", not as a dipole of zero.
`CPMDCDipole` and `cpmdc_last_dipole()` carry the same `valid` flag, and Python `Potential.evaluate()` leaves out the `dipole` key.
The dipole is a Berry-phase quantity: defined only for a neutral periodic cell, and only modulo the lattice quantum.
The same release makes `PotentialResult.gradient` the negated `forces` in the units of the `ForceInput`.
An evaluation of more than 1365 atoms is rejected before it runs, so no wire result carries a gradient for such a system.
