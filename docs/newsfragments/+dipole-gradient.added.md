`cpmdc_last_dipole()` and its session form return a molecular dipole only when
linked OpenCPMD recomputed it for the current geometry. Python `Session` and
`Potential` expose that value, and `CPMDCCalculator` provides ASE's `dipole`
property in e*Angstrom. Reference, isolated, k-point, and unrequested
evaluations report the property as unavailable instead of returning a zero or
stale value.

The dipole is a Berry-phase quantity of the periodic cell. It is defined only
for a neutral cell, and only modulo the lattice quantum: it can jump by an
integer number of electron charges times a lattice vector between geometries.
`cpmdc_session_calculate_au()` has no dipole argument, because its argument
list is fixed for the LAMMPS, GROMACS, and i-PI callers that bind to it; read
`cpmdc_session_last_dipole()` after it returns 0.
`Potential.caps()` reports `dipole` as true only when the evaluator can produce
one: never for the reference evaluator, and for OpenCPMD until an evaluation
comes back without a dipole. `Session.dipole_possible()` gives the same answer.
