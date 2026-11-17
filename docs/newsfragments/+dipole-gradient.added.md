`cpmdc_last_dipole()` and its session form return a molecular dipole only when
linked OpenCPMD recomputed it for the current geometry. Python `Session` and
`Potential` expose that value, and `CPMDCCalculator` provides ASE's `dipole`
property in e*Angstrom. Reference, isolated, k-point, and unrequested
evaluations report the property as unavailable instead of returning a zero or
stale value.

`cpmdc_last_gradient()` and its session form provide a correctly named view of
the nuclear gradient.
