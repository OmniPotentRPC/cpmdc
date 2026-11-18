The header and the C ABI reference now say what `CPMDCPropertySnapshot` holds: `polarizability` is always zeros with a count of 9 (no build computes it), and `hessian` is the nuclear gradient in Hartree/Bohr, cut at 4096 values, so a system of more than 1365 atoms is truncated and `hessian_count` gives the number kept.
They claimed that polarizability was filled when PROP is available and did not mention the limit.
The `cpmdc_energy` comment no longer says that no gradient is allocated; the call evaluates the gradient into a temporary array.
