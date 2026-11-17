`CPMDCPropertySnapshot.hessian` and `hessian_count` are deprecated. The array
contains the nuclear gradient dE/dR, not a Hessian. Existing binaries retain
the field; new callers should use `cpmdc_last_gradient()` or
`cpmdc_session_last_gradient()`.
