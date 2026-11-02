`tools/opencpmd_embed_extrapolate.patch` honours `EXTRAPOLATE WFN` between in-process calls.
OpenCPMD turns the keyword off outside Born-Oppenheimer molecular dynamics.
A host-driven trajectory therefore started each self-consistent field (SCF) from the unchanged orbitals of the last call.
With the patch, `rwfopt` keeps a ring of the converged orbitals of the previous calls and starts the SCF from the `extrapwf` prediction of Car-Parrinello molecular dynamics (CPMD), with the always stable predictor-corrector (ASPC) or the polynomial coefficients.
The in-process route reached an SCF iteration count of 297 on a 40-step water trajectory with ASPC of order 4, against 292 for `cpmd.x` with the same deck and 441 for `cpmd.x` without extrapolation.
`CPMDC_SCF_STEPS=1` prints `cpmdc_wfn_extrapolated K` for a call that started from a prediction over `K` earlier calls.
Decks without the keyword run as before.
