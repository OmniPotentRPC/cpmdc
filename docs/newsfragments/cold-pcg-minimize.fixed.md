A library call with no stored orbitals runs PCG MINIMIZE, whichever
wavefunction optimiser the deck names, and restores that optimiser for the
next warm call. `cpmdc_last_wavefunction_optimiser` and
`cpmdc_session_last_wavefunction_optimiser` report the optimiser that ran.
