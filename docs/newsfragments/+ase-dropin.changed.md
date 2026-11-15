Each `CPMDCCalculator` step runs with `calc.directory` as the working directory, so a relative `CPMDC_PSEUDO_DIR` or `CPMD_PP_LIBRARY_PATH` resolves there, and so does the `RESTART.1` that a deck with `RESTART WAVEFUNCTION` reads when the message sets neither `permanentDir` nor `scratchDir`.
The default directory, `.`, changes nothing.
`charge`, `multiplicity`, `xc`, or any other keyword that is not `params`, `library`, or one of ASE's raises a `TypeError` that says the `CPMDParams` message holds the method, and a calculator without `params` raises `TypeError` with "params is required".
