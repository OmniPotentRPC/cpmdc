GROMACS-generated PBE and BLYP QM/MM inputs can load a matching
`CPMDParams` message from `CPMDC_GROMACS_PARAMS`.

The C ABI exposes `cpmdc_params_method()` and Python exposes
`params_method()` for reading a message's effective functional, charge, and
multiplicity without applying it.
