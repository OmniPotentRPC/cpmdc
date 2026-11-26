OpenCPMD stores Hirshfeld charges, q = Z_v - n, after a converged Gamma-point
SCF.
`Session.charges` and the ASE `charges` property read them in the caller's atom
order, in electrons.
K points, broken symmetry, a species without a weight function, and
`CPMDC_CHARGE_SNAPSHOT=0` leave them unset.
The reference evaluator leaves them unset too.
