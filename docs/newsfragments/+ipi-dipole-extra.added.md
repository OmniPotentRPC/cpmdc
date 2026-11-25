`cpmdc-ipi` and `cpmdc-ipi-c` send `{"dipole":[x,y,z]}` in electron Bohr as the
`FORCEREADY` extra string when the library recomputed a finite dipole for the
step, which i-PI reads as `extras["dipole"]`.
Other steps, including a NaN or infinite dipole, send an empty extra string.
