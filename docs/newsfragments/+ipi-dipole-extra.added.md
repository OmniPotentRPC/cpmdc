`cpmdc-ipi` and `cpmdc-ipi-c` send a recomputed molecular dipole in the i-PI
extra string.
The JSON object is `{"dipole":[x,y,z]}` in electron Bohr.
A step with no dipole still sends an empty extra string.
A non-finite dipole stops the driver before the step is counted.
