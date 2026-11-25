`cpmdc.lammps.FixExternal` and `cpmdc.ase.CPMDCCalculator` record the
composition of a step before refusing a non-finite energy, forces, or virial, so
the next step with another composition opens a new session instead of reusing
the one the refused step ran in.
