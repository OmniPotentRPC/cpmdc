A deck that sets `DIPOLE DYNAMICS` and `KPOINTS` is now rejected when the parameters are installed, with an error that names both keywords, instead of stopping OpenCPMD in `localize2` during the first evaluation.
The check reads the rendered deck, so it also covers parameters given as raw input blocks.
Use a Gamma-point deck to get a dipole.
