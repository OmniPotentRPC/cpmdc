The `cpmdc_warm_start` and `cpmdc_scf_steps` lines of `CPMDC_SCF_STEPS=1` now reach a standard error that is redirected to a file at the end of each call.
The gfortran runtime buffers that unit when it is a regular file, so the lines used to appear only at exit, or not at all after an abort.
