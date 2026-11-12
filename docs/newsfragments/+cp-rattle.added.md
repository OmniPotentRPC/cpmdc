`cpmdc_cp_step` takes an opt-in Car-Parrinello velocity-Verlet step.
`CPMDC_CP_ORTHO_RATTLE` orthonormalises the coefficients with iterative SHAKE and their velocities with RATTLE.
`CPMDC_CP_ORTHO_LOWDIN` diagonalises the overlap instead.
Session force calls are unchanged.
On the three-orbital water model, 2000 steps at a time step of 0.25 atomic units gave a conserved-energy drift of 5.410173338e-2 Ha with RATTLE and 4.314417233e-1 Ha with spectral Löwdin.
