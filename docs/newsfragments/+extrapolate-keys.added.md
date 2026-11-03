`CPMDC_EXTRAP_KEYS=1` keeps one `EXTRAPOLATE WFN` history per orbital key of `cpmdc_session_select_orbitals()`.
A host that alternates two or more trajectories on one calculator then starts each self-consistent field (SCF) from the prediction of that trajectory's own key.
On two 12-step water trajectories under two keys in one process, the in-process route took 88 and 100 SCF iterations.
`cpmd.x` with the always stable predictor-corrector (ASPC) of order 4 took as many on each trajectory alone.
Without the variable, it took 125 and 143, since a change of key empties the history as before.
`tools/opencpmd_embed_extrapolate.patch` adds `embed_extrap_keys` to `embed_ctrl`, and a `libcpmdc` built from this tree needs an archive with the updated patch.
