The new explanation page `explanation/pbe-water-density.org` compares 32 PBE
waters at 350 K and 1 bar with the PBE density of Miceli et al. (2015): the
time per step of each OpenCPMD build from `bench/results/pbe_water_cost.org`,
a longer NPT run at a constant cutoff with its blocking curve and the plateau
test of Lee et al. from `bench/results/pbe_water_npt.org`, and the molecular
pressures of two Troullier-Martins pseudopotential pairs, and of the fixed
basis of Miceli et al., on rigid boxes and on states of that run from
`bench/results/pbe_water_pressure.org` and `pbe_water_states.org`.
`tests/python/test_constant_cutoff_npt_table.py` requires the page to carry
the generated tables.
