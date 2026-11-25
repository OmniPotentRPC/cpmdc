`CPMDCCalculator.name` is now `"cpmdc"`, the name that `cpmdc.ase` registers
with `register_calculator_class`.
It was ASE's default, the lower-cased class name `"cpmdccalculator"`, which
`get_calculator_class` does not know, so a trajectory or an ASE database row
recorded a calculator name that ASE could not load.
Trajectories and rows written before this change keep the old name.
