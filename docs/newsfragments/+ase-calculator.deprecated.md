The `schema` argument of `Session`, `Potential`, and `CPMDCCalculator`, and the `restart` argument of `CPMDCCalculator`, are ignored and warn with `DeprecationWarning`.
cpmdc no longer needs the schema file to encode a cell and reads no restart file.
