`CPMDCCalculator.check_state` no longer lists `cell` while the stored structure
and the new one both have every `pbc` flag false.
A cluster step sends no cell to the library, so a different `atoms.cell`, for
example after `atoms.set_cell`, evaluated the same step again.
A change of `positions`, `numbers`, or `pbc`, or a cell change of a periodic
structure, still evaluates the step.
