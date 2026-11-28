The renderer names the extrapolation scheme explicitly, matching CPMD's own
default.
Empty options with `extrapolateWfnOrder` greater than 0 render
`EXTRAPOLATE WFN ASPC`.
Empty options at order 0, the schema default, render no `EXTRAPOLATE WFN` line.
Options that are only `STORE` render `EXTRAPOLATE WFN STORE ASPC`.
The added token does not change the coefficients.
