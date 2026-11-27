A reference build can evaluate Gaussian MM charges in a periodic cell when
`CPMDC_PERIODIC_EMBEDDING` is `1`.
The added energy is the induction of the reference polarizable density in the
conducting-boundary Ewald field of those charges.
An OpenCPMD build still refuses that cell.
The reference stress snapshot stays invalid on these steps.
