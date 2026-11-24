The GROMACS interface caps a covalent bond between a QM atom and an MM atom
with a link hydrogen, so a QM group no longer has to hold whole molecules.
The library reads each `&LINK` section of the generated input, or of a
`-qmi` file, places a hydrogen at the IMOMM position of CP2K with
`ALPHA_IMOMM` 1.38 unless given, adds it to the QM calculation, and returns
its force to the two atoms of the bond in the shares of the chain rule. The
virial, MPI runs, and mechanical embedding follow, and the MM atom of a cut
bond leaves the QM-MM Coulomb sum.

`cpmdc.gromacs.qmi_input()` and `write_qmi()` take `links` and `link_alpha`.
