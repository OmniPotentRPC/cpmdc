`cpmdc-mdi` is an MDI engine for LAMMPS `fix mdi/qm`.
It takes Bohr coordinates and a cell, subtracts `>CELL_DISPL`, and returns the
cpmdc energy, forces, and stress in atomic units.
`KE_ELEC` is zero because that energy already includes the electronic kinetic
energy.
