`cpmdc-mdi` includes `mdi.h` from a CMake install of MDI that publishes
`mdi_INCLUDE_DIR` without putting the header on the imported target.
`-Dwith_mdi=enabled` stops when that header still cannot be included.
The TCP launch line passes `-hostname localhost`, which an MDI engine requires.
`fix mdi/qm` remaps coordinates into the box before `>COORDS`.
The engine unwraps each atom to the nearest image of its previous position, so a
step across a face keeps the energy and forces of the unwrapped path.
