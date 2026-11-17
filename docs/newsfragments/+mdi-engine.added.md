`cpmdc-mdi` (also `python -m cpmdc.mdi`) is an MDI engine on `cpmdc_session_calculate_au()`, so LAMMPS `fix mdi/qm` and any other MDI driver can request energies, forces, and stress in atomic units from one cpmdc session.
`pymdi` is the optional extra `mdi`, not a dependency of `cpmdc`. The plugin entry is `MDI_Plugin_init_cpmdc`.
Drive cpmdc as an MDI engine (`docs/orgmode/howto/mdi.org`) covers the `@DEFAULT` commands, the units, a Python driver, and LAMMPS.
