GROMACS can run NPT with cpmdc as the QM engine.
`libcpmdc_gromacs` exports `cp2k_get_stress_tensor()`, the strain derivative of the energy it returned: cpmdc's stress with `CELL PERIODIC`, the symmetric part of the QM atoms' `sum x (x) f` with `CELL NONE`, plus the strain derivative of the mechanical-embedding Ewald sum, in cubic and triclinic boxes.
GROMACS' QM/MM module never adds a virial, so `tools/gromacs/patches/qmmm-virial-2024.patch` (GROMACS 2024.6 and 2025.5) and `qmmm-virial-2026.patch` (2026.4) make it add `-0.5 sigma V` in kJ/mol from the main rank; `scripts/ci/build_gromacs_cpmdc.sh` applies them.
The library now reads GROMACS 2026's box, sent one lattice vector per row, as well as the transposed box of 2024 and 2025.
`tests/python/test_gromacs.py` checks the stress against ASE's finite-difference stress, GROMACS' virial and pressure against the stress times the volume and against a strain difference of GROMACS' potential energy, and the virial across domain decompositions; `tests/python/test_gromacs_opencpmd.py` checks the OpenCPMD stress the same way, and `bench/gromacs_opencpmd_npt.py` runs C-rescale or Parrinello-Rahman NPT with an OpenCPMD QM water beside a classical control.
Run GROMACS on cpmdc forces (`docs/orgmode/howto/gromacs.org`) describes the settings.
The C ABI is unchanged.
