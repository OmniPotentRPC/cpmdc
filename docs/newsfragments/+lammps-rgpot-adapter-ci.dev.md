The `lammps` CI job now runs
`tests/python/test_lammps_rgpot_adapter.py` in the same pytest call as `tests/python/test_lammps.py`,
so the session adapter that `cpmdc.lammps` hands to rgpot's LAMMPS frontend is
tested where `rgpot.lammps` is installed.
