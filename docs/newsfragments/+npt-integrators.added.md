`bench/constant_cutoff_npt.py run` takes `--thermostat svr`, stochastic
velocity rescaling of the atoms under the MTK barostat, and `--barostat pr`,
ASE's `MelchionnaNPT` with the cell shape held, a Parrinello-Rahman piston and
a Nose-Hoover thermostat; its state files hold the strain rate, so a run
continues from the state of another integrator.
`--oxygen` and `--hydrogen` name other pseudopotential files for `run` and
`pressure`, `--rs-wfn-keep`, `--aspc-order`, `--odiis-vectors` and
`--hamiltonian-cutoff` set the SCF, and `--keep-every` keeps states along a
run.
`pressure --density` scales the molecules of a configuration rigidly to another
density, and `table` gives the compressibility from the volume fluctuations,
the blocking curve every `--curve-step` levels, and the pseudopotentials and
integrators of the runs.
`cost` gives the wall time and SCF iterations per step of runs that differ
only in their library or SCF options, and the largest difference of their
potential energies.
`tests/python/test_constant_cutoff_npt_integrators.py` checks the integrators
and their states on EMT copper.
