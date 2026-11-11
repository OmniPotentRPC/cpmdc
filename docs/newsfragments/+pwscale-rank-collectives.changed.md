The opt-in scaffold's `cpmdc_pwscale_eval_mpi` runs 2 collectives per
call, down from 12 for the plane-wave split and 11 for the atom split.
`cpmdc_pwscale_mpi_bind_create` checks a problem once, and
`cpmdc_pwscale_eval_mpi_bound` then runs one `MPI_Allgather` per call.
Every rank sums the partials in rank order, so energy and forces keep
their bytes, match on every rank, and repeat from run to run.
`bench/run_pwscale.sh` takes `BIND=core`.
Its default, `none` with `OMP_PROC_BIND=close`, puts every rank on
CPU 0; on an 8-core host 2, 4, and 8 ranks took 0.0960, 0.180, and
0.416 s against 0.0747 s for one thread, and pinned ranks took 0.0380,
0.0211, and 0.0143 s.
At 16 atoms and 64 plane waves an 8-rank strong call went from 20.1 to
7.6 us; at the default size the change stayed within the spread of two
equal builds.
`bench/run_pwscale_abba.sh` compares two builds in ABBA order, and
`tests/test_pwscale_collectives.c` counts collectives through a PMPI
shim.
