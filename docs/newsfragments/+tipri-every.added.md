`CPMDC_TIPRI_EVERY=N` prints CPMD's timing table after every N-th OpenCPMD
evaluation, with the evaluations after the first under
`embed_scf`, and `bench/cpmd_timing.py` turns the last two tables of each output
into the calls and seconds of each routine per evaluation;
`tests/test_embed_tipri_every.py` guards where the timer and the print sit.
Rank 0 reads the variable and broadcasts it before the collective print.
A value longer than 16 characters fails the call before the evaluation, stores
no energy, and leaves the warm basis as it was.
