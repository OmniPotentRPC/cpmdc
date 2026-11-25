Direct cold-deck composition now rejects singular, non-finite, and
volume-underflow cells for isolated decks.
A box whose components are all NaN is rejected instead of falling back to the 12
Angstrom default cell.
A deck whose `&SYSTEM` names `CLUSTER`, `POISSON SOLVER`, `ISOLATED MOLECULE`,
or `MOLECULE ISOLATED` is isolated when the box replaces its cell, as it already
was for periodicity injection and the reference evaluator.
