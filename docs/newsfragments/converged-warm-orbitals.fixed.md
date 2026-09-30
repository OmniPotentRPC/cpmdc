A warm start keeps `c0` only after the SCF converges. An unconverged ODIIS
pass leaves the previous orbitals, and a PCG MINIMIZE continuation restarts
from that copy. A cutoff, cell, or state-count change drops the saved
orbitals before the next SCF.
