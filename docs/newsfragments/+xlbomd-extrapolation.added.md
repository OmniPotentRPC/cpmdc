`extrapolateWfnOptions` set to `XLBOMD` propagates an auxiliary vector with the
extended Lagrangian of Niklasson, Phys. Rev. Lett. 100, 123004 (2008).
Orders 3 through 9 use the dissipative coefficients of Niklasson et al., J.
Chem. Phys. 130, 214109 (2009).
A reference session records each converged vector.
The extrapolated guess drives the surrogate SCF density mix and does not yet
feed the orbital solve's initial guess.
Each session keeps its own samples, and a RESTART image can carry them.
An OpenCPMD build refuses the option at configure time.
ASPC and polynomial decks are unchanged.
