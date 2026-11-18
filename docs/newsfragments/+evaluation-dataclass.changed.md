`cpmdc.Evaluation` is a dataclass, so it has a `repr` that lists the energy, forces and stress, and the keyword form `Evaluation(energy_h=..., forces_h_bohr=..., stress_h_bohr3=...)`.
Instances still compare and hash by identity, as before.
