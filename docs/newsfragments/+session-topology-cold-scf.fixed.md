A session refuses a change of atom count or elemental composition.
Destroy the session and create a new one: the next SCF is cold, and no orbital
store, density, or extrapolation history carries over.
A moved geometry with the same atoms can reuse stored forces, an orbital warm
start, ASPC or POLY extrapolation, and a RESTART wavefunction seed.
