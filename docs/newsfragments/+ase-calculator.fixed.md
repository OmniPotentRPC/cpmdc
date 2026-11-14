`CPMDCCalculator` follows ASE's `Calculator` contract.
`params` and `library` are `parameters` with `default_parameters`, and `todict()` returns them for a trajectory, leaving out a `params` buffer given as bytes.
`set()` rejects unknown keys, opens the new session before it drops the old one, clears the results, and returns the changed keys.
`copy.copy`, `copy.deepcopy`, and `pickle` work, and a copy opens its own session; so does a calculator after `close()`.
A fully periodic structure without three independent cell vectors raises `ValueError` instead of running as a cluster.
Stress on free boundaries or on an isolated deck raises `PropertyNotImplementedError` with the reason, and asking again does not evaluate the step again.
A new composition opens a new session instead of failing; a reordering of the same atoms keeps the session.
A change of `initial_charges` or `initial_magmoms` no longer evaluates again, because `CPMDParams` sets the charge and multiplicity.
`Session` encodes a periodic `ForceInput` in process instead of running `capnp encode` on every step, and a closed `Session` raises `CpmdcError`.
