Wire public **capnp-fortran** into the embed apply path: Fortran decodes
serialized `CPMDParams` (functional / cutOffRy / charge / multiplicity /
cpmdRoot, including system and dft section overrides) via
`cpmdc_embed_apply_params`. Exploded `cpmdc_embed_set_config` is unused on the
serialized-params path; C still renders the input deck and geometry merges
still use `set_deck`.
