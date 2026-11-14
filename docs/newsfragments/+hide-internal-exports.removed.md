`libcpmdc` exports only the 44 functions of `include/cpmdc.h` and `include/cpmdc_features.h`.
The 56 internal `cpmdc_` glue symbols (`cpmdc_embed_*`, `cpmdc_params_*`, `cpmdc_stop_*`, and the rest) and the symbols of the static archives linked into the library (c-capnproto, capnp-fortran, minimage, OpenCPMD) are no longer in its dynamic symbol table.
No public header declared any of them, so no public function changes, and `CPMDC_ABI_VERSION` and the soversion stay as they are.
A program that looked one of them up with `dlsym` no longer finds it.
`scripts/abi_baseline.txt` now lists the 44 public names.
