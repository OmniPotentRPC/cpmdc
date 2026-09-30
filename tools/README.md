# OpenCPMD patches required for live `libcpmdc` embed

Build OpenCPMD’s `libcpmd.a` as usual, then apply these patches to the same
source tree before recompiling the affected objects into `libcpmd.a` and
linking `libcpmdc`.

| Patch | Purpose |
| --- | --- |
| `opencpmd_embed_rwfopt.patch` | Publish `embed_set_warm_orbitals`, `embed_set_need_forces`, and `embed_reset_warm_orbitals`. Leave `fion` allocated only when `embed_need_forces` is true, and save `c0` only when `embed_warm_orbitals` is true and the SCF converged. A shape mismatch on restore calls `stopgm` |
| `opencpmd_converged_state.patch` | Keep the converged `c0` synchronized with the energy and forces computed by `forcedr`. DIIS/PCG/steepest-descent updates run only while the pre-update gradient is unconverged. Steepest descent with `iproj <= 1` raises `gemax` before that check |
| `opencpmd_kpoints_inputfile.patch` | Name the deck CPMD read (`cnts%inputfile`) in the k-point report instead of `argv[1]`, which in an embedding host is the host's own argument and stops CPMD with `STOP 12345` when it is missing or longer than 80 characters |
| `opencpmd_stopgm_return.patch` | Publish `cpmd_stopgm_hook`. While an embed call is armed, cpmdc installs a catch that records the stop code and returns 1, so `stopgm` returns instead of calling `my_stopall`. The wavefunction and forces of that call are undefined, the next call sets CPMD up again, and more than one rank aborts the others. A null hook leaves `cpmd.x` calling `my_stopall` |
| `opencpmd_tistopgm.patch` | `tistopgm` prints the call stack only when `trace_depth` lies inside `trace_names`. Before `tistart` the depth is `HUGE(0)`, and the unguarded walk faults instead of writing `LocalError` |
| PEF stress (no extra OpenCPMD patch) | Embed sets `cntl%tpres` before `wfopts`; snapshots `paiu/omega` (Ha/Bohr^3) into the caller image and `PotentialResult.stress` |

```bash
# from the OpenCPMD/CPMD tree used as -Dcpmd_root=
patch -p1 < /path/to/cpmdc/tools/opencpmd_embed_rwfopt.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_converged_state.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_kpoints_inputfile.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_stopgm_return.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_tistopgm.patch
/path/to/cpmdc/tools/rebuild_opencpmd_embed.sh /path/to/cpmd-root
# rebuild libcpmdc against the updated archive
```

The helper rebuilds the module objects with one Make job, then replaces
those members with `ar r`. OpenCPMD build trees contain a directory named
`lib`, and the generated Makefile may set `AR` to the archiver with no
operation letter, so `make lib` does not refresh `lib/libcpmd.a`.
`FFLAGS` and `CFLAGS` in that Makefile need `-fPIC`: `libcpmdc` is a
shared library, and a non-PIC member fails the link with `R_X86_64_PC32`.

Cold embed path also requires at runtime:

```bash
export CPMDC_PSEUDO_DIR=/path/to/pseudopotentials   # relative *PP basenames
# optional synonym used by OpenCPMD recpnew (trailing slash required):
export CPMD_PP_LIBRARY_PATH=$CPMDC_PSEUDO_DIR/
```

`libcpmdc` cold start calls `cpmdc_prepare_pp_cwd` (chdir into the library)
before OpenCPMD `ratom`/`recpnew`, then `cpmdc_restore_host_cwd` after the first
force. That is required when the host process has `argc>1` (Catch2 filters, eOn
CLI args): stock OpenCPMD `get_pplib` then treats `argv[2]` as the PP library
path and ignores `CPMD_PP_LIBRARY_PATH`, so relative `*PP` basenames only resolve
via the CWD fallback after the chdir.

Without the patches, multi-force may still return energies but nuclear force
buffers can be all zeros, warm re-entry can lose its orbitals, and the saved
orbitals can describe a DIIS update made after the reported result.
