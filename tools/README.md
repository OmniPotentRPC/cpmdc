# OpenCPMD patches required for live `libcpmdc` embed

Build OpenCPMD’s `libcpmd.a` as usual, then apply these patches to the same
source tree before recompiling the affected objects into `libcpmd.a` and
linking `libcpmdc`.

| Patch | Purpose |
| --- | --- |
| `opencpmd_embed_geometry.patch` | Add `embed_ctrl` and `embed_write_files` (default `.TRUE.`). `embed_set_write_files(.FALSE.)` skips the `geofile` `GEOMETRY` write in `initrun` as well as the `zhwwf` and `geofile` writes in `rwfopt`. `wrgeof` still prints coordinates. Apply this before the rwfopt patch |
| `opencpmd_embed_rwfopt.patch` | Publish `embed_set_warm_orbitals`, `embed_set_need_forces`, `embed_reset_warm_orbitals`, and `embed_set_write_files`. Leave `fion` allocated only when `embed_need_forces` is true. An unconverged call leaves `fion` zero. Save `c0` when `embed_warm_orbitals` is true and the SCF converged. An unconverged pass keeps a partial copy only while no converged copy exists, and only a PCG MINIMIZE continuation restores it. `embed_set_warm_orbitals(.FALSE.)` frees the saved copy. A shape mismatch on restore calls `stopgm` |
| `opencpmd_converged_state.patch` | Keep the converged `c0` synchronized with the energy and forces computed by `forcedr`. DIIS/PCG/steepest-descent updates run only while the pre-update gradient is unconverged. Steepest descent with `iproj <= 1` raises `gemax` before that check |
| `opencpmd_kpoints_inputfile.patch` | Name the deck CPMD read (`cnts%inputfile`) in the k-point report instead of `argv[1]`, which in an embedding host is the host's own argument and stops CPMD with `STOP 12345` when it is missing or longer than 80 characters |
| `opencpmd_stopgm_return.patch` | Publish `cpmd_stopgm_hook`. While an embed call is armed, cpmdc installs a catch that records the stop code and returns 1, so `stopgm` returns instead of calling `my_stopall`. The wavefunction and forces of that call are undefined, the next call sets CPMD up again, and more than one rank aborts the others. A null hook leaves `cpmd.x` calling `my_stopall` |
| `opencpmd_tistopgm.patch` | `tistopgm` prints the call stack only when `trace_depth` lies inside `trace_names`. Before `tistart` the depth is `HUGE(0)`, and the unguarded walk faults instead of writing `LocalError` |
| `opencpmd_c_mem_addrs.patch` | `cGetMemAddrs` returns the pointer as `size_t`. `cuda_get_address` stores that integer. GCC 14 rejects the stock return of a `size_t` pointer |
| PEF stress (no extra OpenCPMD patch) | Embed sets `cntl%tpres` before `wfopts`; snapshots `paiu/omega` (Ha/Bohr^3) into the caller image and `PotentialResult.stress` |

```bash
# from the OpenCPMD/CPMD tree used as -Dcpmd_root=
patch -p1 < /path/to/cpmdc/tools/opencpmd_embed_geometry.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_embed_rwfopt.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_converged_state.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_kpoints_inputfile.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_stopgm_return.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_tistopgm.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_c_mem_addrs.patch
/path/to/cpmdc/tools/rebuild_opencpmd_embed.sh /path/to/cpmd-root
# rebuild libcpmdc against the updated archive
```

`header` calls `timetag`. OpenCPMD compiles `src/timetag.F90` to
`obj/timetag.o` outside `libcpmd.a`, and the embed link needs that object:

```bash
make -C /path/to/cpmd-root/obj -f /path/to/cpmd-root/Makefile timetag.o
```

The helper rebuilds the module objects with one Make job, then replaces
those members with `ar r`. OpenCPMD build trees contain a directory named
`lib`, and the generated Makefile may set `AR` to the archiver with no
operation letter, so `make lib` does not refresh `lib/libcpmd.a`.
`FFLAGS` and `CFLAGS` in that Makefile need `-fPIC`: `libcpmdc` is a
shared library, and a non-PIC member fails the link with `R_X86_64_PC32`.
`tools/LINUX-X86_64-GFORTRAN-MPI-PIC` is the configuration the archive
how-to copies into `configure/` before `configure.sh`. It adds `-fPIC`
and links `-lopenblas`.

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
