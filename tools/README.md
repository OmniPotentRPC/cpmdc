# OpenCPMD patches required for live `libcpmdc` embed

Build OpenCPMD’s `libcpmd.a` as usual, then apply these patches to the same
source tree before recompiling the affected objects into `libcpmd.a` and
linking `libcpmdc`.

| Patch | Purpose |
| --- | --- |
| `opencpmd_keep_fion.patch` | Do not `DEALLOCATE(fion)` at end of `rwfopt` so the embed can copy nuclear forces |
| `opencpmd_warm_orbitals.patch` | Store and restore `c0` after initializing each SCF call, reset the store with each applied configuration, and force PEF/BOMD ionic-force evaluation |
| `opencpmd_converged_state.patch` | Keep the converged `c0` synchronized with the energy and forces computed by `forcedr`; DIIS/PCG/steepest-descent updates only run while the pre-update gradient is unconverged |
| PEF stress (no extra OpenCPMD patch) | Embed sets `cntl%tpres` before `wfopts`; snapshots `paiu/omega` (Ha/Bohr^3) into `cpmdc_last_stress` / `PotentialResult.stress` |

```bash
# from the OpenCPMD/CPMD tree used as -Dcpmd_root=
patch -p1 < /path/to/cpmdc/tools/opencpmd_keep_fion.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_warm_orbitals.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_converged_state.patch
/path/to/cpmdc/tools/rebuild_opencpmd_embed.sh /path/to/cpmd-root
# rebuild libcpmdc against the updated archive
```

The helper addresses the generated module objects and archive as explicit Make
targets. OpenCPMD build trees contain a directory named `lib`, so `make lib`
can consider that target satisfied without refreshing `lib/libcpmd.a`.

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
