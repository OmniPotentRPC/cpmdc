# OpenCPMD patches required for live `libcpmdc` embed

Build OpenCPMD’s `libcpmd.a` as usual, then apply these patches to the
**same** `src/rwfopt_utils.mod.F90` tree before recompiling that object into
`libcpmd.a` and linking `libcpmdc`.

| Patch | Purpose |
| --- | --- |
| `opencpmd_keep_fion.patch` | Do not `DEALLOCATE(fion)` at end of `rwfopt` so the embed can copy nuclear forces |
| `opencpmd_warm_orbitals.patch` | `cpmdc_set_warm_orbitals` (store/restore `c0`, skip `initrun` on warm multi-force), `cpmdc_set_need_forces` (OR into `tfor` so PEF/BOMD ionic forces are not zeroed after `forcedr`) |

```bash
# from the OpenCPMD/CPMD tree used as -Dcpmd_root=
patch -p1 < /path/to/cpmdc/tools/opencpmd_keep_fion.patch
patch -p1 < /path/to/cpmdc/tools/opencpmd_warm_orbitals.patch
# rebuild rwfopt object and update lib/libcpmd.a, then rebuild libcpmdc
```

Cold embed path also requires at runtime:

```bash
export CPMDC_PSEUDO_DIR=/path/to/pseudopotentials   # relative *PP basenames
# optional synonym used by OpenCPMD recpnew:
export CPMD_PP_LIBRARY_PATH=$CPMDC_PSEUDO_DIR/
```

Without the patches, multi-force may still return energies but nuclear force
buffers can be all zeros, and warm re-entry may fail to link or skip orbital reuse.
