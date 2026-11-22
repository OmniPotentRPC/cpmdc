Two `DIPOLE DYNAMICS` evaluations in one process no longer fault in `opeigr_para` when the second uses a larger basis, such as 40 Ry followed by 50 Ry.
`localize2` kept the orbital maps and dipole matrices of the first basis; the OpenCPMD teardown patch now resets them like the other routines that keep first-call state.
Rebuild `libcpmd.a` with the updated patch; `tools/rebuild_opencpmd_embed.sh` recompiles `localize_utils`.
