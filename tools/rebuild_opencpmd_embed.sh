#!/usr/bin/env bash
set -euo pipefail

cpmd_root=${1:?usage: rebuild_opencpmd_embed.sh CPMD_ROOT}
make_program=${MAKE:-make}

if [[ ! -f "$cpmd_root/Makefile" ]]; then
  printf 'missing OpenCPMD Makefile: %s\n' "$cpmd_root/Makefile" >&2
  exit 2
fi
if [[ ! -d "$cpmd_root/obj" || ! -d "$cpmd_root/lib" ]]; then
  printf 'incomplete OpenCPMD build tree: %s\n' "$cpmd_root" >&2
  exit 2
fi

# One object per OpenCPMD source file that a patch in tools/ changes, in
# module dependency order:
#   error_handling.mod.o   opencpmd_stopgm_return.patch
#   embed_ctrl.mod.o       opencpmd_embed_geometry.patch (new module),
#                          opencpmd_embed_teardown.patch
#   rinitwf_driver.mod.o   opencpmd_embed_rinitwf.patch
#   initrun_driver.mod.o   opencpmd_embed_geometry.patch
#   rwfopt_utils.mod.o     opencpmd_embed_rwfopt.patch
#   updwf_utils.mod.o      opencpmd_converged_state.patch
#   rkpnt_utils.mod.o      opencpmd_kpoints_inputfile.patch
#   timer.mod.o            opencpmd_tistopgm.patch
#   c_mem_utils.o          opencpmd_c_mem_addrs.patch
#   broyden_utils.mod.o .. opencpmd_embed_teardown.patch, the 32 modules
#   vpsi_utils.mod.o       whose first-call state the teardown resets
#   embed_teardown.mod.o   opencpmd_embed_teardown.patch (new module)
# scex_utils.mod.o carries no patch. gfortran stores a copy of Scex_t in
# rwfopt_utils.mod, so scex_utils.mod is written before rwfopt_utils reads
# it. A parallel rebuild leaves the two module files naming different
# vtable components, and the next USE rwfopt_utils stops.
# rwfopt_utils reads embed_have_orbitals from rinitwf_driver and
# embed_write_files from embed_ctrl; initrun_driver reads both modules.
# The teardown patch adds embed_teardowns to embed_ctrl, which its 32
# modules read, and embed_teardown reads rwfopt_utils, so it is last.
objects=(
  error_handling.mod.o
  scex_utils.mod.o
  embed_ctrl.mod.o
  rinitwf_driver.mod.o
  initrun_driver.mod.o
  rwfopt_utils.mod.o
  updwf_utils.mod.o
  rkpnt_utils.mod.o
  timer.mod.o
  c_mem_utils.o
  broyden_utils.mod.o
  calc_alm_utils.mod.o
  chksym_utils.mod.o
  detdof_utils.mod.o
  drhov_utils.mod.o
  ehpsi_utils.mod.o
  fftnew_utils.mod.o
  forcedr_driver.mod.o
  hfx_utils.mod.o
  initclust_utils.mod.o
  k_odiis_utils.mod.o
  k_pcgrad_utils.mod.o
  mixing_g_utils.mod.o
  mixing_r_utils.mod.o
  moverho_utils.mod.o
  nlccset_utils.mod.o
  numpw_utils.mod.o
  phfac_utils.mod.o
  pw_hfx.mod.o
  qvan2_utils.mod.o
  recpnew_utils.mod.o
  rhodiis_utils.mod.o
  rnlsmd_utils.mod.o
  rpiiint_utils.mod.o
  rwswap_utils.mod.o
  setbasis_utils.mod.o
  symtrz_utils.mod.o
  tauofr_utils.mod.o
  testex_utils.mod.o
  updrho_utils.mod.o
  vdw_utils.mod.o
  vpsi_utils.mod.o
  embed_teardown.mod.o
)
"$make_program" \
  -C "$cpmd_root/obj" \
  -f "$cpmd_root/Makefile" \
  -j1 \
  "${objects[@]}"
# The generated Makefile sets AR to the archiver with no operation
# letter, so its libcpmd.a rule does not replace members.
ar_program=${AR:-ar}
"$ar_program" r "$cpmd_root/lib/libcpmd.a" \
  "${objects[@]/#/$cpmd_root/obj/}"
ranlib "$cpmd_root/lib/libcpmd.a"
