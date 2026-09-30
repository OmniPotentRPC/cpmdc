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

# One object per vendored patch in tools/, in module dependency order:
#   error_handling.mod.o   opencpmd_stopgm_return.patch
#   embed_ctrl.mod.o       opencpmd_embed_geometry.patch (new module)
#   rinitwf_driver.mod.o   opencpmd_embed_rinitwf.patch
#   initrun_driver.mod.o   opencpmd_embed_geometry.patch
#   rwfopt_utils.mod.o     opencpmd_embed_rwfopt.patch
#   updwf_utils.mod.o      opencpmd_converged_state.patch
#   rkpnt_utils.mod.o      opencpmd_kpoints_inputfile.patch
#   timer.mod.o            opencpmd_tistopgm.patch
#   c_mem_utils.o          opencpmd_c_mem_addrs.patch
# scex_utils.mod.o carries no patch. gfortran stores a copy of Scex_t in
# rwfopt_utils.mod, so scex_utils.mod is written before rwfopt_utils reads
# it. A parallel rebuild leaves the two module files naming different
# vtable components, and the next USE rwfopt_utils stops.
# rwfopt_utils reads embed_have_orbitals from rinitwf_driver and
# embed_write_files from embed_ctrl; initrun_driver reads both modules.
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
