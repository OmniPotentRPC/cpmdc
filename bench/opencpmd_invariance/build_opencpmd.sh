#!/bin/bash
# Evidence only, not for upstream as is.
# Build OpenCPMD 062582b with the harness patch and one flag set.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SHA="062582b7cfd832d36f88f504cd08e4ead42eb404"
SRC_BASE="${OPENCPMD_SRC:-/tmp/opencpmd-src}"
WORK="${OPENCPMD_WORK:-/tmp/opencpmd-work}"
DEST="${OPENCPMD_DEST:-/tmp/opencpmd-build}"
VARIANT="${1:-stock}"
BIN_DIR="${ROOT}/bin"

if [[ ! -d "${SRC_BASE}/src" ]]; then
  echo "missing OpenCPMD source at ${SRC_BASE}" >&2
  exit 1
fi
got="$(git -C "${SRC_BASE}" rev-parse HEAD)"
if [[ "${got}" != "${SHA}" ]]; then
  echo "OpenCPMD HEAD ${got} is not ${SHA}" >&2
  exit 1
fi

write_flags() {
  local flags="$1"
  cat > "${WORK}/src/evidence_flags.h" <<EOF
#ifndef CPMDC_EVIDENCE_FLAGS_H
#define CPMDC_EVIDENCE_FLAGS_H
/* Evidence only, not for upstream as is.
   Variant ${VARIANT}. Absent macros keep the stock operation. */
${flags}
#endif
EOF
}

case "${VARIANT}" in
  stock) MACROS="" ;;
  a_gather) MACROS="#define CPMDC_EVIDENCE_MP_GATHER" ;;
  a_exact) MACROS="#define CPMDC_EVIDENCE_MP_EXACT" ;;
  b_gorder) MACROS="#define CPMDC_EVIDENCE_GORDER" ;;
  d_omp) MACROS="#define CPMDC_EVIDENCE_OMP_FIXED" ;;
  ab) MACROS=$'#define CPMDC_EVIDENCE_MP_GATHER\n#define CPMDC_EVIDENCE_GORDER' ;;
  full) MACROS=$'#define CPMDC_EVIDENCE_MP_GATHER\n#define CPMDC_EVIDENCE_GORDER\n#define CPMDC_EVIDENCE_OMP_FIXED' ;;
  fftw) MACROS="#define CPMDC_EVIDENCE_FFTW_ESTIMATE" ;;
  *) echo "unknown variant ${VARIANT}" >&2; exit 1 ;;
esac

if [[ ! -f "${DEST}/Makefile" && "${VARIANT}" != "stock" && "${VARIANT}" != "fftw" ]]; then
  echo "build stock first so ${DEST} exists" >&2
  exit 1
fi

if [[ "${VARIANT}" == "stock" || "${VARIANT}" == "fftw" ]]; then
  if [[ "${OPENCPMD_REUSE_WORK:-}" != "1" ]]; then
    rm -rf "${WORK}"
    cp -a "${SRC_BASE}" "${WORK}"
    python3 "${ROOT}/make_evidence_patch.py" \
      --src "${WORK}/src" \
      --pristine "${SRC_BASE}/src" \
      --patch-out "${ROOT}/patches/evidence.patch"
  fi
  cfg="LINUX-X86_64-GFORTRAN-MPI-EVIDENCE"
  if [[ "${VARIANT}" == "fftw" ]]; then
    cfg="LINUX-X86_64-GFORTRAN-MPI-EVIDENCE-FFTW"
    DEST="${OPENCPMD_DEST:-/tmp/opencpmd-build-fftw}"
  fi
  cp "${ROOT}/configure/${cfg}" "${WORK}/configure/${cfg}"
  write_flags "${MACROS}"
  (cd "${WORK}" && ./configure.sh -omp -DEST="${DEST}" "${cfg}")
fi

write_flags "${MACROS}"
# The Makefile does not depend on the flag header.
while IFS= read -r f; do
  if grep -q -e 'evidence_flags.h' -e 'CPMDC_EVIDENCE_' "$f"; then
    touch "$f"
  fi
done < <(find "${WORK}/src" \( -name '*.F90' -o -name '*.inc' \))

make -C "${DEST}" -j "${CPMDC_BUILD_JOBS:-4}"
mkdir -p "${BIN_DIR}"
cp -a "${DEST}/bin/cpmd.x" "${BIN_DIR}/cpmd.x.${VARIANT}"
echo "built ${BIN_DIR}/cpmd.x.${VARIANT}"
