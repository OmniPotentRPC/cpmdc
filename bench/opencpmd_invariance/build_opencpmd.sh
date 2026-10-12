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
  exact) MACROS="" ;;
  *) echo "unknown variant ${VARIANT}" >&2; exit 1 ;;
esac

if [[ ! -f "${DEST}/Makefile" && "${VARIANT}" != "stock" && "${VARIANT}" != "fftw" && "${VARIANT}" != "exact" ]]; then
  echo "build stock first so ${DEST} exists" >&2
  exit 1
fi

if [[ "${VARIANT}" == "exact" ]]; then
  # Harness plus the exact series. Do not regenerate evidence.patch.
  WORK="${OPENCPMD_WORK:-/tmp/opencpmd-exact-work}"
  DEST="${OPENCPMD_DEST:-/tmp/opencpmd-exact-build}"
  rm -rf "${WORK}"
  cp -a "${SRC_BASE}" "${WORK}"
  side="$(mktemp -d)"
  cp -a "${SRC_BASE}/." "${side}/"
  for p in "${ROOT}/patches/exact/"*.patch; do
    patch -p1 --forward --batch -d "${side}" < "${p}"
  done
  patch -p1 --forward --batch -d "${WORK}" < "${ROOT}/patches/evidence.patch"
  while IFS= read -r rel; do
    rel="${rel#b/}"
    if [[ "${rel}" == "src/SOURCES" ]]; then
      continue
    fi
    cp -a "${side}/${rel}" "${WORK}/${rel}"
  done < <(grep -h '^+++ b/' "${ROOT}/patches/exact/"*.patch | awk '{print $2}' | sort -u)
  python3 - "${WORK}/src/SOURCES" <<'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1])
text = path.read_text()
needle = "mp_interface.mod.F90 "
insert = "mp_interface.mod.F90 exact_acc_utils.mod.F90 "
if "exact_acc_utils.mod.F90" not in text:
    if needle not in text:
        raise SystemExit("SOURCES mp_interface anchor missing")
    text = text.replace(needle, insert, 1)
    path.write_text(text)
PY
  cfg="LINUX-X86_64-GFORTRAN-MPI-EVIDENCE"
  cp "${ROOT}/configure/${cfg}" "${WORK}/configure/${cfg}"
  write_flags ""
  (cd "${WORK}" && ./configure.sh -omp -DEST="${DEST}" "${cfg}")
  # The switch is compile-time and default off. The harness macros stay unset.
  sed -i 's/-D__HAS_FFT_DEFAULT/-D__HAS_FFT_DEFAULT -D__EXACT_ACC/' "${DEST}/Makefile"
  rm -rf "${side}"
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
