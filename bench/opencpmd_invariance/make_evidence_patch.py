#!/usr/bin/env python3
"""Apply the evidence-only OpenCPMD instrumentation and write a unified diff.

Evidence only, not for upstream as is.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HARNESS = Path(__file__).resolve().parent
OMP_FILES = [
    "kin_energy_utils.mod.F90",
    "rhoofr_utils.mod.F90",
    "xcener_utils.mod.F90",
    "vofrhoa_utils.mod.F90",
    "vofrhob_utils.mod.F90",
    "vofrhoh_utils.mod.F90",
    "vofrhos_utils.mod.F90",
    "gcener_utils.mod.F90",
    "stress_utils.mod.F90",
    "vlocst_utils.mod.F90",
    "htrstr_utils.mod.F90",
    "potfor_utils.mod.F90",
    "nlforce_utils.mod.F90",
    "rnlfor_utils.mod.F90",
    "ppener_utils.mod.F90",
]


def read(path: Path) -> str:
    return path.read_text(encoding="latin1")


def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="latin1")


def ensure_include(text: str) -> str:
    if "evidence_flags.h" in text:
        return text
    return '#include "evidence_flags.h"\n' + text


def replace_mpi_sums(text: str) -> tuple[str, int]:
    """Replace real8 and complex8 MPI_SUM reductions. Leave integer and real4."""
    n = 0

    def sub(pattern: str, repl, src: str) -> str:
        nonlocal n
        src2, k = re.subn(pattern, repl, src, flags=re.IGNORECASE)
        n += k
        return src2

    # In-place forms first so MPI_IN_PLACE is not parsed as a buffer.
    text = sub(
        r"call\s+mpi_allreduce\(\s*MPI_IN_PLACE\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_PRECISION\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_allreduce_r8_inplace(\1,\2,\3,\4)",
        text,
    )
    text = sub(
        r"call\s+mpi_allreduce\(\s*MPI_IN_PLACE\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_COMPLEX\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_allreduce_c8_inplace(\1,\2,\3,\4)",
        text,
    )
    text = sub(
        r"call\s+mpi_allreduce\(\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_PRECISION\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_allreduce_r8(\1,\2,\3,\4,\5)",
        text,
    )
    text = sub(
        r"call\s+mpi_allreduce\(\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_COMPLEX\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_allreduce_c8(\1,\2,\3,\4,\5)",
        text,
    )
    text = sub(
        r"call\s+mpi_reduce\(\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_PRECISION\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_reduce_r8(\1,\2,\3,\4,\5,\6)",
        text,
    )
    text = sub(
        r"call\s+mpi_reduce\(\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*"
        r"MPI_DOUBLE_COMPLEX\s*,\s*&\s*MPI_SUM\s*,\s*([^,\n]+)\s*,\s*([^,\n]+)\s*,\s*([^)\n]+)\s*\)",
        r"call evidence_reduce_c8(\1,\2,\3,\4,\5,\6)",
        text,
    )
    text = fix_scalar_calls(text)
    return text, n


def fix_scalar_calls(text: str) -> str:
    """Bare scalars cannot sequence-associate with an assumed-size dummy."""
    text = re.sub(
        r"call evidence_allreduce_r8\(([^,()\n]+),([^,()\n]+),1,([^,()\n]+),([^,()\n]+)\)",
        r"call evidence_allreduce_r8_s(\1,\2,\3,\4)",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"call evidence_allreduce_c8\(([^,()\n]+),([^,()\n]+),1,([^,()\n]+),([^,()\n]+)\)",
        r"call evidence_allreduce_c8_s(\1,\2,\3,\4)",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"call evidence_reduce_r8\(([^,()\n]+),([^,()\n]+),1,([^,()\n]+),([^,()\n]+),([^,()\n]+)\)",
        r"call evidence_reduce_r8_s(\1,\2,\3,\4,\5)",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"call evidence_reduce_c8\(([^,()\n]+),([^,()\n]+),1,([^,()\n]+),([^,()\n]+),([^,()\n]+)\)",
        r"call evidence_reduce_c8_s(\1,\2,\3,\4,\5)",
        text,
        flags=re.IGNORECASE,
    )
    return text


def patch_mp_interface(text: str) -> str:
    needle = "  USE kinds,                           ONLY: int_1,&"
    insert = (
        "  USE evidence_reduce,                 ONLY: evidence_allreduce_c8,&\n"
        "                                             evidence_allreduce_c8_inplace,&\n"
        "                                             evidence_allreduce_r8,&\n"
        "                                             evidence_allreduce_r8_inplace,&\n"
        "                                             evidence_allreduce_r8_s,&\n"
        "                                             evidence_allreduce_c8_s,&\n"
        "                                             evidence_reduce_c8,&\n"
        "                                             evidence_reduce_c8_s,&\n"
        "                                             evidence_reduce_r8,&\n"
        "                                             evidence_reduce_r8_s\n"
    )
    if "USE evidence_reduce" not in text:
        if needle not in text:
            raise SystemExit("mp_interface USE kinds anchor missing")
        text = text.replace(needle, insert + needle, 1)
    return text


def patch_sources(text: str) -> str:
    old = "para_global.mod.F90 \\\n           mp_interface.mod.F90"
    new = "para_global.mod.F90 \\\n           evidence_reduce.mod.F90 mp_interface.mod.F90"
    if "evidence_reduce.mod.F90" not in text:
        if old not in text:
            raise SystemExit("SOURCES mp_interface anchor missing")
        text = text.replace(old, new, 1)
    # rwfopt USEs evidence_dump, so the module must be compiled first.
    text = text.replace("evidence_dump.mod.F90 ", "")
    text = text.replace("evidence_dump.mod.F90", "")
    old = "rwfopt_utils.mod.F90"
    new = "evidence_dump.mod.F90 rwfopt_utils.mod.F90"
    if "evidence_dump.mod.F90 rwfopt_utils.mod.F90" not in text:
        if old not in text:
            raise SystemExit("SOURCES rwfopt anchor missing")
        text = text.replace(old, new, 1)
    return text


def patch_debfor(text: str) -> str:
    if "evidence_write_results" not in text:
        text = text.replace(
            "  USE forcedr_driver,                  ONLY: forcedr\n",
            "  USE forcedr_driver,                  ONLY: forcedr\n"
            "  USE evidence_dump,                   ONLY: evidence_write_results\n"
            "  USE totstr_utils,                    ONLY: totstr\n",
            1,
        )
    old = (
        "    ELSE\n"
        "       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&\n"
        "            nstate,1,.TRUE.,.TRUE.)\n"
        "    ENDIF\n"
        "    CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)\n"
    )
    new = (
        "    ELSE\n"
        "       ropt_mod%calste=cntl%tpres\n"
        "       CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&\n"
        "            nstate,1,.TRUE.,.TRUE.)\n"
        "       IF (cntl%tpres) CALL totstr\n"
        "    ENDIF\n"
        "    CALL evidence_write_results(fion)\n"
        "    CALL dscal(3*maxsys%nax*maxsys%nsx,-1.0_real_8,fion(1,1,1),1)\n"
    )
    if old not in text:
        raise SystemExit("debfor forcedr anchor missing")
    text = text.replace(old, new, 1)
    old = (
        "    ELSE\n"
        "       CALL give_scr_updwf(lupdate,tag,nstate,.FALSE.)\n"
        "       CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.TRUE.)\n"
        "    ENDIF\n"
    )
    new = (
        "    ELSE\n"
        "       CALL give_scr_updwf(lupdate,tag,nstate,.FALSE.)\n"
        "       ev_calste=ropt_mod%calste\n"
        "       ropt_mod%calste=ropt_mod%calste.OR.cntl%tpres\n"
        "       CALL give_scr_forcedr(lforces,tag,nstate,.TRUE.,.TRUE.)\n"
        "       ropt_mod%calste=ev_calste\n"
        "    ENDIF\n"
    )
    if old not in text:
        raise SystemExit("debfor scratch anchor missing")
    text = text.replace(old, new, 1)
    old = (
        "    INTEGER                                  :: lcalc_alm, lcopot, lforces, &\n"
        "                                                linitrun, lortho, lrhoofr, &\n"
        "                                                lrnlsm, lsymmat, ltddft, &\n"
        "                                                lupdate, nstate\n"
    )
    new = (
        "    INTEGER                                  :: lcalc_alm, lcopot, lforces, &\n"
        "                                                linitrun, lortho, lrhoofr, &\n"
        "                                                lrnlsm, lsymmat, ltddft, &\n"
        "                                                lupdate, nstate\n"
        "    LOGICAL                                  :: ev_calste\n"
    )
    if "LOGICAL                                  :: ev_calste" not in text:
        if old not in text:
            raise SystemExit("debfor declaration anchor missing")
        text = text.replace(old, new, 1)
    return text


def patch_rwfopt(text: str) -> str:
    if "USE evidence_dump" not in text:
        needle = "  USE ener\n"
        insert = (
            "  USE ener\n"
            "  USE evidence_dump,                   ONLY: evidence_write_results\n"
        )
        if needle not in text:
            raise SystemExit("rwfopt USE ener anchor missing")
        text = text.replace(needle, insert, 1)
    old = "    CALL zeroing(psi)\n    CALL give_scr_rwfopt(lscr,tag)\n"
    new = (
        "    CALL zeroing(psi)\n"
        "    ropt_mod%calste=cntl%tpres\n"
        "    CALL give_scr_rwfopt(lscr,tag)\n"
    )
    if old not in text and "ropt_mod%calste=cntl%tpres\n    CALL give_scr_rwfopt" not in text:
        raise SystemExit("rwfopt scratch anchor missing")
    text = text.replace(old, new, 1)
    old = "    IF (ropt_mod%convwf.AND.(tfor.OR.cntl%tpres.OR.vdwl%vdwd)) THEN\n"
    new = "    IF (cntl%tpres.OR.(ropt_mod%convwf.AND.(tfor.OR.vdwl%vdwd))) THEN\n"
    if old not in text and new not in text:
        raise SystemExit("rwfopt force guard missing")
    text = text.replace(old, new, 1)
    old = (
        "       CALL prmem('    RWFOPT')\n"
        "       IF (.NOT.cntl%bsymm.OR.ropt_mod%convwf) CALL finalp(tau0,fion,tau0,eigv)\n"
    )
    new = (
        "       CALL prmem('    RWFOPT')\n"
        "       CALL evidence_write_results(fion)\n"
        "       IF (.NOT.cntl%bsymm.OR.ropt_mod%convwf) CALL finalp(tau0,fion,tau0,eigv)\n"
    )
    if "CALL evidence_write_results(fion)" not in text:
        if old not in text:
            raise SystemExit("rwfopt finalp anchor missing")
        text = text.replace(old, new, 1)
    return text


def patch_wrener(text: str) -> str:
    if "USE evidence_dump" not in text:
        needle = "                                             tenergy_ok\n"
        insert = (
            "                                             tenergy_ok\n"
            "  USE evidence_dump,                   ONLY: evidence_log_step\n"
        )
        if needle not in text:
            raise SystemExit("wrener USE ener anchor missing")
        text = text.replace(needle, insert, 1)
    old = (
        "       IF (paral%io_parent) THEN\n"
        "          WRITE(6,'(I4,1PE11.3,1PE12.3,0PF15.6,1PE13.3,0PF10.2) ')&\n"
        "               nfi,gemax,cnorm,etot,detot,tcpU\n"
        "       ENDIF\n"
    )
    new = (
        "       IF (paral%io_parent) THEN\n"
        "          WRITE(6,'(I4,1PE11.3,1PE12.3,0PF15.6,1PE13.3,0PF10.2) ')&\n"
        "               nfi,gemax,cnorm,etot,detot,tcpU\n"
        "          CALL evidence_log_step(tcpu)\n"
        "       ENDIF\n"
    )
    if old not in text and "CALL evidence_log_step(tcpu)" not in text:
        raise SystemExit("wrener TCPU anchor missing")
    text = text.replace(old, new, 1)
    return text


def patch_kin(text: str) -> str:
    text = ensure_include(text)
    text = text.replace(
        "  USE cppt,                            ONLY: hg\n",
        "  USE cppt,                            ONLY: hg,&\n"
        "                                             inyh\n"
        "  USE evidence_reduce,                 ONLY: evidence_gorder_sum\n"
        "  USE parac,                           ONLY: parai\n",
        1,
    )
    anchor = "    ener_com%ekin=xkin*parm%tpiba2\n"
    block = """#if defined(CPMDC_EVIDENCE_GORDER)
    IF (prcp_com%akin.LE.deltakin) THEN
       BLOCK
         INTEGER :: ev_i, ev_ig
         INTEGER(8), ALLOCATABLE :: ev_key(:)
         REAL(real_8), ALLOCATABLE :: ev_val(:)
         REAL(real_8) :: ev_sk, ev_tot
         ALLOCATE(ev_key(ncpw%ngw), ev_val(ncpw%ngw))
         DO ev_ig=1,ncpw%ngw
            ev_key(ev_ig)=INT(inyh(1,ev_ig),8)+ISHFT(INT(inyh(2,ev_ig),8),21)+&
                 ISHFT(INT(inyh(3,ev_ig),8),42)
            ev_sk=0.0_real_8
            DO ev_i=1,nstate
               IF (crge%f(ev_i,1).NE.0._real_8) THEN
                  ev_sk=ev_sk+crge%f(ev_i,1)*REAL(hg(ev_ig)*CONJG(c0(ev_ig,ev_i))*c0(ev_ig,ev_i))
               ENDIF
            ENDDO
            ev_val(ev_ig)=ev_sk
         ENDDO
         CALL evidence_gorder_sum(ncpw%ngw,ev_key,ev_val,parai%allgrp,ev_tot)
         IF (parai%mepos.EQ.0) THEN
            xkin=ev_tot
         ELSE
            xkin=0.0_real_8
         ENDIF
       END BLOCK
    ENDIF
#endif
"""
    if "CPMDC_EVIDENCE_GORDER" not in text:
        if anchor not in text:
            raise SystemExit("kin_energy anchor missing")
        text = text.replace(anchor, block + anchor, 1)
    return text


def patch_dotp(text: str) -> str:
    text = ensure_include(text)
    old = "       dotp=dotp+2.0_real_8*ddot(2*n-2,a(2),1,b(2),1)\n"
    new = (
        "#if defined(CPMDC_EVIDENCE_GORDER)\n"
        "       BLOCK\n"
        "         INTEGER :: ev_k\n"
        "         REAL(real_8) :: ev_d\n"
        "         ev_d=0.0_real_8\n"
        "         DO ev_k=2,n\n"
        "            ev_d=ev_d+REAL(a(ev_k),KIND=real_8)*REAL(b(ev_k),KIND=real_8)+&\n"
        "                 AIMAG(a(ev_k))*AIMAG(b(ev_k))\n"
        "         ENDDO\n"
        "         dotp=dotp+2.0_real_8*ev_d\n"
        "       END BLOCK\n"
        "#else\n"
        "       dotp=dotp+2.0_real_8*ddot(2*n-2,a(2),1,b(2),1)\n"
        "#endif\n"
    )
    if text.count(old) != 1:
        raise SystemExit(f"dotp_c ddot count {text.count(old)}")
    text = text.replace(old, new, 1)
    old = "     dotp=dotp+2.0_real_8*ddot(2*n-2,a(3),1,b(3),1)\n"
    new = (
        "#if defined(CPMDC_EVIDENCE_GORDER)\n"
        "     BLOCK\n"
        "       INTEGER :: ev_k\n"
        "       REAL(real_8) :: ev_d\n"
        "       ev_d=0.0_real_8\n"
        "       DO ev_k=3,2*n\n"
        "          ev_d=ev_d+a(ev_k)*b(ev_k)\n"
        "       ENDDO\n"
        "       dotp=dotp+ev_d+ev_d\n"
        "     END BLOCK\n"
        "#else\n"
        "     dotp=dotp+2.0_real_8*ddot(2*n-2,a(3),1,b(3),1)\n"
        "#endif\n"
    )
    if text.count(old) != 1:
        raise SystemExit(f"dotp ddot count {text.count(old)}")
    text = text.replace(old, new, 1)
    return text


def patch_iraymax(text: str) -> str:
    text = ensure_include(text)
    start = text.find("    thread_buff(:,:) = 0\n")
    end = text.find("    ! ==--------------------------------------------------------------==\n    RETURN\n  END SUBROUTINE iraymax")
    if start < 0 or end < 0:
        raise SystemExit("iraymax anchors missing")
    original = text[start:end]
    serial = (
        "    thread_buff(:,:) = 0\n"
        "#if defined(CPMDC_EVIDENCE_OMP_FIXED)\n"
        "    im=-2**30\n"
        "    i1=1\n"
        "    i2=1\n"
        "    DO j2=1,n2\n"
        "       DO j1=1,n1\n"
        "          IF (ir(j1,j2).GE.im) THEN\n"
        "             im=ir(j1,j2)\n"
        "             i1=j1\n"
        "             i2=j2\n"
        "          ENDIF\n"
        "       ENDDO\n"
        "    ENDDO\n"
        "#else\n"
        + original
        + "#endif\n"
    )
    if "CPMDC_EVIDENCE_OMP_FIXED" not in text:
        text = text[:start] + serial + text[end:]
    return text


def patch_mltfft(text: str) -> str:
    text = ensure_include(text)
    old = (
        "    IF (use_wisdom) THEN\n"
        "       fftw_flags=fftw_me\n"
        "    ELSE\n"
        "       fftw_flags=fftw_es\n"
        "    ENDIF\n"
    )
    new = (
        "#if defined(CPMDC_EVIDENCE_FFTW_ESTIMATE)\n"
        "    fftw_flags=fftw_es\n"
        "#else\n"
        "    IF (use_wisdom) THEN\n"
        "       fftw_flags=fftw_me\n"
        "    ELSE\n"
        "       fftw_flags=fftw_es\n"
        "    ENDIF\n"
        "#endif\n"
    )
    if "CPMDC_EVIDENCE_FFTW_ESTIMATE" not in text:
        if old not in text:
            raise SystemExit("mltfft wisdom anchor missing")
        text = text.replace(old, new, 1)
    return text


def join_continuations(text: str) -> str:
    out = []
    buf = ""
    for line in text.splitlines():
        if buf:
            buf = buf[:-1].rstrip() + " " + line.strip()
        else:
            buf = line
        if buf.rstrip().endswith("&"):
            continue
        out.append(buf)
        buf = ""
    if buf:
        out.append(buf)
    return "\n".join(out)


def declaration_index(src: Path) -> dict[str, tuple[str, int]]:
    """Map lower-case name to (type, rank). type is real, complex, or integer."""
    found: dict[str, tuple[str, int]] = {}
    for path in src.glob("*.F90"):
        joined = join_continuations(read(path))
        for line in joined.splitlines():
            m = re.match(
                r"\s*(REAL|COMPLEX|INTEGER|DOUBLE\s+PRECISION)\b[^:\n]*::(.*)$",
                line,
                flags=re.IGNORECASE,
            )
            if not m:
                continue
            kind = m.group(1).lower()
            if kind.startswith("double"):
                typ = "real"
            elif kind.startswith("complex"):
                typ = "complex"
            elif kind.startswith("integer"):
                typ = "integer"
            else:
                typ = "real"
            for ent in split_top(m.group(2), ","):
                ent = ent.strip()
                if not ent:
                    continue
                name = re.match(r"[A-Za-z_][A-Za-z0-9_]*", ent)
                if not name:
                    continue
                rank = 0
                if "(" in ent:
                    rank = ent.count(",") + 1
                    if ":" in ent and "," not in ent.split("(", 1)[1]:
                        rank = ent.split("(", 1)[1].count(":")
                        if rank == 0:
                            rank = 1
                found[name.group(0).lower()] = (typ, rank)
    return found


def split_top(text: str, sep: str) -> list[str]:
    parts = []
    buf = []
    depth = 0
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == sep and depth == 0:
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    parts.append("".join(buf))
    return parts


def matching_paren(text: str, open_at: int) -> int:
    depth = 0
    for i in range(open_at, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced paren")


def replace_ident(body: str, name: str, replacement: str) -> str:
    return re.sub(
        rf"(?i)(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])",
        replacement,
        body,
    )


def replace_array(body: str, name: str) -> str:
    out = []
    i = 0
    low = body.lower()
    key = name.lower()
    while True:
        j = low.find(key, i)
        if j < 0:
            out.append(body[i:])
            break
        before_ok = j == 0 or not (body[j - 1].isalnum() or body[j - 1] == "_")
        after = j + len(key)
        after_ok = after >= len(body) or not (body[after].isalnum() or body[after] == "_")
        if not (before_ok and after_ok):
            out.append(body[i : j + 1])
            i = j + 1
            continue
        k = after
        while k < len(body) and body[k] in " \t":
            k += 1
        if k >= len(body) or body[k] != "(":
            raise ValueError(f"bare array use of {name}")
        close = matching_paren(body, k)
        out.append(body[i:j])
        out.append(f"ev_{key}{body[k:close]},ev_c)")
        i = close + 1
        low = body.lower()
    return "".join(out)


def parse_do(line: str) -> tuple[str, str, str] | None:
    m = re.match(
        r"\s*DO\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$",
        line,
        flags=re.IGNORECASE,
    )
    if not m:
        return None
    var = m.group(1)
    rest = m.group(2).split("!", 1)[0].strip()
    if rest.lower().startswith("while"):
        return None
    parts = split_top(rest, ",")
    if len(parts) < 2 or len(parts) > 3:
        return None
    if len(parts) == 3 and parts[2].strip() not in {"1", "1_4", "1_8"}:
        return None
    return var, parts[0].strip(), parts[1].strip()


def find_enddo(lines: list[str], do_at: int) -> int | None:
    depth = 1
    for i in range(do_at + 1, len(lines)):
        raw = lines[i]
        code = raw.split("!", 1)[0]
        if re.match(r"\s*DO\b", code, flags=re.IGNORECASE):
            depth += 1
        if re.match(r"\s*END\s*DO\b", code, flags=re.IGNORECASE):
            depth -= 1
            if depth == 0:
                return i
    return None


def directive_text(lines: list[str], start: int) -> tuple[str, int] | None:
    if "!$omp" not in lines[start].lower():
        return None
    chunk = [lines[start]]
    i = start
    while chunk[-1].rstrip().endswith("&") or (
        i + 1 < len(lines) and lines[i + 1].lstrip().lower().startswith("!$omp")
        and "end" not in lines[i + 1].lower()
        and "do" not in lines[start].lower().split("parallel")[-1]
        and i == start
    ):
        # Always consume following !$omp continuations.
        if i + 1 < len(lines) and lines[i + 1].lstrip().lower().startswith("!$omp"):
            i += 1
            chunk.append(lines[i])
            continue
        break
    # The loop above may over-consume a following omp that is not a continuation
    # when the first line has no '&'. Restrict to lines that continue the directive.
    if not lines[start].rstrip().endswith("&"):
        return " ".join(chunk), start
    return " ".join(chunk), i


def collect_directive(lines: list[str], start: int) -> tuple[str, int]:
    parts = [lines[start]]
    i = start
    while parts[-1].rstrip().endswith("&"):
        i += 1
        if i >= len(lines):
            break
        parts.append(lines[i])
    return "\n".join(parts), i


def reduction_names(directive: str) -> list[str] | None:
    low = directive.lower()
    if "reduction" not in low:
        return None
    if re.search(r"reduction\s*\(\s*(max|\*|min|-)", low):
        return None
    names = []
    for m in re.finditer(r"reduction\s*\(\s*\+\s*:\s*([^)]*)\)", directive, flags=re.IGNORECASE):
        for part in m.group(1).split(","):
            name = part.strip()
            if name:
                names.append(name)
    if not names:
        return None
    return names


def private_names(directive: str) -> list[str]:
    names = []
    for m in re.finditer(r"private\s*\(\s*([^)]*)\)", directive, flags=re.IGNORECASE):
        for part in m.group(1).split(","):
            name = part.strip()
            if name:
                names.append(name)
    return names


def region_line(line: str) -> bool:
    stripped = line.strip()
    if stripped == "":
        return True
    low = stripped.lower()
    if low.startswith("!$omp"):
        return True
    if stripped.startswith("#"):
        return True
    if stripped.startswith("!") and not low.startswith("!$"):
        return True
    return False


def pp_balanced(lines: list[str]) -> bool:
    depth = 0
    for line in lines:
        stripped = line.strip().lower()
        if stripped.startswith("#if"):
            depth += 1
        elif stripped.startswith("#endif"):
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def rewrite_omp(text: str, decls: dict[str, tuple[str, int]], log: list[str], fname: str) -> str:
    raw = text.splitlines()
    spans = []
    i = 0
    while i < len(raw):
        parsed = parse_do(raw[i])
        if parsed is None:
            i += 1
            continue
        start = i
        k = i - 1
        depth = 0
        while k >= 0 and region_line(raw[k]):
            stripped = raw[k].strip().lower()
            if stripped.startswith("#endif"):
                depth += 1
            elif stripped.startswith("#if"):
                if depth == 0:
                    break
                depth -= 1
            elif stripped.startswith("#else") or stripped.startswith("#elif"):
                if depth == 0:
                    break
            start = k
            k -= 1
        if depth != 0:
            log.append(f"{fname}:{i+1} skip unbalanced region")
            i += 1
            continue
        if spans and start <= spans[-1][1]:
            i += 1
            continue
        region = raw[start:i]
        if not any("reduction" in ln.lower() for ln in region):
            i += 1
            continue
        if not pp_balanced(region):
            log.append(f"{fname}:{i+1} skip unbalanced preprocessor around reduction")
            i += 1
            continue
        directive = "\n".join(ln for ln in region if "!$omp" in ln.lower())
        names = reduction_names(directive)
        if not names:
            log.append(f"{fname}:{i+1} skip non-plus reduction")
            i += 1
            continue
        uniq = []
        for name in names:
            if name.lower() not in {u.lower() for u in uniq}:
                uniq.append(name)
        names = uniq
        end = find_enddo(raw, i)
        if end is None:
            log.append(f"{fname}:{i+1} skip missing enddo")
            i += 1
            continue
        body = "\n".join(raw[i + 1 : end])
        if "!$omp" in body.lower():
            log.append(f"{fname}:{i+1} skip nested omp")
            i += 1
            continue
        var, lo, hi = parsed
        specs = []
        skip = False
        for name in names:
            info = decls.get(name.lower())
            use_rank = None
            match = re.search(rf"(?i)(?<![A-Za-z0-9_]){name}\s*\(", body)
            if match:
                close = matching_paren(body, match.end() - 1)
                inside = body[match.end() : close]
                use_rank = len(split_top(inside, ","))
            if info is None:
                if use_rank:
                    typ = "complex" if "conjg" in body.lower() else "real"
                    info = (typ, use_rank)
                else:
                    info = ("real", 0)
            elif use_rank is not None:
                info = (info[0], use_rank)
            typ, rank = info
            for line in body.splitlines():
                if "conjg" in line.lower() and re.search(
                    rf"(?i)(?<![A-Za-z0-9_]){name}(?![A-Za-z0-9_])", line
                ):
                    typ = "complex"
                    break
            if typ == "integer":
                skip = True
                break
            specs.append((name, typ, rank))
        if skip:
            log.append(f"{fname}:{i+1} skip integer reduction")
            i += 1
            continue
        try:
            new_body = body
            for name, _typ, rank in specs:
                if rank > 0:
                    new_body = replace_array(new_body, name)
                else:
                    new_body = replace_ident(new_body, name, f"ev_{name.lower()}(ev_c)")
        except ValueError as exc:
            log.append(f"{fname}:{i+1} skip {exc}")
            i += 1
            continue
        priv = private_names(directive)
        seen = {p.lower() for p in priv}
        for extra in (var, "ev_lo", "ev_hi"):
            if extra.lower() not in seen:
                priv.append(extra)
                seen.add(extra.lower())
        spans.append((start, end, specs, priv, var, lo, hi, new_body))
        i = end + 1

    out: list[str] = []
    cursor = 0
    for start, end, specs, priv, var, lo, hi, new_body in spans:
        out.extend(raw[cursor:start])
        deduped = []
        seen_out = set()
        for name in priv:
            if name.lower() in seen_out:
                continue
            seen_out.add(name.lower())
            deduped.append(name)
        priv = deduped
        indent = re.match(r"\s*", raw[start]).group(0)
        decls_b = []
        zeros = []
        combines = []
        for name, typ, rank in specs:
            lname = name.lower()
            if typ == "complex":
                ftype = "COMPLEX(real_8)"
                zero = "(0.0_real_8,0.0_real_8)"
            else:
                ftype = "REAL(real_8)"
                zero = "0.0_real_8"
            if rank == 0:
                decls_b.append(f"{indent}  {ftype}, ALLOCATABLE :: ev_{lname}(:)")
                zeros.append(f"{indent}  ev_{lname}(:)={zero}")
                combines.append(f"{indent}    {name}={name}+ev_{lname}(ev_c)")
            else:
                dims = ", ".join(f"SIZE({name},{d})" for d in range(1, rank + 1))
                colons = ", ".join([":"] * rank)
                decls_b.append(f"{indent}  {ftype}, ALLOCATABLE :: ev_{lname}({colons},:)")
                zeros.append(f"{indent}  ALLOCATE(ev_{lname}({dims},ev_nc))")
                zeros.append(f"{indent}  ev_{lname}={zero}")
                idx = ", ".join([":"] * rank)
                combines.append(f"{indent}    {name}={name}+ev_{lname}({idx},ev_c)")
        scalar_alloc = [name.lower() for name, _t, rank in specs if rank == 0]
        block = ["#if defined(CPMDC_EVIDENCE_OMP_FIXED)", f"{indent}BLOCK"]
        block.append(f"{indent}  INTEGER :: ev_c, ev_nc, ev_lo, ev_hi")
        block.append(f"{indent}  INTEGER, PARAMETER :: ev_cs=64")
        block.extend(decls_b)
        block.append(f"{indent}  ev_nc=({hi})-({lo})+1")
        block.append(f"{indent}  IF (ev_nc.LT.0) ev_nc=0")
        block.append(f"{indent}  ev_nc=(ev_nc+ev_cs-1)/ev_cs")
        block.append(f"{indent}  IF (ev_nc.LT.1) ev_nc=1")
        if scalar_alloc:
            block.append(
                f"{indent}  ALLOCATE(" + ", ".join(f"ev_{n}(ev_nc)" for n in scalar_alloc) + ")"
            )
        block.extend(zeros)
        block.append(f"{indent}!$omp parallel do schedule(static) &")
        groups = []
        acc = []
        for name in priv:
            trial = ", ".join(acc + [name])
            if acc and len(trial) > 80:
                groups.append(acc)
                acc = [name]
            else:
                acc.append(name)
        if acc:
            groups.append(acc)
        for gi, group in enumerate(groups):
            tail = " &" if gi + 1 < len(groups) else ""
            block.append(f"{indent}!$omp  private({', '.join(group)}){tail}")
        block.append(f"{indent}  DO ev_c=1,ev_nc")
        block.append(f"{indent}    ev_lo=({lo})+(ev_c-1)*ev_cs")
        block.append(f"{indent}    ev_hi=MIN(({lo})+ev_c*ev_cs-1,({hi}))")
        block.append(f"{indent}    DO {var}=ev_lo,ev_hi")
        if new_body:
            block.extend(new_body.split("\n"))
        block.append(f"{indent}    ENDDO")
        block.append(f"{indent}  ENDDO")
        block.append(f"{indent}  DO ev_c=1,ev_nc")
        block.extend(combines)
        block.append(f"{indent}  ENDDO")
        block.append(f"{indent}END BLOCK")
        block.append("#else")
        out.extend(block)
        out.extend(raw[start : end + 1])
        out.append("#endif")
        cursor = end + 1
    out.extend(raw[cursor:])
    log.append(f"{fname} rewritten {len(spans)}")
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def refresh_omp(src: Path, pristine: Path, log_path: Path) -> None:
    shutil.copyfile(HARNESS / "src" / "evidence_dump.mod.F90", src / "evidence_dump.mod.F90")
    for name in OMP_FILES:
        shutil.copyfile(pristine / name, src / name)
    write(src / "kin_energy_utils.mod.F90", patch_kin(read(src / "kin_energy_utils.mod.F90")))
    decls = declaration_index(src)
    log: list[str] = []
    for name in OMP_FILES:
        path = src / name
        text = ensure_include(read(path))
        text = rewrite_omp(text, decls, log, name)
        write(path, text)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"omp log {log_path}")

def apply_tree(src: Path, log_path: Path) -> None:
    for name in ("evidence_reduce.mod.F90", "evidence_dump.mod.F90", "evidence_flags.h"):
        shutil.copyfile(HARNESS / "src" / name, src / name)
    sources = src / "SOURCES"
    write(sources, patch_sources(read(sources)))
    mp = src / "mp_interface.mod.F90"
    write(mp, patch_mp_interface(read(mp)))
    total = 0
    for name in ("sum.inc", "sum_root.inc", "sum_scalar_in_place.inc", "sum_in_place.inc"):
        path = src / name
        text, n = replace_mpi_sums(read(path))
        write(path, text)
        total += n
        print(f"{name}: replaced {n}")
    if total < 20:
        raise SystemExit(f"too few MPI replacements: {total}")
    write(src / "debfor_utils.mod.F90", patch_debfor(read(src / "debfor_utils.mod.F90")))
    write(src / "rwfopt_utils.mod.F90", patch_rwfopt(read(src / "rwfopt_utils.mod.F90")))
    write(src / "wrener_utils.mod.F90", patch_wrener(read(src / "wrener_utils.mod.F90")))
    write(src / "kin_energy_utils.mod.F90", patch_kin(read(src / "kin_energy_utils.mod.F90")))
    write(src / "dotp_utils.mod.F90", patch_dotp(read(src / "dotp_utils.mod.F90")))
    write(src / "loadpa_utils.mod.F90", patch_iraymax(read(src / "loadpa_utils.mod.F90")))
    write(src / "mltfft_utils.mod.F90", patch_mltfft(read(src / "mltfft_utils.mod.F90")))
    decls = declaration_index(src)
    log: list[str] = []
    for name in OMP_FILES:
        path = src / name
        text = ensure_include(read(path))
        text = rewrite_omp(text, decls, log, name)
        write(path, text)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(log) + "\n", encoding="utf-8")
    print(f"omp log {log_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", type=Path, required=True, help="OpenCPMD src directory to modify")
    parser.add_argument("--pristine", type=Path, help="Pristine src directory for the unified diff")
    parser.add_argument("--patch-out", type=Path, help="Where to write the unified diff")
    parser.add_argument("--refresh-omp", action="store_true")
    args = parser.parse_args()
    log_path = HARNESS / "patches" / "omp_rewrite.log"
    if args.refresh_omp:
        if not args.pristine:
            raise SystemExit("--refresh-omp needs --pristine")
        refresh_omp(args.src, args.pristine, log_path)
    else:
        apply_tree(args.src, log_path)
    if args.pristine and args.patch_out:
        args.patch_out.parent.mkdir(parents=True, exist_ok=True)
        diff = subprocess.run(
            ["diff", "-ruN", str(args.pristine), str(args.src)],
            check=False,
            capture_output=True,
            text=True,
        )
        header = (
            "Evidence only, not for upstream as is.\n"
            "OpenCPMD 062582b7cfd832d36f88f504cd08e4ead42eb404.\n"
            "The harness selects behavior with src/evidence_flags.h.\n"
            "\n"
        )
        body = diff.stdout.replace(str(args.pristine), "a/src").replace(str(args.src), "b/src")
        args.patch_out.write_text(header + body, encoding="latin1")
        print(f"patch bytes {args.patch_out.stat().st_size}")


if __name__ == "__main__":
    main()
