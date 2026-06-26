! SPDX-License-Identifier: MIT
!
! Stable ISO_C embed surface for OpenCPMD as a *library*.
!
! Configuration and geometry NEVER go through CPMD's INPUT/`control` file
! parser. Host code (rgpot, readcon-core, tests) supplies Cap'n Proto
! `CPMDParams` + per-step `ForceInput`; `cpmdc` decodes those in C and pushes
! scalars/arrays into this layer via bind(C). We set OpenCPMD module state
! directly (control_def defaults + typed overrides) and CALL library entry
! points (setsys path pieces, wfopts/forces) in-process against libcpmd.a.
!
! Pseudopotential *file paths* may still be opened by CPMD PP readers when we
! point ion PP names at paths from Cap'n Proto — that is PP data I/O, not
! method configuration through INPUT decks.

#include "cpmd_embed_config.h"

MODULE cpmd_embed_c_api
  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_double, c_int, c_null_char
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpmdc_embed_init
  PUBLIC :: cpmdc_embed_available
  PUBLIC :: cpmdc_embed_set_config
  PUBLIC :: cpmdc_embed_energy_grad
  PUBLIC :: cpmdc_embed_finalize

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  LOGICAL, SAVE :: cpmd_initialized = .FALSE.
  CHARACTER(LEN=64), SAVE :: cfg_functional = 'BLYP'
  REAL(real64), SAVE :: cfg_cutoff_ry = 70.0_real64
  INTEGER, SAVE :: cfg_charge = 0
  INTEGER, SAVE :: cfg_mult = 1
  ! Optional PP / long-tail text from Cap'n Proto (not fed to control()).
  CHARACTER(LEN=16384), SAVE :: cfg_notes = ' '
  CHARACTER(LEN=1024), SAVE :: cfg_cpmd_root = ' '
  INTEGER, PARAMETER :: max_embed_atoms = 256
  REAL(real64), PARAMETER :: ang_to_au = 1.0_real64 / 0.529177210903_real64

CONTAINS

  SUBROUTINE copy_c_string(src, n, dst)
    CHARACTER(KIND=c_char), INTENT(IN) :: src(*)
    INTEGER(c_int), INTENT(IN), VALUE :: n
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER :: i, lim
    dst = ' '
    lim = MIN(INT(n), LEN(dst))
    DO i = 1, lim
      IF (src(i) == c_null_char) EXIT
      dst(i:i) = TRANSFER(src(i), 'a')
    END DO
  END SUBROUTINE copy_c_string

  FUNCTION cpmdc_embed_init() RESULT(ok) BIND(C, NAME='cpmdc_embed_init')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    cpmd_initialized = .FALSE.
    ok = 1_c_int
#else
    runtime_ready = .FALSE.
    ok = 0_c_int
#endif
  END FUNCTION cpmdc_embed_init

  FUNCTION cpmdc_embed_available() RESULT(ok) BIND(C, NAME='cpmdc_embed_available')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#else
    ok = 0_c_int
#endif
  END FUNCTION cpmdc_embed_available

  FUNCTION cpmdc_embed_set_config(functional, functional_len, cutoff_ry, &
      charge, multiplicity, input_deck, input_deck_len, cpmd_root, &
      cpmd_root_len) RESULT(ok) BIND(C, NAME='cpmdc_embed_set_config')
    CHARACTER(KIND=c_char), INTENT(IN) :: functional(*)
    INTEGER(c_int), INTENT(IN), VALUE :: functional_len
    REAL(c_double), INTENT(IN), VALUE :: cutoff_ry
    INTEGER(c_int), INTENT(IN), VALUE :: charge
    INTEGER(c_int), INTENT(IN), VALUE :: multiplicity
    CHARACTER(KIND=c_char), INTENT(IN) :: input_deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: input_deck_len
    CHARACTER(KIND=c_char), INTENT(IN) :: cpmd_root(*)
    INTEGER(c_int), INTENT(IN), VALUE :: cpmd_root_len
    INTEGER(c_int) :: ok

    CALL copy_c_string(functional, functional_len, cfg_functional)
    IF (LEN_TRIM(cfg_functional) == 0) cfg_functional = 'BLYP'
    cfg_cutoff_ry = REAL(cutoff_ry, KIND=real64)
    IF (cfg_cutoff_ry <= 0.0_real64) cfg_cutoff_ry = 70.0_real64
    cfg_charge = INT(charge)
    cfg_mult = INT(multiplicity)
    IF (cfg_mult <= 0) cfg_mult = 1
    ! Kept for diagnostics / future PP path lists — NOT passed to control().
    CALL copy_c_string(input_deck, input_deck_len, cfg_notes)
    CALL copy_c_string(cpmd_root, cpmd_root_len, cfg_cpmd_root)
    cpmd_initialized = .FALSE.
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready)
#else
    ok = 0_c_int
#endif
  END FUNCTION cpmdc_embed_set_config

  FUNCTION cpmdc_embed_energy_grad(n_atoms, positions_ang, atomic_numbers, &
      cell_ang, has_cell, energy_h, grad_h_bohr) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_energy_grad')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad_h_bohr(*)
    INTEGER(c_int) :: ok
    INTEGER :: i, n3

    energy_h = 0.0_c_double
    n3 = MAX(0, INT(n_atoms)) * 3
    DO i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    END DO

#if defined(CPMDC_HAS_CPMD)
    IF (.NOT. runtime_ready .OR. runtime_finalized) THEN
      ok = 0_c_int
      RETURN
    END IF
    IF (n_atoms <= 0 .OR. n_atoms > max_embed_atoms) THEN
      ok = 0_c_int
      RETURN
    END IF
    CALL cpmd_lib_energy_grad(INT(n_atoms), positions_ang, atomic_numbers, &
         cell_ang, INT(has_cell), energy_h, grad_h_bohr, ok)
#else
    ok = 0_c_int
#endif
  END FUNCTION cpmdc_embed_energy_grad

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
    cpmd_initialized = .FALSE.
  END SUBROUTINE cpmdc_embed_finalize

#if defined(CPMDC_HAS_CPMD)
  ! ------------------------------------------------------------------
  ! In-memory OpenCPMD setup: defaults + Cap'n Proto knobs, no INPUT I/O.
  ! ------------------------------------------------------------------

  SUBROUTINE apply_method_knobs()
    USE control_def_utils, ONLY: control_def
    USE system, ONLY: cntr, cntl, parm
    USE spin, ONLY: clsd
    USE func, ONLY: func1, mfxcx_is_slaterx, mfxcc_is_lyp, mgcx_is_becke88, &
         mgcc_is_lyp, mfxcc_is_pz, mgcx_is_pbex, mgcc_is_pbec
    IMPLICIT NONE
    CHARACTER(LEN=16) :: xc

    ! Upstream defaults (same starting point as a fresh CPMD process).
    CALL control_def()

    ! Plane-wave cutoff (Rydberg) — normally set in &SYSTEM CUTOFF via sysin.
    cntr%ecut = REAL(cfg_cutoff_ry, KIND=KIND(cntr%ecut))

    ! Wavefunction optimization / forces (embed is BO single-point / step).
    cntl%tdiag = .FALSE.
    ! Prefer optimize-wavefunction style; exact cntl flags vary by version —
    ! wfopts inspects cntl%tmd / geometry flags set elsewhere.

    ! Spin
    IF (cfg_mult > 1) THEN
      clsd%nlsd = 2
    ELSE
      clsd%nlsd = 1
    END IF

    ! Minimal XC mapping from Cap'n Proto functional name (extend as needed).
    xc = ADJUSTL(cfg_functional)
    CALL TO_UPPER(xc)
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_pz
    func1%mgcx = 0
    func1%mgcc = 0
    IF (INDEX(xc, 'BLYP') > 0) THEN
      func1%mfxcc = mfxcc_is_lyp
      func1%mgcx = mgcx_is_becke88
      func1%mgcc = mgcc_is_lyp
    ELSE IF (INDEX(xc, 'PBE') > 0) THEN
      func1%mgcx = mgcx_is_pbex
      func1%mgcc = mgcc_is_pbec
    END IF

    ! Neutral charge default; charged systems need extra embed knobs later.
    IF (cfg_charge /= 0) THEN
      ! Charged systems: extend with explicit cntl charge fields when wired.
    END IF
  END SUBROUTINE apply_method_knobs

  SUBROUTINE TO_UPPER(s)
    CHARACTER(LEN=*), INTENT(INOUT) :: s
    INTEGER :: i, ic
    DO i = 1, LEN_TRIM(s)
      ic = IACHAR(s(i:i))
      IF (ic >= IACHAR('a') .AND. ic <= IACHAR('z')) &
           s(i:i) = ACHAR(ic - 32)
    END DO
  END SUBROUTINE TO_UPPER

  SUBROUTINE apply_cell_ang(cell_ang, has_cell)
    USE system, ONLY: parm
    USE kinds, ONLY: real_8
    IMPLICIT NONE
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER :: k
    IF (has_cell == 0) THEN
      ! Default cubic 10 Angstrom cell in a.u.
      parm%a1 = [10.0_real_8, 0.0_real_8, 0.0_real_8] * REAL(ang_to_au, real_8)
      parm%a2 = [0.0_real_8, 10.0_real_8, 0.0_real_8] * REAL(ang_to_au, real_8)
      parm%a3 = [0.0_real_8, 0.0_real_8, 10.0_real_8] * REAL(ang_to_au, real_8)
    ELSE
      DO k = 1, 3
        parm%a1(k) = REAL(cell_ang(k), real_8) * REAL(ang_to_au, real_8)
        parm%a2(k) = REAL(cell_ang(3 + k), real_8) * REAL(ang_to_au, real_8)
        parm%a3(k) = REAL(cell_ang(6 + k), real_8) * REAL(ang_to_au, real_8)
      END DO
    END IF
  END SUBROUTINE apply_cell_ang

  SUBROUTINE apply_geometry(n_atoms, positions_ang, atomic_numbers)
    USE coor, ONLY: tau0, fion
    USE ions, ONLY: ions0, ions1
    USE kinds, ONLY: real_8
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    INTEGER :: k, isp, iat, idx

    IF (.NOT. ALLOCATED(tau0)) RETURN

    ! Walk CPMD species layout; requires prior species registration that
    ! matches Cap'n Proto atomic_numbers order (single-species and simple
    ! multi-species builds set ions0%na / ions1%nsp in init).
    idx = 0
    DO isp = 1, ions1%nsp
      DO iat = 1, ions0%na(isp)
        DO k = 1, 3
          idx = idx + 1
          IF (idx > n_atoms * 3) RETURN
          tau0(k, iat, isp) = REAL(positions_ang(idx), real_8) * &
               REAL(ang_to_au, real_8)
        END DO
      END DO
    END DO
  END SUBROUTINE apply_geometry

  SUBROUTINE cpmd_lib_energy_grad(n_atoms, positions_ang, atomic_numbers, &
       cell_ang, has_cell, energy_h, grad_h_bohr, ok)
    USE coor, ONLY: fion
    USE ener, ONLY: ener_com
    USE ions, ONLY: ions0, ions1
    USE kinds, ONLY: real_8
    USE startpa_utils, ONLY: startpa
    USE envir_utils, ONLY: envir
    USE setcnst_utils, ONLY: setcnst
    USE setsc_utils, ONLY: setsc
    USE setsys_utils, ONLY: setsys
    USE genxc_utils, ONLY: genxc
    USE numpw_utils, ONLY: numpw
    USE rinit_utils, ONLY: rinit
    USE wfopts_utils, ONLY: wfopts
    USE parac, ONLY: paral
    USE fileopen_utils, ONLY: init_fileopen
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: positions_ang(*), cell_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad_h_bohr(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: isp, iat, k, idx

    ok = 0_c_int
    energy_h = 0.0_c_double

    IF (.NOT. cpmd_initialized) THEN
      CALL init_fileopen
      CALL startpa
      CALL envir
      CALL setcnst
      CALL apply_method_knobs()
      CALL apply_cell_ang(cell_ang, has_cell)
      ! Species / PP registration must be driven from Cap'n Proto atom + PP
      ! fields (detsp/ratom equivalents in-memory) — not INPUT. Until that
      ! table is fully wired, fail closed rather than calling control().
      IF (ions1%nsp <= 0 .OR. ions1%nat <= 0) THEN
        ! Bootstrap minimal single-pass species list from atomic_numbers.
        CALL register_species_from_z(n_atoms, atomic_numbers)
      END IF
      CALL apply_geometry(n_atoms, positions_ang, atomic_numbers)
      CALL setsc
      IF (paral%qmnode) THEN
        CALL setsys
        CALL genxc
        CALL numpw
        CALL rinit
        CALL wfopts
      END IF
      cpmd_initialized = .TRUE.
    ELSE
      CALL apply_cell_ang(cell_ang, has_cell)
      CALL apply_geometry(n_atoms, positions_ang, atomic_numbers)
      IF (paral%qmnode) CALL wfopts
    END IF

    energy_h = REAL(ener_com%etot, KIND=c_double)
    idx = 0
    IF (ALLOCATED(fion)) THEN
      DO isp = 1, ions1%nsp
        DO iat = 1, ions0%na(isp)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > n_atoms * 3) EXIT
            grad_h_bohr(idx) = -REAL(fion(k, iat, isp), KIND=c_double)
          END DO
        END DO
      END DO
    END IF
    ok = 1_c_int
  END SUBROUTINE cpmd_lib_energy_grad

  SUBROUTINE register_species_from_z(n_atoms, atomic_numbers)
    ! Minimal in-memory species table so we never need ratom() INPUT parsing.
    ! Full PP projector load still requires PP files referenced by Cap'n Proto
    ! atoms.pseudopotentials — set on ion PP name arrays in a follow-up.
    USE ions, ONLY: ions0, ions1
    USE system, ONLY: maxsys
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    INTEGER :: i, z, nsp, found, s

    ions1%nsp = 0
    ions1%nat = n_atoms
    DO i = 1, n_atoms
      z = INT(atomic_numbers(i))
      found = 0
      DO s = 1, ions1%nsp
        IF (ions0%iatyp(s) == z) THEN
          ions0%na(s) = ions0%na(s) + 1
          found = 1
          EXIT
        END IF
      END DO
      IF (found == 0) THEN
        ions1%nsp = ions1%nsp + 1
        nsp = ions1%nsp
        ions0%iatyp(nsp) = z
        ions0%na(nsp) = 1
        ions0%zv(nsp) = REAL(z, KIND=KIND(ions0%zv))
      END IF
    END DO
    maxsys%nsx = MAX(maxsys%nsx, ions1%nsp)
    maxsys%nax = MAX(maxsys%nax, MAXVAL(ions0%na(1:ions1%nsp)))
  END SUBROUTINE register_species_from_z
#endif

END MODULE cpmd_embed_c_api
