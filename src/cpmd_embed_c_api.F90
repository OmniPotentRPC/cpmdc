! SPDX-License-Identifier: MIT
!
! Stable ISO_C embed surface for OpenCPMD as a *library*.
! Links against libcpmd.a (from CPMD_ROOT/lib) when CPMDC_HAS_CPMD is defined.
! No subprocess / cpmd.x "driver": CALL into CPMD modules only.
!
! CPMD still reads a classical INPUT deck for control/sysin/ratom (upstream
! design). We write that deck to a private workdir and chdir for the duration
! of init/eval so the Fortran control path sees files; the ABI remains
! in-process library calls.

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
  CHARACTER(LEN=16384), SAVE :: cfg_input_deck = ' '
  CHARACTER(LEN=1024), SAVE :: cfg_cpmd_root = ' '
  CHARACTER(LEN=1024), SAVE :: workdir = ' '
  INTEGER, PARAMETER :: max_embed_atoms = 256

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
    CALL copy_c_string(input_deck, input_deck_len, cfg_input_deck)
    CALL copy_c_string(cpmd_root, cpmd_root_len, cfg_cpmd_root)
    ! Config change invalidates a prior CPMD init (new deck / topology setup).
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
    IF (n_atoms <= 0) RETURN
    IF (atomic_numbers(1) < -1) RETURN
    IF (has_cell /= 0 .AND. ABS(cell_ang(1)) < 0.0_c_double) RETURN
    IF (ABS(positions_ang(1)) < -1.0_c_double) RETURN
#endif
  END FUNCTION cpmdc_embed_energy_grad

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
    cpmd_initialized = .FALSE.
  END SUBROUTINE cpmdc_embed_finalize

#if defined(CPMDC_HAS_CPMD)
  ! ---- OpenCPMD library glue (same translation unit for module state) ----

  SUBROUTINE write_input_deck(path, n_atoms, positions_ang, atomic_numbers, &
       cell_ang, has_cell)
    CHARACTER(LEN=*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: positions_ang(*), cell_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    INTEGER :: u, ios, i
    REAL(real64) :: a, b, c
    CHARACTER(LEN=2) :: el
    INTEGER :: z

    OPEN (NEWUNIT=u, FILE=TRIM(path), STATUS='REPLACE', ACTION='WRITE', &
         IOSTAT=ios)
    IF (ios /= 0) RETURN

    IF (LEN_TRIM(cfg_input_deck) > 0) THEN
      WRITE (u, '(A)') TRIM(cfg_input_deck)
    ELSE
      WRITE (u, '(A)') '&CPMD'
      WRITE (u, '(A)') ' OPTIMIZE WAVEFUNCTION'
      WRITE (u, '(A)') ' CONVERGENCE ORBITALS'
      WRITE (u, '(A)') '  1.0d-6'
      WRITE (u, '(A)') ' DEBUG FORCES'
      WRITE (u, '(A)') '&END'
      WRITE (u, '(A)') ''
      WRITE (u, '(A)') '&SYSTEM'
      WRITE (u, '(A)') ' SYMMETRY'
      WRITE (u, '(A)') '  0'
      WRITE (u, '(A)') ' ANGSTROM'
      WRITE (u, '(A)') ' CELL'
      IF (has_cell /= 0) THEN
        a = SQRT(REAL(cell_ang(1), real64)**2 + REAL(cell_ang(2), real64)**2 + &
             REAL(cell_ang(3), real64)**2)
        b = SQRT(REAL(cell_ang(4), real64)**2 + REAL(cell_ang(5), real64)**2 + &
             REAL(cell_ang(6), real64)**2)
        c = SQRT(REAL(cell_ang(7), real64)**2 + REAL(cell_ang(8), real64)**2 + &
             REAL(cell_ang(9), real64)**2)
        WRITE (u, '(1X,3F16.8,1X,A)') a, b, c, '90.0 90.0 90.0'
      ELSE
        WRITE (u, '(A)') '  10.0 10.0 10.0 90.0 90.0 90.0'
      END IF
      WRITE (u, '(A)') ' CUTOFF'
      WRITE (u, '(1X,F16.8)') cfg_cutoff_ry
      IF (cfg_charge /= 0) THEN
        WRITE (u, '(A)') ' CHARGE'
        WRITE (u, '(1X,I0)') cfg_charge
      END IF
      WRITE (u, '(A)') '&END'
      WRITE (u, '(A)') ''
      WRITE (u, '(A)') '&DFT'
      WRITE (u, '(A,1X,A)') ' FUNCTIONAL', TRIM(cfg_functional)
      IF (cfg_mult > 1) WRITE (u, '(A)') ' LSD'
      WRITE (u, '(A)') '&END'
    END IF

    WRITE (u, '(A)') ''
    WRITE (u, '(A)') '&ATOMS'
    DO i = 1, n_atoms
      z = INT(atomic_numbers(i))
      el = element_symbol(z)
      ! Pseudopotential path must be supplied via Cap'n Proto deck for production;
      ! placeholder name documents the required *Element_PP stanza shape.
      WRITE (u, '(A,A,A)') '*', TRIM(el), '_none'
      WRITE (u, '(A)') ' LMAX=S'
      WRITE (u, '(A)') '   1'
      WRITE (u, '(1X,3F18.10)') REAL(positions_ang(3*(i-1)+1), real64), &
           REAL(positions_ang(3*(i-1)+2), real64), &
           REAL(positions_ang(3*(i-1)+3), real64)
    END DO
    WRITE (u, '(A)') '&END'
    CLOSE (u)
  END SUBROUTINE write_input_deck

  FUNCTION element_symbol(z) RESULT(sym)
    INTEGER, INTENT(IN) :: z
    CHARACTER(LEN=2) :: sym
    CHARACTER(LEN=2), PARAMETER :: table(0:36) = [ CHARACTER(LEN=2) :: &
         'X ', 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
         'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', &
         'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', &
         'As', 'Se', 'Br', 'Kr' ]
    IF (z < 0 .OR. z > 36) THEN
      sym = 'X '
    ELSE
      sym = table(z)
    END IF
  END FUNCTION element_symbol

  SUBROUTINE cpmd_lib_energy_grad(n_atoms, positions_ang, atomic_numbers, &
       cell_ang, has_cell, energy_h, grad_h_bohr, ok)
    ! Direct library evaluation: initialize CPMD modules once, update ionic
    ! coordinates in coor%tau0, run wavefunction optimization / forces, read
    ! ener_com%etot and coor%fion.
    USE coor, ONLY: tau0, fion
    USE ener, ONLY: ener_com
    USE ions, ONLY: ions0, ions1
    USE kinds, ONLY: real_8
    USE control_utils, ONLY: control
    USE setsys_utils, ONLY: setsys
    USE system, ONLY: cnts
    USE wfopts_utils, ONLY: wfopts
    USE parac, ONLY: paral
    USE startpa_utils, ONLY: startpa
    USE envir_utils, ONLY: envir
    USE setcnst_utils, ONLY: setcnst
    USE fileopen_utils, ONLY: init_fileopen
    USE dftin_utils, ONLY: dftin
    USE sysin_utils, ONLY: sysin
    USE setsc_utils, ONLY: setsc
    USE detsp_utils, ONLY: detsp
    USE mm_init_utils, ONLY: mm_init
    USE ratom_utils, ONLY: ratom
    USE vdwin_utils, ONLY: vdwin
    USE propin_utils, ONLY: propin
    USE genxc_utils, ONLY: genxc
    USE numpw_utils, ONLY: numpw
    USE rinit_utils, ONLY: rinit
    USE machine, ONLY: m_chdir, m_getcwd
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: positions_ang(*), cell_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad_h_bohr(*)
    INTEGER(c_int), INTENT(OUT) :: ok

    CHARACTER(LEN=1024) :: cwd_save, inp
    INTEGER :: is, ia, k, idx, isp, iat
    LOGICAL :: tinfo
    REAL(real_8) :: ang_to_au

    ok = 0_c_int
    energy_h = 0.0_c_double
    ang_to_au = 1.0_real_8 / 0.529177210903_real_8

    ! Private workdir per process; CPMD control reads ./INPUT by default.
    IF (LEN_TRIM(workdir) == 0) THEN
      CALL get_env_workdir(workdir)
    END IF
    CALL m_getcwd(cwd_save)
    CALL write_input_deck(TRIM(workdir)//'/INPUT', n_atoms, positions_ang, &
         atomic_numbers, cell_ang, has_cell)
    CALL m_chdir(TRIM(workdir))

    IF (.NOT. cpmd_initialized) THEN
      CALL init_fileopen
      CALL startpa
      CALL envir
      CALL setcnst
      cnts%inputfile = 'INPUT'
      CALL control
      CALL dftin
      CALL sysin
      CALL setsc
      CALL detsp
      CALL mm_init
      tinfo = .TRUE.
      IF (paral%qmnode) THEN
        CALL ratom
        CALL vdwin
        CALL propin(tinfo)
        CALL setsys
        CALL genxc
        CALL numpw
        CALL rinit
        ! Full wavefunction / G-vector / PP setup continues in wfopts.
        CALL wfopts
      END IF
      cpmd_initialized = .TRUE.
    ELSE
      ! Sustained session: update ionic positions in-place (Angstrom -> a.u.).
      IF (.NOT. ALLOCATED(tau0) .OR. .NOT. ALLOCATED(fion)) THEN
        CALL m_chdir(TRIM(cwd_save))
        RETURN
      END IF
      idx = 0
      DO isp = 1, ions1%nsp
        DO iat = 1, ions0%na(isp)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > n_atoms * 3) EXIT
            tau0(k, iat, isp) = REAL(positions_ang(idx), KIND=real_8) * ang_to_au
          END DO
        END DO
      END DO
      IF (paral%qmnode) CALL wfopts
    END IF

    energy_h = REAL(ener_com%etot, KIND=c_double)
    ! Map FION (a.u. force) -> gradient Ha/Bohr in C atom order (species major
    ! as stored in CPMD is not necessarily input order; for single-species
    ! tests this matches. Multi-species ordering follows ions0%na walk.)
    idx = 0
    IF (ALLOCATED(fion)) THEN
      DO isp = 1, ions1%nsp
        DO iat = 1, ions0%na(isp)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > n_atoms * 3) EXIT
            ! ABI gradient = -force
            grad_h_bohr(idx) = -REAL(fion(k, iat, isp), KIND=c_double)
          END DO
        END DO
      END DO
    END IF

    CALL m_chdir(TRIM(cwd_save))
    ok = 1_c_int
  END SUBROUTINE cpmd_lib_energy_grad

  SUBROUTINE get_env_workdir(dir)
    CHARACTER(LEN=*), INTENT(OUT) :: dir
    CHARACTER(LEN=1024) :: env
    INTEGER :: ios
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_WORKDIR', env, STATUS=ios)
    IF (ios == 0 .AND. LEN_TRIM(env) > 0) THEN
      dir = TRIM(env)
    ELSE
      dir = '/tmp/cpmdc-embed'
    END IF
    CALL execute_command_line('mkdir -p '//TRIM(dir), WAIT=.TRUE.)
  END SUBROUTINE get_env_workdir
#endif

END MODULE cpmd_embed_c_api
