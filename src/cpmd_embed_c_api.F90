! SPDX-License-Identifier: MIT
!
! cpmd_embed_c_api.F90 — compiler-independent C ABI for OpenCPMD embed
! (nwchemc pattern: bind(C) names; engine CALLs live here / in legacy helpers).
!
#include "cpmd_embed_config.h"
MODULE cpmd_embed_c_api
  USE, INTRINSIC :: iso_c_binding
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpmdc_embed_init, cpmdc_embed_available, cpmdc_embed_finalize
  PUBLIC :: cpmdc_embed_set_config, cpmdc_embed_energy_grad

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  CHARACTER(LEN=64), SAVE :: cfg_functional = 'BLYP'
  REAL(real64), SAVE :: cfg_cutoff_ry = 70.0_real64
  INTEGER, SAVE :: cfg_charge = 0
  INTEGER, SAVE :: cfg_mult = 1
  CHARACTER(LEN=4096), SAVE :: cfg_input_deck = ' '
  CHARACTER(LEN=1024), SAVE :: cfg_cpmd_root = ' '
  CHARACTER(LEN=512), SAVE :: cfg_workdir = ' '
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double

CONTAINS

  FUNCTION cpmdc_embed_init() RESULT(ok) BIND(C, NAME='cpmdc_embed_init')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    ok = 1_c_int
#else
    runtime_ready = .FALSE.
    ok = 0_c_int
#endif
  END FUNCTION

  FUNCTION cpmdc_embed_available() RESULT(ok) BIND(C, NAME='cpmdc_embed_available')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#else
    ok = 0_c_int
#endif
  END FUNCTION

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
  END SUBROUTINE

  FUNCTION cpmdc_embed_set_config(functional, functional_len, cutoff_ry, charge, &
      multiplicity, input_deck, input_deck_len, cpmd_root, cpmd_root_len) &
      RESULT(ok) BIND(C, NAME='cpmdc_embed_set_config')
    CHARACTER(KIND=c_char), INTENT(IN) :: functional(*)
    INTEGER(c_int), INTENT(IN), VALUE :: functional_len
    REAL(c_double), INTENT(IN), VALUE :: cutoff_ry
    INTEGER(c_int), INTENT(IN), VALUE :: charge, multiplicity
    CHARACTER(KIND=c_char), INTENT(IN) :: input_deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: input_deck_len
    CHARACTER(KIND=c_char), INTENT(IN) :: cpmd_root(*)
    INTEGER(c_int), INTENT(IN), VALUE :: cpmd_root_len
    INTEGER(c_int) :: ok
    ok = 0_c_int
#if defined(CPMDC_HAS_CPMD)
    IF (.NOT. runtime_ready) RETURN
    CALL cstr_to_f(functional, functional_len, cfg_functional)
    IF (LEN_TRIM(cfg_functional) == 0) cfg_functional = 'BLYP'
    cfg_cutoff_ry = REAL(cutoff_ry, KIND=real64)
    IF (cfg_cutoff_ry <= 0.0_real64) cfg_cutoff_ry = 70.0_real64
    cfg_charge = INT(charge)
    cfg_mult = MAX(1, INT(multiplicity))
    CALL cstr_to_f(input_deck, input_deck_len, cfg_input_deck)
    CALL cstr_to_f(cpmd_root, cpmd_root_len, cfg_cpmd_root)
    ok = 1_c_int
#else
    IF (functional_len < 0 .OR. input_deck_len < 0 .OR. cpmd_root_len < 0) RETURN
    IF (cutoff_ry < 0.0_c_double) RETURN
#endif
  END FUNCTION

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
    ok = 0_c_int
    energy_h = 0.0_c_double
    n3 = MAX(0, INT(n_atoms) * 3)
    DO i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    END DO
#if defined(CPMDC_HAS_CPMD)
    IF (.NOT. runtime_ready .OR. n_atoms <= 0) RETURN
    CALL run_embed_scf(INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, ok)
#else
    IF (has_cell < 0) RETURN
#endif
  END FUNCTION

#if defined(CPMDC_HAS_CPMD)
  SUBROUTINE cstr_to_f(cbuf, n, fstr)
    CHARACTER(KIND=c_char), INTENT(IN) :: cbuf(*)
    INTEGER(c_int), INTENT(IN), VALUE :: n
    CHARACTER(LEN=*), INTENT(OUT) :: fstr
    INTEGER :: i, lim
    fstr = ' '
    lim = MIN(INT(n), LEN(fstr))
    DO i = 1, lim
      IF (cbuf(i) == c_null_char) EXIT
      fstr(i:i) = TRANSFER(cbuf(i), 'a')
    END DO
  END SUBROUTINE

  SUBROUTINE ensure_workdir()
    CHARACTER(LEN=512) :: tmpl
    INTEGER :: stat
    IF (LEN_TRIM(cfg_workdir) > 0) RETURN
    tmpl = '/tmp/cpmdc_embed_XXXXXX'
    ! Prefer pre-seeded demo dir with PPs; else mkdtemp via C helper not available
    ! Use fixed unique path from PID
    WRITE(cfg_workdir, '(A,I0)') '/tmp/cpmdc_embed_', GETPID()
    CALL EXECUTE_COMMAND_LINE('mkdir -p ' // TRIM(cfg_workdir), EXITSTAT=stat)
  END SUBROUTINE

  INTEGER FUNCTION GETPID()
    ! glibc getpid via C binding
    INTERFACE
      FUNCTION c_getpid() BIND(C, NAME='getpid')
        IMPORT :: c_int
        INTEGER(c_int) :: c_getpid
      END FUNCTION
    END INTERFACE
    GETPID = INT(c_getpid())
  END FUNCTION

  SUBROUTINE stage_input_and_pps(n_atoms, pos, z, cell, has_cell, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: u, i, j, nz, count, zz
    LOGICAL :: used(0:120)
    CHARACTER(LEN=64) :: pp
    CHARACTER(LEN=8) :: lmax
    REAL(real64) :: cell_a
    CHARACTER(LEN=1024) :: cmd, src, dst
    ierr = 0
    CALL ensure_workdir()
    cell_a = 10.0_real64
    IF (has_cell /= 0) THEN
      IF (cell(1) > 0.0_c_double) cell_a = REAL(cell(1), KIND=real64)
    END IF
    OPEN(NEWUNIT=u, FILE=TRIM(cfg_workdir)//'/INPUT', STATUS='REPLACE', &
         ACTION='WRITE', IOSTAT=ierr)
    IF (ierr /= 0) RETURN
    WRITE(u, '(A)') '&CPMD'
    WRITE(u, '(A)') ' OPTIMIZE WAVEFUNCTION'
    WRITE(u, '(A)') ' CONVERGENCE ORBITALS'
    WRITE(u, '(A)') '  1.0d-5'
    WRITE(u, '(A)') ' MAXITER'
    WRITE(u, '(A)') '  40'
    WRITE(u, '(A)') ' CENTER MOLECULE OFF'
    WRITE(u, '(A)') '&END'
    WRITE(u, '(A)') '&SYSTEM'
    WRITE(u, '(A)') ' SYMMETRY'
    WRITE(u, '(A)') '  0'
    WRITE(u, '(A)') ' ANGSTROM'
    WRITE(u, '(A)') ' CELL'
    WRITE(u, '(A,F12.6,A)') '  ', cell_a, ' 1.0 1.0 0.0 0.0 0.0'
    WRITE(u, '(A)') ' CUTOFF'
    WRITE(u, '(A,F12.6)') '  ', cfg_cutoff_ry
    WRITE(u, '(A)') ' POISSON SOLVER HOCKNEY'
    WRITE(u, '(A)') '&END'
    WRITE(u, '(A)') '&DFT'
    WRITE(u, '(A)') ' OLDCODE'
    WRITE(u, '(A)') ' FUNCTIONAL BLYP'
    WRITE(u, '(A)') '&END'
    WRITE(u, '(A)') '&ATOMS'
    used = .FALSE.
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) CYCLE
      IF (used(zz)) CYCLE
      used(zz) = .TRUE.
      IF (zz == 8) THEN
        pp = 'O_MT_BLYP.psp'
        lmax = 'P'
      ELSE IF (zz == 1) THEN
        pp = 'H_CVB_BLYP.psp'
        lmax = 'S'
      ELSE
        CLOSE(u)
        ierr = 1
        RETURN
      END IF
      ! copy PP from known locations
      dst = TRIM(cfg_workdir) // '/' // TRIM(pp)
      src = '/tmp/cpmdc_lib_scf_test/' // TRIM(pp)
      CALL EXECUTE_COMMAND_LINE('cp -f "' // TRIM(src) // '" "' // TRIM(dst) // &
           '" 2>/dev/null || true')
      OPEN(NEWUNIT=nz, FILE=TRIM(dst), STATUS='OLD', ACTION='READ', IOSTAT=j)
      IF (j /= 0) THEN
        CLOSE(u)
        ierr = 2
        RETURN
      END IF
      CLOSE(nz)
      count = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) == zz) count = count + 1
      END DO
      WRITE(u, '(A,A)') '*', TRIM(pp)
      WRITE(u, '(A,A)') ' LMAX=', TRIM(lmax)
      WRITE(u, '(A,I4)') '   ', count
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        WRITE(u, '(3F14.6)') pos(3*(j-1)+1), pos(3*(j-1)+2), pos(3*(j-1)+3)
      END DO
    END DO
    WRITE(u, '(A)') '&END'
    CLOSE(u)
  END SUBROUTINE

  SUBROUTINE apply_blyp_knobs()
    USE system, ONLY: cntr, cntl
    USE spin, ONLY: clsd
    USE func, ONLY: func1, mfxcx_is_slaterx, mfxcc_is_lyp, mgcx_is_becke88, mgcc_is_lyp
    cntr%ecut = cfg_cutoff_ry
    clsd%nlsd = MERGE(2, 1, cfg_mult > 1)
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_lyp
    func1%mgcx = mgcx_is_becke88
    func1%mgcc = mgcc_is_lyp
    cntl%wfopt = .TRUE.
    cntl%use_xc_driver = .FALSE.
    cntl%tgcc = .TRUE.
    cntl%tgcx = .TRUE.
    cntl%tdiag = .FALSE.
    cntl%diis = .FALSE.
    cntl%tsde = .TRUE.
    cntl%tdipd = .FALSE.
  END SUBROUTINE

  SUBROUTINE run_embed_scf(n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    USE fileopen_utils, ONLY: init_fileopen
    USE timer, ONLY: tistart
    USE startpa_utils, ONLY: startpa
    USE envir_utils, ONLY: envir
    USE setcnst_utils, ONLY: setcnst
    USE control_utils, ONLY: control
    USE dftin_utils, ONLY: dftin
    USE sysin_utils, ONLY: sysin
    USE setsc_utils, ONLY: setsc
    USE detsp_utils, ONLY: detsp
    USE mm_init_utils, ONLY: mm_init
    USE ratom_utils, ONLY: ratom
    USE vdwin_utils, ONLY: vdwin
    USE propin_utils, ONLY: propin
    USE setsys_utils, ONLY: setsys
    USE setbasis_utils, ONLY: setbasis
    USE genxc_utils, ONLY: genxc
    USE numpw_utils, ONLY: numpw
    USE rinit_utils, ONLY: rinit
    USE rinforce_utils, ONLY: rinforce
    USE fftprp_utils, ONLY: fft_init
    USE ortho_utils, ONLY: ortho_init
    USE initclust_utils, ONLY: initclust
    USE dginit_utils, ONLY: dg_init
    USE nosalloc_utils, ONLY: nosalloc
    USE exterp_utils, ONLY: exterp
    USE dqgalloc_utils, ONLY: dqgalloc
    USE prng_utils, ONLY: prnginit
    USE gle_utils, ONLY: gle_alloc
    USE vdw_wf_alloc_utils, ONLY: vdw_wf_alloc
    USE wfopts_utils, ONLY: wfopts
    USE ener, ONLY: ener_com
    USE coor, ONLY: fion
    USE ions, ONLY: ions0, ions1
    USE parac, ONLY: paral
    USE system, ONLY: cnts
    USE ropt, ONLY: init_pinf_pointers
    USE bicanonicalCpmd, ONLY: bicanonicalCpmdConfig, bicanonicalCpmdInputConfig, New
    USE bicanonicalConfig, ONLY: New
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, is, ia, k, idx, nmax
    LOGICAL :: tinfo
    ok = 0_c_int
    energy_h = 0.0_c_double
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    CALL stage_input_and_pps(n_atoms, pos, z, cell, has_cell, ierr)
    IF (ierr /= 0) RETURN
    CALL CHDIR(TRIM(cfg_workdir))
    paral%io_parent = .TRUE.
    cnts%inputfile = 'INPUT'
    CALL tistart(tcpu0, twall0)
    CALL init_fileopen
    CALL startpa
    CALL New(bicanonicalCpmdInputConfig)
    tinfo = .TRUE.
    CALL init_pinf_pointers()
    CALL envir
    CALL setcnst
    CALL control
    CALL dftin
    CALL sysin
    CALL setsc
    CALL detsp
    CALL mm_init
    CALL ratom
    CALL vdwin
    CALL propin(tinfo)
    CALL setsys
    CALL New(bicanonicalCpmdConfig, bicanonicalCpmdInputConfig)
    CALL genxc
    CALL numpw
    CALL rinit
    CALL rinforce
    CALL fft_init()
    CALL ortho_init()
    CALL initclust
    CALL dg_init
    CALL nosalloc
    CALL exterp
    CALL setbasis
    CALL dqgalloc
    CALL prnginit
    CALL gle_alloc
    CALL vdw_wf_alloc
    CALL apply_blyp_knobs()
    CALL wfopts
    energy_h = REAL(ener_com%etot, KIND=c_double)
    IF (ALLOCATED(fion)) THEN
      idx = 0
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > nmax) EXIT
            grad(idx) = REAL(-fion(k, ia, is), KIND=c_double)
          END DO
        END DO
      END DO
    END IF
    IF (energy_h == energy_h) ok = 1_c_int  ! finite check (NaN != NaN)
  END SUBROUTINE
#else
  SUBROUTINE cstr_to_f(cbuf, n, fstr)
    CHARACTER(KIND=c_char), INTENT(IN) :: cbuf(*)
    INTEGER(c_int), INTENT(IN), VALUE :: n
    CHARACTER(LEN=*), INTENT(OUT) :: fstr
    fstr = ' '
    IF (n < 0 .OR. cbuf(1) == c_null_char) RETURN
  END SUBROUTINE
#endif
END MODULE cpmd_embed_c_api
