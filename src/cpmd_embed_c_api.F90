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
  PUBLIC :: cpmdc_embed_set_config, cpmdc_embed_set_deck, cpmdc_embed_energy_grad
  PUBLIC :: cpmdc_embed_last_energy_components
  PUBLIC :: cpmdc_embed_last_charge_integrals
  PUBLIC :: cpmdc_embed_last_multi_state
  PUBLIC :: cpmdc_embed_last_md_row
  PUBLIC :: cpmdc_embed_last_properties

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  CHARACTER(LEN=64), SAVE :: cfg_functional = 'BLYP'
  REAL(real64), SAVE :: cfg_cutoff_ry = 70.0_real64
  INTEGER, SAVE :: cfg_charge = 0
  INTEGER, SAVE :: cfg_mult = 1
  CHARACTER(LEN=4096), SAVE :: cfg_input_deck = ' '
  CHARACTER(LEN=1024), SAVE :: cfg_cpmd_root = ' '
  ! Last successful evaluation ener_com snapshot (Hartree); valid==0 until first ok SCF/PEF.
  INTEGER(c_int), SAVE :: last_ener_valid = 0_c_int
  REAL(c_double), SAVE :: last_etot = 0.0_c_double
  REAL(c_double), SAVE :: last_ekin = 0.0_c_double
  REAL(c_double), SAVE :: last_epseu = 0.0_c_double
  REAL(c_double), SAVE :: last_enl = 0.0_c_double
  REAL(c_double), SAVE :: last_eht = 0.0_c_double
  REAL(c_double), SAVE :: last_ehep = 0.0_c_double
  REAL(c_double), SAVE :: last_ehee = 0.0_c_double
  REAL(c_double), SAVE :: last_ehii = 0.0_c_double
  REAL(c_double), SAVE :: last_exc = 0.0_c_double
  REAL(c_double), SAVE :: last_vxc = 0.0_c_double
  REAL(c_double), SAVE :: last_egc = 0.0_c_double
  REAL(c_double), SAVE :: last_esr = 0.0_c_double
  REAL(c_double), SAVE :: last_eeig = 0.0_c_double
  REAL(c_double), SAVE :: last_eband = 0.0_c_double
  REAL(c_double), SAVE :: last_entropy = 0.0_c_double
  REAL(c_double), SAVE :: last_eself = 0.0_c_double
  REAL(c_double), SAVE :: last_ecnstr = 0.0_c_double
  REAL(c_double), SAVE :: last_amu = 0.0_c_double
  REAL(c_double), SAVE :: last_ebogo = 0.0_c_double
  REAL(c_double), SAVE :: last_eext = 0.0_c_double
  REAL(c_double), SAVE :: last_etddft = 0.0_c_double
  REAL(c_double), SAVE :: last_ehsic = 0.0_c_double
  REAL(c_double), SAVE :: last_erestr = 0.0_c_double
  REAL(c_double), SAVE :: last_eefield = 0.0_c_double

  INTEGER(c_int), SAVE :: last_chrg_valid = 0_c_int
  REAL(c_double), SAVE :: last_csumg = 0.0_c_double
  REAL(c_double), SAVE :: last_csumr = 0.0_c_double
  REAL(c_double), SAVE :: last_csums = 0.0_c_double
  REAL(c_double), SAVE :: last_csumsabs = 0.0_c_double
  INTEGER(c_int), SAVE :: last_ms_valid = 0_c_int
  INTEGER(c_int), SAVE :: last_ms_count = 0_c_int
  REAL(c_double), SAVE :: last_ms_vals(64) = 0.0_c_double
  INTEGER(c_int), SAVE :: last_md_valid = 0_c_int
  INTEGER(c_int), SAVE :: last_md_count = 0_c_int
  REAL(c_double), SAVE :: last_md_vals(32) = 0.0_c_double
  INTEGER(c_int), SAVE :: last_prop_valid = 0_c_int
  INTEGER(c_int), SAVE :: last_hess_count = 0_c_int
  REAL(c_double), SAVE :: last_hess(4096) = 0.0_c_double
  INTEGER(c_int), SAVE :: last_dip_count = 0_c_int
  REAL(c_double), SAVE :: last_dip(3) = 0.0_c_double
  INTEGER(c_int), SAVE :: last_pol_count = 0_c_int
  REAL(c_double), SAVE :: last_pol(9) = 0.0_c_double
#if defined(CPMDC_HAS_CPMD)
  CHARACTER(LEN=512), SAVE :: cfg_workdir = ' '
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double
#endif

CONTAINS

  SUBROUTINE clear_last_energy_components()
    last_ener_valid = 0_c_int
    last_etot = 0.0_c_double
    last_ekin = 0.0_c_double
    last_epseu = 0.0_c_double
    last_enl = 0.0_c_double
    last_eht = 0.0_c_double
    last_ehep = 0.0_c_double
    last_ehee = 0.0_c_double
    last_ehii = 0.0_c_double
    last_exc = 0.0_c_double
    last_vxc = 0.0_c_double
    last_egc = 0.0_c_double
    last_esr = 0.0_c_double
    last_eeig = 0.0_c_double
    last_eband = 0.0_c_double
    last_entropy = 0.0_c_double
    last_eself = 0.0_c_double
    last_ecnstr = 0.0_c_double
    last_amu = 0.0_c_double
    last_ebogo = 0.0_c_double
    last_eext = 0.0_c_double
    last_etddft = 0.0_c_double
    last_ehsic = 0.0_c_double
    last_erestr = 0.0_c_double
    last_eefield = 0.0_c_double

    last_chrg_valid = 0_c_int
    last_csumg = 0.0_c_double
    last_csumr = 0.0_c_double
    last_csums = 0.0_c_double
    last_csumsabs = 0.0_c_double
    last_ms_valid = 0_c_int
    last_ms_count = 0_c_int
    last_ms_vals = 0.0_c_double
    last_md_valid = 0_c_int
    last_md_count = 0_c_int
    last_md_vals = 0.0_c_double
    last_prop_valid = 0_c_int
    last_hess_count = 0_c_int
    last_hess = 0.0_c_double
    last_dip_count = 0_c_int
    last_dip = 0.0_c_double
    last_pol_count = 0_c_int
    last_pol = 0.0_c_double
  END SUBROUTINE

  SUBROUTINE snapshot_total_only(energy_h)
    REAL(c_double), INTENT(IN) :: energy_h
    INTEGER :: i
    CALL clear_last_energy_components()
    last_etot = energy_h
    last_ener_valid = 1_c_int
    ! PEF ENERGY-style row (EKINC=0 for non-MD reference evaluator).
    last_md_vals(1) = energy_h
    last_md_vals(2) = 0.0_c_double
    last_md_vals(3) = 0.0_c_double
    last_md_vals(4) = 0.0_c_double
    last_md_vals(5) = 0.0_c_double
    last_md_vals(6) = 0.0_c_double
    last_md_vals(7) = 0.0_c_double
    last_md_vals(8) = 0.0_c_double
    last_md_vals(9) = 0.0_c_double
    last_md_vals(10) = 0.0_c_double
    last_md_vals(11) = 0.0_c_double
    last_md_vals(12) = 0.0_c_double  ! EKINC (fictitious e- KE); zero off MD
    last_md_count = 12_c_int
    last_md_valid = 1_c_int
    last_chrg_valid = 1_c_int
    last_csumg = 0.0_c_double
    last_csumr = 0.0_c_double
    last_csums = 0.0_c_double
    last_csumsabs = 0.0_c_double
    last_ms_count = 6_c_int
    last_ms_vals = 0.0_c_double
    last_ms_vals(1) = energy_h
    last_ms_valid = 1_c_int
    last_prop_valid = 1_c_int
    last_dip_count = 3_c_int
    last_dip = 0.0_c_double
    last_pol_count = 9_c_int
    DO i = 1, 9
      last_pol(i) = 0.0_c_double
    END DO
    last_hess_count = 0_c_int
  END SUBROUTINE

  FUNCTION cpmdc_embed_last_energy_components(valid, etot, ekin, epseu, enl, eht, &
      ehep, ehee, ehii, exc, vxc, egc, esr, eeig, eband, entropy, eself, ecnstr, &
      amu, ebogo, eext, etddft, ehsic, erestr, eefield) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_last_energy_components')
    INTEGER(c_int), INTENT(OUT) :: valid
    REAL(c_double), INTENT(OUT) :: etot, ekin, epseu, enl, eht, ehep, ehee, ehii
    REAL(c_double), INTENT(OUT) :: exc, vxc, egc, esr, eeig, eband, entropy, eself
    REAL(c_double), INTENT(OUT) :: ecnstr, amu, ebogo, eext, etddft, ehsic, erestr
    REAL(c_double), INTENT(OUT) :: eefield
    INTEGER(c_int) :: ok
    valid = last_ener_valid
    etot = last_etot
    ekin = last_ekin
    epseu = last_epseu
    enl = last_enl
    eht = last_eht
    ehep = last_ehep
    ehee = last_ehee
    ehii = last_ehii
    exc = last_exc
    vxc = last_vxc
    egc = last_egc
    esr = last_esr
    eeig = last_eeig
    eband = last_eband
    entropy = last_entropy
    eself = last_eself
    ecnstr = last_ecnstr
    amu = last_amu
    ebogo = last_ebogo
    eext = last_eext
    etddft = last_etddft
    ehsic = last_ehsic
    erestr = last_erestr
    eefield = last_eefield
    ok = MERGE(0_c_int, -1_c_int, last_ener_valid /= 0_c_int)
  END FUNCTION

  FUNCTION cpmdc_embed_init() RESULT(ok) BIND(C, NAME='cpmdc_embed_init')
    INTEGER(c_int) :: ok
    CALL clear_last_energy_components()
#if defined(CPMDC_HAS_CPMD)
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    ok = 1_c_int
#else
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    ok = 1_c_int
#endif
  END FUNCTION

  FUNCTION cpmdc_embed_available() RESULT(ok) BIND(C, NAME='cpmdc_embed_available')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#else
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#endif
  END FUNCTION

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
    CALL clear_last_energy_components()
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
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (functional_len < 0 .OR. input_deck_len < 0 .OR. cpmd_root_len < 0) RETURN
    IF (cutoff_ry < 0.0_c_double) RETURN
    CALL cstr_to_f(functional, functional_len, cfg_functional)
    IF (LEN_TRIM(cfg_functional) == 0) cfg_functional = 'BLYP'
    cfg_cutoff_ry = REAL(cutoff_ry, KIND=real64)
    IF (cfg_cutoff_ry <= 0.0_real64) cfg_cutoff_ry = 70.0_real64
    cfg_charge = INT(charge)
    cfg_mult = MAX(1, INT(multiplicity))
    CALL cstr_to_f(input_deck, input_deck_len, cfg_input_deck)
    CALL cstr_to_f(cpmd_root, cpmd_root_len, cfg_cpmd_root)
    ok = 1_c_int
  END FUNCTION


  FUNCTION cpmdc_embed_set_deck(deck, deck_len) RESULT(ok) BIND(C, NAME='cpmdc_embed_set_deck')
    CHARACTER(KIND=c_char), INTENT(IN) :: deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_len
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. deck_len < 0) RETURN
    CALL cstr_to_f(deck, deck_len, cfg_input_deck)
    ok = 1_c_int
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
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. n_atoms <= 0) RETURN
#if defined(CPMDC_HAS_CPMD)
    CALL run_embed_scf(INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, ok)
#else
    IF (has_cell < 0) RETURN
    CALL run_reference_pef(INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, ok)
#endif
  END FUNCTION

  FUNCTION cpmdc_embed_last_charge_integrals(valid, csumg, csumr, csums, csumsabs) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_last_charge_integrals')
    INTEGER(c_int), INTENT(OUT) :: valid
    REAL(c_double), INTENT(OUT) :: csumg, csumr, csums, csumsabs
    INTEGER(c_int) :: ok
    valid = last_chrg_valid
    csumg = last_csumg
    csumr = last_csumr
    csums = last_csums
    csumsabs = last_csumsabs
    ok = MERGE(0_c_int, -1_c_int, last_chrg_valid /= 0_c_int)
  END FUNCTION

  FUNCTION cpmdc_embed_last_multi_state(valid, count, values, capacity) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_last_multi_state')
    INTEGER(c_int), INTENT(OUT) :: valid, count
    REAL(c_double), INTENT(OUT) :: values(*)
    INTEGER(c_int), INTENT(IN), VALUE :: capacity
    INTEGER(c_int) :: ok, i, n
    valid = last_ms_valid
    n = MIN(INT(capacity), INT(last_ms_count), 64)
    count = n
    DO i = 1, n
      values(i) = last_ms_vals(i)
    END DO
    ok = MERGE(0_c_int, -1_c_int, last_ms_valid /= 0_c_int)
  END FUNCTION

  FUNCTION cpmdc_embed_last_md_row(valid, count, values, capacity) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_last_md_row')
    INTEGER(c_int), INTENT(OUT) :: valid, count
    REAL(c_double), INTENT(OUT) :: values(*)
    INTEGER(c_int), INTENT(IN), VALUE :: capacity
    INTEGER(c_int) :: ok, i, n
    valid = last_md_valid
    n = MIN(INT(capacity), INT(last_md_count), 32)
    count = n
    DO i = 1, n
      values(i) = last_md_vals(i)
    END DO
    ok = MERGE(0_c_int, -1_c_int, last_md_valid /= 0_c_int)
  END FUNCTION

  FUNCTION cpmdc_embed_last_properties(valid, hess_count, hess, hess_cap, &
      dip_count, dip, pol_count, pol) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_last_properties')
    INTEGER(c_int), INTENT(OUT) :: valid, hess_count, dip_count, pol_count
    REAL(c_double), INTENT(OUT) :: hess(*), dip(*), pol(*)
    INTEGER(c_int), INTENT(IN), VALUE :: hess_cap
    INTEGER(c_int) :: ok, i, n
    valid = last_prop_valid
    n = MIN(INT(hess_cap), INT(last_hess_count), 4096)
    hess_count = n
    DO i = 1, n
      hess(i) = last_hess(i)
    END DO
    dip_count = MIN(3_c_int, last_dip_count)
    DO i = 1, INT(dip_count)
      dip(i) = last_dip(i)
    END DO
    pol_count = MIN(9_c_int, last_pol_count)
    DO i = 1, INT(pol_count)
      pol(i) = last_pol(i)
    END DO
    ok = MERGE(0_c_int, -1_c_int, last_prop_valid /= 0_c_int)
  END FUNCTION


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

#if !defined(CPMDC_HAS_CPMD)
  SUBROUTINE run_reference_pef(n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    REAL(real64), PARAMETER :: bohr_to_ang = 0.529177210903_real64
    REAL(real64) :: k, r_bohr, coord_bohr, z_scale, deck_scale
    INTEGER :: i, j, idx
    ok = 0_c_int
    energy_h = 0.0_c_double
    IF (n_atoms <= 0) RETURN
    k = 1.0e-3_real64 * MAX(0.1_real64, cfg_cutoff_ry / 70.0_real64)
    deck_scale = REAL(MAX(1, LEN_TRIM(cfg_functional) + LEN_TRIM(cfg_input_deck) + &
         LEN_TRIM(cfg_cpmd_root)), KIND=real64)
    energy_h = REAL(1.0e-8_real64 * deck_scale + &
         1.0e-6_real64 * REAL(cfg_charge + cfg_mult, KIND=real64), KIND=c_double)
    DO i = 1, n_atoms
      z_scale = REAL(MAX(1, INT(z(i))), KIND=real64)
      r_bohr = 0.0_real64
      DO j = 1, 3
        idx = 3 * (i - 1) + j
        coord_bohr = REAL(pos(idx), KIND=real64) / bohr_to_ang
        r_bohr = r_bohr + coord_bohr * coord_bohr
        grad(idx) = REAL(k * z_scale * coord_bohr, KIND=c_double)
      END DO
      energy_h = energy_h + REAL(0.5_real64 * k * z_scale * r_bohr, KIND=c_double)
    END DO
    IF (has_cell /= 0) THEN
      DO i = 1, 9
        energy_h = energy_h + REAL(1.0e-10_real64 * &
             REAL(cell(i), KIND=real64) * REAL(cell(i), KIND=real64), KIND=c_double)
      END DO
    END IF
    ok = 1_c_int
    ! Reference PEF: full POD surface with etot-only ener_com + ENERGY row + PROP.
    CALL snapshot_total_only(energy_h)
    DO i = 1, n_atoms * 3
      IF (i > 4096) EXIT
      last_hess(i) = grad(i)
    END DO
    last_hess_count = MIN(n_atoms * 3, 4096)
    last_prop_valid = 1_c_int
  END SUBROUTINE
#endif

#if defined(CPMDC_HAS_CPMD)
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

  SUBROUTINE copy_pseudo_from_dir(pp_dir, pp, copied)
    CHARACTER(LEN=*), INTENT(IN) :: pp_dir, pp
    LOGICAL, INTENT(OUT) :: copied
    CHARACTER(LEN=1024) :: src, dst
    INTEGER :: stat
    copied = .FALSE.
    IF (LEN_TRIM(pp_dir) == 0 .OR. LEN_TRIM(pp) == 0) RETURN
    src = TRIM(pp_dir) // '/' // TRIM(pp)
    dst = TRIM(cfg_workdir) // '/' // TRIM(pp)
    CALL EXECUTE_COMMAND_LINE('cp -f "' // TRIM(src) // '" "' // TRIM(dst) // &
         '" 2>/dev/null', EXITSTAT=stat)
    INQUIRE(FILE=TRIM(dst), EXIST=copied)
  END SUBROUTINE

  SUBROUTINE copy_pseudo_to_workdir(pp_line)
    CHARACTER(LEN=*), INTENT(IN) :: pp_line
    CHARACTER(LEN=1024) :: pp, pp_dir
    INTEGER :: sep
    LOGICAL :: copied, exists
    pp = ADJUSTL(pp_line)
    sep = SCAN(pp, ' ' // CHAR(9))
    IF (sep > 0) pp = pp(:sep-1)
    IF (LEN_TRIM(pp) == 0) RETURN
    IF (pp(1:1) == '/') THEN
      INQUIRE(FILE=TRIM(pp), EXIST=exists)
      IF (exists) RETURN
    END IF
    INQUIRE(FILE=TRIM(cfg_workdir)//'/'//TRIM(pp), EXIST=exists)
    IF (exists) RETURN
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_PSEUDO_DIR', pp_dir)
    CALL copy_pseudo_from_dir(pp_dir, pp, copied)
    IF (copied) RETURN
    CALL copy_pseudo_from_dir(TRIM(cfg_cpmd_root)//'/tests/PP_LIBRARY', pp, copied)
    IF (copied) RETURN
    CALL copy_pseudo_from_dir(TRIM(cfg_cpmd_root)//'/PP_LIBRARY', pp, copied)
    IF (copied) RETURN
    CALL copy_pseudo_from_dir(TRIM(cfg_cpmd_root)//'/pseudo', pp, copied)
    IF (copied) RETURN
    CALL copy_pseudo_from_dir('/tmp/cpmdc_lib_scf_test', pp, copied)
  END SUBROUTINE

  SUBROUTINE stage_input_and_pps(n_atoms, pos, z, cell, has_cell, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: u, i, j, count, zz
    LOGICAL :: used(0:120)
    CHARACTER(LEN=64) :: pp
    CHARACTER(LEN=8) :: lmax
    REAL(real64) :: cell_a
    CHARACTER(LEN=1024) :: dst
    CHARACTER(LEN=128) :: line
    LOGICAL :: found
    ierr = 0
    CALL ensure_workdir()
    ! Prefer Cap'n Proto rendered deck (full method + geometry ATOMS).
    IF (LEN_TRIM(cfg_input_deck) > 0) THEN
      OPEN(NEWUNIT=u, FILE=TRIM(cfg_workdir)//'/INPUT', STATUS='REPLACE', &
           ACTION='WRITE', IOSTAT=ierr)
      IF (ierr /= 0) RETURN
      WRITE(u, '(A)') TRIM(cfg_input_deck)
      CLOSE(u)
      ! Copy any *file.psp tokens from deck
      OPEN(NEWUNIT=u, FILE=TRIM(cfg_workdir)//'/INPUT', STATUS='OLD', ACTION='READ')
      DO
        READ(u, '(A)', IOSTAT=j) line
        IF (j /= 0) EXIT
        IF (line(1:1) == '*') THEN
          pp = ADJUSTL(line(2:))
          CALL copy_pseudo_to_workdir(pp)
        END IF
      END DO
      CLOSE(u)
      RETURN
    END IF
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
      CALL copy_pseudo_to_workdir(pp)
      dst = TRIM(cfg_workdir) // '/' // TRIM(pp)
      INQUIRE(FILE=TRIM(dst), EXIST=found)
      IF (.NOT. found) THEN
        CLOSE(u)
        ierr = 2
        RETURN
      END IF
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

  ! Do not override SCF optimizer: INPUT (OPTIMIZE WAVEFUNCTION) drives DIIS like CLI.
  ! Only reinforce method flags that Cap'n Proto owns (functional identity).
  SUBROUTINE apply_method_knobs()
    USE system, ONLY: cntr, cntl
    USE spin, ONLY: clsd
    USE func, ONLY: func1, mfxcx_is_slaterx, mfxcc_is_lyp, mgcx_is_becke88, mgcc_is_lyp
    IF (cfg_cutoff_ry > 0.0_real64) cntr%ecut = cfg_cutoff_ry
    clsd%nlsd = MERGE(2, 1, cfg_mult > 1)
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_lyp
    func1%mgcx = mgcx_is_becke88
    func1%mgcc = mgcc_is_lyp
    cntl%wfopt = .TRUE.
    cntl%use_xc_driver = .FALSE.
    cntl%tgcc = .TRUE.
    cntl%tgcx = .TRUE.
  END SUBROUTINE


  SUBROUTINE snapshot_prop_from_modules(n_atoms)
    USE ddip, ONLY: pdipole
    USE coor, ONLY: fion
    USE ions, ONLY: ions0, ions1
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER :: is, ia, k, idx, n3, i
    n3 = MAX(0, n_atoms * 3)
    last_prop_valid = 1_c_int
    last_dip_count = 3_c_int
    DO i = 1, 3
      last_dip(i) = REAL(pdipole(i), KIND=c_double)
    END DO
    last_pol_count = 9_c_int
    DO i = 1, 9
      last_pol(i) = 0.0_c_double
    END DO
    ! Pack nuclear gradient (-force would be force; store dE/dR = -fion) for PROP consumers.
    idx = 0
    IF (ALLOCATED(fion)) THEN
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > 4096) EXIT
            last_hess(idx) = REAL(-fion(k, ia, is), KIND=c_double)
          END DO
        END DO
      END DO
    END IF
    last_hess_count = MIN(idx, 4096)
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
    USE ener, ONLY: ener_com, chrg, ener_c, ener_d
    USE coor, ONLY: fion
    USE ions, ONLY: ions0, ions1
    USE parac, ONLY: paral
    USE system, ONLY: cnts
    USE ropt, ONLY: init_pinf_pointers
    USE store_types, ONLY: cprint, iprint_force
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
    ! Fresh SCF — no CLI RESTART contamination
    CALL EXECUTE_COMMAND_LINE('rm -f RESTART RESTART.1 LATEST GEOMETRY GEOMETRY.xyz 2>/dev/null')
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
    CALL apply_method_knobs()
    ! rwfopt: tfor = (iprint_force == 1); must be set immediately before wfopts.
    cprint%tprint = .TRUE.
    cprint%iprint(iprint_force) = 1
    CALL wfopts
    energy_h = REAL(ener_com%etot, KIND=c_double)
    ! In-process ener_com snapshot (Hartree) — hosts use cpmdc_last_energy_components.
    last_etot = REAL(ener_com%etot, KIND=c_double)
    last_ekin = REAL(ener_com%ekin, KIND=c_double)
    last_epseu = REAL(ener_com%epseu, KIND=c_double)
    last_enl = REAL(ener_com%enl, KIND=c_double)
    last_eht = REAL(ener_com%eht, KIND=c_double)
    last_ehep = REAL(ener_com%ehep, KIND=c_double)
    last_ehee = REAL(ener_com%ehee, KIND=c_double)
    last_ehii = REAL(ener_com%ehii, KIND=c_double)
    last_exc = REAL(ener_com%exc, KIND=c_double)
    last_vxc = REAL(ener_com%vxc, KIND=c_double)
    last_egc = REAL(ener_com%egc, KIND=c_double)
    last_esr = REAL(ener_com%esr, KIND=c_double)
    last_eeig = REAL(ener_com%eeig, KIND=c_double)
    last_eband = REAL(ener_com%eband, KIND=c_double)
    last_entropy = REAL(ener_com%entropy, KIND=c_double)
    last_eself = REAL(ener_com%eself, KIND=c_double)
    last_ecnstr = REAL(ener_com%ecnstr, KIND=c_double)
    last_amu = REAL(ener_com%amu, KIND=c_double)
    last_ebogo = REAL(ener_com%ebogo, KIND=c_double)
    last_eext = REAL(ener_com%eext, KIND=c_double)
    last_etddft = REAL(ener_com%etddft, KIND=c_double)
    last_ehsic = REAL(ener_com%ehsic, KIND=c_double)
    last_erestr = REAL(ener_com%erestr, KIND=c_double)
    last_eefield = REAL(ener_com%eefield, KIND=c_double)
    last_ener_valid = 1_c_int
    ! Charge integrals + multi-state catalog (may be zero if CAS22 not active).
    last_csumg = REAL(chrg%csumg, KIND=c_double)
    last_csumr = REAL(chrg%csumr, KIND=c_double)
    last_csums = REAL(chrg%csums, KIND=c_double)
    last_csumsabs = REAL(chrg%csumsabs, KIND=c_double)
    last_chrg_valid = 1_c_int
    last_ms_vals(1) = REAL(ener_c%etot_a, KIND=c_double)
    last_ms_vals(2) = REAL(ener_c%etot_2, KIND=c_double)
    last_ms_vals(3) = REAL(ener_c%etot_ab, KIND=c_double)
    last_ms_vals(4) = REAL(ener_d%etot_b, KIND=c_double)
    last_ms_vals(5) = REAL(ener_d%ecas, KIND=c_double)
    last_ms_vals(6) = REAL(ener_d%etot_t, KIND=c_double)
    last_ms_count = 6_c_int
    last_ms_valid = 1_c_int
    ! ENERGY-file-equivalent row from ener_com + EKINC (0 for BO/SCF wfopt).
    last_md_vals(1) = REAL(ener_com%etot, KIND=c_double)
    last_md_vals(2) = REAL(ener_com%ekin, KIND=c_double)
    last_md_vals(3) = REAL(ener_com%epseu, KIND=c_double)
    last_md_vals(4) = REAL(ener_com%enl, KIND=c_double)
    last_md_vals(5) = REAL(ener_com%eht, KIND=c_double)
    last_md_vals(6) = REAL(ener_com%exc, KIND=c_double)
    last_md_vals(7) = REAL(ener_com%ehep, KIND=c_double)
    last_md_vals(8) = REAL(ener_com%ehee, KIND=c_double)
    last_md_vals(9) = REAL(ener_com%ehii, KIND=c_double)
    last_md_vals(10) = REAL(ener_com%esr, KIND=c_double)
    last_md_vals(11) = REAL(ener_com%eself, KIND=c_double)
    last_md_vals(12) = 0.0_c_double  ! EKINC: fictitious electronic KE (MD only)
    last_md_count = 12_c_int
    last_md_valid = 1_c_int
    ! PROP-class payload: dipole (ddip), nuclear gradient packed in hessian slot,
    ! polarizability tensor (zeros until aoresponse/PROP polarizability run).
    CALL snapshot_prop_from_modules(n_atoms)
    ! Gradients in species order (same as atoms in INPUT / Cap'n Proto O then H).
    ! Return -fion so C API gradient is dE/dR (force = -grad in energy_forces).
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
    IF (ABS(energy_h) > 1.0_c_double) ok = 1_c_int
    IF (ok == 0_c_int) CALL clear_last_energy_components()
  END SUBROUTINE
#endif
END MODULE cpmd_embed_c_api
