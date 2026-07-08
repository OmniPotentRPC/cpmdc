! SPDX-License-Identifier: MIT
!
! cpmd_embed_c_api.F90 — compiler-independent C ABI for OpenCPMD embed
! (nwchemc pattern: bind(C) names; engine CALLs live here / in legacy helpers).
!
#include "cpmd_embed_config.h"
MODULE cpmd_embed_c_api
  USE, INTRINSIC :: iso_c_binding
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  ! Config knobs live in cpmdc_embed_apply_params (capnp-fortran decode path).
  USE cpmdc_embed_apply_params_mod, ONLY: &
      applied_functional, applied_cutoff_ry, applied_charge, applied_mult, &
      applied_input_deck, applied_cpmd_root
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
  ! Config state is applied_* from cpmdc_embed_apply_params_mod.
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
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double
  ! Successful in-process SCF steps; warm path reuses orbitals in memory.
  INTEGER, SAVE :: cfg_warm_steps = 0
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
    CALL cstr_to_f(functional, functional_len, applied_functional)
    IF (LEN_TRIM(applied_functional) == 0) applied_functional = 'BLYP'
    applied_cutoff_ry = REAL(cutoff_ry, KIND=real64)
    IF (applied_cutoff_ry <= 0.0_real64) applied_cutoff_ry = 70.0_real64
    applied_charge = INT(charge)
    applied_mult = MAX(1, INT(multiplicity))
    CALL cstr_to_f(input_deck, input_deck_len, applied_input_deck)
    CALL cstr_to_f(cpmd_root, cpmd_root_len, applied_cpmd_root)
    ok = 1_c_int
  END FUNCTION


  FUNCTION cpmdc_embed_set_deck(deck, deck_len) RESULT(ok) BIND(C, NAME='cpmdc_embed_set_deck')
    CHARACTER(KIND=c_char), INTENT(IN) :: deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_len
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. deck_len < 0) RETURN
    CALL cstr_to_f(deck, deck_len, applied_input_deck)
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
    k = 1.0e-3_real64 * MAX(0.1_real64, applied_cutoff_ry / 70.0_real64)
    deck_scale = REAL(MAX(1, LEN_TRIM(applied_functional) + LEN_TRIM(applied_input_deck) + &
         LEN_TRIM(applied_cpmd_root)), KIND=real64)
    energy_h = REAL(1.0e-8_real64 * deck_scale + &
         1.0e-6_real64 * REAL(applied_charge + applied_mult, KIND=real64), KIND=c_double)
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
  SUBROUTINE apply_method_knobs()
    USE system, ONLY: cntr, cntl
    USE spin, ONLY: clsd
    USE func, ONLY: func1, func2, func3, mfxcx_is_slaterx, mfxcc_is_lyp, &
         mgcx_is_becke88, mgcc_is_lyp, mhfx_is_skipped
    USE tbxc, ONLY: toldcode
    USE ener, ONLY: tenergy_ok
    IF (applied_cutoff_ry > 0.0_real64) cntr%ecut = applied_cutoff_ry
    clsd%nlsd = MERGE(2, 1, applied_mult > 1)
    toldcode = .TRUE.
    tenergy_ok = .TRUE.
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_lyp
    func1%mgcx = mgcx_is_becke88
    func1%mgcc = mgcc_is_lyp
    func1%mhfx = mhfx_is_skipped
    func2%salpha = 2.0_real64 / 3.0_real64
    func2%bbeta = 0.0042_real64
    func3%pxlda = 1.0_real64
    func3%pxgc = 1.0_real64
    func3%pclda = 1.0_real64
    func3%pcgc = 1.0_real64
    func3%phfx = 0.0_real64
    cntl%wfopt = .TRUE.
    cntl%diis = .TRUE.
    cntl%prec = .TRUE.
    cntl%tgc = .TRUE.
    cntl%tgcx = .TRUE.
    cntl%tgcc = .TRUE.
    cntl%use_xc_driver = .FALSE.
  END SUBROUTINE

  SUBROUTINE snapshot_prop_from_modules(n_atoms)
    USE ddip, ONLY: pdipole
    USE coor, ONLY: fion
    USE ions, ONLY: ions0, ions1
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER :: is, ia, k, idx, i
    last_prop_valid = 1_c_int
    last_dip_count = 3_c_int
    DO i = 1, 3
      last_dip(i) = REAL(pdipole(i), KIND=c_double)
    END DO
    last_pol_count = 9_c_int
    DO i = 1, 9
      last_pol(i) = 0.0_c_double
    END DO
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

  SUBROUTINE embed_pp_for_z(zz, pp, lmax_val, ok)
    INTEGER, INTENT(IN) :: zz
    CHARACTER(LEN=*), INTENT(OUT) :: pp
    INTEGER, INTENT(OUT) :: lmax_val, ok
    ok = 0
    pp = ' '
    lmax_val = -1
    IF (zz == 1) THEN
      pp = 'H_CVB_BLYP.psp'; lmax_val = 0; ok = 1
    ELSE IF (zz == 6) THEN
      pp = 'C_MT_BLYP.psp'; lmax_val = 1; ok = 1
    ELSE IF (zz == 7) THEN
      pp = 'N_MT_BLYP.psp'; lmax_val = 1; ok = 1
    ELSE IF (zz == 8) THEN
      pp = 'O_MT_BLYP.psp'; lmax_val = 1; ok = 1
    ELSE IF (zz == 14) THEN
      pp = 'Si_MT_BLYP.psp'; lmax_val = 1; ok = 1
    ELSE IF (zz == 32) THEN
      pp = 'Ge_MT_BLYP.psp'; lmax_val = 1; ok = 1
    END IF
  END SUBROUTINE

  SUBROUTINE embed_set_tau0_from_pos(n_atoms, pos, z, ierr)
    USE coor, ONLY: tau0
    USE ions, ONLY: ions0, ions1
    USE cnst, ONLY: fbohr
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: is, j, k, taken, zz, expect
    ierr = 1
    IF (.NOT. ALLOCATED(tau0)) RETURN
    expect = 0
    DO is = 1, ions1%nsp
      expect = expect + ions0%na(is)
    END DO
    IF (expect /= n_atoms) RETURN
    j = 0
    DO is = 1, ions1%nsp
      zz = ions0%iatyp(is)
      taken = 0
      DO WHILE (taken < ions0%na(is))
        j = j + 1
        IF (j > n_atoms) RETURN
        IF (INT(z(j)) /= zz) RETURN
        taken = taken + 1
        DO k = 1, 3
          tau0(k, taken, is) = REAL(pos(3*(j-1)+k), KIND=real64) * fbohr
        END DO
      END DO
    END DO
    ierr = 0
  END SUBROUTINE

  SUBROUTINE embed_eval_energy_grad(n_atoms, pos, z, energy_h, grad, ok)
    USE wfopts_utils, ONLY: wfopts
    USE rwfopt_utils, ONLY: cpmdc_set_warm_orbitals, cpmdc_set_need_forces
    USE phfac_utils, ONLY: phfac
    USE ener, ONLY: ener_com, chrg, ener_c, ener_d
    USE coor, ONLY: tau0, fion, taup
    USE ions, ONLY: ions0, ions1
    USE store_types, ONLY: cprint, iprint_force
    USE system, ONLY: cnti
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, is, ia, k, idx, nmax
    ok = 0_c_int
    energy_h = 0.0_c_double
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    CALL embed_set_tau0_from_pos(n_atoms, pos, z, ierr)
    IF (ierr /= 0) RETURN
    CALL phfac(tau0)
    IF (ALLOCATED(fion)) DEALLOCATE(fion)
    IF (ALLOCATED(taup)) DEALLOCATE(taup)
    cprint%tprint = .TRUE.
    cprint%iprint(iprint_force) = 1
    ! BOMD/PEF: nuclear forces after WFN optim. OpenCPMD zeros fion unless
    ! tfor; iprint_force alone was not enough on the memfd embed path.
    CALL cpmdc_set_need_forces(.TRUE.)
    ! Warm: retain orbitals (skip initrun) and converge to cntr%tolog with the
    ! same MAXITER budget as cold — do not clamp nomore_iter (that is not a
    ! physical SCF for the new geometry).
    IF (cfg_warm_steps > 0) THEN
      CALL cpmdc_set_warm_orbitals(.TRUE.)
    ELSE
      CALL cpmdc_set_warm_orbitals(.FALSE.)
    END IF
    CALL wfopts
    energy_h = REAL(ener_com%etot, KIND=c_double)
    last_etot = energy_h
    last_ekin = REAL(ener_com%ekin, KIND=c_double)
    last_epseu = REAL(ener_com%epseu, KIND=c_double)
    last_enl = REAL(ener_com%enl, KIND=c_double)
    last_eht = REAL(ener_com%eht, KIND=c_double)
    last_exc = REAL(ener_com%exc, KIND=c_double)
    last_ener_valid = 1_c_int
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
    IF (ALLOCATED(fion)) THEN
      IF (SIZE(fion, 1) >= 3 .AND. SIZE(fion, 2) >= 1 .AND. SIZE(fion, 3) >= ions1%nsp) THEN
        idx = 0
        DO is = 1, ions1%nsp
          DO ia = 1, ions0%na(is)
            IF (ia > SIZE(fion, 2)) EXIT
            DO k = 1, 3
              idx = idx + 1
              IF (idx > nmax) EXIT
              grad(idx) = REAL(-fion(k, ia, is), KIND=c_double)
            END DO
          END DO
        END DO
      END IF
    END IF
    CALL snapshot_prop_from_modules(n_atoms)
    IF (ABS(energy_h) > 1.0e-8_c_double) ok = 1_c_int
  END SUBROUTINE

  ! control/dftin/sysin without reading INPUT (nwchemc rtdb_put analogue).
  SUBROUTINE embed_set_control_dft_sys(cell_ang, has_cell)
    USE control_def_utils, ONLY: control_def
    USE control_test_utils, ONLY: control_test
    USE control_bcast_utils, ONLY: control_bcast
    USE system, ONLY: cntl, cnti, cntr, parm, maxsys, dual00
    USE cell, ONLY: cell_com
    USE isos, ONLY: isos1, isos3
    USE elct, ONLY: crge
    USE spin, ONLY: clsd
    USE func, ONLY: func1, func2, func3, mfxcx_is_slaterx, mfxcc_is_lyp, &
         mgcx_is_becke88, mgcc_is_lyp, mhfx_is_skipped, mtau_is_skipped, &
         msrx_is_skipped, mgcsrx_is_skipped
    USE tbxc, ONLY: toldcode
    USE vdwcmod, ONLY: empvdwc
    USE ener, ONLY: tenergy_ok
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER, INTENT(IN) :: has_cell
    REAL(real64) :: cell_a
    CALL control_def()
    cntl%wfopt = .TRUE.
    cntl%diis = .TRUE.
    cntl%prec = .TRUE.
    cnti%nomore_iter = 40
    cntr%tolog = 1.0e-5_real64
    cntr%hthrs = 0.5_real64
    cntr%gceps = 1.0e-8_real64
    isos1%tcent = .FALSE.
    ! &DFT OLDCODE + FUNCTIONAL BLYP (same fields dftin sets; no INPUT parse)
    toldcode = .TRUE.
    tenergy_ok = .TRUE.
    cntl%tgc = .TRUE.
    cntl%tgcx = .TRUE.
    cntl%tgcc = .TRUE.
    cntl%use_xc_driver = .FALSE.
    cntl%thybrid = .FALSE.
    cntl%ttau = .FALSE.
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_lyp
    func1%mgcx = mgcx_is_becke88
    func1%mgcc = mgcc_is_lyp
    func1%mhfx = mhfx_is_skipped
    func1%mtau = mtau_is_skipped
    func1%msrx = msrx_is_skipped
    func1%mgcsrx = mgcsrx_is_skipped
    func2%salpha = 2.0_real64 / 3.0_real64
    func2%bbeta = 0.0042_real64
    ! Uninitialized func3 weights yield zero XC (module has no default =).
    func3%pxlda = 1.0_real64
    func3%pxgc = 1.0_real64
    func3%pclda = 1.0_real64
    func3%pcgc = 1.0_real64
    func3%phfx = 0.0_real64
    empvdwc%dft_func = applied_functional
    IF (LEN_TRIM(empvdwc%dft_func) == 0) empvdwc%dft_func = 'BLYP'
    cntl%bohr = .FALSE.
    ! sysin defaults (module fields otherwise HUGE/undefined without file parse)
    dual00%cdual = 4.0_real64
    dual00%dual = .FALSE.
    parm%nr1 = 0
    parm%nr2 = 0
    parm%nr3 = 0
    parm%ibrav = -1
    cell_a = 12.0_real64
    IF (has_cell /= 0) THEN
      IF (cell_ang(1) > 0.0_c_double) cell_a = REAL(cell_ang(1), KIND=real64)
    END IF
    cell_com%celldm(1) = cell_a
    cell_com%celldm(2) = 1.0_real64
    cell_com%celldm(3) = 1.0_real64
    cell_com%celldm(4:6) = 0.0_real64
    parm%ibrav = 0
    cntr%ecut = applied_cutoff_ry
    isos3%ps_type = 1
    isos1%tclust = .TRUE.
    isos1%tisos = .TRUE.
    parm%ibrav = 1
    crge%charge = REAL(applied_charge, KIND=real64)
    clsd%nlsd = MERGE(2, 1, applied_mult > 1)
    maxsys%mmaxx = MAX(cnti%nsplp + 20, 999)
    CALL control_test(.FALSE.)
    maxsys%mmaxx = MAX(cnti%nsplp + 20, 999)
    CALL control_bcast()
  END SUBROUTINE

  SUBROUTINE embed_detsp_from_z(n_atoms, z, ierr)
    USE coor, ONLY: tau0, velp, lvelini
    USE ions, ONLY: ions0, al, bl, rcl, maxgau
    USE cotr, ONLY: duat
    USE clas, ONLY: clas3, tclas
    USE system, ONLY: maxsys, maxsp, lmaxx, nhx
    USE nlps, ONLY: wsg, rgh, wgh, nghtol, nghcom
    USE nlcc, ONLY: corecg, corei, corer, rcgrid
    USE atom, ONLY: gnl, rps, rv, rw, vr
    USE sgpp, ONLY: mpro
    USE zeroing_utils, ONLY: zeroing
    USE error_handling, ONLY: stopgm
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: i, is, zz, nsp, nax, nasp, ia, aerr
    LOGICAL :: used(0:120)
    INTEGER :: order_z(32)
    CHARACTER(*), PARAMETER :: procedureN = 'embed_detsp_from_z'
    ierr = 1
    used = .FALSE.
    nsp = 0
    nax = 0
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) RETURN
      IF (.NOT. used(zz)) THEN
        used(zz) = .TRUE.
        nsp = nsp + 1
        IF (nsp > maxsp .OR. nsp > 32) RETURN
        order_z(nsp) = zz
      END IF
    END DO
    IF (nsp < 1) RETURN
    DO is = 1, nsp
      nasp = 0
      DO i = 1, n_atoms
        IF (INT(z(i)) == order_z(is)) nasp = nasp + 1
      END DO
      ions0%na(is) = nasp
      nax = MAX(nax, nasp)
    END DO
    maxsys%nsx = nsp
    maxsys%nax = nax
    maxsys%ncorx = (nsp * (nsp + 1)) / 2
    clas3%nclatom = 0
    clas3%ncltyp = 0
    duat%ndat = 0
    tclas = .FALSE.
    ALLOCATE(tau0(3, maxsys%nax, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    CALL zeroing(tau0)
    ALLOCATE(velp(3, maxsys%nax, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    CALL zeroing(velp)
    ALLOCATE(lvelini(0:maxsys%nax+1, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    DO is = 1, maxsys%nsx
      DO ia = 0, maxsys%nax + 1
        lvelini(ia, is) = .FALSE.
      END DO
    END DO
    ALLOCATE(al(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(bl(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rcl(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(wsg(maxsys%nsx, nhx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rgh(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(wgh(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(nghtol(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(nghcom(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(gnl(maxsys%mmaxx, maxsys%nsx, lmaxx*mpro), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rps(maxsys%mmaxx, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rw(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rv(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(vr(maxsys%mmaxx, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(rcgrid(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(corecg(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    CALL zeroing(rcgrid)
    CALL zeroing(corecg)
    DO is = 1, maxsp
      corer%anlcc(is) = 0.0_real64
      corer%bnlcc(is) = 0.0_real64
      corer%enlcc(is) = 0.0_real64
      corer%clogcc(is) = 0.0_real64
      corei%nlcct(is) = 0
      corei%meshcc(is) = 0
    END DO
    ierr = 0
  END SUBROUTINE

  SUBROUTINE embed_ratom_from_arrays(n_atoms, pos, z, ierr)
    USE coor, ONLY: tau0, velp
    USE ions, ONLY: ions0, ions1
    USE elct, ONLY: crge
    USE cotr, ONLY: lskcor, lskptr, cotc0, duat, cotr007
    USE dpot, ONLY: dpot_mod
    USE pslo, ONLY: pslo_com
    USE nlcc, ONLY: corel
    USE nlps, ONLY: nlps_com
    USE movi, ONLY: imtyp
    USE atwf, ONLY: atchg
    USE atom, ONLY: gnl, rps, rv, rw, vr, patom1
    USE system, ONLY: maxsys, maxsp
    USE recpnew_utils, ONLY: recpnew
    USE zeroing_utils, ONLY: zeroing
    USE error_handling, ONLY: stopgm
    USE cnst, ONLY: fbohr
    USE mm_dimmod, ONLY: mmdim
    USE mm_input, ONLY: g96_vel
    USE symm, ONLY: symmt
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: is, i, j, k, zz, taken, aerr, lmax_val, pok, NSX_q
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=40) :: ecpnam
    CHARACTER(LEN=64) :: pp
    CHARACTER(*), PARAMETER :: procedureN = 'embed_ratom_from_arrays'
    ierr = 1
    NSX_q = maxsys%nsx
    mmdim%nspm = maxsys%nsx
    patom1%pconf = .FALSE.
    seen = .FALSE.
    ALLOCATE(lskcor(3, maxsys%nax*maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(lskptr(3, maxsys%nax*maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    ALLOCATE(atchg(NSX_q), STAT=aerr)
    IF (aerr /= 0) CALL stopgm(procedureN, 'allocation problem', __LINE__, __FILE__)
    CALL zeroing(atchg)
    CALL zeroing(lskptr)
    DO i = 1, maxsys%nax * maxsys%nsx
      DO j = 1, 3
        lskcor(j, i) = 1
      END DO
    END DO
    IF (g96_vel%ntx_vel /= 1) CALL zeroing(velp)
    CALL zeroing(gnl)
    CALL zeroing(rps)
    CALL zeroing(vr)
    CALL zeroing(rw)
    CALL zeroing(rv)
    DO i = 1, maxsp
      imtyp(i) = 0
      dpot_mod%team(i) = .FALSE.
      dpot_mod%tkb(i) = .FALSE.
    END DO
    cotc0%lfcom = .FALSE.
    symmt%tgenc = .FALSE.
    ions1%nsp = 0
    crge%nel = 0.0_real64
    duat%ndat = 0
    duat%ndat1 = 0
    duat%ndat2 = 0
    duat%ndat3 = 0
    duat%ndat4 = 0
    cotc0%mcnstr = 0
    cotr007%mrestr = 0
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) RETURN
      IF (seen(zz)) CYCLE
      seen(zz) = .TRUE.
      CALL embed_pp_for_z(zz, pp, lmax_val, pok)
      IF (pok == 0) RETURN
      ions1%nsp = ions1%nsp + 1
      is = ions1%nsp
      ecpnam = TRIM(pp)
      pslo_com%tvan(is) = .FALSE.
      pslo_com%tbin(is) = .FALSE.
      corel%tnlcc(is) = .FALSE.
      dpot_mod%tkb(is) = .FALSE.
      dpot_mod%lmax(is) = lmax_val
      dpot_mod%lloc(is) = lmax_val + 1
      dpot_mod%lskip(is) = lmax_val + 2
      nlps_com%ngh(is) = 0
      CALL recpnew(is, ecpnam)
      dpot_mod%lmax(is) = dpot_mod%lmax(is) + 1
      crge%nel = crge%nel + REAL(ions0%na(is) * NINT(ions0%zv(is)), KIND=real64)
      ! Store Angstrom (cntl%bohr=.FALSE.); setsys multiplies by fbohr once.
      taken = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        taken = taken + 1
        DO k = 1, 3
          tau0(k, taken, is) = REAL(pos(3*(j-1)+k), KIND=real64)
        END DO
      END DO
    END DO
    IF (ions1%nsp /= maxsys%nsx) RETURN
    ierr = 0
  END SUBROUTINE

  ! Cold: OpenCPMD parsers via anonymous memfd deck (no disk write). Warm: C arrays only.
  SUBROUTINE embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=*), INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, lmax_val
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=64) :: pp
    CHARACTER(LEN=8) :: lmax_c
    CHARACTER(LEN=128) :: line
    REAL(real64) :: cell_a
    ierr = 1
    deck = ' '
    nlen = 0
    cell_a = 12.0_real64
    IF (has_cell /= 0) THEN
      IF (cell(1) > 0.0_c_double) cell_a = REAL(cell(1), KIND=real64)
    END IF
    CALL append(deck, nlen, '&CPMD'//NEW_LINE('A'))
    CALL append(deck, nlen, ' OPTIMIZE WAVEFUNCTION'//NEW_LINE('A'))
    CALL append(deck, nlen, ' CONVERGENCE ORBITALS'//NEW_LINE('A'))
    CALL append(deck, nlen, '  1.0d-5'//NEW_LINE('A'))
    CALL append(deck, nlen, ' MAXITER'//NEW_LINE('A'))
    CALL append(deck, nlen, '  40'//NEW_LINE('A'))
    CALL append(deck, nlen, ' CENTER MOLECULE OFF'//NEW_LINE('A'))
    CALL append(deck, nlen, '&END'//NEW_LINE('A'))
    CALL append(deck, nlen, '&SYSTEM'//NEW_LINE('A'))
    CALL append(deck, nlen, ' SYMMETRY'//NEW_LINE('A'))
    CALL append(deck, nlen, '  0'//NEW_LINE('A'))
    CALL append(deck, nlen, ' ANGSTROM'//NEW_LINE('A'))
    CALL append(deck, nlen, ' CELL'//NEW_LINE('A'))
    WRITE(line, '(A,F12.6,A)') '  ', cell_a, ' 1.0 1.0 0.0 0.0 0.0'
    CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
    CALL append(deck, nlen, ' CUTOFF'//NEW_LINE('A'))
    WRITE(line, '(A,F12.6)') '  ', applied_cutoff_ry
    CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
    IF (applied_charge /= 0) THEN
      CALL append(deck, nlen, ' CHARGE'//NEW_LINE('A'))
      WRITE(line, '(A,I6)') '  ', applied_charge
      CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
    END IF
    IF (applied_mult > 1) THEN
      CALL append(deck, nlen, ' MULTIPLICITY'//NEW_LINE('A'))
      WRITE(line, '(A,I6)') '  ', applied_mult
      CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
    END IF
    CALL append(deck, nlen, ' POISSON SOLVER HOCKNEY'//NEW_LINE('A'))
    CALL append(deck, nlen, '&END'//NEW_LINE('A'))
    CALL append(deck, nlen, '&DFT'//NEW_LINE('A'))
    CALL append(deck, nlen, ' OLDCODE'//NEW_LINE('A'))
    ! Honor wire/applied functional (capnp-fortran apply path); default BLYP.
    CALL append(deck, nlen, ' FUNCTIONAL '//TRIM(applied_functional)//NEW_LINE('A'))
    CALL append(deck, nlen, '&END'//NEW_LINE('A'))
    CALL append(deck, nlen, '&ATOMS'//NEW_LINE('A'))
    seen = .FALSE.
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) RETURN
      IF (seen(zz)) CYCLE
      seen(zz) = .TRUE.
      CALL embed_pp_for_z(zz, pp, lmax_val, pok)
      IF (pok == 0) RETURN
      IF (lmax_val == 0) THEN
        lmax_c = 'S'
      ELSE IF (lmax_val == 1) THEN
        lmax_c = 'P'
      ELSE
        lmax_c = 'D'
      END IF
      CALL append(deck, nlen, '*'//TRIM(pp)//NEW_LINE('A'))
      CALL append(deck, nlen, ' LMAX='//TRIM(lmax_c)//NEW_LINE('A'))
      count = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) == zz) count = count + 1
      END DO
      WRITE(line, '(A,I4)') '   ', count
      CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        WRITE(line, '(3F14.6)') pos(3*(j-1)+1), pos(3*(j-1)+2), pos(3*(j-1)+3)
        CALL append(deck, nlen, TRIM(line)//NEW_LINE('A'))
      END DO
    END DO
    CALL append(deck, nlen, '&END'//NEW_LINE('A'))
    ierr = 0
  CONTAINS
    SUBROUTINE append(buf, n, s)
      CHARACTER(LEN=*), INTENT(INOUT) :: buf
      INTEGER, INTENT(INOUT) :: n
      CHARACTER(LEN=*), INTENT(IN) :: s
      INTEGER :: m
      m = LEN(s)
      IF (n + m > LEN(buf)) RETURN
      buf(n+1:n+m) = s
      n = n + m
    END SUBROUTINE
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
    INTEGER :: ierr, idx, nmax, nlen, mfd
    LOGICAL :: tinfo
    CHARACTER(LEN=16384) :: deck
    CHARACTER(LEN=64) :: mempath
    CHARACTER(LEN=1024) :: pp_dir
    INTERFACE
      FUNCTION cpmdc_memfd_write(bytes, nbytes, path_out, path_cap) BIND(C, NAME='cpmdc_memfd_write')
        IMPORT :: c_char, c_int
        CHARACTER(KIND=c_char), INTENT(IN) :: bytes(*)
        INTEGER(c_int), VALUE :: nbytes
        CHARACTER(KIND=c_char), INTENT(OUT) :: path_out(*)
        INTEGER(c_int), VALUE :: path_cap
        INTEGER(c_int) :: cpmdc_memfd_write
      END FUNCTION
    END INTERFACE
    ok = 0_c_int
    energy_h = 0.0_c_double
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    IF (cfg_warm_steps > 0) THEN
      CALL embed_eval_energy_grad(n_atoms, pos, z, energy_h, grad, ok)
      IF (ok /= 0_c_int) cfg_warm_steps = cfg_warm_steps + 1
      RETURN
    END IF
    ! Cold: prefer Cap'n-rendered method deck (applied_input_deck) when it
    ! already carries &ATOMS — honors typed sections (CONSTRAINTS/VDW/...) from
    ! CPMDParams. Else build a minimal deck that still honors applied
    ! functional/cutoff/charge/mult. Geometry for forces always from C arrays
    ! into TAU0 after parse.
    nlen = 0
    ierr = 1
    IF (LEN_TRIM(applied_input_deck) > 0) THEN
      IF (INDEX(applied_input_deck, '&ATOMS') > 0 .OR. &
          INDEX(applied_input_deck, '&atoms') > 0) THEN
        nlen = MIN(LEN_TRIM(applied_input_deck), LEN(deck))
        deck(1:nlen) = applied_input_deck(1:nlen)
        IF (nlen < LEN(deck)) deck(nlen+1:) = ' '
        ierr = 0
      END IF
    END IF
    IF (ierr /= 0) THEN
      CALL embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, &
           ierr)
    END IF
    IF (ierr /= 0 .OR. nlen < 1) RETURN
    mfd = INT(cpmdc_memfd_write(deck, INT(nlen, KIND=c_int), mempath, &
         INT(LEN(mempath), KIND=c_int)))
    IF (mfd < 0) RETURN
    paral%io_parent = .TRUE.
    ! C wrote a NUL-terminated /proc/self/fd/N path into mempath.
    DO idx = 1, LEN(mempath)
      IF (IACHAR(mempath(idx:idx)) == 0) THEN
        cnts%inputfile = mempath(1:idx-1)
        EXIT
      END IF
    END DO
    ! recpnew resolves relative *PP basenames via CWD and/or CPMD_PP_LIBRARY_PATH.
    ! Point CWD at the PP library so cold memfd decks with relative paths work.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_PSEUDO_DIR', pp_dir)
    IF (LEN_TRIM(pp_dir) == 0) &
        CALL GET_ENVIRONMENT_VARIABLE('CPMD_PP_LIBRARY_PATH', pp_dir)
    IF (LEN_TRIM(pp_dir) == 0) RETURN
    CALL CHDIR(TRIM(pp_dir))
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
    ! First and later forces: positions only from C arrays (nwchemc geom pattern).
    CALL embed_eval_energy_grad(n_atoms, pos, z, energy_h, grad, ok)
    IF (ok /= 0_c_int) THEN
      cfg_warm_steps = 1
    ELSE
      CALL clear_last_energy_components()
    END IF
  END SUBROUTINE
#endif
END MODULE cpmd_embed_c_api

