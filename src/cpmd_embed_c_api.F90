! ISO_C embed: Cap'n Proto knobs -> OpenCPMD modules (libcpmd.a). No INPUT/control I/O.
#include "cpmd_embed_config.h"
MODULE cpmd_embed_c_api
  USE, INTRINSIC :: iso_c_binding
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cpmdc_embed_init, cpmdc_embed_available, cpmdc_embed_set_config
  PUBLIC :: cpmdc_embed_energy_grad, cpmdc_embed_finalize
  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  CHARACTER(LEN=64), SAVE :: cfg_functional = 'BLYP'
  REAL(real64), SAVE :: cfg_cutoff_ry = 70.0_real64
  INTEGER, SAVE :: cfg_charge = 0, cfg_mult = 1
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
  END SUBROUTINE
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
  FUNCTION cpmdc_embed_set_config(functional, functional_len, cutoff_ry, charge, &
      multiplicity, input_deck, input_deck_len, cpmd_root, cpmd_root_len) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_set_config')
    CHARACTER(KIND=c_char), INTENT(IN) :: functional(*), input_deck(*), cpmd_root(*)
    INTEGER(c_int), INTENT(IN), VALUE :: functional_len, charge, multiplicity
    INTEGER(c_int), INTENT(IN), VALUE :: input_deck_len, cpmd_root_len
    REAL(c_double), INTENT(IN), VALUE :: cutoff_ry
    INTEGER(c_int) :: ok
    CALL copy_c_string(functional, functional_len, cfg_functional)
    IF (LEN_TRIM(cfg_functional) == 0) cfg_functional = 'BLYP'
    cfg_cutoff_ry = REAL(cutoff_ry, real64)
    IF (cfg_cutoff_ry <= 0.0_real64) cfg_cutoff_ry = 70.0_real64
    cfg_charge = INT(charge)
    cfg_mult = MAX(1, INT(multiplicity))
#if defined(CPMDC_HAS_CPMD)
    CALL apply_method_knobs()
    ok = MERGE(1_c_int, 0_c_int, runtime_ready)
#else
    ok = 0_c_int
#endif
  END FUNCTION
  FUNCTION cpmdc_embed_energy_grad(n_atoms, positions_ang, atomic_numbers, cell_ang, &
      has_cell, energy_h, grad_h_bohr) RESULT(ok) BIND(C, NAME='cpmdc_embed_energy_grad')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: positions_ang(*), cell_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(OUT) :: energy_h, grad_h_bohr(*)
    INTEGER(c_int) :: ok
    INTEGER :: i
    energy_h = 0.0_c_double
    DO i = 1, MAX(0, INT(n_atoms)) * 3
      grad_h_bohr(i) = 0.0_c_double
    END DO
#if defined(CPMDC_HAS_CPMD)
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. n_atoms <= 0) THEN
      ok = 0_c_int
      RETURN
    END IF
    ! Full wfopts/setsys SCF requires complete libcpmd object set. When those
    ! modules are linked, extend here. Until then: apply knobs (control_def)
    ! and report failure so we never invent energies.
    CALL apply_method_knobs()
    ok = 0_c_int
#else
    ok = 0_c_int
#endif
  END FUNCTION
  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
  END SUBROUTINE
#if defined(CPMDC_HAS_CPMD)
  SUBROUTINE apply_method_knobs()
    USE control_def_utils, ONLY: control_def
    USE system, ONLY: cntr
    USE spin, ONLY: clsd
    USE func, ONLY: func1, mfxcx_is_slaterx, mfxcc_is_lyp, mgcx_is_becke88, mgcc_is_lyp
    IMPLICIT NONE
    CALL control_def()
    cntr%ecut = REAL(cfg_cutoff_ry, KIND=KIND(cntr%ecut))
    IF (cfg_mult > 1) THEN
      clsd%nlsd = 2
    ELSE
      clsd%nlsd = 1
    END IF
    ! BLYP defaults for the water demo functional name.
    func1%mfxcx = mfxcx_is_slaterx
    func1%mfxcc = mfxcc_is_lyp
    func1%mgcx = mgcx_is_becke88
    func1%mgcc = mgcc_is_lyp
  END SUBROUTINE
#endif
END MODULE cpmd_embed_c_api
