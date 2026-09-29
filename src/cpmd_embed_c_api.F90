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
  PUBLIC :: cpmdc_embed_bind_calculator
  PUBLIC :: cpmdc_embed_reset_state
  PUBLIC :: cpmdc_embed_set_config, cpmdc_embed_set_deck, cpmdc_embed_energy_grad
  ! cpmdc_embed_compose_cold_deck: BIND(C) in HAS_CPMD / stub branches (not listed
  ! in PUBLIC — gfortran rejects forward PUBLIC when the body is ifdef-gated).

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  ! Config state is applied_* from cpmdc_embed_apply_params_mod.
  ! Results and the warm cell live in the caller image (CPMDCEmbedImage).
  ! runtime_ready is process-wide. tcpu0 is the timer origin of one SCF.

  TYPE, BIND(C) :: cpmdc_energy_components
    INTEGER(c_int) :: valid
    REAL(c_double) :: etot, ekin, epseu, enl, eht, ehep, ehee, ehii
    REAL(c_double) :: exc, vxc, egc, esr, eeig, eband, entropy, eself
    REAL(c_double) :: ecnstr, amu, ebogo, eext, etddft, ehsic, erestr, eefield
  END TYPE
  TYPE, BIND(C) :: cpmdc_charge_integrals
    INTEGER(c_int) :: valid
    REAL(c_double) :: csumg, csumr, csums, csumsabs
  END TYPE
  TYPE, BIND(C) :: cpmdc_multi_state
    INTEGER(c_int) :: valid
    INTEGER(c_size_t) :: count
    REAL(c_double) :: values(64)
  END TYPE
  TYPE, BIND(C) :: cpmdc_md_row
    INTEGER(c_int) :: valid
    INTEGER(c_size_t) :: count
    REAL(c_double) :: values(32)
  END TYPE
  TYPE, BIND(C) :: cpmdc_property_snapshot
    INTEGER(c_int) :: valid
    INTEGER(c_size_t) :: hessian_count
    REAL(c_double) :: hessian(4096)
    INTEGER(c_size_t) :: dipole_count
    REAL(c_double) :: dipole(3)
    INTEGER(c_size_t) :: polarizability_count
    REAL(c_double) :: polarizability(9)
  END TYPE
  TYPE, BIND(C) :: cpmdc_stress_tensor
    INTEGER(c_int) :: valid
    REAL(c_double) :: values(9)
  END TYPE
  TYPE, BIND(C) :: cpmdc_embed_image
    TYPE(cpmdc_energy_components) :: energy
    TYPE(cpmdc_charge_integrals) :: charge
    TYPE(cpmdc_multi_state) :: multi
    TYPE(cpmdc_md_row) :: md
    TYPE(cpmdc_property_snapshot) :: prop
    TYPE(cpmdc_stress_tensor) :: stress
    REAL(c_double) :: warm_cell(9)
    INTEGER(c_int) :: warm_cell_set
    INTEGER(c_int) :: warm_has_cell
    INTEGER(c_int) :: cfg_warm_steps
  END TYPE
#if defined(CPMDC_HAS_CPMD)
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double
#endif

CONTAINS

  SUBROUTINE clear_last_energy_components(image)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    image%energy%valid = 0_c_int
    image%energy%etot = 0.0_c_double
    image%energy%ekin = 0.0_c_double
    image%energy%epseu = 0.0_c_double
    image%energy%enl = 0.0_c_double
    image%energy%eht = 0.0_c_double
    image%energy%ehep = 0.0_c_double
    image%energy%ehee = 0.0_c_double
    image%energy%ehii = 0.0_c_double
    image%energy%exc = 0.0_c_double
    image%energy%vxc = 0.0_c_double
    image%energy%egc = 0.0_c_double
    image%energy%esr = 0.0_c_double
    image%energy%eeig = 0.0_c_double
    image%energy%eband = 0.0_c_double
    image%energy%entropy = 0.0_c_double
    image%energy%eself = 0.0_c_double
    image%energy%ecnstr = 0.0_c_double
    image%energy%amu = 0.0_c_double
    image%energy%ebogo = 0.0_c_double
    image%energy%eext = 0.0_c_double
    image%energy%etddft = 0.0_c_double
    image%energy%ehsic = 0.0_c_double
    image%energy%erestr = 0.0_c_double
    image%energy%eefield = 0.0_c_double

    image%charge%valid = 0_c_int
    image%charge%csumg = 0.0_c_double
    image%charge%csumr = 0.0_c_double
    image%charge%csums = 0.0_c_double
    image%charge%csumsabs = 0.0_c_double
    image%multi%valid = 0_c_int
    image%multi%count = 0_c_size_t
    image%multi%values = 0.0_c_double
    image%md%valid = 0_c_int
    image%md%count = 0_c_size_t
    image%md%values = 0.0_c_double
    image%prop%valid = 0_c_int
    image%prop%hessian_count = 0_c_size_t
    image%prop%hessian = 0.0_c_double
    image%prop%dipole_count = 0_c_size_t
    image%prop%dipole = 0.0_c_double
    image%prop%polarizability_count = 0_c_size_t
    image%prop%polarizability = 0.0_c_double
    image%stress%valid = 0_c_int
    image%stress%values = 0.0_c_double
  END SUBROUTINE

  SUBROUTINE snapshot_total_only(image, energy_h)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    REAL(c_double), INTENT(IN) :: energy_h
    INTEGER :: i
    CALL clear_last_energy_components(image)
    image%energy%etot = energy_h
    image%energy%valid = 1_c_int
    ! PEF ENERGY-style row (EKINC=0 for non-MD reference evaluator).
    image%md%values(1) = energy_h
    image%md%values(2) = 0.0_c_double
    image%md%values(3) = 0.0_c_double
    image%md%values(4) = 0.0_c_double
    image%md%values(5) = 0.0_c_double
    image%md%values(6) = 0.0_c_double
    image%md%values(7) = 0.0_c_double
    image%md%values(8) = 0.0_c_double
    image%md%values(9) = 0.0_c_double
    image%md%values(10) = 0.0_c_double
    image%md%values(11) = 0.0_c_double
    image%md%values(12) = 0.0_c_double  ! EKINC (fictitious e- KE); zero off MD
    image%md%count = 12_c_size_t
    image%md%valid = 1_c_int
    image%charge%valid = 1_c_int
    image%charge%csumg = 0.0_c_double
    image%charge%csumr = 0.0_c_double
    image%charge%csums = 0.0_c_double
    image%charge%csumsabs = 0.0_c_double
    image%multi%count = 6_c_size_t
    image%multi%values = 0.0_c_double
    image%multi%values(1) = energy_h
    image%multi%valid = 1_c_int
    image%prop%valid = 1_c_int
    image%prop%dipole_count = 3_c_size_t
    image%prop%dipole = 0.0_c_double
    image%prop%polarizability_count = 9_c_size_t
    DO i = 1, 9
      image%prop%polarizability(i) = 0.0_c_double
    END DO
    image%prop%hessian_count = 0_c_size_t
    image%stress%valid = 0_c_int
    image%stress%values = 0.0_c_double
  END SUBROUTINE


  FUNCTION cpmdc_embed_bind_calculator(ranks_per_calc) RESULT(calc) &
      BIND(C, NAME='cpmdc_embed_bind_calculator')
    INTEGER(c_int), INTENT(IN), VALUE :: ranks_per_calc
    INTEGER(c_int) :: calc
    calc = -1_c_int
#if defined(CPMDC_HAS_CPMD)
    BLOCK
      USE mpi
      USE mp_interface, ONLY: mp_comm_set, mp_comm_world
      INTEGER :: ierr, rank, npe, rpc, color, key, comm
      LOGICAL :: inited
      rpc = INT(ranks_per_calc)
      CALL MPI_INITIALIZED(inited, ierr)
      IF (.NOT. inited) CALL MPI_Init(ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD, npe, ierr)
      IF (rpc <= 0) rpc = npe
      IF (rpc > npe .OR. MOD(npe, rpc) /= 0) RETURN
      calc = INT(rank / rpc, c_int)
      IF (mp_comm_set) RETURN
      color = rank / rpc
      key = MOD(rank, rpc)
      CALL MPI_Comm_split(MPI_COMM_WORLD, color, key, comm, ierr)
      IF (ierr /= 0) THEN
        calc = -1_c_int
        RETURN
      END IF
      mp_comm_world = comm
      mp_comm_set = .TRUE.
    END BLOCK
#endif
  END FUNCTION cpmdc_embed_bind_calculator

  FUNCTION cpmdc_embed_init() RESULT(ok) BIND(C, NAME='cpmdc_embed_init')
    INTEGER(c_int) :: ok
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    ok = 1_c_int
  END FUNCTION

  FUNCTION cpmdc_embed_available() RESULT(ok) BIND(C, NAME='cpmdc_embed_available')
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#else
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#endif
  END FUNCTION

  FUNCTION cpmdc_embed_reset_state(image_c) RESULT(ok) BIND(C, NAME='cpmdc_embed_reset_state')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: cpmdc_reset_warm_orbitals
#endif
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
    CALL clear_last_energy_components(image)
    image%cfg_warm_steps = 0_c_int
    image%warm_cell_set = 0_c_int
    image%warm_has_cell = 0_c_int
    image%warm_cell = 0.0_c_double
#if defined(CPMDC_HAS_CPMD)
    CALL cpmdc_reset_warm_orbitals()
#endif
    ok = 1_c_int
  END FUNCTION

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
  END SUBROUTINE

  FUNCTION cpmdc_embed_set_config(functional, functional_len, cutoff_ry, charge, &
      multiplicity, input_deck, input_deck_len, cpmd_root, cpmd_root_len, image_c) &
      RESULT(ok) BIND(C, NAME='cpmdc_embed_set_config')
    CHARACTER(KIND=c_char), INTENT(IN) :: functional(*)
    INTEGER(c_int), INTENT(IN), VALUE :: functional_len
    REAL(c_double), INTENT(IN), VALUE :: cutoff_ry
    INTEGER(c_int), INTENT(IN), VALUE :: charge, multiplicity
    CHARACTER(KIND=c_char), INTENT(IN) :: input_deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: input_deck_len
    CHARACTER(KIND=c_char), INTENT(IN) :: cpmd_root(*)
    INTEGER(c_int), INTENT(IN), VALUE :: cpmd_root_len
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (functional_len < 0 .OR. input_deck_len < 0 .OR. cpmd_root_len < 0) RETURN
    IF (cutoff_ry < 0.0_c_double) RETURN
    IF (cpmdc_embed_reset_state(image_c) == 0_c_int) RETURN
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
      cell_ang, has_cell, energy_h, grad_h_bohr, image_c) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_energy_grad')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad_h_bohr(*)
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    INTEGER :: i, n3
    ok = 0_c_int
    energy_h = 0.0_c_double
    n3 = MAX(0, INT(n_atoms) * 3)
    DO i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    END DO
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. n_atoms <= 0) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
#if defined(CPMDC_HAS_CPMD)
    CALL run_embed_scf(image, INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, ok)
#else
    IF (has_cell < 0) RETURN
    CALL run_reference_pef(image, INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, ok)
#endif
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
  SUBROUTINE run_reference_pef(image, n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
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
    CALL snapshot_total_only(image, energy_h)
    DO i = 1, n_atoms * 3
      IF (i > 4096) EXIT
      image%prop%hessian(i) = grad(i)
    END DO
    image%prop%hessian_count = INT(MIN(n_atoms * 3, 4096), KIND=c_size_t)
    image%prop%valid = 1_c_int
    ! Toy isotropic stress (Ha/Bohr^3) so PotentialResult.stress is exercised
    ! without OpenCPMD: sigma_ii ~ energy / (cell volume in Bohr^3) when cell
    ! present, else zeros with valid set.
    image%stress%values = 0.0_c_double
    IF (has_cell /= 0) THEN
      CALL reference_pef_fill_stress(image, cell, energy_h)
    END IF
    image%stress%valid = 1_c_int
  END SUBROUTINE

  SUBROUTINE reference_pef_fill_stress(image, cell, energy_h)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    REAL(c_double), INTENT(IN) :: cell(*)
    REAL(c_double), INTENT(IN) :: energy_h
    REAL(real64), PARAMETER :: bohr_to_ang = 0.529177210903_real64
    REAL(real64) :: a(3), b(3), c(3), vol, inv_b, sig
    INTEGER :: i
    inv_b = 1.0_real64 / bohr_to_ang
    DO i = 1, 3
      a(i) = REAL(cell(i), KIND=real64) * inv_b
      b(i) = REAL(cell(3 + i), KIND=real64) * inv_b
      c(i) = REAL(cell(6 + i), KIND=real64) * inv_b
    END DO
    vol = ABS(a(1) * (b(2) * c(3) - b(3) * c(2)) + &
         a(2) * (b(3) * c(1) - b(1) * c(3)) + &
         a(3) * (b(1) * c(2) - b(2) * c(1)))
    IF (vol <= 1.0e-12_real64) RETURN
    sig = REAL(energy_h, KIND=real64) / vol
    image%stress%values(1) = REAL(sig, KIND=c_double)
    image%stress%values(5) = REAL(sig, KIND=c_double)
    image%stress%values(9) = REAL(sig, KIND=c_double)
  END SUBROUTINE

  ! Cold-deck helpers available without linking OpenCPMD (cmocka stub path).
  SUBROUTINE embed_pp_for_z_local(zz, pp, lmax_val, ok)
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
      pp = 'Si_MT_BLYP.psp'; lmax_val = 2; ok = 1
    ELSE IF (zz == 32) THEN
      pp = 'Ge_MT_BLYP.psp'; lmax_val = 1; ok = 1
    END IF
  END SUBROUTINE

  ! True only when &ATOMS has PP stars AND at least one coordinate triple
  ! before its &END. Params often render *PP.psp stubs without coords; those
  ! must take the method+geometry merge path so ForceInput positions inject.
  LOGICAL FUNCTION deck_has_real_atoms_local(d)
    CHARACTER(LEN=*), INTENT(IN) :: d
    INTEGER :: ia, star, n, iend, ls, le
    deck_has_real_atoms_local = .FALSE.
    ia = INDEX(d, '&ATOMS')
    IF (ia <= 0) ia = INDEX(d, '&atoms')
    IF (ia <= 0) RETURN
    star = INDEX(d(ia:), '*')
    IF (star <= 0) RETURN
    iend = INDEX(d(ia:), '&END')
    IF (iend <= 0) iend = INDEX(d(ia:), '&end')
    IF (iend > 0) THEN
      n = ia + iend - 2
    ELSE
      n = LEN_TRIM(d)
    END IF
    ! A coordinate line has three or more numeric tokens. '*file' lines
    ! and option lines (LMAX=, LOC=, ...) never count.
    ls = ia
    DO WHILE (ls <= n)
      le = INDEX(d(ls:n), NEW_LINE('A'))
      IF (le == 0) THEN
        le = n
      ELSE
        le = ls + le - 2
      END IF
      IF (coord_line(d(ls:le))) THEN
        deck_has_real_atoms_local = .TRUE.
        RETURN
      END IF
      ls = le + 2
    END DO
  CONTAINS
    LOGICAL FUNCTION coord_line(line)
      CHARACTER(LEN=*), INTENT(IN) :: line
      INTEGER :: j, ntok
      LOGICAL :: in_tok, tok_ok, tok_digit
      CHARACTER(LEN=1) :: ch
      coord_line = .FALSE.
      IF (INDEX(line, '*') > 0 .OR. INDEX(line, '=') > 0) RETURN
      ntok = 0
      in_tok = .FALSE.
      tok_ok = .TRUE.
      tok_digit = .FALSE.
      DO j = 1, LEN(line) + 1
        IF (j <= LEN(line)) THEN
          ch = line(j:j)
        ELSE
          ch = ' '
        END IF
        IF (ch == ' ' .OR. ch == CHAR(9) .OR. ch == ',') THEN
          IF (in_tok .AND. tok_ok .AND. tok_digit) ntok = ntok + 1
          in_tok = .FALSE.
          tok_ok = .TRUE.
          tok_digit = .FALSE.
        ELSE
          in_tok = .TRUE.
          IF (ch >= '0' .AND. ch <= '9') THEN
            tok_digit = .TRUE.
          ELSE IF (INDEX('.+-eEdD', ch) == 0) THEN
            tok_ok = .FALSE.
          END IF
        END IF
      END DO
      coord_line = ntok >= 3
    END FUNCTION
  END FUNCTION

  LOGICAL FUNCTION deck_has_method_sections_local(d)
    CHARACTER(LEN=*), INTENT(IN) :: d
    deck_has_method_sections_local = &
        INDEX(d, '&DFT') > 0 .OR. INDEX(d, '&dft') > 0 .OR. &
        INDEX(d, '&SYSTEM') > 0 .OR. INDEX(d, '&system') > 0 .OR. &
        INDEX(d, '&CPMD') > 0 .OR. INDEX(d, '&cpmd') > 0
  END FUNCTION

  SUBROUTINE strip_atoms_sections_local(src, dst, nlen)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER, INTENT(OUT) :: nlen
    INTEGER :: i, n, end_at, j
    CHARACTER(LEN=16) :: tag
    dst = ' '
    nlen = 0
    n = LEN_TRIM(src)
    IF (n < 1) RETURN
    i = 1
    DO WHILE (i <= n)
      IF (i + 5 <= n) THEN
        tag = src(i:MIN(i + 5, n))
        IF (tag == '&ATOMS' .OR. tag == '&atoms') THEN
          end_at = 0
          j = i + 6
          DO WHILE (j + 3 <= n)
            IF (src(j:j+3) == '&END' .OR. src(j:j+3) == '&end') THEN
              end_at = j + 3
              EXIT
            END IF
            j = j + 1
          END DO
          IF (end_at > 0) THEN
            i = end_at + 1
            DO WHILE (i <= n .AND. (src(i:i) == NEW_LINE('A') .OR. src(i:i) == ' '))
              i = i + 1
            END DO
            CYCLE
          END IF
        END IF
      END IF
      IF (nlen < LEN(dst)) THEN
        nlen = nlen + 1
        dst(nlen:nlen) = src(i:i)
      END IF
      i = i + 1
    END DO
  END SUBROUTINE

  SUBROUTINE embed_method_deck_plus_atoms_local(n_atoms, pos, z, cell, has_cell, &
      deck, nlen, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=*), INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, lmax_val, base, mlen
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=64) :: pp
    CHARACTER(LEN=8) :: lmax_c
    CHARACTER(LEN=128) :: line
    CHARACTER(LEN=4096) :: method
    ierr = 1
    deck = ' '
    nlen = 0
    CALL strip_atoms_sections_local(applied_input_deck, method, mlen)
    base = MIN(mlen, LEN(deck) - 64)
    IF (base < 1) RETURN
    deck(1:base) = method(1:base)
    nlen = base
    IF (deck(nlen:nlen) /= NEW_LINE('A')) THEN
      IF (nlen < LEN(deck)) THEN
        nlen = nlen + 1
        deck(nlen:nlen) = NEW_LINE('A')
      END IF
    END IF
    CALL append_local(deck, nlen, '&ATOMS'//NEW_LINE('A'))
    seen = .FALSE.
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) RETURN
      IF (seen(zz)) CYCLE
      seen(zz) = .TRUE.
      CALL embed_pp_for_z_local(zz, pp, lmax_val, pok)
      IF (pok == 0) RETURN
      IF (lmax_val == 0) THEN
        lmax_c = 'S'
      ELSE IF (lmax_val == 1) THEN
        lmax_c = 'P'
      ELSE
        lmax_c = 'D'
      END IF
      IF (INDEX(applied_input_deck, 'KLEINMAN-BYLANDER') > 0) THEN
        CALL append_local(deck, nlen, '*'//TRIM(pp)//' KLEINMAN-BYLANDER'// &
             NEW_LINE('A'))
      ELSE
        CALL append_local(deck, nlen, '*'//TRIM(pp)//NEW_LINE('A'))
      END IF
      IF (zz == 14) THEN
        CALL append_local(deck, nlen, ' LMAX='//TRIM(lmax_c)//' LOC=D'//NEW_LINE('A'))
      ELSE
        CALL append_local(deck, nlen, ' LMAX='//TRIM(lmax_c)//NEW_LINE('A'))
      END IF
      count = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) == zz) count = count + 1
      END DO
      WRITE(line, '(A,I4)') '   ', count
      CALL append_local(deck, nlen, TRIM(line)//NEW_LINE('A'))
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        WRITE(line, '(3F14.6)') pos(3*(j-1)+1), pos(3*(j-1)+2), pos(3*(j-1)+3)
        CALL append_local(deck, nlen, TRIM(line)//NEW_LINE('A'))
      END DO
    END DO
    CALL append_local(deck, nlen, '&END'//NEW_LINE('A'))
    IF (has_cell < 0) RETURN
    IF (cell(1) < -1.0e300_c_double) RETURN
    ierr = 0
  CONTAINS
    SUBROUTINE append_local(buf, n, s)
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

  SUBROUTINE embed_compose_cold_deck_local(n_atoms, pos, z, cell, has_cell, &
      deck, nlen, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=*), INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    nlen = 0
    ierr = 1
    deck = ' '
    IF (LEN_TRIM(applied_input_deck) > 0) THEN
      IF (deck_has_real_atoms_local(applied_input_deck)) THEN
        nlen = MIN(LEN_TRIM(applied_input_deck), LEN(deck))
        deck(1:nlen) = applied_input_deck(1:nlen)
        IF (nlen < LEN(deck)) deck(nlen+1:) = ' '
        ierr = 0
      ELSE IF (deck_has_method_sections_local(applied_input_deck)) THEN
        CALL embed_method_deck_plus_atoms_local(n_atoms, pos, z, cell, has_cell, &
             deck, nlen, ierr)
      END IF
    END IF
  END SUBROUTINE

  FUNCTION cpmdc_embed_compose_cold_deck(n_atoms, positions_ang, atomic_numbers, &
      cell_ang, has_cell, deck_out, deck_cap, deck_len) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_compose_cold_deck')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    CHARACTER(KIND=c_char), INTENT(OUT) :: deck_out(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_cap
    INTEGER(c_int), INTENT(OUT) :: deck_len
    INTEGER(c_int) :: ok
    CHARACTER(LEN=16384) :: deck
    INTEGER :: nlen, ierr, i, ncopy
    ok = 0_c_int
    deck_len = 0_c_int
    IF (deck_cap < 2) RETURN
    CALL embed_compose_cold_deck_local(INT(n_atoms), positions_ang, &
         atomic_numbers, cell_ang, INT(has_cell), deck, nlen, ierr)
    IF (ierr /= 0 .OR. nlen < 1) RETURN
    ncopy = MIN(nlen, INT(deck_cap) - 1)
    DO i = 1, ncopy
      deck_out(i) = deck(i:i)
    END DO
    deck_out(ncopy + 1) = c_null_char
    deck_len = INT(ncopy, KIND=c_int)
    ok = 1_c_int
  END FUNCTION
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

  SUBROUTINE snapshot_prop_from_modules(image, n_atoms)
    USE ddip, ONLY: pdipole
    USE coor, ONLY: fion
    USE ions, ONLY: ions0, ions1
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER :: is, ia, k, idx, i
    image%prop%valid = 1_c_int
    image%prop%dipole_count = 3_c_size_t
    DO i = 1, 3
      image%prop%dipole(i) = REAL(pdipole(i), KIND=c_double)
    END DO
    image%prop%polarizability_count = 9_c_size_t
    DO i = 1, 9
      image%prop%polarizability(i) = 0.0_c_double
    END DO
    idx = 0
    IF (ALLOCATED(fion)) THEN
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          DO k = 1, 3
            idx = idx + 1
            IF (idx > 4096) EXIT
            image%prop%hessian(idx) = REAL(-fion(k, ia, is), KIND=c_double)
          END DO
        END DO
      END DO
    END IF
    image%prop%hessian_count = INT(MIN(idx, 4096), KIND=c_size_t)
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
      pp = 'Si_MT_BLYP.psp'; lmax_val = 2; ok = 1
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

  SUBROUTINE embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, ok)
    USE wfopts_utils, ONLY: wfopts
    USE rwfopt_utils, ONLY: cpmdc_set_warm_orbitals, cpmdc_set_need_forces
    USE phfac_utils, ONLY: phfac
    USE ener, ONLY: ener_com, chrg, ener_c, ener_d
    USE coor, ONLY: tau0, fion, taup
    USE ions, ONLY: ions0, ions1
    USE store_types, ONLY: cprint, iprint_force, restart1
    USE system, ONLY: cnti, cntl, parm
    USE strs, ONLY: paiu
    USE isos, ONLY: isos1
    USE ropt, ONLY: ropt_mod
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, is, ia, k, idx, nmax, i, j
    REAL(real64) :: omega
    ok = 0_c_int
    energy_h = 0.0_c_double
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    image%stress%valid = 0_c_int
    image%stress%values = 0.0_c_double
    CALL embed_set_tau0_from_pos(n_atoms, pos, z, ierr)
    IF (ierr /= 0) RETURN
    ! The caller (eOn, rgmin, rgsaddle) owns the ionic geometry and the
    ! ionic velocities. initrun -> zhrwf would replace tau0 from RESTART
    ! section 3 when RESTART COORDINATES or RESTART ALL is in the deck.
    ! The wavefunction section stays under the deck: RESTART WAVEFUNCTION
    ! still fills c0 on a cold entry. Plane-wave velocities stay unread
    ! for a wavefunction optimisation, which is setirec's BOMD rule.
    restart1%rco = .FALSE.
    restart1%rvel = .FALSE.
    restart1%rgeo = .FALSE.
    CALL phfac(tau0)
    IF (ALLOCATED(fion)) DEALLOCATE(fion)
    IF (ALLOCATED(taup)) DEALLOCATE(taup)
    cprint%tprint = .TRUE.
    cprint%iprint(iprint_force) = 1
    ! BOMD/PEF: nuclear forces after WFN optim. OpenCPMD zeros fion unless
    ! tfor; iprint_force alone was not enough on the memfd embed path.
    CALL cpmdc_set_need_forces(.TRUE.)
    ! PEF stress: totstr fills paiu when cntl%tpres. Skip isolated/Hockney
    ! (tclust): rinitwf does tpres→newcell→gf_periodic but scg is only
    ! allocated for periodic cells in initclust — SEGV on cluster decks.
    IF (.NOT. isos1%tclust .AND. embed_stress_wanted()) cntl%tpres = .TRUE.
    ! Warm: retain orbitals (skip initrun) and converge to cntr%tolog with the
    ! same MAXITER budget as cold — do not clamp nomore_iter (that is not a
    ! physical SCF for the new geometry).
    IF (image%cfg_warm_steps > 0) THEN
      CALL cpmdc_set_warm_orbitals(.TRUE.)
    ELSE
      CALL cpmdc_set_warm_orbitals(.FALSE.)
    END IF
    CALL wfopts
    energy_h = REAL(ener_com%etot, KIND=c_double)
    image%energy%etot = energy_h
    image%energy%ekin = REAL(ener_com%ekin, KIND=c_double)
    image%energy%epseu = REAL(ener_com%epseu, KIND=c_double)
    image%energy%enl = REAL(ener_com%enl, KIND=c_double)
    image%energy%eht = REAL(ener_com%eht, KIND=c_double)
    image%energy%exc = REAL(ener_com%exc, KIND=c_double)
    image%energy%valid = 1_c_int
    image%charge%csumg = REAL(chrg%csumg, KIND=c_double)
    image%charge%csumr = REAL(chrg%csumr, KIND=c_double)
    image%charge%csums = REAL(chrg%csums, KIND=c_double)
    image%charge%csumsabs = REAL(chrg%csumsabs, KIND=c_double)
    image%charge%valid = 1_c_int
    image%multi%values(1) = REAL(ener_c%etot_a, KIND=c_double)
    image%multi%values(2) = REAL(ener_c%etot_2, KIND=c_double)
    image%multi%values(3) = REAL(ener_c%etot_ab, KIND=c_double)
    image%multi%values(4) = REAL(ener_d%etot_b, KIND=c_double)
    image%multi%values(5) = REAL(ener_d%ecas, KIND=c_double)
    image%multi%values(6) = REAL(ener_d%etot_t, KIND=c_double)
    image%multi%count = 6_c_size_t
    image%multi%valid = 1_c_int
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
    ! Cartesian stress Ha/Bohr^3: OpenCPMD stores virial in paiu (energy),
    ! true stress is paiu/omega (see totstr/wrstress). Row-major for Cap'n Proto.
    omega = parm%omega
    IF (omega > 1.0e-30_real64) THEN
      DO i = 1, 3
        DO j = 1, 3
          image%stress%values(3 * (i - 1) + j) = REAL(paiu(i, j) / omega, KIND=c_double)
        END DO
      END DO
      image%stress%valid = 1_c_int
    END IF
    CALL snapshot_prop_from_modules(image, n_atoms)
    ! rwfopt computes the ionic forces only for converged orbitals: an SCF
    ! that ran out of MAXITER leaves fion zero, which a caller would read as
    ! a stationary point. Report it as a failure instead.
    IF (ropt_mod%convwf .AND. ABS(energy_h) > 1.0e-8_c_double) ok = 1_c_int
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
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: i, is, zz, nsp, nax, nasp, ia, aerr
    LOGICAL :: used(0:120)
    INTEGER :: order_z(32)
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
    IF (aerr /= 0) RETURN
    CALL zeroing(tau0)
    ALLOCATE(velp(3, maxsys%nax, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    CALL zeroing(velp)
    ALLOCATE(lvelini(0:maxsys%nax+1, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    DO is = 1, maxsys%nsx
      DO ia = 0, maxsys%nax + 1
        lvelini(ia, is) = .FALSE.
      END DO
    END DO
    ALLOCATE(al(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(bl(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rcl(maxgau, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(wsg(maxsys%nsx, nhx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rgh(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(wgh(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(nghtol(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(nghcom(nhx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(gnl(maxsys%mmaxx, maxsys%nsx, lmaxx*mpro), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rps(maxsys%mmaxx, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rw(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rv(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(vr(maxsys%mmaxx, maxsys%nsx, lmaxx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(rcgrid(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(corecg(maxsys%mmaxx, maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
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
    ierr = 1
    NSX_q = maxsys%nsx
    mmdim%nspm = maxsys%nsx
    patom1%pconf = .FALSE.
    seen = .FALSE.
    ALLOCATE(lskcor(3, maxsys%nax*maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(lskptr(3, maxsys%nax*maxsys%nsx), STAT=aerr)
    IF (aerr /= 0) RETURN
    ALLOCATE(atchg(NSX_q), STAT=aerr)
    IF (aerr /= 0) RETURN
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
    CHARACTER(LEN=400) :: celltxt
    INTEGER :: celln
    ierr = 1
    deck = ' '
    nlen = 0
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
    CALL format_cell_lines(cell, has_cell, celltxt, celln)
    CALL append(deck, nlen, celltxt(1:celln))
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
      IF (zz == 14) THEN
        CALL append(deck, nlen, ' LMAX='//TRIM(lmax_c)//' LOC=D'//NEW_LINE('A'))
      ELSE
        CALL append(deck, nlen, ' LMAX='//TRIM(lmax_c)//NEW_LINE('A'))
      END IF
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

  ! True only when &ATOMS has PP stars AND at least one coordinate triple
  ! before its &END. C render *PP.psp stubs without coords must NOT block
  ! method+geometry merge from ForceInput positions.
  ! CPMDC_STRESS=0 skips the stress tensor on periodic cells. A caller that
  ! takes only energy and forces does not need totstr.
  LOGICAL FUNCTION embed_stress_wanted()
    CHARACTER(LEN=16) :: v
    INTEGER :: st
    embed_stress_wanted = .TRUE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_STRESS', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '0') embed_stress_wanted = .FALSE.
  END FUNCTION

  LOGICAL FUNCTION deck_has_real_atoms(d)
    CHARACTER(LEN=*), INTENT(IN) :: d
    INTEGER :: ia, star, n, iend, ls, le
    deck_has_real_atoms = .FALSE.
    ia = INDEX(d, '&ATOMS')
    IF (ia <= 0) ia = INDEX(d, '&atoms')
    IF (ia <= 0) RETURN
    star = INDEX(d(ia:), '*')
    IF (star <= 0) RETURN
    iend = INDEX(d(ia:), '&END')
    IF (iend <= 0) iend = INDEX(d(ia:), '&end')
    IF (iend > 0) THEN
      n = ia + iend - 2
    ELSE
      n = LEN_TRIM(d)
    END IF
    ! A coordinate line has three or more numeric tokens. '*file' lines
    ! and option lines (LMAX=, LOC=, ...) never count.
    ls = ia
    DO WHILE (ls <= n)
      le = INDEX(d(ls:n), NEW_LINE('A'))
      IF (le == 0) THEN
        le = n
      ELSE
        le = ls + le - 2
      END IF
      IF (coord_line(d(ls:le))) THEN
        deck_has_real_atoms = .TRUE.
        RETURN
      END IF
      ls = le + 2
    END DO
  CONTAINS
    LOGICAL FUNCTION coord_line(line)
      CHARACTER(LEN=*), INTENT(IN) :: line
      INTEGER :: j, ntok
      LOGICAL :: in_tok, tok_ok, tok_digit
      CHARACTER(LEN=1) :: ch
      coord_line = .FALSE.
      IF (INDEX(line, '*') > 0 .OR. INDEX(line, '=') > 0) RETURN
      ntok = 0
      in_tok = .FALSE.
      tok_ok = .TRUE.
      tok_digit = .FALSE.
      DO j = 1, LEN(line) + 1
        IF (j <= LEN(line)) THEN
          ch = line(j:j)
        ELSE
          ch = ' '
        END IF
        IF (ch == ' ' .OR. ch == CHAR(9) .OR. ch == ',') THEN
          IF (in_tok .AND. tok_ok .AND. tok_digit) ntok = ntok + 1
          in_tok = .FALSE.
          tok_ok = .TRUE.
          tok_digit = .FALSE.
        ELSE
          in_tok = .TRUE.
          IF (ch >= '0' .AND. ch <= '9') THEN
            tok_digit = .TRUE.
          ELSE IF (INDEX('.+-eEdD', ch) == 0) THEN
            tok_ok = .FALSE.
          END IF
        END IF
      END DO
      coord_line = ntok >= 3
    END FUNCTION
  END FUNCTION

  LOGICAL FUNCTION deck_has_method_sections(d)
    CHARACTER(LEN=*), INTENT(IN) :: d
    deck_has_method_sections = &
        INDEX(d, '&DFT') > 0 .OR. INDEX(d, '&dft') > 0 .OR. &
        INDEX(d, '&SYSTEM') > 0 .OR. INDEX(d, '&system') > 0 .OR. &
        INDEX(d, '&CPMD') > 0 .OR. INDEX(d, '&cpmd') > 0
  END FUNCTION

  ! Drop &ATOMS ... &END blocks so geometry can be re-appended from C arrays.
  SUBROUTINE strip_atoms_sections(src, dst, nlen)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER, INTENT(OUT) :: nlen
    INTEGER :: i, n, start_at, end_at, j
    CHARACTER(LEN=16) :: tag
    dst = ' '
    nlen = 0
    n = LEN_TRIM(src)
    IF (n < 1) RETURN
    i = 1
    DO WHILE (i <= n)
      IF (i + 5 <= n) THEN
        tag = src(i:MIN(i + 5, n))
        IF (tag == '&ATOMS' .OR. tag == '&atoms') THEN
          start_at = i
          end_at = 0
          j = i + 6
          DO WHILE (j + 3 <= n)
            IF (src(j:j+3) == '&END' .OR. src(j:j+3) == '&end') THEN
              end_at = j + 3
              EXIT
            END IF
            j = j + 1
          END DO
          IF (end_at > 0) THEN
            i = end_at + 1
            ! skip trailing newlines after &END
            DO WHILE (i <= n .AND. (src(i:i) == NEW_LINE('A') .OR. src(i:i) == ' '))
              i = i + 1
            END DO
            CYCLE
          END IF
        END IF
      END IF
      IF (nlen < LEN(dst)) THEN
        nlen = nlen + 1
        dst(nlen:nlen) = src(i:i)
      END IF
      i = i + 1
    END DO
  END SUBROUTINE

  ! Cap'n method deck (empty/missing &ATOMS) + geometry atoms from C arrays.
  SUBROUTINE embed_method_deck_plus_atoms(n_atoms, pos, z, cell, has_cell, &
      deck, nlen, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=*), INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, lmax_val, base
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=64) :: pp
    CHARACTER(LEN=8) :: lmax_c
    CHARACTER(LEN=128) :: line
    CHARACTER(LEN=4096) :: method
    INTEGER :: mlen
    ierr = 1
    deck = ' '
    nlen = 0
    CALL strip_atoms_sections(applied_input_deck, method, mlen)
    base = MIN(mlen, LEN(deck) - 64)
    IF (base < 1) RETURN
    deck(1:base) = method(1:base)
    nlen = base
    IF (deck(nlen:nlen) /= NEW_LINE('A')) THEN
      IF (nlen < LEN(deck)) THEN
        nlen = nlen + 1
        deck(nlen:nlen) = NEW_LINE('A')
      END IF
    END IF
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
      IF (INDEX(applied_input_deck, 'KLEINMAN-BYLANDER') > 0) THEN
        CALL append(deck, nlen, '*'//TRIM(pp)//' KLEINMAN-BYLANDER'//NEW_LINE('A'))
      ELSE
        CALL append(deck, nlen, '*'//TRIM(pp)//NEW_LINE('A'))
      END IF
      IF (zz == 14) THEN
        CALL append(deck, nlen, ' LMAX='//TRIM(lmax_c)//' LOC=D'//NEW_LINE('A'))
      ELSE
        CALL append(deck, nlen, ' LMAX='//TRIM(lmax_c)//NEW_LINE('A'))
      END IF
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
    ! Lattice is injected in embed_compose_cold_deck via inject_cell_if_missing.
    IF (has_cell < 0) RETURN
    IF (cell(1) < -1.0e300_c_double) RETURN
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

  ! Cap'n C render of bare scalars often omits MAXITER; OpenCPMD then uses
  ! 10000 SC steps (~minutes per force). Clamp cold force decks to a finite
  ! bound unless the wire already set one.
  SUBROUTINE inject_maxiter_if_missing(deck, nlen)
    CHARACTER(LEN=*), INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    CHARACTER(LEN=48) :: insert
    CHARACTER(LEN=16384) :: tmp
    INTEGER :: icpmd, ins_at, m, k
    IF (nlen < 1) RETURN
    IF (INDEX(deck(1:nlen), 'MAXITER') > 0 .OR. &
        INDEX(deck(1:nlen), 'maxiter') > 0) RETURN
    insert = ' MAXITER'//NEW_LINE('A')//'  40'//NEW_LINE('A')
    m = LEN_TRIM(insert)
    IF (nlen + m > LEN(deck)) RETURN
    icpmd = INDEX(deck(1:nlen), '&CPMD')
    IF (icpmd == 0) icpmd = INDEX(deck(1:nlen), '&cpmd')
    IF (icpmd <= 0) RETURN
    ins_at = icpmd
    DO WHILE (ins_at <= nlen .AND. deck(ins_at:ins_at) /= NEW_LINE('A'))
      ins_at = ins_at + 1
    END DO
    IF (ins_at <= nlen) ins_at = ins_at + 1
    tmp = deck(1:nlen)
    deck = ' '
    IF (ins_at > 1) deck(1:ins_at-1) = tmp(1:ins_at-1)
    deck(ins_at:ins_at+m-1) = insert(1:m)
    IF (ins_at <= nlen) THEN
      k = nlen - ins_at + 1
      deck(ins_at+m:ins_at+m+k-1) = tmp(ins_at:nlen)
    END IF
    nlen = nlen + m
  END SUBROUTINE

  SUBROUTINE format_cell_lines(cell, has_cell, text, ntext)
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    CHARACTER(LEN=*), INTENT(OUT) :: text
    INTEGER, INTENT(OUT) :: ntext
    REAL(real64) :: cell_a, b_over_a, c_over_a
    LOGICAL :: diagonal
    text = ' '
    ntext = 0
    cell_a = 12.0_real64
    b_over_a = 1.0_real64
    c_over_a = 1.0_real64
    diagonal = .TRUE.
    IF (has_cell /= 0) THEN
      IF (cell(1) > 0.0_c_double) cell_a = REAL(cell(1), KIND=real64)
      IF (cell(5) > 0.0_c_double .AND. cell_a > 0.0_real64) &
          b_over_a = REAL(cell(5), KIND=real64) / cell_a
      IF (cell(9) > 0.0_c_double .AND. cell_a > 0.0_real64) &
          c_over_a = REAL(cell(9), KIND=real64) / cell_a
      IF (ABS(cell(2)) + ABS(cell(3)) + ABS(cell(4)) + ABS(cell(6)) + &
          ABS(cell(7)) + ABS(cell(8)) > 1.0e-8_c_double) diagonal = .FALSE.
    END IF
    IF (.NOT. diagonal) THEN
      WRITE(text, '(A,/,3F16.8,/,3F16.8,/,3F16.8,A)') ' CELL VECTORS', &
           cell(1), cell(2), cell(3), cell(4), cell(5), cell(6), &
           cell(7), cell(8), cell(9), NEW_LINE('A')
    ELSE
      WRITE(text, '(A,3F12.6,A)') ' CELL'//NEW_LINE('A')//'  ', &
           cell_a, b_over_a, c_over_a, ' 0.0 0.0 0.0'//NEW_LINE('A')
    END IF
    ntext = LEN_TRIM(text)
  END SUBROUTINE

  ! Cap'n C render often emits &SYSTEM without CELL (lattice comes from the
  ! ForceInput box). OpenCPMD sysin stopgms if the lattice constant is zero.
  SUBROUTINE inject_cell_if_missing(deck, nlen, cell, has_cell)
    CHARACTER(LEN=*), INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    CHARACTER(LEN=400) :: insert
    CHARACTER(LEN=16384) :: tmp
    INTEGER :: isys, iang, ins_at, m, k
    IF (nlen < 1) RETURN
    ! Already has a CELL keyword (do not second-guess explicit decks).
    IF (INDEX(deck(1:nlen), 'CELL') > 0 .OR. INDEX(deck(1:nlen), 'cell') > 0) &
        RETURN
    CALL format_cell_lines(cell, has_cell, insert, m)
    IF (nlen + m > LEN(deck)) RETURN
    isys = INDEX(deck(1:nlen), '&SYSTEM')
    IF (isys == 0) isys = INDEX(deck(1:nlen), '&system')
    IF (isys <= 0) RETURN
    iang = INDEX(deck(isys:nlen), 'ANGSTROM')
    IF (iang == 0) iang = INDEX(deck(isys:nlen), 'angstrom')
    IF (iang > 0) THEN
      ins_at = isys + iang - 1
    ELSE
      ins_at = isys
    END IF
    DO WHILE (ins_at <= nlen .AND. deck(ins_at:ins_at) /= NEW_LINE('A'))
      ins_at = ins_at + 1
    END DO
    IF (ins_at <= nlen) ins_at = ins_at + 1
    tmp = deck(1:nlen)
    deck = ' '
    IF (ins_at > 1) deck(1:ins_at-1) = tmp(1:ins_at-1)
    deck(ins_at:ins_at+m-1) = insert(1:m)
    IF (ins_at <= nlen) THEN
      k = nlen - ins_at + 1
      deck(ins_at+m:ins_at+m+k-1) = tmp(ins_at:nlen)
    END IF
    nlen = nlen + m
  END SUBROUTINE

  ! Shared cold-deck assembly used by SCF and by compose preview for tests.
  SUBROUTINE embed_compose_cold_deck(n_atoms, pos, z, cell, has_cell, deck, &
      nlen, ierr)
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=*), INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    nlen = 0
    ierr = 1
    deck = ' '
    IF (LEN_TRIM(applied_input_deck) > 0) THEN
      IF (deck_has_real_atoms(applied_input_deck)) THEN
        nlen = MIN(LEN_TRIM(applied_input_deck), LEN(deck))
        deck(1:nlen) = applied_input_deck(1:nlen)
        IF (nlen < LEN(deck)) deck(nlen+1:) = ' '
        ierr = 0
      ELSE IF (deck_has_method_sections(applied_input_deck)) THEN
        CALL embed_method_deck_plus_atoms(n_atoms, pos, z, cell, has_cell, &
             deck, nlen, ierr)
      END IF
    END IF
    IF (ierr /= 0) THEN
      CALL embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, &
           ierr)
    END IF
    IF (ierr == 0 .AND. nlen > 0) THEN
      CALL inject_cell_if_missing(deck, nlen, cell, has_cell)
      CALL inject_maxiter_if_missing(deck, nlen)
    END IF
  END SUBROUTINE

  LOGICAL FUNCTION warm_cell_matches(image, cell, has_cell)
    TYPE(cpmdc_embed_image), INTENT(IN) :: image
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER :: i
    warm_cell_matches = .FALSE.
    IF (image%warm_cell_set == 0) RETURN
    IF (has_cell /= image%warm_has_cell) RETURN
    IF (has_cell == 0) THEN
      warm_cell_matches = .TRUE.
      RETURN
    END IF
    warm_cell_matches = .TRUE.
    DO i = 1, 9
      IF (ABS(cell(i) - image%warm_cell(i)) > 1.0e-8_c_double) &
          warm_cell_matches = .FALSE.
    END DO
  END FUNCTION

  SUBROUTINE latch_warm_cell(image, cell, has_cell)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER :: i
    image%warm_has_cell = has_cell
    image%warm_cell_set = 1
    IF (has_cell == 0) RETURN
    DO i = 1, 9
      image%warm_cell(i) = cell(i)
    END DO
  END SUBROUTINE

  SUBROUTINE run_embed_scf(image, n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    USE rwfopt_utils, ONLY: cpmdc_reset_warm_orbitals
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
    USE system, ONLY: cnts, cntl
    USE isos, ONLY: isos1
    USE ropt, ONLY: init_pinf_pointers
    USE bicanonicalCpmd, ONLY: bicanonicalCpmdConfig, bicanonicalCpmdInputConfig, New
    USE bicanonicalConfig, ONLY: New
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
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
      FUNCTION cpmdc_prepare_pp_cwd(pseudo_dir) BIND(C, NAME='cpmdc_prepare_pp_cwd')
        IMPORT :: c_char, c_int
        CHARACTER(KIND=c_char), INTENT(IN) :: pseudo_dir(*)
        INTEGER(c_int) :: cpmdc_prepare_pp_cwd
      END FUNCTION
      FUNCTION cpmdc_restore_host_cwd() BIND(C, NAME='cpmdc_restore_host_cwd')
        IMPORT :: c_int
        INTEGER(c_int) :: cpmdc_restore_host_cwd
      END FUNCTION
    END INTERFACE
    ok = 0_c_int
    energy_h = 0.0_c_double
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    IF (image%cfg_warm_steps > 0 .AND. warm_cell_matches(image, cell, has_cell)) THEN
      CALL embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, ok)
      IF (ok /= 0_c_int) image%cfg_warm_steps = image%cfg_warm_steps + 1
      RETURN
    END IF
    IF (image%cfg_warm_steps > 0) THEN
      image%cfg_warm_steps = 0
      CALL cpmdc_reset_warm_orbitals()
    END IF
    ! Cold: honor Cap'n-rendered applied_input_deck for method sections.
    ! 1) Deck with real &ATOMS PP lines (*...) → use as-is.
    ! 2) Method sections with empty/missing &ATOMS (C render placeholder) →
    !    keep method text and append geometry &ATOMS from C arrays.
    ! 3) Else minimal deck with applied functional/cutoff/charge/mult.
    ! Geometry for forces always from C arrays into TAU0 after parse.
    CALL embed_compose_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, &
         ierr)
    IF (ierr /= 0 .OR. nlen < 1) RETURN
    mfd = INT(cpmdc_memfd_write(deck, INT(nlen, KIND=c_int), mempath, &
         INT(LEN(mempath), KIND=c_int)))
    IF (mfd < 0) RETURN
    paral%io_parent = .TRUE.
    ! (compose path above honors method-only Cap'n decks.)
    ! C wrote a NUL-terminated /proc/self/fd/N path into mempath.
    DO idx = 1, LEN(mempath)
      IF (IACHAR(mempath(idx:idx)) == 0) THEN
        cnts%inputfile = mempath(1:idx-1)
        EXIT
      END IF
    END DO
    ! OpenCPMD get_pplib uses argv[2] as the PP library whenever argc>1,
    ! ignoring CPMD_PP_LIBRARY_PATH. Hosts like Catch2/eOn always pass filters
    ! so argc>1; relative *PP basenames then only resolve via recpnew's
    ! second-chance CWD lookup. cpmdc_prepare_pp_cwd chdirs to the library
    ! and also exports CPMD_PP_LIBRARY_PATH (trailing slash) for argc==1 hosts.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_PSEUDO_DIR', pp_dir)
    IF (LEN_TRIM(pp_dir) == 0) &
        CALL GET_ENVIRONMENT_VARIABLE('CPMD_PP_LIBRARY_PATH', pp_dir)
    IF (LEN_TRIM(pp_dir) == 0) RETURN
    IF (cpmdc_prepare_pp_cwd(TRIM(pp_dir)//c_null_char) /= 0_c_int) RETURN
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
    ! PEF stress needs cntl%tpres before dqgalloc. Isolated/Hockney (tclust)
    ! must not set tpres: rinitwf→newcell→gf_periodic needs scg, which
    ! initclust only allocates for periodic cells.
    IF (.NOT. isos1%tclust .AND. embed_stress_wanted()) cntl%tpres = .TRUE.
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
    CALL embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, ok)
    IF (ok /= 0_c_int) THEN
      image%cfg_warm_steps = 1
      CALL latch_warm_cell(image, cell, has_cell)
    ELSE
      CALL clear_last_energy_components(image)
    END IF
    ! Hand CWD back to the host (eOn workdir); PP loads already finished.
    IF (cpmdc_restore_host_cwd() /= 0_c_int) THEN
      ! Keep ok from SCF; lost host CWD is non-fatal for the energy itself.
    END IF
  END SUBROUTINE

  FUNCTION cpmdc_embed_compose_cold_deck(n_atoms, positions_ang, atomic_numbers, &
      cell_ang, has_cell, deck_out, deck_cap, deck_len) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_compose_cold_deck')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    CHARACTER(KIND=c_char), INTENT(OUT) :: deck_out(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_cap
    INTEGER(c_int), INTENT(OUT) :: deck_len
    INTEGER(c_int) :: ok
    CHARACTER(LEN=16384) :: deck
    INTEGER :: nlen, ierr, i, ncopy
    ok = 0_c_int
    deck_len = 0_c_int
    IF (deck_cap < 2) RETURN
    CALL embed_compose_cold_deck(INT(n_atoms), positions_ang, atomic_numbers, &
         cell_ang, INT(has_cell), deck, nlen, ierr)
    IF (ierr /= 0 .OR. nlen < 1) RETURN
    ncopy = MIN(nlen, INT(deck_cap) - 1)
    DO i = 1, ncopy
      deck_out(i) = deck(i:i)
    END DO
    deck_out(ncopy + 1) = c_null_char
    deck_len = INT(ncopy, KIND=c_int)
    ok = 1_c_int
  END FUNCTION
#endif
END MODULE cpmd_embed_c_api
