! SPDX-License-Identifier: MIT
!
! cpmd_embed_c_api.F90 — compiler-independent C ABI for OpenCPMD embed
! (nwchemc pattern: bind(C) names; engine CALLs live here / in legacy helpers).
!
#include "cpmd_embed_config.h"

MODULE cpmdc_embed_host_iface
  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_int, c_double, c_size_t
  IMPLICIT NONE
  INTERFACE
    FUNCTION cpmdc_species_order_map(n_atoms, atomic_numbers, n_species, &
        species_z, species_count, map_out) BIND(C, NAME='cpmdc_species_order_map')
      IMPORT :: c_int
      INTEGER(c_int), INTENT(IN), VALUE :: n_atoms, n_species
      INTEGER(c_int), INTENT(IN) :: atomic_numbers(*), species_z(*), species_count(*)
      INTEGER(c_int), INTENT(OUT) :: map_out(*)
      INTEGER(c_int) :: cpmdc_species_order_map
    END FUNCTION
    SUBROUTINE cpmdc_scatter_species_gradient(n_atoms, map, species_grad, grad) &
        BIND(C, NAME='cpmdc_scatter_species_gradient')
      IMPORT :: c_int, c_double
      INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
      INTEGER(c_int), INTENT(IN) :: map(*)
      REAL(c_double), INTENT(IN) :: species_grad(*)
      REAL(c_double), INTENT(OUT) :: grad(*)
    END SUBROUTINE
    SUBROUTINE cpmdc_note_embed_failure(msg) BIND(C, NAME='cpmdc_note_embed_failure')
      IMPORT :: c_char
      CHARACTER(KIND=c_char), INTENT(IN) :: msg(*)
    END SUBROUTINE
    FUNCTION cpmdc_enter_output_cwd(output_dir) BIND(C, NAME='cpmdc_enter_output_cwd')
      IMPORT :: c_char, c_int
      CHARACTER(KIND=c_char), INTENT(IN) :: output_dir(*)
      INTEGER(c_int) :: cpmdc_enter_output_cwd
    END FUNCTION
    FUNCTION cpmdc_pseudopotential_directory(buf, cap) &
        BIND(C, NAME='cpmdc_pseudopotential_directory')
      IMPORT :: c_char, c_int, c_size_t
      CHARACTER(KIND=c_char), INTENT(OUT) :: buf(*)
      INTEGER(c_size_t), INTENT(IN), VALUE :: cap
      INTEGER(c_int) :: cpmdc_pseudopotential_directory
    END FUNCTION
    FUNCTION cpmdc_memfd_write(bytes, nbytes, path_out, path_cap) &
        BIND(C, NAME='cpmdc_memfd_write')
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
END MODULE

MODULE cpmd_embed_c_api
  USE, INTRINSIC :: iso_c_binding
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64, error_unit
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpmdc_embed_init, cpmdc_embed_available, cpmdc_embed_finalize
  PUBLIC :: cpmdc_embed_bind_calculator
  PUBLIC :: cpmdc_embed_reset_state
  PUBLIC :: cpmdc_embed_set_config, cpmdc_embed_set_deck, cpmdc_embed_energy_grad
  PUBLIC :: cpmdc_embed_compose_cold_deck

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  ! Results, the warm cell, and the method knobs live in the caller image.
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
    INTEGER(c_int) :: cfg_set
    CHARACTER(KIND=c_char) :: functional(64)
    REAL(c_double) :: cutoff_ry
    INTEGER(c_int) :: cfg_charge
    INTEGER(c_int) :: multiplicity
    TYPE(c_ptr) :: input_deck
    INTEGER(c_int) :: input_deck_len
    CHARACTER(KIND=c_char) :: cpmd_root(1024)
    CHARACTER(KIND=c_char) :: output_dir(1024)
  END TYPE
#if defined(CPMDC_HAS_CPMD)
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double
  ! mp_start assigns mp_comm_world only when it calls MPI_Init.
  ! The flag records that this process already split a calculator.
  LOGICAL, SAVE :: embed_calculator_bound = .FALSE.
  ! The saved c0 matches one cutoff, cell, charge, multiplicity,
  ! functional, deck, and elemental composition. A later call with a
  ! different basis clears it before rwfopt can restore it.
  LOGICAL, SAVE :: embed_basis_latched = .FALSE.
  REAL(c_double), SAVE :: embed_saved_cutoff = -1.0_c_double
  INTEGER, SAVE :: embed_saved_charge = -999
  INTEGER, SAVE :: embed_saved_mult = -999
  INTEGER, SAVE :: embed_saved_natoms = -1
  INTEGER, ALLOCATABLE, SAVE :: embed_saved_z(:)
  CHARACTER(LEN=64), SAVE :: embed_saved_functional = ''
  CHARACTER(LEN=:), ALLOCATABLE, SAVE :: embed_saved_deck
  INTEGER, SAVE :: embed_saved_has_cell = 0
  INTEGER, SAVE :: embed_saved_cell_set = 0
  REAL(c_double), SAVE :: embed_saved_cell(9) = 0.0_c_double
#endif

  TYPE :: embed_knobs
    CHARACTER(LEN=64) :: functional = 'BLYP'
    REAL(real64) :: cutoff_ry = 70.0_real64
    INTEGER :: charge = 0
    INTEGER :: mult = 1
    CHARACTER(LEN=:), ALLOCATABLE :: input_deck
    CHARACTER(LEN=1024) :: cpmd_root = ' '
  END TYPE

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
      USE mp_interface, ONLY: mp_comm_world
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
      IF (embed_calculator_bound) RETURN
      color = rank / rpc
      key = MOD(rank, rpc)
      CALL MPI_Comm_split(MPI_COMM_WORLD, color, key, comm, ierr)
      IF (ierr /= 0) THEN
        calc = -1_c_int
        RETURN
      END IF
      mp_comm_world = comm
      embed_calculator_bound = .TRUE.
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
    ok = 0_c_int
#if defined(CPMDC_HAS_CPMD)
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
#endif
  END FUNCTION

  FUNCTION cpmdc_embed_reset_state(image_c) RESULT(ok) BIND(C, NAME='cpmdc_embed_reset_state')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_reset_warm_orbitals
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
    CALL embed_reset_warm_orbitals()
    CALL clear_embed_basis_latch()
#endif
    ok = 1_c_int
  END FUNCTION

  ! A new image has no warm counter of its own. The process copy of c0
  ! stays until a force call sees a different basis.
  FUNCTION cpmdc_embed_detach_image(image_c) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_detach_image')
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
    image%cfg_warm_steps = 0_c_int
    image%warm_cell_set = 0_c_int
    image%warm_has_cell = 0_c_int
    image%warm_cell = 0.0_c_double
    ok = 1_c_int
  END FUNCTION

  SUBROUTINE cpmdc_embed_finalize() BIND(C, NAME='cpmdc_embed_finalize')
    runtime_ready = .FALSE.
    runtime_finalized = .TRUE.
  END SUBROUTINE

  SUBROUTINE copy_f_to_cchars(src, dst, n)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(KIND=c_char), INTENT(OUT) :: dst(*)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, m
    m = MIN(LEN_TRIM(src), n)
    DO i = 1, m
      dst(i) = src(i:i)
    END DO
    IF (m < n) dst(m + 1) = c_null_char
  END SUBROUTINE

  SUBROUTINE copy_cchars_to_f(src, n, dst)
    CHARACTER(KIND=c_char), INTENT(IN) :: src(*)
    INTEGER, INTENT(IN) :: n
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER :: i, m
    dst = ' '
    m = MIN(n, LEN(dst))
    DO i = 1, m
      IF (src(i) == c_null_char) EXIT
      dst(i:i) = src(i)
    END DO
  END SUBROUTINE

  FUNCTION knobs_of(image) RESULT(k)
    TYPE(cpmdc_embed_image), INTENT(IN) :: image
    TYPE(embed_knobs) :: k
    k%input_deck = ''
    IF (image%cfg_set == 0_c_int) RETURN
    CALL copy_cchars_to_f(image%functional, 64, k%functional)
    IF (LEN_TRIM(k%functional) == 0) k%functional = 'BLYP'
    k%cutoff_ry = REAL(image%cutoff_ry, KIND=real64)
    IF (k%cutoff_ry <= 0.0_real64) k%cutoff_ry = 70.0_real64
    k%charge = INT(image%cfg_charge)
    k%mult = MAX(1, INT(image%multiplicity))
    IF (C_ASSOCIATED(image%input_deck) .AND. image%input_deck_len > 0) &
        CALL cptr_to_alloc(image%input_deck, INT(image%input_deck_len), k%input_deck)
    CALL copy_cchars_to_f(image%cpmd_root, 1024, k%cpmd_root)
  END FUNCTION

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
    TYPE(cpmdc_embed_image), POINTER :: image
    CHARACTER(LEN=64) :: functional_l
    CHARACTER(LEN=1024) :: root_l
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (functional_len < 0 .OR. input_deck_len < 0 .OR. cpmd_root_len < 0) RETURN
    IF (cutoff_ry < 0.0_c_double) RETURN
    IF (cpmdc_embed_reset_state(image_c) == 0_c_int) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
    CALL cstr_to_f(functional, functional_len, functional_l)
    IF (LEN_TRIM(functional_l) == 0) functional_l = 'BLYP'
    CALL cstr_to_f(cpmd_root, cpmd_root_len, root_l)
    CALL copy_f_to_cchars(functional_l, image%functional, 64)
    image%cutoff_ry = REAL(cutoff_ry, KIND=c_double)
    IF (image%cutoff_ry <= 0.0_c_double) image%cutoff_ry = 70.0_c_double
    image%cfg_charge = INT(charge, KIND=c_int)
    image%multiplicity = MAX(1_c_int, INT(multiplicity, KIND=c_int))
    IF (store_image_deck(image_c, input_deck, input_deck_len) /= 0_c_int) THEN
      image%cfg_set = 0_c_int
      RETURN
    END IF
    CALL copy_f_to_cchars(root_l, image%cpmd_root, 1024)
    image%cfg_set = 1_c_int
    ok = 1_c_int
  END FUNCTION


  FUNCTION cpmdc_embed_set_deck(deck, deck_len, image_c) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_set_deck')
    CHARACTER(KIND=c_char), INTENT(IN) :: deck(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_len
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized .OR. deck_len < 0) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    IF (store_image_deck(image_c, deck, deck_len) /= 0_c_int) RETURN
    CALL C_F_POINTER(image_c, image)
    image%cfg_set = 1_c_int
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
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(n3) / 3) RETURN
    n3 = INT(n_atoms) * 3
    DO i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    END DO
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
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

  FUNCTION store_image_deck(image_c, text, text_len) RESULT(rc)
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    CHARACTER(KIND=c_char), INTENT(IN) :: text(*)
    INTEGER(c_int), INTENT(IN), VALUE :: text_len
    INTEGER(c_int) :: rc
    INTERFACE
      FUNCTION cpmdc_embed_image_store_deck(image, text, text_len) &
          BIND(C, NAME='cpmdc_embed_image_store_deck')
        IMPORT :: c_ptr, c_char, c_int
        TYPE(c_ptr), VALUE :: image
        CHARACTER(KIND=c_char), INTENT(IN) :: text(*)
        INTEGER(c_int), VALUE :: text_len
        INTEGER(c_int) :: cpmdc_embed_image_store_deck
      END FUNCTION
    END INTERFACE
    rc = cpmdc_embed_image_store_deck(image_c, text, text_len)
  END FUNCTION

  SUBROUTINE cptr_to_alloc(src, n, dst)
    TYPE(c_ptr), INTENT(IN) :: src
    INTEGER, INTENT(IN) :: n
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: dst
    CHARACTER(KIND=c_char), POINTER :: bytes(:)
    INTEGER :: i, stat
    IF (n <= 0 .OR. .NOT. C_ASSOCIATED(src)) THEN
      dst = ''
      RETURN
    END IF
    ALLOCATE(CHARACTER(LEN=n) :: dst, STAT=stat)
    IF (stat /= 0) THEN
      dst = ''
      RETURN
    END IF
    CALL C_F_POINTER(src, bytes, [n])
    DO i = 1, n
      dst(i:i) = bytes(i)
    END DO
  END SUBROUTINE

  ! Grow buf so the next m characters fit. A short buf is an error, not a cut.
  SUBROUTINE append_grow(buf, n, s, ierr)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: buf
    INTEGER, INTENT(INOUT) :: n
    CHARACTER(LEN=*), INTENT(IN) :: s
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=:), ALLOCATABLE :: grown
    INTEGER :: m, need, stat
    ierr = 1
    m = LEN(s)
    IF (n < 0 .OR. m < 0) RETURN
    IF (m > HUGE(need) - MAX(n, 0)) RETURN
    need = n + m
    IF (.NOT. ALLOCATED(buf) .OR. need > LEN(buf)) THEN
      ALLOCATE(CHARACTER(LEN=need) :: grown, STAT=stat)
      IF (stat /= 0) RETURN
      IF (n > 0 .AND. ALLOCATED(buf)) grown(1:n) = buf(1:n)
      CALL MOVE_ALLOC(grown, buf)
    END IF
    IF (m > 0) buf(n + 1:n + m) = s(1:m)
    n = need
    ierr = 0
  END SUBROUTINE

  ! Copy src without &ATOMS blocks. dst is sized from src, so nothing is dropped.
  SUBROUTINE strip_atoms_alloc(src, dst, nlen, ierr)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: dst
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, n, end_at, j, stat
    CHARACTER(LEN=16) :: tag
    ierr = 0
    nlen = 0
    n = LEN_TRIM(src)
    IF (n < 1) THEN
      dst = ''
      RETURN
    END IF
    ALLOCATE(CHARACTER(LEN=n) :: dst, STAT=stat)
    IF (stat /= 0) THEN
      dst = ''
      ierr = 1
      RETURN
    END IF
    i = 1
    DO WHILE (i <= n)
      IF (i + 5 <= n) THEN
        tag = src(i:MIN(i + 5, n))
        IF (tag == '&ATOMS' .OR. tag == '&atoms') THEN
          end_at = 0
          j = i + 6
          DO WHILE (j + 3 <= n)
            IF (src(j:j + 3) == '&END' .OR. src(j:j + 3) == '&end') THEN
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
      IF (nlen >= LEN(dst)) THEN
        ierr = 1
        RETURN
      END IF
      nlen = nlen + 1
      dst(nlen:nlen) = src(i:i)
      i = i + 1
    END DO
  END SUBROUTINE


  ! Rendered &ATOMS stubs tag each pseudopotential with '!SPECIES SYM'.
  ! KLEINMAN-BYLANDER is read from that species' *file line only.
  FUNCTION embed_z_of_symbol(sym) RESULT(zz)
    CHARACTER(LEN=*), INTENT(IN) :: sym
    INTEGER :: zz
    CHARACTER(LEN=*), PARAMETER :: k_sym = &
      'H ' // 'HE' // 'LI' // 'BE' // 'B ' // 'C ' // 'N ' // 'O ' // 'F ' // 'NE' // &
      'NA' // 'MG' // 'AL' // 'SI' // 'P ' // 'S ' // 'CL' // 'AR' // 'K ' // 'CA' // &
      'SC' // 'TI' // 'V ' // 'CR' // 'MN' // 'FE' // 'CO' // 'NI' // 'CU' // 'ZN' // &
      'GA' // 'GE' // 'AS' // 'SE' // 'BR' // 'KR' // 'RB' // 'SR' // 'Y ' // 'ZR' // &
      'NB' // 'MO' // 'TC' // 'RU' // 'RH' // 'PD' // 'AG' // 'CD' // 'IN' // 'SN' // &
      'SB' // 'TE' // 'I ' // 'XE' // 'CS' // 'BA' // 'LA' // 'CE' // 'PR' // 'ND' // &
      'PM' // 'SM' // 'EU' // 'GD' // 'TB' // 'DY' // 'HO' // 'ER' // 'TM' // 'YB' // &
      'LU' // 'HF' // 'TA' // 'W ' // 'RE' // 'OS' // 'IR' // 'PT' // 'AU' // 'HG' // &
      'TL' // 'PB' // 'BI' // 'PO' // 'AT' // 'RN' // 'FR' // 'RA' // 'AC' // 'TH' // &
      'PA' // 'U ' // 'NP' // 'PU' // 'AM' // 'CM' // 'BK' // 'CF' // 'ES' // 'FM' // &
      'MD' // 'NO' // 'LR' // 'RF' // 'DB' // 'SG' // 'BH' // 'HS' // 'MT' // 'DS' // &
      'RG' // 'CN' // 'NH' // 'FL' // 'MC' // 'LV' // 'TS' // 'OG'
    CHARACTER(LEN=2) :: key
    INTEGER :: i, n
    CHARACTER(LEN=1) :: ch
    zz = -1
    key = ' '
    n = 0
    DO i = 1, LEN_TRIM(sym)
      ch = sym(i:i)
      IF (ch == ' ' .OR. ch == ACHAR(9)) EXIT
      IF (ch >= 'a' .AND. ch <= 'z') ch = ACHAR(IACHAR(ch) - 32)
      IF (ch < 'A' .OR. ch > 'Z') RETURN
      n = n + 1
      IF (n > 2) RETURN
      key(n:n) = ch
    END DO
    IF (n < 1) RETURN
    DO i = 1, LEN(k_sym) / 2
      IF (k_sym(2 * i - 1:2 * i) == key) THEN
        zz = i
        RETURN
      END IF
    END DO
  END FUNCTION

  SUBROUTINE embed_ltrim(src, dst)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER :: i, n, k
    dst = ' '
    n = LEN_TRIM(src)
    i = 1
    DO WHILE (i <= n .AND. (src(i:i) == ' ' .OR. src(i:i) == ACHAR(9)))
      i = i + 1
    END DO
    k = n - i + 1
    IF (k < 1) RETURN
    IF (k > LEN(dst)) k = LEN(dst)
    dst(1:k) = src(i:i + k - 1)
  END SUBROUTINE

  SUBROUTINE embed_upper_copy(src, dst)
    CHARACTER(LEN=*), INTENT(IN) :: src
    CHARACTER(LEN=*), INTENT(OUT) :: dst
    INTEGER :: i, n
    CHARACTER(LEN=1) :: ch
    dst = ' '
    n = MIN(LEN_TRIM(src), LEN(dst))
    DO i = 1, n
      ch = src(i:i)
      IF (ch >= 'a' .AND. ch <= 'z') ch = ACHAR(IACHAR(ch) - 32)
      dst(i:i) = ch
    END DO
  END SUBROUTINE

  SUBROUTINE embed_next_tok(line, pos, tok)
    CHARACTER(LEN=*), INTENT(IN) :: line
    INTEGER, INTENT(INOUT) :: pos
    CHARACTER(LEN=*), INTENT(OUT) :: tok
    INTEGER :: n, i, k
    tok = ' '
    n = LEN_TRIM(line)
    DO WHILE (pos <= n .AND. (line(pos:pos) == ' ' .OR. line(pos:pos) == ACHAR(9)))
      pos = pos + 1
    END DO
    IF (pos > n) RETURN
    i = pos
    DO WHILE (i <= n .AND. line(i:i) /= ' ' .AND. line(i:i) /= ACHAR(9))
      i = i + 1
    END DO
    k = MIN(i - pos, LEN(tok))
    IF (k > 0) tok(1:k) = line(pos:pos + k - 1)
    pos = i
  END SUBROUTINE

  SUBROUTINE embed_take_line(src, zend, pos, line)
    CHARACTER(LEN=*), INTENT(IN) :: src
    INTEGER, INTENT(IN) :: zend
    INTEGER, INTENT(INOUT) :: pos
    CHARACTER(LEN=*), INTENT(OUT) :: line
    INTEGER :: i, k
    line = ' '
    IF (pos < 1) pos = 1
    IF (pos > zend) RETURN
    i = pos
    DO WHILE (i <= zend .AND. src(i:i) /= NEW_LINE('A') .AND. src(i:i) /= ACHAR(13))
      i = i + 1
    END DO
    k = MIN(i - pos, LEN(line))
    IF (k > 0) line(1:k) = src(pos:pos + k - 1)
    IF (i <= zend .AND. src(i:i) == ACHAR(13)) i = i + 1
    IF (i <= zend .AND. src(i:i) == NEW_LINE('A')) i = i + 1
    pos = i
  END SUBROUTINE

  LOGICAL FUNCTION embed_is_int_token(tok)
    CHARACTER(LEN=*), INTENT(IN) :: tok
    INTEGER :: i, n, a
    embed_is_int_token = .FALSE.
    n = LEN_TRIM(tok)
    IF (n < 1) RETURN
    a = 1
    IF (tok(1:1) == '+' .OR. tok(1:1) == '-') THEN
      IF (n == 1) RETURN
      a = 2
    END IF
    DO i = a, n
      IF (tok(i:i) < '0' .OR. tok(i:i) > '9') RETURN
    END DO
    embed_is_int_token = .TRUE.
  END FUNCTION

  LOGICAL FUNCTION embed_is_count_line(work)
    CHARACTER(LEN=*), INTENT(IN) :: work
    INTEGER :: tpos, ntok
    CHARACTER(LEN=64) :: tok
    embed_is_count_line = .FALSE.
    tpos = 1
    ntok = 0
    DO
      CALL embed_next_tok(work, tpos, tok)
      IF (LEN_TRIM(tok) == 0) EXIT
      ntok = ntok + 1
      IF (ntok > 1) RETURN
      IF (.NOT. embed_is_int_token(tok)) RETURN
    END DO
    embed_is_count_line = ntok == 1
  END FUNCTION

  LOGICAL FUNCTION embed_is_species_tag(up)
    CHARACTER(LEN=*), INTENT(IN) :: up
    embed_is_species_tag = .FALSE.
    IF (LEN_TRIM(up) < 8) RETURN
    IF (up(1:8) /= '!SPECIES') RETURN
    IF (LEN_TRIM(up) == 8) THEN
      embed_is_species_tag = .TRUE.
      RETURN
    END IF
    embed_is_species_tag = up(9:9) == ' ' .OR. up(9:9) == ACHAR(9)
  END FUNCTION

  SUBROUTINE embed_absorb_option_line(work, lmaxc, locc, extra)
    CHARACTER(LEN=*), INTENT(IN) :: work
    CHARACTER(LEN=*), INTENT(INOUT) :: lmaxc, locc, extra
    INTEGER :: tpos
    CHARACTER(LEN=192) :: tok
    CHARACTER(LEN=192) :: up
    tpos = 1
    DO
      CALL embed_next_tok(work, tpos, tok)
      IF (LEN_TRIM(tok) == 0) EXIT
      CALL embed_upper_copy(tok, up)
      IF (LEN_TRIM(up) >= 6 .AND. up(1:5) == 'LMAX=') THEN
        lmaxc = up(6:6)
      ELSE IF (LEN_TRIM(up) >= 5 .AND. up(1:4) == 'LOC=') THEN
        locc = up(5:5)
      ELSE IF (TRIM(up) /= 'KLEINMAN-BYLANDER') THEN
        IF (LEN_TRIM(extra) == 0) THEN
          extra = TRIM(tok)
        ELSE
          extra = TRIM(extra) // ' ' // TRIM(tok)
        END IF
      END IF
    END DO
  END SUBROUTINE

  SUBROUTINE embed_parse_star(work, path, kb)
    CHARACTER(LEN=*), INTENT(IN) :: work
    CHARACTER(LEN=*), INTENT(OUT) :: path
    LOGICAL, INTENT(OUT) :: kb
    INTEGER :: tpos
    CHARACTER(LEN=512) :: tok
    CHARACTER(LEN=512) :: up
    path = ' '
    kb = .FALSE.
    tpos = 2
    CALL embed_next_tok(work, tpos, path)
    DO
      CALL embed_next_tok(work, tpos, tok)
      IF (LEN_TRIM(tok) == 0) EXIT
      CALL embed_upper_copy(tok, up)
      IF (TRIM(up) == 'KLEINMAN-BYLANDER') kb = .TRUE.
    END DO
  END SUBROUTINE

  SUBROUTINE embed_parse_following_options(src, zend, pos, lmaxc, locc, extra)
    CHARACTER(LEN=*), INTENT(IN) :: src
    INTEGER, INTENT(IN) :: zend
    INTEGER, INTENT(INOUT) :: pos
    CHARACTER(LEN=*), INTENT(INOUT) :: lmaxc, locc, extra
    INTEGER :: save
    CHARACTER(LEN=512) :: line, work, up
    DO WHILE (pos <= zend)
      save = pos
      CALL embed_take_line(src, zend, pos, line)
      CALL embed_ltrim(line, work)
      IF (LEN_TRIM(work) == 0) CYCLE
      IF (work(1:1) == '*' .OR. work(1:1) == '&') THEN
        pos = save
        RETURN
      END IF
      CALL embed_upper_copy(work, up)
      IF (embed_is_species_tag(up)) THEN
        pos = save
        RETURN
      END IF
      IF (embed_is_count_line(work)) RETURN
      CALL embed_absorb_option_line(work, lmaxc, locc, extra)
    END DO
  END SUBROUTINE

  SUBROUTINE embed_collect_message_pp(src, have, kbflag, paths, lmaxs, locs, extras)
    CHARACTER(LEN=*), INTENT(IN) :: src
    LOGICAL, INTENT(OUT) :: have(0:120), kbflag(0:120)
    CHARACTER(LEN=512), INTENT(OUT) :: paths(0:120)
    CHARACTER(LEN=8), INTENT(OUT) :: lmaxs(0:120), locs(0:120)
    CHARACTER(LEN=192), INTENT(OUT) :: extras(0:120)
    INTEGER :: nsrc, pos, rel, ia, iend, zend, pending, tpos
    CHARACTER(LEN=512) :: line, work, up, path
    CHARACTER(LEN=64) :: sym
    CHARACTER(LEN=8) :: lmaxc, locc
    CHARACTER(LEN=192) :: extra
    LOGICAL :: kb
    have = .FALSE.
    kbflag = .FALSE.
    paths = ' '
    lmaxs = ' '
    locs = ' '
    extras = ' '
    nsrc = LEN_TRIM(src)
    pos = 1
    DO WHILE (pos <= nsrc)
      rel = INDEX(src(pos:nsrc), '&ATOMS')
      IF (rel == 0) rel = INDEX(src(pos:nsrc), '&atoms')
      IF (rel <= 0) EXIT
      ia = pos + rel - 1
      iend = INDEX(src(ia:nsrc), '&END')
      IF (iend == 0) iend = INDEX(src(ia:nsrc), '&end')
      IF (iend <= 0) EXIT
      zend = ia + iend - 2
      pos = ia
      pending = -1
      DO WHILE (pos <= zend)
        CALL embed_take_line(src, zend, pos, line)
        CALL embed_ltrim(line, work)
        IF (LEN_TRIM(work) == 0) CYCLE
        CALL embed_upper_copy(work, up)
        IF (embed_is_species_tag(up)) THEN
          tpos = 9
          CALL embed_next_tok(work, tpos, sym)
          pending = embed_z_of_symbol(sym)
        ELSE IF (work(1:1) == '*') THEN
          CALL embed_parse_star(work, path, kb)
          lmaxc = 'S'
          locc = ' '
          extra = ' '
          CALL embed_parse_following_options(src, zend, pos, lmaxc, locc, extra)
          IF (pending >= 1 .AND. pending <= 120 .AND. LEN_TRIM(path) > 0) THEN
            have(pending) = .TRUE.
            kbflag(pending) = kb
            paths(pending) = path
            lmaxs(pending) = lmaxc
            locs(pending) = locc
            extras(pending) = extra
          END IF
          pending = -1
        END IF
      END DO
      pos = ia + iend + 3
    END DO
  END SUBROUTINE

  SUBROUTINE embed_table_pp(zz, pp, lmax_val, ok)
    INTEGER, INTENT(IN) :: zz
    CHARACTER(LEN=*), INTENT(OUT) :: pp
    INTEGER, INTENT(OUT) :: lmax_val, ok
    ok = 0
    pp = ' '
    lmax_val = -1
    IF (zz == 1) THEN
      pp = 'H_CVB_BLYP.psp'
      lmax_val = 0
      ok = 1
    ELSE IF (zz == 6) THEN
      pp = 'C_MT_BLYP.psp'
      lmax_val = 1
      ok = 1
    ELSE IF (zz == 7) THEN
      pp = 'N_MT_BLYP.psp'
      lmax_val = 1
      ok = 1
    ELSE IF (zz == 8) THEN
      pp = 'O_MT_BLYP.psp'
      lmax_val = 1
      ok = 1
    ELSE IF (zz == 14) THEN
      pp = 'Si_MT_BLYP.psp'
      lmax_val = 2
      ok = 1
    ELSE IF (zz == 32) THEN
      pp = 'Ge_MT_BLYP.psp'
      lmax_val = 1
      ok = 1
    END IF
  END SUBROUTINE

  SUBROUTINE embed_pp_lines_for_z(src, zz, star, opt, ok)
    CHARACTER(LEN=*), INTENT(IN) :: src
    INTEGER, INTENT(IN) :: zz
    CHARACTER(LEN=*), INTENT(OUT) :: star, opt
    INTEGER, INTENT(OUT) :: ok
    LOGICAL, ALLOCATABLE :: have(:), kbflag(:)
    CHARACTER(LEN=512), ALLOCATABLE :: paths(:)
    CHARACTER(LEN=8), ALLOCATABLE :: lmaxs(:), locs(:)
    CHARACTER(LEN=192), ALLOCATABLE :: extras(:)
    CHARACTER(LEN=64) :: tpath
    CHARACTER(LEN=1) :: lc
    INTEGER :: astat, lmax_val, pok
    ok = 0
    star = ' '
    opt = ' '
    IF (zz < 1 .OR. zz > 120) THEN
      WRITE(error_unit, '(A,I0,A)') 'cpmdc: atomic number ', zz, ' is outside 1..120'
      FLUSH(error_unit)
      RETURN
    END IF
    ALLOCATE(have(0:120), kbflag(0:120), paths(0:120), lmaxs(0:120), locs(0:120), &
         extras(0:120), STAT=astat)
    IF (astat /= 0) THEN
      WRITE(error_unit, '(A)') 'cpmdc: failed to allocate pseudopotential table'
      FLUSH(error_unit)
      RETURN
    END IF
    CALL embed_collect_message_pp(src, have, kbflag, paths, lmaxs, locs, extras)
    IF (have(zz)) THEN
      star = '*' // TRIM(paths(zz))
      IF (kbflag(zz)) star = TRIM(star) // ' KLEINMAN-BYLANDER'
      opt = ' LMAX=' // TRIM(lmaxs(zz))
      IF (LEN_TRIM(locs(zz)) > 0) opt = TRIM(opt) // ' LOC=' // TRIM(locs(zz))
      IF (LEN_TRIM(extras(zz)) > 0) opt = TRIM(opt) // ' ' // TRIM(extras(zz))
      ok = 1
      RETURN
    END IF
    CALL embed_table_pp(zz, tpath, lmax_val, pok)
    IF (pok == 0) THEN
      WRITE(error_unit, '(A,I0,A)') &
           'cpmdc: no pseudopotential for atomic number ', zz, &
           ' (absent from the message and from the built-in table H C N O Si Ge)'
      FLUSH(error_unit)
      RETURN
    END IF
    IF (lmax_val <= 0) THEN
      lc = 'S'
    ELSE IF (lmax_val == 1) THEN
      lc = 'P'
    ELSE
      lc = 'D'
    END IF
    star = '*' // TRIM(tpath)
    opt = ' LMAX=' // lc
    IF (zz == 14) opt = TRIM(opt) // ' LOC=D'
    ok = 1
  END SUBROUTINE

  ! One cold-deck procedure for the stub and the engine.
  ! Cold: OpenCPMD parsers via anonymous memfd deck (no disk write). Warm: C arrays only.
  SUBROUTINE embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, ierr, knobs)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, failed
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=600) :: star
    CHARACTER(LEN=256) :: opt
    CHARACTER(LEN=128) :: line
    CHARACTER(LEN=400) :: celltxt
    INTEGER :: celln
    ierr = 1
    nlen = 0
    failed = 0
    deck = ''
    CALL append('&CPMD'//NEW_LINE('A'))
    CALL append(' OPTIMIZE WAVEFUNCTION'//NEW_LINE('A'))
    CALL append(' CONVERGENCE ORBITALS'//NEW_LINE('A'))
    CALL append('  1.0d-5'//NEW_LINE('A'))
    CALL append(' MAXITER'//NEW_LINE('A'))
    CALL append('  40'//NEW_LINE('A'))
    CALL append(' CENTER MOLECULE OFF'//NEW_LINE('A'))
    CALL append('&END'//NEW_LINE('A'))
    CALL append('&SYSTEM'//NEW_LINE('A'))
    CALL append(' SYMMETRY'//NEW_LINE('A'))
    CALL append('  0'//NEW_LINE('A'))
    CALL append(' ANGSTROM'//NEW_LINE('A'))
    CALL format_cell_lines(cell, has_cell, celltxt, celln)
    CALL append(celltxt(1:celln))
    CALL append(' CUTOFF'//NEW_LINE('A'))
    WRITE(line, '(A,F12.6)') '  ', knobs%cutoff_ry
    CALL append(TRIM(line)//NEW_LINE('A'))
    IF (knobs%charge /= 0) THEN
      CALL append(' CHARGE'//NEW_LINE('A'))
      WRITE(line, '(A,I6)') '  ', knobs%charge
      CALL append(TRIM(line)//NEW_LINE('A'))
    END IF
    IF (knobs%mult > 1) THEN
      CALL append(' MULTIPLICITY'//NEW_LINE('A'))
      WRITE(line, '(A,I6)') '  ', knobs%mult
      CALL append(TRIM(line)//NEW_LINE('A'))
    END IF
    CALL append(' POISSON SOLVER HOCKNEY'//NEW_LINE('A'))
    CALL append('&END'//NEW_LINE('A'))
    CALL append('&DFT'//NEW_LINE('A'))
    CALL append(' OLDCODE'//NEW_LINE('A'))
    ! Honor wire/applied functional (capnp-fortran apply path); default BLYP.
    CALL append(' FUNCTIONAL '//TRIM(knobs%functional)//NEW_LINE('A'))
    CALL append('&END'//NEW_LINE('A'))
    CALL append('&ATOMS'//NEW_LINE('A'))
    seen = .FALSE.
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) THEN
        CALL embed_pp_lines_for_z(knobs%input_deck, zz, star, opt, pok)
        ierr = 2
        RETURN
      END IF
      IF (seen(zz)) CYCLE
      seen(zz) = .TRUE.
      CALL embed_pp_lines_for_z(knobs%input_deck, zz, star, opt, pok)
      IF (pok == 0) THEN
        ierr = 2
        RETURN
      END IF
      CALL append(TRIM(star)//NEW_LINE('A'))
      CALL append(TRIM(opt)//NEW_LINE('A'))
      count = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) == zz) count = count + 1
      END DO
      WRITE(line, '(A,I4)') '   ', count
      CALL append(TRIM(line)//NEW_LINE('A'))
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        WRITE(line, '(3F14.6)') pos(3*(j-1)+1), pos(3*(j-1)+2), pos(3*(j-1)+3)
        CALL append(TRIM(line)//NEW_LINE('A'))
      END DO
    END DO
    CALL append('&END'//NEW_LINE('A'))
    IF (failed /= 0) RETURN
    ierr = 0
  CONTAINS
    SUBROUTINE append(s)
      CHARACTER(LEN=*), INTENT(IN) :: s
      INTEGER :: st
      IF (failed /= 0) RETURN
      CALL append_grow(deck, nlen, s, st)
      IF (st /= 0) failed = 1
    END SUBROUTINE
  END SUBROUTINE

  ! True only when &ATOMS has a pseudopotential star and one coordinate triple
  ! before its &END. A star with no coordinates is not a geometry.
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

  ! Cap'n method deck (empty/missing &ATOMS) + geometry atoms from C arrays.
  SUBROUTINE embed_method_deck_plus_atoms(n_atoms, pos, z, cell, has_cell, &
      deck, nlen, ierr, knobs)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, mlen, failed, stat
    LOGICAL :: seen(0:120)
    CHARACTER(LEN=600) :: star
    CHARACTER(LEN=256) :: opt
    CHARACTER(LEN=128) :: line
    CHARACTER(LEN=:), ALLOCATABLE :: method
    ierr = 1
    nlen = 0
    failed = 0
    deck = ''
    IF (.NOT. ALLOCATED(knobs%input_deck)) RETURN
    CALL strip_atoms_alloc(knobs%input_deck, method, mlen, stat)
    IF (stat /= 0 .OR. mlen < 1) RETURN
    deck = method(1:mlen)
    nlen = mlen
    IF (deck(nlen:nlen) /= NEW_LINE('A')) CALL append(NEW_LINE('A'))
    CALL append('&ATOMS'//NEW_LINE('A'))
    seen = .FALSE.
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 0 .OR. zz > 120) THEN
        CALL embed_pp_lines_for_z(knobs%input_deck, zz, star, opt, pok)
        ierr = 2
        RETURN
      END IF
      IF (seen(zz)) CYCLE
      seen(zz) = .TRUE.
      CALL embed_pp_lines_for_z(knobs%input_deck, zz, star, opt, pok)
      IF (pok == 0) THEN
        ierr = 2
        RETURN
      END IF
      CALL append(TRIM(star)//NEW_LINE('A'))
      CALL append(TRIM(opt)//NEW_LINE('A'))
      count = 0
      DO j = 1, n_atoms
        IF (INT(z(j)) == zz) count = count + 1
      END DO
      WRITE(line, '(A,I4)') '   ', count
      CALL append(TRIM(line)//NEW_LINE('A'))
      DO j = 1, n_atoms
        IF (INT(z(j)) /= zz) CYCLE
        WRITE(line, '(3F14.6)') pos(3*(j-1)+1), pos(3*(j-1)+2), pos(3*(j-1)+3)
        CALL append(TRIM(line)//NEW_LINE('A'))
      END DO
    END DO
    CALL append('&END'//NEW_LINE('A'))
    ! Lattice is injected in embed_compose_cold_deck via inject_cell_if_missing.
    IF (failed /= 0) RETURN
    IF (has_cell < 0) RETURN
    IF (cell(1) < -1.0e300_c_double) RETURN
    ierr = 0
  CONTAINS
    SUBROUTINE append(s)
      CHARACTER(LEN=*), INTENT(IN) :: s
      INTEGER :: st
      IF (failed /= 0) RETURN
      CALL append_grow(deck, nlen, s, st)
      IF (st /= 0) failed = 1
    END SUBROUTINE
  END SUBROUTINE

  ! Cap'n C render of bare scalars often omits MAXITER; OpenCPMD then uses
  ! 10000 SC steps (~minutes per force). Clamp cold force decks to a finite
  ! bound unless the wire already set one.
  SUBROUTINE inject_maxiter_if_missing(deck, nlen, ierr)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=:), ALLOCATABLE :: tmp
    CHARACTER(LEN=48) :: insert
    INTEGER :: icpmd, ins_at, m, k, stat, newn
    ierr = 0
    IF (.NOT. ALLOCATED(deck) .OR. nlen < 1) RETURN
    IF (INDEX(deck(1:nlen), 'MAXITER') > 0 .OR. &
        INDEX(deck(1:nlen), 'maxiter') > 0) RETURN
    insert = ' MAXITER'//NEW_LINE('A')//'  40'//NEW_LINE('A')
    m = LEN_TRIM(insert)
    icpmd = INDEX(deck(1:nlen), '&CPMD')
    IF (icpmd == 0) icpmd = INDEX(deck(1:nlen), '&cpmd')
    IF (icpmd <= 0) RETURN
    ins_at = icpmd
    DO WHILE (ins_at <= nlen .AND. deck(ins_at:ins_at) /= NEW_LINE('A'))
      ins_at = ins_at + 1
    END DO
    IF (ins_at <= nlen) ins_at = ins_at + 1
    IF (m > HUGE(newn) - nlen) THEN
      ierr = 1
      RETURN
    END IF
    newn = nlen + m
    ALLOCATE(CHARACTER(LEN=newn) :: tmp, STAT=stat)
    IF (stat /= 0) THEN
      ierr = 1
      RETURN
    END IF
    IF (ins_at > 1) tmp(1:ins_at-1) = deck(1:ins_at-1)
    tmp(ins_at:ins_at+m-1) = insert(1:m)
    IF (ins_at <= nlen) THEN
      k = nlen - ins_at + 1
      tmp(ins_at+m:ins_at+m+k-1) = deck(ins_at:nlen)
    END IF
    CALL MOVE_ALLOC(tmp, deck)
    nlen = newn
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
  SUBROUTINE inject_cell_if_missing(deck, nlen, cell, has_cell, ierr)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=400) :: insert
    CHARACTER(LEN=:), ALLOCATABLE :: tmp
    INTEGER :: isys, iang, ins_at, m, k, stat, newn
    ierr = 0
    IF (.NOT. ALLOCATED(deck) .OR. nlen < 1) RETURN
    ! Already has a CELL keyword (do not second-guess explicit decks).
    IF (INDEX(deck(1:nlen), 'CELL') > 0 .OR. INDEX(deck(1:nlen), 'cell') > 0) &
        RETURN
    CALL format_cell_lines(cell, has_cell, insert, m)
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
    IF (m < 0 .OR. m > HUGE(newn) - nlen) THEN
      ierr = 1
      RETURN
    END IF
    newn = nlen + m
    ALLOCATE(CHARACTER(LEN=newn) :: tmp, STAT=stat)
    IF (stat /= 0) THEN
      ierr = 1
      RETURN
    END IF
    IF (ins_at > 1) tmp(1:ins_at-1) = deck(1:ins_at-1)
    IF (m > 0) tmp(ins_at:ins_at+m-1) = insert(1:m)
    IF (ins_at <= nlen) THEN
      k = nlen - ins_at + 1
      tmp(ins_at+m:ins_at+m+k-1) = deck(ins_at:nlen)
    END IF
    CALL MOVE_ALLOC(tmp, deck)
    nlen = newn
  END SUBROUTINE

  ! Shared cold-deck assembly used by SCF and by compose preview for tests.
  SUBROUTINE embed_compose_cold_deck(n_atoms, pos, z, cell, has_cell, deck, &
      nlen, ierr, knobs)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    nlen = 0
    ierr = 1
    deck = ''
    IF (ALLOCATED(knobs%input_deck)) THEN
      IF (LEN_TRIM(knobs%input_deck) > 0) THEN
        IF (deck_has_real_atoms(knobs%input_deck)) THEN
          nlen = LEN_TRIM(knobs%input_deck)
          deck = knobs%input_deck(1:nlen)
          ierr = 0
        ELSE IF (deck_has_method_sections(knobs%input_deck)) THEN
          CALL embed_method_deck_plus_atoms(n_atoms, pos, z, cell, has_cell, &
               deck, nlen, ierr, knobs)
        END IF
      END IF
    END IF
    ! ierr 2 is a missing pseudopotential. The minimal deck cannot supply one.
    IF (ierr /= 0 .AND. ierr /= 2) THEN
      CALL embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, &
           ierr, knobs)
    END IF
    IF (ierr == 0 .AND. nlen > 0) THEN
      CALL inject_cell_if_missing(deck, nlen, cell, has_cell, ierr)
      IF (ierr == 0) CALL inject_maxiter_if_missing(deck, nlen, ierr)
    END IF
  END SUBROUTINE

  FUNCTION cpmdc_embed_compose_cold_deck(n_atoms, positions_ang, atomic_numbers, &
      cell_ang, has_cell, deck_out, deck_cap, deck_len, image_c) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_compose_cold_deck_image')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    CHARACTER(KIND=c_char), INTENT(OUT) :: deck_out(*)
    INTEGER(c_int), INTENT(IN), VALUE :: deck_cap
    INTEGER(c_int), INTENT(OUT) :: deck_len
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    CHARACTER(LEN=:), ALLOCATABLE :: deck
    INTEGER :: nlen, ierr, i
    TYPE(embed_knobs) :: knobs
    ok = 0_c_int
    deck_len = 0_c_int
    IF (deck_cap < 2) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
    knobs = knobs_of(image)
    CALL embed_compose_cold_deck(INT(n_atoms), positions_ang, atomic_numbers, &
         cell_ang, INT(has_cell), deck, nlen, ierr, knobs)
    IF (ierr /= 0 .OR. nlen < 1 .OR. .NOT. ALLOCATED(deck)) RETURN
    ! One byte stays for the trailing NUL. A short buffer is an error.
    IF (nlen >= INT(deck_cap)) RETURN
    DO i = 1, nlen
      deck_out(i) = deck(i:i)
    END DO
    deck_out(nlen + 1) = c_null_char
    deck_len = INT(nlen, KIND=c_int)
    ok = 1_c_int
  END FUNCTION

#if !defined(CPMDC_HAS_CPMD)
  SUBROUTINE run_reference_pef(image, n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    TYPE(embed_knobs) :: knobs
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
    knobs = knobs_of(image)
    ok = 0_c_int
    energy_h = 0.0_c_double
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(i) / 3) RETURN
    k = 1.0e-3_real64 * MAX(0.1_real64, knobs%cutoff_ry / 70.0_real64)
    deck_scale = REAL(MAX(1, LEN_TRIM(knobs%functional) + LEN_TRIM(knobs%input_deck) + &
         LEN_TRIM(knobs%cpmd_root)), KIND=real64)
    energy_h = REAL(1.0e-8_real64 * deck_scale + &
         1.0e-6_real64 * REAL(knobs%charge + knobs%mult, KIND=real64), KIND=c_double)
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
    ! Toy isotropic stress (Ha/Bohr^3) for a periodic cell with a positive
    ! volume, so PotentialResult.stress is exercised without OpenCPMD.
    ! An isolated deck has a box volume and no tensor. CPMDC_STRESS=0 skips
    ! the tensor on a periodic cell as well.
    image%stress%values = 0.0_c_double
    image%stress%valid = 0_c_int
    IF (reference_stress_wanted() .AND. has_cell /= 0 .AND. &
         .NOT. reference_deck_isolated(knobs%input_deck)) THEN
      IF (reference_pef_fill_stress(image, cell, energy_h) /= 0) &
           image%stress%valid = 1_c_int
    END IF
  END SUBROUTINE

  LOGICAL FUNCTION reference_stress_wanted()
    CHARACTER(LEN=32) :: v
    INTEGER :: st
    reference_stress_wanted = .TRUE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_STRESS', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '0') reference_stress_wanted = .FALSE.
  END FUNCTION

  ! Symmetry 0, a missing symmetry line, CLUSTER, or an isolated-molecule
  ! keyword is a cluster. CHECK SYMMETRY is not the symmetry code. A positive
  ! symmetry code is a periodic cell.
  LOGICAL FUNCTION reference_deck_isolated(deck)
    CHARACTER(LEN=*), INTENT(IN) :: deck
    INTEGER :: p, n, sym, saw
    reference_deck_isolated = .TRUE.
    n = LEN_TRIM(deck)
    IF (INDEX(deck, ' ISOLATED MOLECULE') > 0) RETURN
    IF (INDEX(deck, ' MOLECULE ISOLATED') > 0) RETURN
    IF (INDEX(deck, ' CLUSTER') > 0) RETURN
    p = 1
    DO WHILE (p <= n - 7)
      IF (deck(p:p+7) == 'SYMMETRY') THEN
        IF (p >= 7) THEN
          IF (deck(p-6:p-1) == 'CHECK ') THEN
            p = p + 1
            CYCLE
          END IF
        END IF
        p = p + 8
        saw = 0
        sym = 0
        DO WHILE (p <= n)
          IF (deck(p:p) >= '0' .AND. deck(p:p) <= '9') THEN
            saw = 1
            sym = sym * 10 + IACHAR(deck(p:p)) - IACHAR('0')
            p = p + 1
          ELSE IF (saw == 1) THEN
            EXIT
          ELSE
            p = p + 1
          END IF
        END DO
        IF (saw == 1 .AND. sym > 0) reference_deck_isolated = .FALSE.
        RETURN
      END IF
      p = p + 1
    END DO
  END FUNCTION

  INTEGER FUNCTION reference_pef_fill_stress(image, cell, energy_h)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    REAL(c_double), INTENT(IN) :: cell(*)
    REAL(c_double), INTENT(IN) :: energy_h
    REAL(real64), PARAMETER :: bohr_to_ang = 0.529177210903_real64
    REAL(real64) :: a(3), b(3), c(3), vol, inv_b, sig
    INTEGER :: i
    reference_pef_fill_stress = 0
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
    reference_pef_fill_stress = 1
  END FUNCTION

#endif

#if defined(CPMDC_HAS_CPMD)
  SUBROUTINE snapshot_prop_from_modules(image, n_atoms, grad)
    USE ddip, ONLY: pdipole
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: grad(*)
    INTEGER :: i, ncopy
    image%prop%valid = 1_c_int
    image%prop%dipole_count = 3_c_size_t
    DO i = 1, 3
      image%prop%dipole(i) = REAL(pdipole(i), KIND=c_double)
    END DO
    image%prop%polarizability_count = 9_c_size_t
    DO i = 1, 9
      image%prop%polarizability(i) = 0.0_c_double
    END DO
    ncopy = n_atoms * 3
    IF (ncopy > 4096) ncopy = 4096
    IF (ncopy < 0) ncopy = 0
    DO i = 1, ncopy
      image%prop%hessian(i) = grad(i)
    END DO
    image%prop%hessian_count = INT(ncopy, KIND=c_size_t)
  END SUBROUTINE

  SUBROUTINE embed_set_tau0_from_pos(n_atoms, pos, z, origin, ierr)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_species_order_map, cpmdc_note_embed_failure
    USE coor, ONLY: tau0
    USE ions, ONLY: ions0, ions1
    USE cnst, ONLY: fbohr
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER(c_int), INTENT(OUT) :: origin(*)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: is, ia, k, slot, j, astat
    INTEGER(c_int), ALLOCATABLE :: species_z(:), species_n(:)
    ierr = 1
    IF (.NOT. ALLOCATED(tau0) .OR. ions1%nsp < 1) THEN
      CALL cpmdc_note_embed_failure( &
           'CPMD species tables are not initialized'//c_null_char)
      RETURN
    END IF
    ALLOCATE(species_z(ions1%nsp), species_n(ions1%nsp), STAT=astat)
    IF (astat /= 0) THEN
      CALL cpmdc_note_embed_failure('out of memory'//c_null_char)
      RETURN
    END IF
    DO is = 1, ions1%nsp
      species_z(is) = INT(ions0%iatyp(is), KIND=c_int)
      species_n(is) = INT(ions0%na(is), KIND=c_int)
    END DO
    ! CPMD stores one contiguous block per species. ForceInput may list the
    ! same element in more than one run; origin maps each block slot back.
    IF (cpmdc_species_order_map(INT(n_atoms, KIND=c_int), z, &
        INT(ions1%nsp, KIND=c_int), species_z, species_n, origin) /= 0_c_int) THEN
      CALL cpmdc_note_embed_failure( &
           'atomic numbers do not match the CPMD species blocks'//c_null_char)
      RETURN
    END IF
    slot = 0
    DO is = 1, ions1%nsp
      DO ia = 1, ions0%na(is)
        slot = slot + 1
        j = INT(origin(slot)) + 1
        DO k = 1, 3
          tau0(k, ia, is) = REAL(pos(3*(j-1)+k), KIND=real64) * fbohr
        END DO
      END DO
    END DO
    ierr = 0
  END SUBROUTINE

  SUBROUTINE embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, ok)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_scatter_species_gradient, &
        cpmdc_note_embed_failure
    USE wfopts_utils, ONLY: wfopts
    USE rwfopt_utils, ONLY: embed_set_warm_orbitals, embed_set_need_forces, &
        embed_set_write_files
    USE phfac_utils, ONLY: phfac
    USE ener, ONLY: ener_com, chrg, ener_c, ener_d
    USE coor, ONLY: tau0, fion, taup
    USE ions, ONLY: ions0, ions1
    USE store_types, ONLY: cprint, iprint_force, restart1
    USE system, ONLY: cnti, cntl, parm
    USE benc, ONLY: ibench
    USE strs, ONLY: paiu
    USE isos, ONLY: isos1
    USE ropt, ONLY: ropt_mod, iteropt
    USE parac, ONLY: paral
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, is, ia, k, idx, nmax, i, j, astat, slot
    INTEGER(c_int), ALLOCATABLE :: origin(:)
    REAL(c_double), ALLOCATABLE :: species_grad(:)
    REAL(real64) :: omega
    LOGICAL :: stress_computed, was_diis, was_pcg, was_pcgmin, was_prec
    INTEGER :: inwfun_deck, ibench_deck
    INTERFACE
      FUNCTION cpmdc_stop_code() BIND(C, NAME='cpmdc_stop_code') RESULT(code)
        IMPORT :: c_int
        INTEGER(c_int) :: code
      END FUNCTION cpmdc_stop_code
    END INTERFACE
    ok = 0_c_int
    energy_h = 0.0_c_double
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(nmax) / 3) RETURN
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    image%stress%valid = 0_c_int
    image%stress%values = 0.0_c_double
    ALLOCATE(origin(n_atoms), species_grad(nmax), STAT=astat)
    IF (astat /= 0) THEN
      CALL cpmdc_note_embed_failure('out of memory'//c_null_char)
      RETURN
    END IF
    CALL embed_set_tau0_from_pos(n_atoms, pos, z, origin, ierr)
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
    CALL embed_set_need_forces(.TRUE.)
    ! Orbitals for the next call stay in embed_c0_store. rwfopt then skips
    ! zhwwf and geofile, so this call writes no RESTART, LATEST, or GEOMETRY.
    ! cpmd.x leaves embed_write_files at its default, which is true.
    CALL embed_set_write_files(.FALSE.)
    ! PEF stress: totstr fills paiu when cntl%tpres. Skip isolated/Hockney
    ! (tclust): rinitwf does tpres then newcell then gf_periodic but scg is only
    ! allocated for periodic cells in initclust — SEGV on cluster decks.
    stress_computed = .NOT. isos1%tclust .AND. embed_stress_wanted()
    IF (stress_computed) cntl%tpres = .TRUE.
    ! Store a converged c0 for the next call. Restore does nothing until
    ! that store exists, and an unconverged SCF leaves the previous
    ! converged copy in place, so a later SCF is still a warm start.
    ! Converge to cntr%tolog with the deck MAXITER. Do not clamp nomore_iter.
    ! initrun builds starting orbitals before the stored c0 replaces them.
    ! Once a converged copy exists (a warm call), the simple atomic
    ! superposition (inwfun 3: loadc, one orthogonalisation, no force
    ! evaluation) stands in for the Lanczos guess and for the random start,
    ! both of which run a full SCF step on orbitals that are then discarded.
    inwfun_deck = cnti%inwfun
    ibench_deck = ibench(1)
    ! The cheap guess is for a process that already holds c0. A new
    ! session's counter is 0; the latch says the basis still matches.
    IF (image%cfg_warm_steps > 0 .OR. embed_basis_latched) THEN
      cnti%inwfun = 3
      ibench(1) = 1
    END IF
    CALL embed_set_warm_orbitals(.TRUE.)
    CALL wfopts
    CALL note_scf_steps(iteropt%nfi)
    ! ODIIS can exhaust MAXITER short of the orbital threshold. That pass
    ! does not replace the stored c0. Continue once with PCG MINIMIZE from
    ! the previous converged copy, which still counts as a warm start.
    ! A stopgm in the first pass leaves the CPMD state undefined, so it
    ! does not start this continuation.
    IF (.NOT. ropt_mod%convwf .AND. cntl%diis .AND. &
        cpmdc_stop_code() == 0_c_int) THEN
      was_diis = cntl%diis
      was_pcg = cntl%pcg
      was_pcgmin = cntl%pcgmin
      was_prec = cntl%prec
      cntl%diis = .FALSE.
      cntl%pcg = .TRUE.
      cntl%pcgmin = .TRUE.
      cntl%prec = .TRUE.
      ropt_mod%spcg = .TRUE.
      IF (paral%io_parent) WRITE(6, '(A)') &
          ' cpmdc: ODIIS did not converge; continuing with PCG MINIMIZE'
      CALL embed_set_warm_orbitals(.TRUE.)
      IF (ALLOCATED(fion)) DEALLOCATE(fion)
      IF (ALLOCATED(taup)) DEALLOCATE(taup)
      CALL wfopts
      CALL note_scf_steps(iteropt%nfi)
      cntl%diis = was_diis
      cntl%pcg = was_pcg
      cntl%pcgmin = was_pcgmin
      cntl%prec = was_prec
    END IF
    cnti%inwfun = inwfun_deck
    ibench(1) = ibench_deck
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
    ! ENERGY-file row. EKINC stays zero off an MD step.
    image%md%values = 0.0_c_double
    image%md%values(1) = energy_h
    image%md%values(2) = image%energy%ekin
    image%md%values(3) = image%energy%epseu
    image%md%values(4) = image%energy%enl
    image%md%values(5) = image%energy%eht
    image%md%values(6) = image%energy%exc
    image%md%count = 12_c_size_t
    image%md%valid = 1_c_int
    species_grad = 0.0_c_double
    IF (ALLOCATED(fion)) THEN
      IF (SIZE(fion, 1) >= 3 .AND. SIZE(fion, 2) >= 1 .AND. SIZE(fion, 3) >= ions1%nsp) THEN
        slot = 0
        DO is = 1, ions1%nsp
          DO ia = 1, ions0%na(is)
            IF (ia > SIZE(fion, 2)) EXIT
            slot = slot + 1
            IF (slot > n_atoms) EXIT
            DO k = 1, 3
              species_grad(3 * (slot - 1) + k) = REAL(-fion(k, ia, is), KIND=c_double)
            END DO
          END DO
        END DO
        CALL cpmdc_scatter_species_gradient(INT(n_atoms, KIND=c_int), origin, &
             species_grad, grad)
      END IF
    END IF
    ! Cartesian stress Ha/Bohr^3: OpenCPMD stores virial in paiu (energy),
    ! true stress is paiu/omega (see totstr/wrstress). Row-major for Cap'n Proto.
    ! totstr runs only when cntl%tpres was set above. A positive omega is the
    ! cell volume, including an isolated box whose paiu was not computed.
    omega = parm%omega
    IF (stress_computed .AND. omega > 1.0e-30_real64) THEN
      DO i = 1, 3
        DO j = 1, 3
          image%stress%values(3 * (i - 1) + j) = REAL(paiu(i, j) / omega, KIND=c_double)
        END DO
      END DO
      image%stress%valid = 1_c_int
    END IF
    CALL snapshot_prop_from_modules(image, n_atoms, grad)
    ! rwfopt computes the ionic forces only for converged orbitals: an SCF
    ! that ran out of MAXITER leaves fion zero, which a caller would read as
    ! a stationary point. Report it as a failure instead.
    IF (ropt_mod%convwf .AND. ABS(energy_h) > 1.0e-8_c_double) ok = 1_c_int
  END SUBROUTINE


  ! CPMDC_SCF_STEPS=1 reports the iteration count of the SCF that just
  ! returned. The error unit stays attached when unit 6 is a file.
  SUBROUTINE note_scf_steps(nfi)
    USE ISO_FORTRAN_ENV, ONLY: error_unit
    USE parac, ONLY: paral
    INTEGER, INTENT(IN) :: nfi
    CHARACTER(LEN=8) :: v
    INTEGER :: st
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_SCF_STEPS', v, STATUS=st)
    IF (st /= 0) RETURN
    IF (TRIM(v) /= '1') RETURN
    IF (.NOT. paral%io_parent) RETURN
    WRITE(error_unit, '(A,I0)') 'cpmdc_scf_steps ', nfi
  END SUBROUTINE

  ! CPMDC_STRESS=0 skips the stress tensor on periodic cells. A caller that
  ! takes only energy and forces does not need totstr.
  LOGICAL FUNCTION embed_stress_wanted()
    CHARACTER(LEN=16) :: v
    INTEGER :: st
    embed_stress_wanted = .TRUE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_STRESS', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '0') embed_stress_wanted = .FALSE.
  END FUNCTION

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
    embed_saved_has_cell = has_cell
    embed_saved_cell_set = 1
    embed_saved_cell = 0.0_c_double
    IF (has_cell == 0) RETURN
    DO i = 1, 9
      image%warm_cell(i) = cell(i)
      embed_saved_cell(i) = cell(i)
    END DO
  END SUBROUTINE

  LOGICAL FUNCTION embed_use_output_dir(image)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_enter_output_cwd, cpmdc_note_embed_failure
    TYPE(cpmdc_embed_image), INTENT(IN) :: image
    embed_use_output_dir = .FALSE.
    IF (cpmdc_enter_output_cwd(image%output_dir) /= 0_c_int) THEN
      CALL cpmdc_note_embed_failure( &
           'CPMD output directory is not a directory'//c_null_char)
      RETURN
    END IF
    embed_use_output_dir = .TRUE.
  END FUNCTION

  SUBROUTINE clear_embed_basis_latch()
    embed_basis_latched = .FALSE.
    embed_saved_cutoff = -1.0_c_double
    embed_saved_charge = -999
    embed_saved_mult = -999
    embed_saved_natoms = -1
    embed_saved_functional = ''
    embed_saved_has_cell = 0
    embed_saved_cell_set = 0
    embed_saved_cell = 0.0_c_double
    IF (ALLOCATED(embed_saved_z)) DEALLOCATE(embed_saved_z)
    IF (ALLOCATED(embed_saved_deck)) DEALLOCATE(embed_saved_deck)
  END SUBROUTINE

  SUBROUTINE latch_embed_basis(n_atoms, z, knobs)
    USE rwfopt_utils, ONLY: embed_reset_warm_orbitals
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: z(*)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER :: i, astat
    IF (ALLOCATED(embed_saved_z)) DEALLOCATE(embed_saved_z)
    ALLOCATE(embed_saved_z(n_atoms), STAT=astat)
    IF (astat /= 0) THEN
      CALL embed_reset_warm_orbitals()
      CALL clear_embed_basis_latch()
      RETURN
    END IF
    DO i = 1, n_atoms
      embed_saved_z(i) = INT(z(i))
    END DO
    embed_saved_cutoff = REAL(knobs%cutoff_ry, KIND=c_double)
    embed_saved_charge = knobs%charge
    embed_saved_mult = knobs%mult
    embed_saved_natoms = n_atoms
    embed_saved_functional = knobs%functional
    IF (ALLOCATED(embed_saved_deck)) DEALLOCATE(embed_saved_deck)
    embed_saved_deck = ''
    IF (ALLOCATED(knobs%input_deck)) embed_saved_deck = knobs%input_deck
    embed_basis_latched = .TRUE.
  END SUBROUTINE

  LOGICAL FUNCTION same_element_counts(n_atoms, z)
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER :: counts(118), i, zz
    same_element_counts = .FALSE.
    IF (.NOT. ALLOCATED(embed_saved_z)) RETURN
    IF (SIZE(embed_saved_z) /= n_atoms) RETURN
    counts = 0
    DO i = 1, n_atoms
      zz = embed_saved_z(i)
      IF (zz < 1 .OR. zz > 118) RETURN
      counts(zz) = counts(zz) + 1
    END DO
    DO i = 1, n_atoms
      zz = INT(z(i))
      IF (zz < 1 .OR. zz > 118) RETURN
      counts(zz) = counts(zz) - 1
    END DO
    DO i = 1, 118
      IF (counts(i) /= 0) RETURN
    END DO
    same_element_counts = .TRUE.
  END FUNCTION

  LOGICAL FUNCTION embed_basis_changed(image, n_atoms, z, cell, has_cell, knobs)
    TYPE(cpmdc_embed_image), INTENT(IN) :: image
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(IN) :: cell(*)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER :: i
    REAL(c_double) :: cutoff
    embed_basis_changed = .FALSE.
    IF (.NOT. embed_basis_latched) RETURN
    cutoff = REAL(knobs%cutoff_ry, KIND=c_double)
    ! A detached image has no cell of its own. The cell stored with the
    ! orbitals is the one a new session has to match. c0 is plane-wave
    ! coefficients, so the count of each element is the composition.
    IF (image%warm_cell_set /= 0) THEN
      IF (.NOT. warm_cell_matches(image, cell, has_cell)) &
          embed_basis_changed = .TRUE.
    ELSE IF (embed_saved_cell_set == 0) THEN
      embed_basis_changed = .TRUE.
    ELSE IF (has_cell /= embed_saved_has_cell) THEN
      embed_basis_changed = .TRUE.
    ELSE IF (has_cell /= 0) THEN
      DO i = 1, 9
        IF (ABS(cell(i) - embed_saved_cell(i)) > 1.0e-8_c_double) &
            embed_basis_changed = .TRUE.
      END DO
    END IF
    IF (ABS(cutoff - embed_saved_cutoff) > 1.0e-8_c_double) &
        embed_basis_changed = .TRUE.
    IF (knobs%charge /= embed_saved_charge) embed_basis_changed = .TRUE.
    IF (knobs%mult /= embed_saved_mult) embed_basis_changed = .TRUE.
    IF (TRIM(knobs%functional) /= TRIM(embed_saved_functional)) &
        embed_basis_changed = .TRUE.
    IF (.NOT. ALLOCATED(embed_saved_deck)) THEN
      IF (ALLOCATED(knobs%input_deck)) THEN
        IF (LEN_TRIM(knobs%input_deck) > 0) embed_basis_changed = .TRUE.
      END IF
    ELSE IF (.NOT. ALLOCATED(knobs%input_deck)) THEN
      IF (LEN_TRIM(embed_saved_deck) > 0) embed_basis_changed = .TRUE.
    ELSE IF (knobs%input_deck /= embed_saved_deck) THEN
      embed_basis_changed = .TRUE.
    END IF
    IF (n_atoms /= embed_saved_natoms) embed_basis_changed = .TRUE.
    IF (.NOT. same_element_counts(n_atoms, z)) embed_basis_changed = .TRUE.
  END FUNCTION

  SUBROUTINE run_embed_scf(image, n_atoms, pos, z, cell, has_cell, energy_h, grad, ok)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_memfd_write, &
        cpmdc_pseudopotential_directory, cpmdc_prepare_pp_cwd, &
        cpmdc_restore_host_cwd, cpmdc_note_embed_failure
    USE rwfopt_utils, ONLY: embed_reset_warm_orbitals
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
    TYPE(embed_knobs) :: knobs
    LOGICAL :: tinfo
    CHARACTER(LEN=:), ALLOCATABLE :: deck
    CHARACTER(LEN=64) :: mempath
    CHARACTER(KIND=c_char) :: pp_c(1024)
    ok = 0_c_int
    energy_h = 0.0_c_double
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(nmax) / 3) RETURN
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    ! Saved orbitals match one method, cell, and elemental composition.
    ! A shape mismatch in embed_restore_orbitals calls stopgm, so drop
    ! the copy before rwfopt when any of those change. A new session and
    ! a reordering of the same atomic numbers are not a new basis.
    knobs = knobs_of(image)
    IF (embed_basis_latched .AND. &
        .NOT. embed_basis_changed(image, n_atoms, z, cell, has_cell, knobs)) THEN
      IF (.NOT. embed_use_output_dir(image)) THEN
        ierr = cpmdc_restore_host_cwd()
        RETURN
      END IF
      CALL embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, ok)
      ierr = cpmdc_restore_host_cwd()
      ! A failed SCF leaves the counter and the previous converged c0.
      ! The next call on this basis is still warm.
      IF (ok /= 0_c_int) THEN
        image%cfg_warm_steps = image%cfg_warm_steps + 1
        CALL latch_warm_cell(image, cell, has_cell)
        CALL latch_embed_basis(n_atoms, z, knobs)
        IF (.NOT. embed_basis_latched) image%cfg_warm_steps = 0
      END IF
      RETURN
    END IF
    IF (embed_basis_latched) THEN
      CALL embed_reset_warm_orbitals()
      CALL clear_embed_basis_latch()
      image%cfg_warm_steps = 0
    END IF
    ! Cold: honor the rendered method deck.
    ! 1) &ATOMS that already has coordinates is kept.
    ! 2) Method sections keep their text. &ATOMS is rebuilt from the step.
    !    A '!SPECIES' entry supplies that element's file, LMAX, LOC, and
    !    KLEINMAN-BYLANDER. Other elements use the built-in table.
    ! 3) Else a minimal deck with the applied functional, cutoff, charge, and
    !    multiplicity.
    ! Geometry for forces always comes from the C arrays into TAU0 after parse.
    CALL embed_compose_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, &
         ierr, knobs)
    IF (ierr /= 0 .OR. nlen < 1 .OR. .NOT. ALLOCATED(deck)) RETURN
    IF (nlen > HUGE(0_c_int)) RETURN
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
    ! The process leaves that directory again once the files are read.
    IF (cpmdc_pseudopotential_directory(pp_c, INT(1024, KIND=c_size_t)) /= 0_c_int) THEN
      CALL cpmdc_note_embed_failure(pp_c)
      RETURN
    END IF
    IF (cpmdc_prepare_pp_cwd(pp_c) /= 0_c_int) THEN
      CALL cpmdc_note_embed_failure( &
           'cannot enter the pseudopotential directory'//c_null_char)
      RETURN
    END IF
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
    ! Pseudopotential files are in memory. Leave that directory before the
    ! SCF so a stray write cannot land on the library.
    IF (.NOT. embed_use_output_dir(image)) THEN
      ierr = cpmdc_restore_host_cwd()
      RETURN
    END IF
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
      CALL latch_embed_basis(n_atoms, z, knobs)
      IF (.NOT. embed_basis_latched) image%cfg_warm_steps = 0
    ELSE
      CALL clear_last_energy_components(image)
    END IF
    ! The force call writes no RESTART, LATEST, or GEOMETRY.
    IF (cpmdc_restore_host_cwd() /= 0_c_int) THEN
      ! Keep ok from SCF; lost host CWD is non-fatal for the energy itself.
    END IF
  END SUBROUTINE

#endif
END MODULE cpmd_embed_c_api
