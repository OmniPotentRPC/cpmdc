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
    FUNCTION cpmdc_embed_failure_pending() &
        BIND(C, NAME='cpmdc_embed_failure_pending') RESULT(pending)
      IMPORT :: c_int
      INTEGER(c_int) :: pending
    END FUNCTION
    SUBROUTINE cpmdc_note_scf_unconverged(minimizer, iterations, maxiter, gemax, &
        tolog, first_minimizer, first_iterations) &
        BIND(C, NAME='cpmdc_note_scf_unconverged')
      IMPORT :: c_int, c_double
      INTEGER(c_int), VALUE :: minimizer, iterations, maxiter, first_minimizer, &
          first_iterations
      REAL(c_double), VALUE :: gemax, tolog
    END SUBROUTINE
    SUBROUTINE cpmdc_note_deck_external_potential() &
        BIND(C, NAME='cpmdc_note_deck_external_potential')
    END SUBROUTINE
    SUBROUTINE cpmdc_note_periodic_embedding() &
        BIND(C, NAME='cpmdc_note_periodic_embedding')
    END SUBROUTINE
    FUNCTION cpmdc_deck_defines_external_potential(deck, deck_len) &
        BIND(C, NAME='cpmdc_deck_defines_external_potential')
      IMPORT :: c_char, c_int, c_size_t
      CHARACTER(KIND=c_char), INTENT(IN) :: deck(*)
      INTEGER(c_size_t), INTENT(IN), VALUE :: deck_len
      INTEGER(c_int) :: cpmdc_deck_defines_external_potential
    END FUNCTION
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
    FUNCTION cpmdc_memfd_close(fd) BIND(C, NAME='cpmdc_memfd_close')
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: fd
      INTEGER(c_int) :: cpmdc_memfd_close
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
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cpmdc_embed_init, cpmdc_embed_available, cpmdc_embed_finalize
  PUBLIC :: cpmdc_embed_bind_calculator
  PUBLIC :: cpmdc_embed_store_comm, cpmdc_embed_adopted_fcomm
  PUBLIC :: cpmdc_embed_reset_state, cpmdc_embed_teardown
  PUBLIC :: cpmdc_embed_set_config, cpmdc_embed_set_deck, cpmdc_embed_energy_grad
  PUBLIC :: cpmdc_embed_compose_cold_deck
  PUBLIC :: cpmdc_embed_select_orbitals, cpmdc_embed_orbital_slots
  PUBLIC :: cpmdc_embed_clear_orbital_key, cpmdc_embed_extrapolation_state
  PUBLIC :: cpmdc_embed_prepare_orbitals, cpmdc_embed_abort_orbitals
  PUBLIC :: cpmdc_embed_commit_orbitals
  PUBLIC :: cpmdc_embed_drop_stale_orbitals
  PUBLIC :: cpmdc_embed_copy_orbitals
  PUBLIC :: cpmdc_embed_last_warm_source

  LOGICAL, SAVE :: runtime_ready = .FALSE.
  LOGICAL, SAVE :: runtime_finalized = .FALSE.
  ! Where the last SCF started; see note_warm_start for the values.
  INTEGER, SAVE :: last_warm_source = -2
  REAL(c_double), SAVE :: embed_cell_in_place_limit = 0.0_c_double
  ! OpenCPMD's tistart runs once per process: a second call prints
  ! CALLED ONCE and executes STOP, which exits the host with status 0.
  LOGICAL, SAVE :: timer_started = .FALSE.
  ! True after a cold setup list finishes without a stopgm, until
  ! cpmdc_embed_teardown frees what that setup allocated. A second cold
  ! setup tears the first down.
  LOGICAL, SAVE :: cpmd_set_up = .FALSE.
  ! True from the first call of a cold setup list until that list
  ! finishes or teardown frees the partial state. cpmd_set_up stays
  ! false while this is set.
  LOGICAL, SAVE :: cpmd_setup_started = .FALSE.
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
  ! Beside the public snapshot, same layout as CPMDCEmbedCharges.
  TYPE, BIND(C) :: cpmdc_embed_charges
    INTEGER(c_int) :: valid
    INTEGER(c_size_t) :: count
    REAL(c_double) :: values(1365)
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
    TYPE(cpmdc_embed_charges) :: charges
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
  ! The rigid shift of one call on an isolated deck; see embed_recentre_cluster.
  TYPE :: embed_recentring
    LOGICAL :: applied = .FALSE.
    REAL(real64) :: shift_ang(3) = 0.0_real64
    REAL(real64) :: proj(3, 3) = 0.0_real64
    REAL(real64), ALLOCATABLE :: weights(:)
  END TYPE
#if defined(CPMDC_HAS_CPMD)
  REAL(c_double), SAVE :: tcpu0 = 0.0_c_double, twall0 = 0.0_c_double
  ! mp_start assigns mp_comm_world only when it calls MPI_Init.
  ! The flag records that this process has stored a calculator communicator.
  ! cpmdc_embed_bind_calculator splits once.
  ! cpmdc_embed_store_comm may replace that communicator. The C caller
  ! refuses the replacement while a session that has evaluated is live.
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
  REAL(real64), ALLOCATABLE, SAVE :: embed_external_pos(:,:), embed_external_q(:)
  REAL(real64), ALLOCATABLE, SAVE :: embed_external_sigma(:), embed_external_force(:,:)
  REAL(real64), ALLOCATABLE, TARGET, SAVE :: embed_external_potential(:)
  INTEGER, SAVE :: embed_external_count = 0
  ! Lattice rows, Bohr, of the cell whose G sphere the setup built. The
  ! stored orbitals are coefficients on those G vectors wherever an
  ! in-place cell move takes them; see embed_cell_moves_in_place.
  LOGICAL, SAVE :: embed_setup_ht_set = .FALSE.
  REAL(real64), SAVE :: embed_setup_ht(3, 3) = 0.0_real64
#endif

  TYPE :: embed_knobs
    CHARACTER(LEN=64) :: functional = 'BLYP'
    REAL(real64) :: cutoff_ry = 70.0_real64
    INTEGER :: charge = 0
    INTEGER :: mult = 1
    CHARACTER(LEN=:), ALLOCATABLE :: input_deck
    CHARACTER(LEN=1024) :: cpmd_root = ' '
  END TYPE

  INTERFACE
    FUNCTION cpmdc_cell_rank3_state(cell) BIND(C, NAME='cpmdc_cell_rank3_state')
      IMPORT :: c_double, c_int
      REAL(c_double), INTENT(IN) :: cell(9)
      INTEGER(c_int) :: cpmdc_cell_rank3_state
    END FUNCTION
    FUNCTION cpmdc_cell_volume3(cell, volume) BIND(C, NAME='cpmdc_cell_volume3')
      IMPORT :: c_double, c_int
      REAL(c_double), INTENT(IN) :: cell(9)
      REAL(c_double), INTENT(OUT) :: volume
      INTEGER(c_int) :: cpmdc_cell_volume3
    END FUNCTION
  END INTERFACE

  ! At r/sigma = 1/8 the closed force coefficient loses at most about
  ! 64*EPSILON to cancellation. Below it, the series through (r/sigma)^14
  ! has a smaller truncation error than double precision at the boundary.
  REAL(real64), PARAMETER :: embed_gaussian_series_switch = &
      1.0_real64 / 8.0_real64

CONTAINS

  SUBROUTINE embed_gaussian_kernel(rr, sigma, potential, force_coefficient, &
      force_derivative)
    REAL(real64), INTENT(IN) :: rr, sigma
    REAL(real64), INTENT(OUT) :: potential, force_coefficient, force_derivative
    REAL(real64) :: root, t, t2, erfv, expo
    root = SQRT(2.0_real64 / ACOS(-1.0_real64))
    t = rr / sigma
    t2 = t * t
    IF (t <= embed_gaussian_series_switch) THEN
      potential = root / sigma * (1.0_real64 - t2 / 6.0_real64 + &
          t2**2 / 40.0_real64 - t2**3 / 336.0_real64 + &
          t2**4 / 3456.0_real64 - t2**5 / 42240.0_real64 + &
          t2**6 / 599040.0_real64 - t2**7 / 9676800.0_real64)
      force_coefficient = root / sigma**3 * (1.0_real64 / 3.0_real64 - &
          t2 / 10.0_real64 + t2**2 / 56.0_real64 - &
          t2**3 / 432.0_real64 + t2**4 / 4224.0_real64 - &
          t2**5 / 49920.0_real64 + t2**6 / 691200.0_real64 - &
          t2**7 / 10967040.0_real64)
      force_derivative = root / sigma**4 * t * (-1.0_real64 / 5.0_real64 + &
          t2 / 14.0_real64 - t2**2 / 72.0_real64 + &
          t2**3 / 528.0_real64 - t2**4 / 4992.0_real64 + &
          t2**5 / 57600.0_real64 - t2**6 / 783360.0_real64)
    ELSE
      erfv = ERF(t / SQRT(2.0_real64))
      expo = EXP(-0.5_real64 * t2)
      potential = erfv / rr
      force_coefficient = erfv / rr**3 - root * expo / (sigma * rr**2)
      force_derivative = 3.0_real64 * root * expo / (sigma * rr**3) + &
          root * expo / (sigma**3 * rr) - 3.0_real64 * erfv / rr**4
    END IF
  END SUBROUTINE

#if defined(CPMDC_HAS_CPMD)
  SUBROUTINE cpmdc_external_charge_density_force(rhoe) &
      BIND(C, NAME='cpmdc_external_charge_density_force')
    USE system, ONLY: fpar, parap, parm, spar
    USE parac, ONLY: parai
    USE mp_interface, ONLY: mp_sum
    REAL(c_double), INTENT(IN) :: rhoe(*)
    REAL(real64) :: point(3), d(3), rr, sigma, potential, coef, derivative, dv
    INTEGER :: ix, iy, iz, lx, idx, m
    IF (embed_external_count <= 0) RETURN
    embed_external_force = 0.0_real64
    dv = parm%omega / REAL(spar%nr1s * spar%nr2s * spar%nr3s, KIND=real64)
    DO iz = 1, spar%nr3s
      DO iy = 1, spar%nr2s
        DO ix = parap%nrxpl(parai%mepos,1), parap%nrxpl(parai%mepos,2)
          lx = ix - parap%nrxpl(parai%mepos,1) + 1
          idx = lx + fpar%kr1 * ((iy - 1) + fpar%kr2s * (iz - 1))
          point = REAL(ix - 1, real64) / REAL(spar%nr1s, real64) * parm%a1 + &
               REAL(iy - 1, real64) / REAL(spar%nr2s, real64) * parm%a2 + &
               REAL(iz - 1, real64) / REAL(spar%nr3s, real64) * parm%a3
          DO m = 1, embed_external_count
            d = point - embed_external_pos(:,m)
            rr = NORM2(d)
            sigma = embed_external_sigma(m)
            CALL embed_gaussian_kernel(rr, sigma, potential, coef, derivative)
            coef = embed_external_q(m) * coef
            embed_external_force(:,m) = embed_external_force(:,m) + &
                 REAL(rhoe(idx), real64) * coef * d * dv
          END DO
        END DO
      END DO
    END DO
    CALL mp_sum(embed_external_force, 3 * embed_external_count, parai%allgrp)
  END SUBROUTINE

  SUBROUTINE prepare_external_charges(n_external, external_pos, external_q, &
      external_width, shift_ang, ierr)
    USE cnst, ONLY: fbohr
    USE efld, ONLY: extf, textfld, set_cpmdc_external_force_hook
    USE isos, ONLY: isos1
    USE system, ONLY: fpar, parap, parm, spar
    USE parac, ONLY: parai
    INTEGER, INTENT(IN) :: n_external
    REAL(c_double), INTENT(IN) :: external_pos(*), external_q(*), external_width(*)
    REAL(real64), INTENT(IN) :: shift_ang(3)
    INTEGER, INTENT(OUT) :: ierr
    REAL(real64) :: point(3), d(3), rr, value, potential, coef, derivative
    INTEGER :: ix, iy, iz, lx, idx, m, astat
    LOGICAL :: own_potential
    ierr = 0
    embed_external_count = 0
    IF (ALLOCATED(embed_external_pos)) DEALLOCATE(embed_external_pos)
    IF (ALLOCATED(embed_external_q)) DEALLOCATE(embed_external_q)
    IF (ALLOCATED(embed_external_sigma)) DEALLOCATE(embed_external_sigma)
    IF (ALLOCATED(embed_external_force)) DEALLOCATE(embed_external_force)
    own_potential = .FALSE.
    IF (ALLOCATED(embed_external_potential)) THEN
      own_potential = ASSOCIATED(extf, embed_external_potential)
      IF (own_potential) NULLIFY(extf)
      DEALLOCATE(embed_external_potential)
      IF (own_potential) textfld = .FALSE.
    END IF
    IF (n_external <= 0) THEN
      CALL set_cpmdc_external_force_hook()
      RETURN
    END IF
    ! The kernel has no image sum, so only an isolated system is valid.
    IF (.NOT. isos1%tisos) THEN
      ierr = 3
      CALL set_cpmdc_external_force_hook()
      RETURN
    END IF
    IF (textfld) THEN
      ierr = 2
      CALL set_cpmdc_external_force_hook()
      RETURN
    END IF
    ALLOCATE(embed_external_pos(3,n_external), embed_external_q(n_external), &
        embed_external_sigma(n_external), embed_external_force(3,n_external), &
        embed_external_potential(fpar%kr1*fpar%kr2s*fpar%kr3s), STAT=astat)
    IF (astat /= 0) THEN
      ierr = 1
      RETURN
    END IF
    DO m = 1, n_external
      embed_external_pos(:,m) = &
          (REAL(external_pos(3*m-2:3*m), real64) + shift_ang) * fbohr
      embed_external_q(m) = REAL(external_q(m), real64)
      embed_external_sigma(m) = REAL(external_width(m), real64) * fbohr
    END DO
    embed_external_force = 0.0_real64
    embed_external_potential = 0.0_real64
    DO iz = 1, spar%nr3s
      DO iy = 1, spar%nr2s
        DO ix = parap%nrxpl(parai%mepos,1), parap%nrxpl(parai%mepos,2)
          lx = ix - parap%nrxpl(parai%mepos,1) + 1
          idx = lx + fpar%kr1 * ((iy - 1) + fpar%kr2s * (iz - 1))
          point = REAL(ix - 1, real64) / REAL(spar%nr1s, real64) * parm%a1 + &
               REAL(iy - 1, real64) / REAL(spar%nr2s, real64) * parm%a2 + &
               REAL(iz - 1, real64) / REAL(spar%nr3s, real64) * parm%a3
          value = 0.0_real64
          DO m = 1, n_external
            d = point - embed_external_pos(:,m)
            rr = NORM2(d)
            CALL embed_gaussian_kernel(rr, embed_external_sigma(m), potential, &
                 coef, derivative)
            value = value - embed_external_q(m) * potential
          END DO
          embed_external_potential(idx) = value
        END DO
      END DO
    END DO
    embed_external_count = n_external
    extf => embed_external_potential
    textfld = .TRUE.
    CALL set_cpmdc_external_force_hook(cpmdc_external_charge_density_force)
  END SUBROUTINE

  SUBROUTINE add_external_ion_terms(energy)
    USE ions, ONLY: ions0, ions1
    USE coor, ONLY: tau0, fion
    REAL(real64), INTENT(INOUT) :: energy
    REAL(real64) :: d(3), rr, sigma, potential, coef, derivative, pair_force(3)
    INTEGER :: is, ia, m
    IF (embed_external_count <= 0) RETURN
    DO is = 1, ions1%nsp
      DO ia = 1, ions0%na(is)
        DO m = 1, embed_external_count
          d = tau0(:,ia,is) - embed_external_pos(:,m)
          rr = NORM2(d)
          sigma = embed_external_sigma(m)
          CALL embed_gaussian_kernel(rr, sigma, potential, coef, derivative)
          energy = energy + ions0%zv(is) * embed_external_q(m) * potential
          coef = ions0%zv(is) * embed_external_q(m) * coef
          pair_force = coef * d
          fion(:,ia,is) = fion(:,ia,is) + pair_force
          embed_external_force(:,m) = embed_external_force(:,m) - pair_force
        END DO
      END DO
    END DO
  END SUBROUTINE
#endif

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
    image%charges%valid = 0_c_int
    image%charges%count = 0_c_size_t
    image%charges%values = 0.0_c_double
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
    image%prop%dipole_count = 0_c_size_t
    image%prop%dipole = 0.0_c_double
    image%charges%valid = 0_c_int
    image%charges%count = 0_c_size_t
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

  FUNCTION cpmdc_embed_store_comm(fcomm) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_store_comm')
    INTEGER(c_int), INTENT(IN), VALUE :: fcomm
    INTEGER(c_int) :: ok
    ok = -1_c_int
#if defined(CPMDC_HAS_CPMD)
    BLOCK
      USE mpi
      USE mp_interface, ONLY: mp_comm_world
      INTEGER :: ierr
      LOGICAL :: inited
      ! Rank, size, and the calculator index are decided in C from this
      ! communicator. This routine only installs it, including a replacement
      ! when the C caller has already allowed a rebind.
      ! Open MPI's Fortran handle for MPI_COMM_WORLD is 0, and MPI_COMM_NULL
      ! is a different value, so 0 is a communicator the host can adopt.
      IF (INT(fcomm) == MPI_COMM_NULL) RETURN
      CALL MPI_INITIALIZED(inited, ierr)
      IF (.NOT. inited) RETURN
      mp_comm_world = INT(fcomm)
      embed_calculator_bound = .TRUE.
      ok = 0_c_int
    END BLOCK
#endif
  END FUNCTION cpmdc_embed_store_comm

  FUNCTION cpmdc_embed_adopted_fcomm() RESULT(fcomm) &
      BIND(C, NAME='cpmdc_embed_adopted_fcomm')
    INTEGER(c_int) :: fcomm
    ! -1 means nothing is stored. 0 is a valid Fortran handle: Open MPI uses
    ! it for MPI_COMM_WORLD.
    fcomm = -1_c_int
#if defined(CPMDC_HAS_CPMD)
    BLOCK
      USE mp_interface, ONLY: mp_comm_world
      IF (.NOT. embed_calculator_bound) RETURN
      fcomm = INT(mp_comm_world, c_int)
    END BLOCK
#endif
  END FUNCTION cpmdc_embed_adopted_fcomm

  FUNCTION cpmdc_embed_init() RESULT(ok) BIND(C, NAME='cpmdc_embed_init')
    INTEGER(c_int) :: ok
    runtime_ready = .TRUE.
    runtime_finalized = .FALSE.
    ok = 1_c_int
  END FUNCTION

  FUNCTION cpmdc_embed_available() RESULT(ok) BIND(C, NAME='cpmdc_embed_available')
    INTEGER(c_int) :: ok
    ok = MERGE(1_c_int, 0_c_int, runtime_ready .AND. .NOT. runtime_finalized)
  END FUNCTION

  FUNCTION cpmdc_embed_reset_state(image_c) RESULT(ok) BIND(C, NAME='cpmdc_embed_reset_state')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_reset_warm_orbitals
#endif
    INTERFACE
      SUBROUTINE cpmdc_orbital_ledger_reset() BIND(C, NAME='cpmdc_orbital_ledger_reset')
      END SUBROUTINE cpmdc_orbital_ledger_reset
    END INTERFACE
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    CALL cpmdc_orbital_ledger_reset()
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
    ok = cpmdc_embed_teardown()
  END FUNCTION

  ! Frees the module arrays, communicators and first-call state of the
  ! OpenCPMD setup in this process (tools/opencpmd_embed_teardown.patch),
  ! so the next force call sets CPMD up from its deck again, for any
  ! basis. Without a setup it does nothing. The timers stay: tistart
  ! runs once per process.
  FUNCTION cpmdc_embed_teardown() RESULT(ok) BIND(C, NAME='cpmdc_embed_teardown')
#if defined(CPMDC_HAS_CPMD)
    USE embed_teardown, ONLY: embed_teardown_state
    USE rmas, ONLY: rmass, rmass_t
    USE system, ONLY: parm
    USE efld, ONLY: extf, textfld, set_cpmdc_external_force_hook
#endif
    INTEGER(c_int) :: ok
#if defined(CPMDC_HAS_CPMD)
    LOGICAL :: own_potential
#endif
    ok = 1_c_int
#if defined(CPMDC_HAS_CPMD)
    embed_external_count = 0
    IF (ALLOCATED(embed_external_pos)) DEALLOCATE(embed_external_pos)
    IF (ALLOCATED(embed_external_q)) DEALLOCATE(embed_external_q)
    IF (ALLOCATED(embed_external_sigma)) DEALLOCATE(embed_external_sigma)
    IF (ALLOCATED(embed_external_force)) DEALLOCATE(embed_external_force)
    own_potential = .FALSE.
    IF (ALLOCATED(embed_external_potential)) THEN
      own_potential = ASSOCIATED(extf, embed_external_potential)
      IF (own_potential) NULLIFY(extf)
      DEALLOCATE(embed_external_potential)
      IF (own_potential) textfld = .FALSE.
    END IF
    CALL set_cpmdc_external_force_hook()
    IF (.NOT. cpmd_set_up .AND. .NOT. cpmd_setup_started) RETURN
    CALL embed_teardown_state()
    ! The pseudopotential readers set a species mass only while it is zero,
    ! so the next deck's species would keep this deck's masses, and the
    ! centre of mass of embed_recentre_cluster would be wrong.
    rmass = rmass_t()
    ! setsc writes the minimum-image half lengths for a CELL deck only.
    ! A CELL VECTORS deck after a SYMMETRY 8 one would keep them, and pbc3
    ! would fold its triclinic cell as a monoclinic one of the old lengths
    ! instead of through htm1.
    parm%apbc = 0.0_real64
    CALL clear_embed_basis_latch()
    embed_setup_ht_set = .FALSE.
    cpmd_set_up = .FALSE.
    cpmd_setup_started = .FALSE.
#endif
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

  ! Drops every key's stored orbitals when the next call is on a new
  ! basis, the same test run_embed_scf applies. A keyed evaluation calls
  ! this before it selects its key, so a saturated slot limit cannot
  ! refuse a park on a basis whose old copies are about to be freed.
  ! Returns 1, or 0 when the image is missing or the embed is not ready.
  FUNCTION cpmdc_embed_drop_stale_orbitals(n_atoms, atomic_numbers, cell_ang, &
      has_cell, image_c) RESULT(ok) BIND(C, NAME='cpmdc_embed_drop_stale_orbitals')
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
#if defined(CPMDC_HAS_CPMD)
    TYPE(embed_knobs) :: knobs
#endif
    INTEGER(c_int) :: ok
    ok = 0_c_int
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
#if defined(CPMDC_HAS_CPMD)
    IF (embed_basis_latched .AND. n_atoms > 0) THEN
      knobs = knobs_of(image)
      IF (embed_basis_changed(image, INT(n_atoms), atomic_numbers, cell_ang, &
          INT(has_cell), knobs)) THEN
        IF (.NOT. embed_cell_moves_in_place(INT(n_atoms), atomic_numbers, &
            cell_ang, INT(has_cell), knobs)) CALL drop_stale_basis(image)
      END IF
    END IF
#else
    IF (n_atoms < 0 .OR. has_cell < 0) CONTINUE
#endif
    ok = 1_c_int
  END FUNCTION

  ! The stored c0 follows the key the caller names: the next SCF starts
  ! from the orbitals the last converged call under that key left. A key
  ! with none starts from the orbitals already stored. Every rank of the
  ! calculator selects the same key before the same call.
  SUBROUTINE cpmdc_embed_select_orbitals(key) &
      BIND(C, NAME='cpmdc_embed_select_orbitals')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_select_orbital_slot
#endif
    INTEGER(c_long_long), INTENT(IN), VALUE :: key
#if defined(CPMDC_HAS_CPMD)
    CALL embed_select_orbital_slot(INT(key, KIND=SELECTED_INT_KIND(18)))
#else
    IF (key == 0_c_long_long) CONTINUE
#endif
  END SUBROUTINE

  ! drop_history nonzero empties the extrapolation ring even when no key
  ! is selected, for a caller whose trajectory ends.
  SUBROUTINE cpmdc_embed_clear_orbital_key(drop_history) &
      BIND(C, NAME='cpmdc_embed_clear_orbital_key')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_clear_orbital_key
#endif
    INTEGER(c_int), INTENT(IN), VALUE :: drop_history
#if defined(CPMDC_HAS_CPMD)
    CALL embed_clear_orbital_key(drop_history /= 0_c_int)
#else
    IF (drop_history == 0_c_int) CONTINUE
#endif
  END SUBROUTINE

  ! The converged c0 in the selected extrapolation ring and the count the
  ! last SCF combined. Both are 0 on the reference evaluator.
  SUBROUTINE cpmdc_embed_extrapolation_state(count, used) &
      BIND(C, NAME='cpmdc_embed_extrapolation_state')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_extrapolation_count, &
        embed_extrapolation_used
#endif
    INTEGER(c_int), INTENT(OUT) :: count, used
#if defined(CPMDC_HAS_CPMD)
    count = INT(embed_extrapolation_count(), c_int)
    used = INT(embed_extrapolation_used(), c_int)
#else
    count = 0_c_int
    used = 0_c_int
#endif
  END SUBROUTINE

  ! Reserve a key selection without changing orbital state. C agrees the
  ! result across ranks before calling the matching commit entry point.
  FUNCTION cpmdc_embed_prepare_orbitals(key, fail_allocation) RESULT(rc) &
      BIND(C, NAME='cpmdc_embed_prepare_orbitals')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_prepare_orbital_slot
#endif
    INTEGER(c_long_long), INTENT(IN), VALUE :: key
    INTEGER(c_int), INTENT(IN), VALUE :: fail_allocation
    INTEGER(c_int) :: rc
#if defined(CPMDC_HAS_CPMD)
    rc = INT(embed_prepare_orbital_slot( &
         INT(key, KIND=SELECTED_INT_KIND(18)), fail_allocation /= 0_c_int), &
         c_int)
#else
    rc = 0_c_int
    IF (key == 0_c_long_long .OR. fail_allocation == 0_c_int) CONTINUE
#endif
  END FUNCTION

  SUBROUTINE cpmdc_embed_abort_orbitals() &
      BIND(C, NAME='cpmdc_embed_abort_orbitals')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_abort_orbital_slot
    CALL embed_abort_orbital_slot()
#endif
  END SUBROUTINE

  SUBROUTINE cpmdc_embed_commit_orbitals() &
      BIND(C, NAME='cpmdc_embed_commit_orbitals')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_commit_orbital_slot
    CALL embed_commit_orbital_slot()
#endif
  END SUBROUTINE

  ! Duplicate from_key's converged orbitals into to_key. 0 on success, -1
  ! when from_key has no converged orbitals. Without OpenCPMD this link
  ! stub returns -1, including a same key copy. The reference build's
  ! session copy does not call it and returns 0.
  FUNCTION cpmdc_embed_copy_orbitals(from_key, to_key) RESULT(rc) &
      BIND(C, NAME='cpmdc_embed_copy_orbitals')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_copy_orbital_slot
#endif
    INTEGER(c_long_long), INTENT(IN), VALUE :: from_key, to_key
    INTEGER(c_int) :: rc
#if defined(CPMDC_HAS_CPMD)
    rc = MERGE(0_c_int, -1_c_int, embed_copy_orbital_slot( &
         INT(from_key, KIND=SELECTED_INT_KIND(18)), &
         INT(to_key, KIND=SELECTED_INT_KIND(18))))
#else
    rc = -1_c_int
    IF (from_key == to_key) rc = -1_c_int
#endif
  END FUNCTION

  ! limit > 0 is the most parked orbital copies plus populated
  ! extrapolation history entries, the selected key's entries included.
  ! -1 leaves them unbounded. -2 makes the next park read
  ! CPMDC_ORBITAL_SLOT_LIMIT.
  SUBROUTINE cpmdc_embed_set_orbital_slot_limit(limit) &
      BIND(C, NAME='cpmdc_embed_set_orbital_slot_limit')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_set_orbital_slot_limit
#endif
    INTEGER(c_int), INTENT(IN), VALUE :: limit
#if defined(CPMDC_HAS_CPMD)
    CALL embed_set_orbital_slot_limit(INT(limit))
#else
    IF (limit == 0_c_int) CONTINUE
#endif
  END SUBROUTINE

  ! The C entry gate reads CPMDC_CELL_IN_PLACE on calculator rank zero and
  ! broadcasts it before setting this identical value on every rank.
  SUBROUTINE cpmdc_embed_set_cell_in_place(limit) &
      BIND(C, NAME='cpmdc_embed_set_cell_in_place')
    REAL(c_double), INTENT(IN), VALUE :: limit
    embed_cell_in_place_limit = MAX(0.0_c_double, limit)
  END SUBROUTINE

  ! 0 the last park was allowed, 1 CPMDC_ORBITAL_SLOT_LIMIT is not a
  ! positive integer, 2 the park would have passed that limit and the
  ! stored copies were left unchanged.
  FUNCTION cpmdc_embed_slot_limit_code() RESULT(code) &
      BIND(C, NAME='cpmdc_embed_slot_limit_code')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_orbital_slot_limit_code
#endif
    INTEGER(c_int) :: code
    code = 0_c_int
#if defined(CPMDC_HAS_CPMD)
    code = INT(embed_orbital_slot_limit_code(), KIND=c_int)
#endif
  END FUNCTION

  ! Keys with parked orbitals in this process, the selected key excluded.
  FUNCTION cpmdc_embed_orbital_slots() RESULT(n) &
      BIND(C, NAME='cpmdc_embed_orbital_slots')
#if defined(CPMDC_HAS_CPMD)
    USE rwfopt_utils, ONLY: embed_orbital_slot_count
#endif
    INTEGER(c_int) :: n
    n = 0_c_int
#if defined(CPMDC_HAS_CPMD)
    n = INT(embed_orbital_slot_count(), KIND=c_int)
#endif
  END FUNCTION

  ! Where the last SCF of this process started: 0 the selected key's
  ! stored orbitals, 1 another key's parked orbitals, -1 the guess because
  ! every stored geometry lay beyond the CPMDC_WARM_RMSD limit, -2 the
  ! guess because nothing was stored.
  FUNCTION cpmdc_embed_last_warm_source() RESULT(source) &
      BIND(C, NAME='cpmdc_embed_last_warm_source')
    INTEGER(c_int) :: source
    source = INT(last_warm_source, KIND=c_int)
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
      cell_ang, has_cell, energy_h, grad_h_bohr, n_external, &
      external_positions_ang, external_charges_e, external_widths_ang, &
      external_forces_h_bohr, image_c) RESULT(ok) &
      BIND(C, NAME='cpmdc_embed_energy_grad')
    USE cpmdc_embed_host_iface, ONLY: cpmdc_note_embed_failure
    INTEGER(c_int), INTENT(IN), VALUE :: n_atoms
    REAL(c_double), INTENT(IN) :: positions_ang(*)
    INTEGER(c_int), INTENT(IN) :: atomic_numbers(*)
    REAL(c_double), INTENT(IN) :: cell_ang(*)
    INTEGER(c_int), INTENT(IN), VALUE :: has_cell
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad_h_bohr(*)
    INTEGER(c_int), INTENT(IN), VALUE :: n_external
    REAL(c_double), INTENT(IN) :: external_positions_ang(*), external_charges_e(*)
    REAL(c_double), INTENT(IN) :: external_widths_ang(*)
    REAL(c_double), INTENT(OUT) :: external_forces_h_bohr(*)
    TYPE(c_ptr), INTENT(IN), VALUE :: image_c
    TYPE(cpmdc_embed_image), POINTER :: image
    INTEGER(c_int) :: ok
    INTEGER :: i, n3
    ok = 0_c_int
    energy_h = 0.0_c_double
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(n3) / 3) RETURN
    IF (.NOT. C_ASSOCIATED(image_c)) RETURN
    CALL C_F_POINTER(image_c, image)
    n3 = INT(n_atoms) * 3
    IF (n3 > SIZE(image%prop%hessian)) THEN
      CALL cpmdc_note_embed_failure( &
          'property gradient snapshot supports at most 1365 atoms'//c_null_char)
      RETURN
    END IF
    DO i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    END DO
    DO i = 1, 3 * INT(n_external)
      external_forces_h_bohr(i) = 0.0_c_double
    END DO
    IF (.NOT. runtime_ready .OR. runtime_finalized) RETURN
#if defined(CPMDC_HAS_CPMD)
    CALL run_embed_scf(image, INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, INT(n_external), external_positions_ang, &
         external_charges_e, external_widths_ang, external_forces_h_bohr, ok)
#else
    IF (has_cell < 0) RETURN
    CALL run_reference_pef(image, INT(n_atoms), positions_ang, atomic_numbers, cell_ang, &
         INT(has_cell), energy_h, grad_h_bohr, INT(n_external), external_positions_ang, &
         external_charges_e, external_widths_ang, external_forces_h_bohr, ok)
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
    INTEGER :: m, need, cur_len, stat
    ierr = 1
    m = LEN(s)
    IF (n < 0 .OR. m < 0) RETURN
    IF (m > HUGE(need) - MAX(n, 0)) RETURN
    need = n + m
    cur_len = 0
    IF (ALLOCATED(buf)) cur_len = LEN(buf)
    IF (.NOT. ALLOCATED(buf) .OR. need > cur_len) THEN
      ALLOCATE(CHARACTER(LEN=need) :: grown, STAT=stat)
      IF (stat /= 0) RETURN
      IF (ALLOCATED(buf)) THEN
        IF (n > 0) grown(1:n) = buf(1:n)
      END IF
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
            DO WHILE (i <= n)
              IF (src(i:i) /= NEW_LINE('A') .AND. src(i:i) /= ' ') EXIT
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
    DO WHILE (i <= n)
      IF (src(i:i) /= ' ' .AND. src(i:i) /= ACHAR(9)) EXIT
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
    DO WHILE (pos <= n)
      IF (line(pos:pos) /= ' ' .AND. line(pos:pos) /= ACHAR(9)) EXIT
      pos = pos + 1
    END DO
    IF (pos > n) RETURN
    i = pos
    DO WHILE (i <= n)
      IF (line(i:i) == ' ' .OR. line(i:i) == ACHAR(9)) EXIT
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
    DO WHILE (i <= zend)
      IF (src(i:i) == NEW_LINE('A') .OR. src(i:i) == ACHAR(13)) EXIT
      i = i + 1
    END DO
    k = MIN(i - pos, LEN(line))
    IF (k > 0) line(1:k) = src(pos:pos + k - 1)
    IF (i <= zend) THEN
      IF (src(i:i) == ACHAR(13)) i = i + 1
    END IF
    IF (i <= zend) THEN
      IF (src(i:i) == NEW_LINE('A')) i = i + 1
    END IF
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

  ! 0 is an isolated box. 1 is cubic, 8 orthorhombic.
  ! -2 is a tilted cell: CELL VECTORS and no SYMMETRY line.
  ! The cubic test matches OpenCPMD's 1e-5 check on b/a and c/a.
  INTEGER FUNCTION periodic_symmetry_code(cell, has_cell)
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER(c_int) :: state
    REAL(real64) :: b_over_a, c_over_a, off
    periodic_symmetry_code = 0
    IF (has_cell == 0) RETURN
    state = cpmdc_cell_rank3_state(cell(1:9))
    IF (state == 0_c_int) RETURN
    IF (state /= 1_c_int) THEN
      periodic_symmetry_code = -3
      RETURN
    END IF
    off = ABS(cell(2)) + ABS(cell(3)) + ABS(cell(4)) + ABS(cell(6)) + &
         ABS(cell(7)) + ABS(cell(8))
    IF (off > 0.0_real64 .OR. cell(1) == 0.0_c_double) THEN
      periodic_symmetry_code = -2
      RETURN
    END IF
    b_over_a = REAL(cell(5), KIND=real64) / REAL(cell(1), KIND=real64)
    c_over_a = REAL(cell(9), KIND=real64) / REAL(cell(1), KIND=real64)
    IF (ABS(b_over_a - 1.0_real64) + ABS(c_over_a - 1.0_real64) < &
        1.0e-5_real64) THEN
      periodic_symmetry_code = 1
    ELSE
      periodic_symmetry_code = 8
    END IF
  END FUNCTION

  ! One cold-deck procedure for the stub and the engine.
  ! Cold: OpenCPMD parsers via anonymous memfd deck (no disk write). Warm: C arrays only.
  SUBROUTINE embed_build_cold_deck(n_atoms, pos, z, cell, has_cell, deck, nlen, ierr, knobs)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: deck
    INTEGER, INTENT(OUT) :: nlen, ierr
    INTEGER :: i, j, zz, count, pok, failed, icode
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
    ! The built-in fallback is an isolated Hockney cluster.
    icode = 0
    CALL append('&SYSTEM'//NEW_LINE('A'))
    IF (icode == 0 .OR. icode == 1 .OR. icode == 8) THEN
      CALL append(' SYMMETRY'//NEW_LINE('A'))
      IF (icode == 0) CALL append('  0'//NEW_LINE('A'))
      IF (icode == 1) CALL append('  1'//NEW_LINE('A'))
      IF (icode == 8) CALL append('  8'//NEW_LINE('A'))
    END IF
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
    IF (icode == 0) CALL append(' POISSON SOLVER HOCKNEY'//NEW_LINE('A'))
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
    DO WHILE (ins_at <= nlen)
      IF (deck(ins_at:ins_at) == NEW_LINE('A')) EXIT
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
          ABS(cell(7)) + ABS(cell(8)) > 0.0_c_double) diagonal = .FALSE.
    END IF
    IF (.NOT. diagonal) THEN
      ! An internal write to one CHARACTER variable has one record; a '/'
      ! edit descriptor would end the write. The rows are joined by
      ! NEW_LINE instead.
      WRITE(text, '(2A,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A,' // &
           'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A,' // &
           'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A)') ' CELL VECTORS', &
           NEW_LINE('A'), cell(1), cell(2), cell(3), NEW_LINE('A'), &
           cell(4), cell(5), cell(6), NEW_LINE('A'), &
           cell(7), cell(8), cell(9), NEW_LINE('A')
    ELSE
      WRITE(text, '(A,ES23.16E3,1X,ES23.16E3,1X,ES23.16E3,A)') &
           ' CELL'//NEW_LINE('A')//'  ', &
           cell_a, b_over_a, c_over_a, ' 0 0 0'//NEW_LINE('A')
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
    DO WHILE (ins_at <= nlen)
      IF (deck(ins_at:ins_at) == NEW_LINE('A')) EXIT
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

  ! True when a SYMMETRY keyword is present and it is not CHECK SYMMETRY.
  LOGICAL FUNCTION deck_names_symmetry(deck, nlen)
    CHARACTER(LEN=*), INTENT(IN) :: deck
    INTEGER, INTENT(IN) :: nlen
    INTEGER :: p
    deck_names_symmetry = .FALSE.
    IF (nlen < 8) RETURN
    p = 1
    DO WHILE (p <= nlen - 7)
      IF (deck(p:p+7) == 'SYMMETRY') THEN
        IF (p >= 7) THEN
          IF (deck(p-6:p-1) == 'CHECK ') THEN
            p = p + 1
            CYCLE
          END IF
        END IF
        deck_names_symmetry = .TRUE.
        RETURN
      END IF
      p = p + 1
    END DO
  END FUNCTION

  ! True when a line of &SYSTEM has CLUSTER, in any case, as its first
  ! token. CHECK CLUSTER is not CLUSTER, as CHECK SYMMETRY is not SYMMETRY,
  ! and the comment text around the sections and other sections do not
  ! count.
  LOGICAL FUNCTION deck_names_cluster(deck, nlen)
    CHARACTER(LEN=*), INTENT(IN) :: deck
    INTEGER, INTENT(IN) :: nlen
    CHARACTER(LEN=80) :: head, tok
    INTEGER :: n, ls, le, pos
    LOGICAL :: in_sys
    deck_names_cluster = .FALSE.
    n = MIN(nlen, LEN(deck))
    in_sys = .FALSE.
    ls = 1
    DO WHILE (ls <= n)
      le = ls
      DO WHILE (le <= n)
        IF (deck(le:le) == NEW_LINE('A')) EXIT
        le = le + 1
      END DO
      CALL embed_upper_copy(deck(ls:le - 1), head)
      pos = 1
      CALL embed_next_tok(head, pos, tok)
      IF (tok(1:7) == '&SYSTEM') THEN
        in_sys = .TRUE.
      ELSE IF (in_sys .AND. tok(1:4) == '&END') THEN
        in_sys = .FALSE.
      ELSE IF (in_sys .AND. TRIM(tok) == 'CLUSTER') THEN
        deck_names_cluster = .TRUE.
        RETURN
      END IF
      ls = le + 1
    END DO
  END FUNCTION

  ! A deck with no symmetry line takes cubic, orthorhombic, or isolated
  ! symmetry from the box. An explicit symmetry, Poisson solver, cluster,
  ! or isolated-molecule keyword stays as written. CELL VECTORS is the
  ! tilted cell and gets no SYMMETRY line.
  SUBROUTINE inject_periodicity_if_unset(deck, nlen, cell, has_cell, ierr)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_note_embed_failure
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: isys, ins_at, code, iend
    CHARACTER(LEN=32) :: symtxt
    CHARACTER(LEN=32) :: hock
    ierr = 0
    IF (.NOT. ALLOCATED(deck) .OR. nlen < 1) RETURN
    IF (INDEX(deck(1:nlen), 'POISSON SOLVER') > 0) RETURN
    IF (deck_names_cluster(deck, nlen)) RETURN
    IF (INDEX(deck(1:nlen), ' ISOLATED MOLECULE') > 0) RETURN
    IF (INDEX(deck(1:nlen), ' MOLECULE ISOLATED') > 0) RETURN
    IF (deck_names_symmetry(deck, nlen)) RETURN
    IF (INDEX(deck(1:nlen), 'CELL VECTORS') > 0) RETURN
    code = periodic_symmetry_code(cell, has_cell)
    IF (code == -3) THEN
      CALL cpmdc_note_embed_failure('the cell does not have a representable ' // &
           'rank-3 volume' // c_null_char)
      ierr = 3
      RETURN
    END IF
    IF (code < 0) RETURN
    isys = INDEX(deck(1:nlen), '&SYSTEM')
    IF (isys == 0) isys = INDEX(deck(1:nlen), '&system')
    IF (isys <= 0) RETURN
    ins_at = isys
    DO WHILE (ins_at <= nlen)
      IF (deck(ins_at:ins_at) == NEW_LINE('A')) EXIT
      ins_at = ins_at + 1
    END DO
    IF (ins_at <= nlen) ins_at = ins_at + 1
    IF (code == 0) THEN
      symtxt = ' SYMMETRY'//NEW_LINE('A')//'  0'//NEW_LINE('A')
    ELSE IF (code == 1) THEN
      symtxt = ' SYMMETRY'//NEW_LINE('A')//'  1'//NEW_LINE('A')
    ELSE IF (code == 8) THEN
      symtxt = ' SYMMETRY'//NEW_LINE('A')//'  8'//NEW_LINE('A')
    ELSE
      RETURN
    END IF
    CALL insert_at(ins_at, symtxt)
    IF (ierr /= 0 .OR. code /= 0) RETURN
    iend = INDEX(deck(isys:nlen), '&END')
    IF (iend == 0) iend = INDEX(deck(isys:nlen), '&end')
    IF (iend <= 0) RETURN
    iend = isys + iend - 1
    hock = ' POISSON SOLVER HOCKNEY'//NEW_LINE('A')
    CALL insert_at(iend, hock)
  CONTAINS
    SUBROUTINE insert_at(at, text)
      INTEGER, INTENT(IN) :: at
      CHARACTER(LEN=*), INTENT(IN) :: text
      CHARACTER(LEN=:), ALLOCATABLE :: tmp
      INTEGER :: m, k, stat, newn
      m = LEN_TRIM(text)
      IF (m < 1) RETURN
      IF (at < 1 .OR. at > nlen + 1) THEN
        ierr = 1
        RETURN
      END IF
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
      IF (at > 1) tmp(1:at-1) = deck(1:at-1)
      tmp(at:at+m-1) = text(1:m)
      IF (at <= nlen) THEN
        k = nlen - at + 1
        tmp(at+m:at+m+k-1) = deck(at:nlen)
      END IF
      CALL MOVE_ALLOC(tmp, deck)
      nlen = newn
    END SUBROUTINE
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
    LOGICAL :: use_box
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
      ! The caller's box is the cell. A deck CELL is replaced by it. An
      ! all-zero box is a system without a cell (eOn's non-periodic
      ! convention), so the deck's CELL, or the default, stays. A NaN
      ! component selects impose_box_cell so that routine rejects it.
      use_box = .FALSE.
      IF (has_cell /= 0) use_box = ANY(ABS(cell(1:9)) > 0.0_c_double) .OR. &
          ANY(IEEE_IS_NAN(cell(1:9)))
      IF (use_box) THEN
        CALL impose_box_cell(deck, nlen, cell, ierr)
      ELSE
        CALL inject_cell_if_missing(deck, nlen, cell, has_cell, ierr)
      END IF
      IF (ierr == 0) CALL inject_periodicity_if_unset(deck, nlen, cell, &
           has_cell, ierr)
      IF (ierr == 0) CALL inject_maxiter_if_missing(deck, nlen, ierr)
    END IF
  END SUBROUTINE

  ! The box of the call (Angstrom, rows are lattice vectors) is the cell
  ! CPMD uses, whatever CELL the deck states. Every CELL line of &SYSTEM
  ! and its data lines go, and so does SYMMETRY with its value line; then
  ! the box goes back in before the &END of &SYSTEM:
  !   isolated deck (SYMMETRY 0 or ISOLATED, CLUSTER, POISSON SOLVER,
  !     ISOLATED MOLECULE, or MOLECULE ISOLATED in &SYSTEM): SYMMETRY 0 and
  !     CELL, for a diagonal box; CLUSTER and CELL VECTORS for another box,
  !     because CPMD permits a non-orthorhombic isolated cell only with
  !     CLUSTER. A deck that already has CLUSTER keeps its one line.
  !   diagonal box: SYMMETRY 1 when the deck said cubic and a = b = c,
  !     otherwise SYMMETRY 8, with CELL a b/a c/a.
  !   other box: CELL VECTORS and no SYMMETRY, which CPMD refuses beside it.
  ! Without ANGSTROM in &SYSTEM the deck reads lengths in Bohr, so the box
  ! is written in Bohr (CPMD's fbohr, 1/0.529177210859).
  SUBROUTINE impose_box_cell(deck, nlen, cell, ierr)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_note_embed_failure
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: deck
    INTEGER, INTENT(INOUT) :: nlen
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(OUT) :: ierr
    REAL(real64), PARAMETER :: ang_per_bohr = 0.529177210859_real64
    CHARACTER(LEN=:), ALLOCATABLE :: out, symmetry_indent
    CHARACTER(LEN=80) :: head, tok, symtok
    CHARACTER(LEN=600) :: text
    REAL(real64) :: c(9), scale
    INTEGER :: ls, le, nxt, skip, pos, ntext, sys_end, osym, k, cell_state
    LOGICAL :: in_sys, angstrom, cubic_deck, isolated_deck, diagonal, done
    LOGICAL :: has_cluster, isolated_keyword
    ierr = 0
    IF (.NOT. ALLOCATED(deck) .OR. nlen < 1) RETURN
    ! First pass: units and the deck's symmetry, inside &SYSTEM only.
    in_sys = .FALSE.
    angstrom = .FALSE.
    cubic_deck = .FALSE.
    isolated_deck = .FALSE.
    has_cluster = deck_names_cluster(deck, nlen)
    isolated_keyword = .FALSE.
    skip = 0
    osym = 0
    ls = 1
    DO WHILE (ls <= nlen)
      le = ls
      DO WHILE (le <= nlen)
        IF (deck(le:le) == NEW_LINE('A')) EXIT
        le = le + 1
      END DO
      CALL embed_upper_copy(deck(ls:MIN(le - 1, nlen)), head)
      pos = 1
      CALL embed_next_tok(head, pos, tok)
      IF (osym == 1) THEN
        ! The value line of SYMMETRY.
        CALL embed_ltrim(head, symtok)
        isolated_deck = symtok(1:1) == '0' .OR. INDEX(symtok, 'ISOLATED') > 0
        cubic_deck = (symtok(1:1) == '1' .AND. LEN_TRIM(symtok) == 1) .OR. &
            (INDEX(symtok, 'CUBIC') > 0 .AND. INDEX(symtok, 'FACE') == 0 &
             .AND. INDEX(symtok, 'BODY') == 0)
        osym = 2
      ELSE IF (tok(1:7) == '&SYSTEM') THEN
        in_sys = .TRUE.
      ELSE IF (in_sys .AND. tok(1:4) == '&END') THEN
        in_sys = .FALSE.
      ELSE IF (in_sys .AND. TRIM(tok) == 'ANGSTROM') THEN
        angstrom = .TRUE.
      ELSE IF (in_sys .AND. TRIM(tok) == 'SYMMETRY' .AND. osym == 0) THEN
        osym = 1
      ELSE IF (in_sys .AND. (INDEX(head, 'POISSON SOLVER') > 0 .OR. &
               INDEX(head, 'ISOLATED MOLECULE') > 0 .OR. &
               INDEX(head, 'MOLECULE ISOLATED') > 0)) THEN
        isolated_keyword = .TRUE.
      END IF
      ls = le + 1
    END DO
    IF (has_cluster .OR. isolated_keyword) isolated_deck = .TRUE.
    DO k = 1, 9
      c(k) = REAL(cell(k), KIND=real64)
    END DO
    cell_state = cpmdc_cell_rank3_state(cell(1:9))
    IF (cell_state /= 1_c_int) THEN
      SELECT CASE (cell_state)
      CASE (-3_c_int)
        CALL cpmdc_note_embed_failure( &
             'cell volume is too small for binary64' // c_null_char)
      CASE (-2_c_int)
        CALL cpmdc_note_embed_failure( &
             'cell contains a non-finite component' // c_null_char)
      CASE DEFAULT
        CALL cpmdc_note_embed_failure( &
             'cell does not have three independent vectors' // c_null_char)
      END SELECT
      ierr = 3
      RETURN
    END IF
    diagonal = ABS(c(2)) + ABS(c(3)) + ABS(c(4)) + ABS(c(6)) + ABS(c(7)) + &
        ABS(c(8)) == 0.0_real64
    scale = 1.0_real64
    IF (.NOT. angstrom) scale = 1.0_real64 / ang_per_bohr
    c = c * scale
    IF (diagonal) THEN
      IF (c(1) <= 0.0_real64 .OR. c(5) <= 0.0_real64 .OR. c(9) <= 0.0_real64) THEN
        CALL cpmdc_note_embed_failure('the ForceInput box has a ' // &
             'non-positive diagonal' // c_null_char)
        ierr = 3
        RETURN
      END IF
      IF (isolated_deck) THEN
        k = 0
      ELSE IF ((cubic_deck .OR. osym == 0) .AND. ABS(c(5) - c(1)) <= 1.0e-10_real64 * c(1) &
               .AND. ABS(c(9) - c(1)) <= 1.0e-10_real64 * c(1)) THEN
        k = 1
      ELSE
        k = 8
      END IF
      symmetry_indent = ''
      IF (osym == 0) symmetry_indent = '  '
      WRITE(text, '(A,A,A,I0,A,A,A,ES23.16E3,1X,ES23.16E3,1X,ES23.16E3,A,A)') &
           ' SYMMETRY', NEW_LINE('A'), &
           symmetry_indent, k, NEW_LINE('A'), ' CELL', NEW_LINE('A'), &
           c(1), c(5) / c(1), &
           c(9) / c(1), ' 0 0 0', NEW_LINE('A')
    ELSE
      IF (isolated_deck .AND. .NOT. has_cluster) THEN
        text = ' CLUSTER' // NEW_LINE('A')
      ELSE
        text = ''
      END IF
      ntext = LEN_TRIM(text)
      WRITE(text(ntext + 1:), '(2A,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A,' // &
           'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A,' // &
           'ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,A)') ' CELL VECTORS', &
           NEW_LINE('A'), c(1), c(2), c(3), NEW_LINE('A'), &
           c(4), c(5), c(6), NEW_LINE('A'), c(7), c(8), c(9), NEW_LINE('A')
    END IF
    ntext = LEN_TRIM(text)
    IF (text(ntext:ntext) /= NEW_LINE('A')) THEN
      ntext = ntext + 1
      text(ntext:ntext) = NEW_LINE('A')
    END IF
    ! Second pass: copy every line but CELL and SYMMETRY of &SYSTEM, and
    ! put the box in before its &END.
    out = ''
    in_sys = .FALSE.
    done = .FALSE.
    skip = 0
    sys_end = 0
    ls = 1
    DO WHILE (ls <= nlen)
      le = ls
      DO WHILE (le <= nlen)
        IF (deck(le:le) == NEW_LINE('A')) EXIT
        le = le + 1
      END DO
      nxt = MIN(le, nlen)
      CALL embed_upper_copy(deck(ls:MIN(le - 1, nlen)), head)
      pos = 1
      CALL embed_next_tok(head, pos, tok)
      IF (skip > 0) THEN
        skip = skip - 1
      ELSE IF (tok(1:7) == '&SYSTEM') THEN
        in_sys = .TRUE.
        out = out // deck(ls:nxt)
      ELSE IF (in_sys .AND. tok(1:4) == '&END') THEN
        in_sys = .FALSE.
        out = out // text(1:ntext) // deck(ls:nxt)
        done = .TRUE.
      ELSE IF (in_sys .AND. TRIM(tok) == 'CELL') THEN
        skip = 1
        IF (INDEX(head, 'VECTORS') > 0) skip = 3
      ELSE IF (in_sys .AND. TRIM(tok) == 'SYMMETRY') THEN
        skip = 1
      ELSE
        out = out // deck(ls:nxt)
      END IF
      ls = le + 1
    END DO
    IF (.NOT. done) THEN
      ! No &SYSTEM section to carry the cell: leave the deck as it was.
      RETURN
    END IF
    sys_end = LEN(out)
    CALL MOVE_ALLOC(out, deck)
    nlen = sys_end
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
  SUBROUTINE run_reference_pef(image, n_atoms, pos, z, cell, has_cell, &
      energy_h, grad, n_external, external_pos, external_q, &
      external_width, external_force, ok)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_note_embed_failure, &
        cpmdc_note_periodic_embedding
    TYPE(embed_knobs) :: knobs
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*)
    INTEGER, INTENT(IN) :: n_external
    REAL(c_double), INTENT(IN) :: external_pos(*), external_q(*), &
        external_width(*)
    REAL(c_double), INTENT(OUT) :: external_force(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    REAL(real64), PARAMETER :: bohr_to_ang = 0.529177210903_real64
    REAL(real64) :: k, r_bohr, coord_bohr, z_scale, deck_scale
    REAL(real64) :: d(3), field(3), jac(3,3), force(3)
    REAL(real64) :: sigma, q, rr, potential, a, aprime, alpha
    INTEGER :: i, j, idx, m, u, v, stress_rc, induction_rc
    LOGICAL :: periodic_induction
    INTERFACE
      FUNCTION cpmdc_reference_periodic_induction(n_atoms, pos, znum, n_ext, &
          ext_pos, ext_q, ext_width, cell_in, energy, grad_out, ext_force) &
          BIND(C, NAME='cpmdc_reference_periodic_induction') RESULT(rc)
        IMPORT :: c_int, c_double
        INTEGER(c_int), VALUE :: n_atoms, n_ext
        REAL(c_double), INTENT(IN) :: pos(*), ext_pos(*), ext_q(*), ext_width(*)
        REAL(c_double), INTENT(IN) :: cell_in(*)
        INTEGER(c_int), INTENT(IN) :: znum(*)
        REAL(c_double), INTENT(INOUT) :: energy, grad_out(*), ext_force(*)
        INTEGER(c_int) :: rc
      END FUNCTION
    END INTERFACE
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
    ! The reference evaluator carries a small isotropic polarizable density
    ! at each QM centre. It is not a chemical method; it gives the public
    ! electrostatic-embedding ABI an energy-conserving analytic oracle.
    periodic_induction = .FALSE.
    IF (n_external > 0 .AND. &
        .NOT. reference_deck_isolated(knobs%input_deck, cell, has_cell)) THEN
      IF (.NOT. reference_periodic_embedding() .OR. has_cell == 0) THEN
        CALL cpmdc_note_periodic_embedding()
        RETURN
      END IF
      induction_rc = cpmdc_reference_periodic_induction(INT(n_atoms, c_int), &
           pos, z, INT(n_external, c_int), external_pos, external_q, &
           external_width, cell, energy_h, grad, external_force)
      IF (induction_rc == 2_c_int) THEN
        CALL cpmdc_note_embed_failure( &
            'periodic embedding did not converge'//c_null_char)
        RETURN
      END IF
      IF (induction_rc /= 0_c_int) THEN
        CALL cpmdc_note_embed_failure( &
            'periodic embedding cell is not usable'//c_null_char)
        RETURN
      END IF
      periodic_induction = .TRUE.
    END IF
    IF (.NOT. periodic_induction) THEN
    DO i = 1, n_atoms
      field = 0.0_real64
      alpha = 0.20_real64 * REAL(MAX(1, INT(z(i))), KIND=real64)
      DO m = 1, n_external
        d = REAL(pos(3*i-2:3*i) - external_pos(3*m-2:3*m), KIND=real64) / &
            bohr_to_ang
        sigma = REAL(external_width(m), KIND=real64) / bohr_to_ang
        rr = NORM2(d)
        CALL embed_gaussian_kernel(rr, sigma, potential, a, aprime)
        q = REAL(external_q(m), KIND=real64)
        field = field + q * a * d
      END DO
      energy_h = energy_h - 0.5_real64 * alpha * DOT_PRODUCT(field, field)
      DO m = 1, n_external
        d = REAL(pos(3*i-2:3*i) - external_pos(3*m-2:3*m), KIND=real64) / &
            bohr_to_ang
        sigma = REAL(external_width(m), KIND=real64) / bohr_to_ang
        rr = NORM2(d)
        CALL embed_gaussian_kernel(rr, sigma, potential, a, aprime)
        q = REAL(external_q(m), KIND=real64)
        jac = 0.0_real64
        DO u = 1, 3
          jac(u,u) = q * a
          IF (rr > 0.0_real64) THEN
            DO v = 1, 3
              jac(u,v) = jac(u,v) + q * aprime * d(u) * d(v) / rr
            END DO
          END IF
        END DO
        force = alpha * MATMUL(TRANSPOSE(jac), field)
        grad(3*i-2:3*i) = grad(3*i-2:3*i) - REAL(force, KIND=c_double)
        external_force(3*m-2:3*m) = external_force(3*m-2:3*m) - &
             REAL(force, KIND=c_double)
      END DO
    END DO
    END IF
    ok = 1_c_int
    ! Reference PEF: full POD surface with etot-only ener_com + ENERGY row + PROP.
    CALL snapshot_total_only(image, energy_h)
    DO i = 1, n_atoms * 3
      image%prop%hessian(i) = grad(i)
    END DO
    image%prop%hessian_count = INT(n_atoms * 3, KIND=c_size_t)
    image%prop%valid = 1_c_int
    ! Strain derivative of this reference energy, Hartree/Bohr^3, with the
    ! OpenCPMD sign: volume * stress = -dE/d(eps). An isolated deck has a
    ! box volume and no tensor. CPMDC_STRESS=0 skips the tensor too.
    ! A periodic induction step leaves it invalid: this fill is the harmonic
    ! reference, not the Ewald field.
    image%stress%values = 0.0_c_double
    image%stress%valid = 0_c_int
    IF (reference_stress_wanted() .AND. has_cell /= 0 .AND. &
        .NOT. periodic_induction) THEN
      IF (ANY(ABS(cell(1:9)) > 0.0_c_double) .AND. &
           .NOT. reference_deck_isolated(knobs%input_deck, cell, has_cell)) THEN
        stress_rc = reference_pef_fill_stress(image, n_atoms, pos, z, cell, k)
        IF (stress_rc > 0) THEN
          image%stress%valid = 1_c_int
        ELSE
          ok = 0_c_int
          energy_h = 0.0_c_double
          grad(1:3 * n_atoms) = 0.0_c_double
          CALL clear_last_energy_components(image)
          IF (stress_rc < 0) THEN
            CALL cpmdc_note_embed_failure('cell volume is too large for ' // &
                 'binary64' // c_null_char)
          ELSE
            CALL cpmdc_note_embed_failure('cell volume is not ' // &
                 'representable for stress' // c_null_char)
          END IF
        END IF
      END IF
    END IF
  END SUBROUTINE

  LOGICAL FUNCTION reference_stress_wanted()
    CHARACTER(LEN=32) :: v
    INTEGER :: st
    reference_stress_wanted = .TRUE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_STRESS', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '0') reference_stress_wanted = .FALSE.
  END FUNCTION

  LOGICAL FUNCTION reference_periodic_embedding()
    CHARACTER(LEN=32) :: v
    INTEGER :: st
    reference_periodic_embedding = .FALSE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_PERIODIC_EMBEDDING', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '1') reference_periodic_embedding = .TRUE.
  END FUNCTION

  ! Symmetry 0 or ISOLATED, a Poisson solver, CLUSTER in &SYSTEM, or an
  ! isolated-molecule keyword is isolated. CHECK CLUSTER is not CLUSTER
  ! and CHECK SYMMETRY is not the symmetry code. Positive SYMMETRY,
  ! SURFACE, POLYMER, or CELL VECTORS without an isolated keyword is
  ! periodic. With periodicity unset, a positive-volume cell is periodic
  ! and a zero or absent box stays isolated.
  LOGICAL FUNCTION reference_deck_isolated(deck, cell, has_cell)
    CHARACTER(LEN=*), INTENT(IN) :: deck
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER, INTENT(IN) :: has_cell
    INTEGER :: ls, le, n, pos
    CHARACTER(LEN=160) :: head, trimmed, tok, symtok
    LOGICAL :: in_sys, expect_sym, sym_isolated, sym_periodic
    LOGICAL :: has_isolated_kw, has_periodic_kw, has_cell_vectors
    n = LEN_TRIM(deck)
    in_sys = .FALSE.
    expect_sym = .FALSE.
    sym_isolated = .FALSE.
    sym_periodic = .FALSE.
    has_isolated_kw = deck_names_cluster(deck, n)
    has_periodic_kw = .FALSE.
    has_cell_vectors = .FALSE.
    ls = 1
    DO WHILE (ls <= n)
      le = ls
      DO WHILE (le <= n)
        IF (deck(le:le) == NEW_LINE('A')) EXIT
        le = le + 1
      END DO
      CALL embed_upper_copy(deck(ls:MIN(le - 1, n)), head)
      CALL embed_ltrim(head, trimmed)
      IF (LEN_TRIM(trimmed) > 0 .AND. trimmed(1:1) /= '!') THEN
        pos = 1
        CALL embed_next_tok(trimmed, pos, tok)
        IF (TRIM(tok) == '&SYSTEM') THEN
          in_sys = .TRUE.
          expect_sym = .FALSE.
        ELSE IF (in_sys .AND. TRIM(tok) == '&END') THEN
          in_sys = .FALSE.
          expect_sym = .FALSE.
        ELSE IF (in_sys .AND. expect_sym) THEN
          IF (TRIM(tok) == '0' .OR. tok(1:8) == 'ISOLATED') THEN
            sym_isolated = .TRUE.
          ELSE
            sym_periodic = .TRUE.
          END IF
          expect_sym = .FALSE.
        ELSE IF (in_sys) THEN
          IF (TRIM(tok) == 'SYMMETRY') THEN
            CALL embed_next_tok(trimmed, pos, symtok)
            IF (LEN_TRIM(symtok) > 0 .AND. symtok(1:1) /= '!') THEN
              IF (TRIM(symtok) == '0' .OR. symtok(1:8) == 'ISOLATED') THEN
                sym_isolated = .TRUE.
              ELSE
                sym_periodic = .TRUE.
              END IF
            ELSE
              expect_sym = .TRUE.
            END IF
          ELSE IF (INDEX(trimmed, 'POISSON SOLVER') == 1 .OR. &
                   INDEX(trimmed, 'ISOLATED MOLECULE') == 1 .OR. &
                   INDEX(trimmed, 'MOLECULE ISOLATED') == 1) THEN
            has_isolated_kw = .TRUE.
          ELSE IF (TRIM(tok) == 'SURFACE' .OR. TRIM(tok) == 'POLYMER') THEN
            has_periodic_kw = .TRUE.
          ELSE IF (INDEX(trimmed, 'CELL VECTORS') == 1) THEN
            has_cell_vectors = .TRUE.
          END IF
        END IF
      END IF
      ls = le + 1
    END DO
    IF (sym_periodic .OR. has_periodic_kw) THEN
      reference_deck_isolated = .FALSE.
    ELSE IF (sym_isolated .OR. has_isolated_kw) THEN
      reference_deck_isolated = .TRUE.
    ELSE IF (has_cell_vectors) THEN
      reference_deck_isolated = .FALSE.
    ELSE
      reference_deck_isolated = (periodic_symmetry_code(cell, has_cell) == 0)
    END IF
  END FUNCTION

  INTEGER FUNCTION reference_pef_fill_stress(image, n_atoms, pos, z, cell, k)
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*), cell(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(real64), INTENT(IN) :: k
    REAL(real64), PARAMETER :: bohr_to_ang = 0.529177210903_real64
    REAL(real64) :: h(3, 3), r(3), ded(3, 3)
    REAL(real64) :: vol, inv_b, z_scale
    REAL(c_double) :: cell_bohr(9), vol_c
    INTEGER(c_int) :: vol_status
    INTEGER :: i, j, ia, vec, slot
    reference_pef_fill_stress = 0
    inv_b = 1.0_real64 / bohr_to_ang
    DO i = 1, 3
      h(i, 1) = REAL(cell(i), KIND=real64)
      h(i, 2) = REAL(cell(3 + i), KIND=real64)
      h(i, 3) = REAL(cell(6 + i), KIND=real64)
    END DO
    DO i = 1, 9
      cell_bohr(i) = REAL(REAL(cell(i), KIND=real64) * inv_b, KIND=c_double)
    END DO
    vol_status = cpmdc_cell_volume3(cell_bohr, vol_c)
    IF (vol_status /= 0_c_int) THEN
      reference_pef_fill_stress = -1
      RETURN
    END IF
    vol = REAL(vol_c, KIND=real64)
    IF (.NOT. (vol > 0.0_real64)) RETURN
    ded = 0.0_real64
    DO ia = 1, n_atoms
      z_scale = REAL(MAX(1, INT(z(ia))), KIND=real64)
      DO i = 1, 3
        r(i) = REAL(pos(3 * (ia - 1) + i), KIND=real64)
      END DO
      DO i = 1, 3
        DO j = 1, 3
          ded(i, j) = ded(i, j) + k * z_scale * r(i) * r(j) * inv_b * inv_b
        END DO
      END DO
    END DO
    DO vec = 1, 3
      DO i = 1, 3
        DO j = 1, 3
          ded(i, j) = ded(i, j) + 2.0e-10_real64 * h(i, vec) * h(j, vec)
        END DO
      END DO
    END DO
    slot = 0
    DO i = 1, 3
      DO j = 1, 3
        slot = slot + 1
        image%stress%values(slot) = REAL(-ded(i, j) / vol, KIND=c_double)
      END DO
    END DO
    reference_pef_fill_stress = 1
  END FUNCTION

#endif

#if defined(CPMDC_HAS_CPMD)
  SUBROUTINE snapshot_prop_from_modules(image, n_atoms, grad, origin)
    USE wfopts_utils, ONLY: embed_get_dipole, embed_get_charges
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: grad(*)
    INTEGER(c_int), INTENT(IN) :: origin(*)
    INTEGER :: i, n, slot, j, astat
    LOGICAL :: dipole_valid, charges_valid
    LOGICAL, ALLOCATABLE :: seen(:)
    REAL(real64) :: dipole(3)
    REAL(real64), ALLOCATABLE :: species_q(:)
    image%prop%valid = 1_c_int
    image%prop%dipole_count = 0_c_size_t
    image%prop%dipole = 0.0_c_double
    CALL embed_get_dipole(dipole, dipole_valid)
    IF (dipole_valid) THEN
      image%prop%dipole_count = 3_c_size_t
      DO i = 1, 3
        image%prop%dipole(i) = REAL(dipole(i), KIND=c_double)
      END DO
    END IF
    image%prop%polarizability_count = 9_c_size_t
    DO i = 1, 9
      image%prop%polarizability(i) = 0.0_c_double
    END DO
    DO i = 1, n_atoms * 3
      image%prop%hessian(i) = grad(i)
    END DO
    image%prop%hessian_count = INT(n_atoms * 3, KIND=c_size_t)
    image%charges%valid = 0_c_int
    image%charges%count = 0_c_size_t
    image%charges%values = 0.0_c_double
    IF (n_atoms < 1 .OR. n_atoms > 1365) RETURN
    ALLOCATE(species_q(n_atoms), seen(n_atoms), STAT=astat)
    IF (astat /= 0) RETURN
    seen = .FALSE.
    CALL embed_get_charges(species_q, n, charges_valid, n_atoms)
    IF (.NOT. charges_valid .OR. n /= n_atoms) THEN
      DEALLOCATE(species_q, seen)
      RETURN
    END IF
    DO slot = 1, n_atoms
      j = INT(origin(slot)) + 1
      IF (j < 1 .OR. j > n_atoms .OR. seen(j)) THEN
        image%charges%values = 0.0_c_double
        DEALLOCATE(species_q, seen)
        RETURN
      END IF
      seen(j) = .TRUE.
      image%charges%values(j) = REAL(species_q(slot), KIND=c_double)
    END DO
    image%charges%valid = 1_c_int
    image%charges%count = INT(n_atoms, KIND=c_size_t)
    DEALLOCATE(species_q, seen)
  END SUBROUTINE

  SUBROUTINE embed_set_tau0_from_pos(n_atoms, pos, z, origin, recentre, ierr)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_species_order_map, cpmdc_note_embed_failure
    USE coor, ONLY: tau0
    USE ions, ONLY: ions0, ions1
    USE cnst, ONLY: fbohr
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    INTEGER(c_int), INTENT(OUT) :: origin(*)
    TYPE(embed_recentring), INTENT(OUT) :: recentre
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
    CALL embed_recentre_cluster(n_atoms, origin, recentre, ierr)
  END SUBROUTINE

  ! An isolated deck (SYMMETRY 0, a Poisson solver, CLUSTER, SURFACE or
  ! POLYMER) describes its density on a box of fixed size with its origin
  ! at a corner. A region that drifts toward a face of that box puts
  ! density across it, and the energy, forces, and any virial from them
  ! go wrong: one BLYP water whose outer hydrogen is 0.5 Angstrom from a
  ! face of the 10 Angstrom Hockney box is off by 5e-4 Hartree, and by
  ! 8e-3 Hartree with a net force of 0.17 Hartree/Bohr when it is 0.5
  ! Angstrom past it.
  ! Each call therefore moves the atoms rigidly so that their centre of
  ! mass (CPMDC_RECENTER=mass, the default) or of geometry
  ! (CPMDC_RECENTER=geometry) is at the centre of the box along every
  ! non-periodic lattice vector; CPMDC_RECENTER=off keeps the caller's
  ! coordinates. The energy is then a function of the shifted coordinates
  ! y_i = x_i + s(x), and its gradient with respect to the caller's x_i is
  ! g_i - w_i P sum_j g_j, with w the weights of that centre and P the
  ! projector on the shifted directions; embed_project_gradient applies it.
  ! That is the residual net force of the grid removed, so energy, forces,
  ! and sum_i x_i (x) f_i do not change when the caller translates the
  ! atoms. CPMDC_RECENTER_MARGIN (Angstrom, default 0) is the least
  ! distance an atom may keep from a face of the box after the shift; a
  ! closer atom fails the call. A negative margin turns the check off.
  SUBROUTINE embed_recentre_cluster(n_atoms, origin, recentre, ierr)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_note_embed_failure
    USE coor, ONLY: tau0
    USE ions, ONLY: ions0, ions1
    USE isos, ONLY: isos1, isos3
    USE metr, ONLY: metr_com
    USE rmas, ONLY: rmass
    USE cnst, ONLY: fbohr
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: origin(*)
    TYPE(embed_recentring), INTENT(OUT) :: recentre
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: policy, is, ia, slot, j, k, m, astat
    LOGICAL :: free_axis(3)
    REAL(real64) :: margin, total, centre(3), frac(3), shift(3), height, gap, closest
    CHARACTER(LEN=256) :: msg
    CHARACTER(LEN=12) :: num(2)
    ierr = 1
    recentre%applied = .FALSE.
    recentre%shift_ang = 0.0_real64
    recentre%proj = 0.0_real64
    policy = recentre_policy()
    margin = recentre_margin()
    IF (policy < 0 .OR. margin /= margin) THEN
      CALL cpmdc_note_embed_failure('CPMDC_RECENTER is mass, geometry, or off, and '// &
           'CPMDC_RECENTER_MARGIN a number of Angstrom'//c_null_char)
      RETURN
    END IF
    ierr = 0
    IF (.NOT. isos1%tclust) RETURN
    free_axis = .TRUE.
    IF (isos1%toned) free_axis(1) = .FALSE.
    IF (isos1%ttwod) THEN
      free_axis = .FALSE.
      free_axis(isos3%snormal) = .TRUE.
    END IF
    IF (policy > 0) THEN
      ALLOCATE(recentre%weights(n_atoms), STAT=astat)
      IF (astat /= 0) THEN
        CALL cpmdc_note_embed_failure('out of memory'//c_null_char)
        ierr = 1
        RETURN
      END IF
      total = 0.0_real64
      centre = 0.0_real64
      slot = 0
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          slot = slot + 1
          j = INT(origin(slot)) + 1
          recentre%weights(j) = 1.0_real64
          IF (policy == 2) recentre%weights(j) = rmass%pma0(is)
          total = total + recentre%weights(j)
          centre = centre + recentre%weights(j) * tau0(:, ia, is)
        END DO
      END DO
      IF (.NOT. total > 0.0_real64) THEN
        CALL cpmdc_note_embed_failure('CPMDC_RECENTER=mass found no atomic masses'//c_null_char)
        ierr = 1
        RETURN
      END IF
      recentre%weights = recentre%weights / total
      centre = centre / total
      ! ht(k,:) is lattice vector k and htm1 its inverse, so x = s ht and
      ! s = x htm1 in rows.
      DO k = 1, 3
        frac(k) = DOT_PRODUCT(centre, metr_com%htm1(:, k))
      END DO
      shift = 0.0_real64
      DO k = 1, 3
        IF (free_axis(k)) shift = shift + (0.5_real64 - frac(k)) * metr_com%ht(k, :)
      END DO
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          tau0(:, ia, is) = tau0(:, ia, is) + shift
        END DO
      END DO
      DO m = 1, 3
        DO j = 1, 3
          DO k = 1, 3
            IF (free_axis(k)) recentre%proj(m, j) = recentre%proj(m, j) + &
                 metr_com%htm1(m, k) * metr_com%ht(k, j)
          END DO
        END DO
      END DO
      recentre%shift_ang = shift / fbohr
      recentre%applied = .TRUE.
    END IF
    IF (margin < 0.0_real64) RETURN
    closest = HUGE(1.0_real64)
    DO k = 1, 3
      IF (.NOT. free_axis(k)) CYCLE
      height = 1.0_real64 / NORM2(metr_com%htm1(:, k))
      DO is = 1, ions1%nsp
        DO ia = 1, ions0%na(is)
          frac(k) = DOT_PRODUCT(tau0(:, ia, is), metr_com%htm1(:, k))
          gap = MIN(frac(k), 1.0_real64 - frac(k)) * height / fbohr
          closest = MIN(closest, gap)
        END DO
      END DO
    END DO
    IF (closest < margin) THEN
      WRITE(num(1), '(F12.3)') closest
      WRITE(num(2), '(F12.3)') margin
      WRITE(msg, '(5A)') 'an atom is ', TRIM(ADJUSTL(num(1))), &
           ' Angstrom from a face of the isolated cell (negative outside it), closer than '// &
           'CPMDC_RECENTER_MARGIN ', TRIM(ADJUSTL(num(2))), ': enlarge the cell or set CPMDC_RECENTER'
      CALL cpmdc_note_embed_failure(TRIM(msg)//c_null_char)
      ierr = 1
    END IF
  END SUBROUTINE

  ! The gradient of E(x + s(x), Rmm + s(x)) at shifted positions. External
  ! charges share the rigid shift, so their gradient contributes to ds/dx.
  SUBROUTINE embed_project_gradient(n_atoms, recentre, grad, n_external, &
      external_force)
    INTEGER, INTENT(IN) :: n_atoms, n_external
    TYPE(embed_recentring), INTENT(IN) :: recentre
    REAL(c_double), INTENT(INOUT) :: grad(*)
    REAL(c_double), INTENT(IN) :: external_force(*)
    REAL(real64) :: total(3), along(3)
    INTEGER :: i
    IF (.NOT. recentre%applied) RETURN
    total = 0.0_real64
    DO i = 1, n_atoms
      total = total + REAL(grad(3*i-2:3*i), KIND=real64)
    END DO
    DO i = 1, n_external
      total = total - REAL(external_force(3*i-2:3*i), KIND=real64)
    END DO
    along = MATMUL(recentre%proj, total)
    DO i = 1, n_atoms
      grad(3*i-2:3*i) = REAL(REAL(grad(3*i-2:3*i), KIND=real64) - &
           recentre%weights(i) * along, KIND=c_double)
    END DO
  END SUBROUTINE

  ! 2 mass, 1 geometry, 0 off, -1 unknown.
  INTEGER FUNCTION recentre_policy() RESULT(policy)
    CHARACTER(LEN=32) :: v
    INTEGER :: st, n
    policy = 2
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_RECENTER', v, LENGTH=n, STATUS=st)
    IF (st /= 0 .OR. n == 0) RETURN
    SELECT CASE (TRIM(ADJUSTL(v)))
    CASE ('mass', 'MASS')
      policy = 2
    CASE ('geometry', 'GEOMETRY')
      policy = 1
    CASE ('off', 'OFF', '0')
      policy = 0
    CASE DEFAULT
      policy = -1
    END SELECT
  END FUNCTION

  ! NaN when the value does not parse.
  REAL(real64) FUNCTION recentre_margin() RESULT(margin)
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan
    CHARACTER(LEN=64) :: v
    INTEGER :: st, n, ios
    margin = 0.0_real64
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_RECENTER_MARGIN', v, LENGTH=n, STATUS=st)
    IF (st /= 0 .OR. n == 0) RETURN
    READ(v, *, IOSTAT=ios) margin
    IF (ios /= 0) margin = ieee_value(margin, ieee_quiet_nan)
  END FUNCTION

  SUBROUTINE embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, &
      n_external, external_pos, external_q, external_width, external_force, ok)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_scatter_species_gradient, &
        cpmdc_note_embed_failure, cpmdc_note_scf_unconverged, &
        cpmdc_note_deck_external_potential, &
        cpmdc_note_periodic_embedding
    USE norm, ONLY: scf_gemax => gemax
    USE wfopts_utils, ONLY: wfopts, embed_set_dipole_snapshot, &
        embed_set_charges_snapshot
    USE rwfopt_utils, ONLY: embed_set_warm_orbitals, embed_set_need_forces, &
        embed_set_write_files, embed_choose_warm_start_metric_v2
    USE phfac_utils, ONLY: phfac
    USE metr, ONLY: metr_com
    USE cnst, ONLY: fbohr
    USE ener, ONLY: ener_com, chrg, ener_c, ener_d
    USE coor, ONLY: tau0, fion, taup
    USE ions, ONLY: ions0, ions1
    USE store_types, ONLY: cprint, iprint_force, restart1
    USE system, ONLY: cnti, cntl, cntr, parm
    USE benc, ONLY: ibench
    USE strs, ONLY: paiu
    USE isos, ONLY: isos1
    USE ropt, ONLY: ropt_mod, iteropt
    USE parac, ONLY: paral
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms
    REAL(c_double), INTENT(IN) :: pos(*), external_pos(*), external_q(*)
    REAL(c_double), INTENT(IN) :: external_width(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*), external_force(*)
    INTEGER, INTENT(IN) :: n_external
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, is, ia, k, idx, nmax, i, j, astat, slot
    INTEGER(c_int), ALLOCATABLE :: origin(:)
    REAL(c_double), ALLOCATABLE :: species_grad(:), warm_geo(:)
    TYPE(embed_recentring) :: recentre
    REAL(real64) :: omega
    REAL(real64) :: external_ion_energy
    LOGICAL :: stress_computed, was_diis, was_pcg, was_pcgmin, was_prec
    INTEGER :: inwfun_deck, ibench_deck, warm_source
    INTEGER(c_int) :: scf_minimizer, scf_steps, first_minimizer, first_steps
    REAL(c_double) :: scf_gemax_final
    REAL(real64) :: warm_dist, warm_cell(9)
    INTEGER(SELECTED_INT_KIND(18)) :: warm_donor
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
    ALLOCATE(origin(n_atoms), species_grad(nmax), warm_geo(nmax), STAT=astat)
    IF (astat /= 0) THEN
      CALL cpmdc_note_embed_failure('out of memory'//c_null_char)
      RETURN
    END IF
    CALL embed_set_tau0_from_pos(n_atoms, pos, z, origin, recentre, ierr)
    IF (ierr /= 0) RETURN
    CALL prepare_external_charges(n_external, external_pos, external_q, &
         external_width, recentre%shift_ang, ierr)
    IF (ierr /= 0) THEN
      IF (ierr == 2) THEN
        CALL cpmdc_note_deck_external_potential()
      ELSE IF (ierr == 3) THEN
        CALL cpmdc_note_periodic_embedding()
      ELSE
        CALL cpmdc_note_embed_failure('cannot allocate external-charge grid'//c_null_char)
      END IF
      RETURN
    END IF
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
    ! When that copy is stored, rinitwf returns after the phase factors
    ! and does not build a guess. inwfun 3 is the start when the copy is
    ! absent: loadc, one orthogonalisation, and no force evaluation.
    ! That stands in for the Lanczos guess and for the random start.
    ! Stored orbitals start this SCF only when they converged at a
    ! geometry within warm_rmsd_limit() of this one; farther, the guess
    ! below is the cheaper start.
    ! Orbitals and their geometry use the engine's species order, independent
    ! of the interleaving of species in a caller's coordinate array.
    DO slot = 1, n_atoms
      j = INT(origin(slot)) + 1
      warm_geo(3*slot-2:3*slot) = pos(3*j-2:3*j) + recentre%shift_ang
    END DO
    ! metr_com%ht contains lattice rows in Bohr, including deck-only cells.
    ! The orbital history belongs to the setup's G vectors, so a cell
    ! moved in place keeps it.
    DO i = 1, 3
      IF (embed_setup_ht_set) THEN
        warm_cell(3*i-2:3*i) = embed_setup_ht(i,1:3) / fbohr
      ELSE
        warm_cell(3*i-2:3*i) = metr_com%ht(i,1:3) / fbohr
      END IF
    END DO
    CALL embed_choose_warm_start_metric_v2(warm_geo, warm_cell, .NOT. isos1%tclust, &
        warm_rmsd_limit(), warm_dist, warm_source, warm_donor, warm_geometry_distance)
    CALL note_warm_start(warm_source, warm_dist, warm_donor, .NOT. isos1%tclust)
    last_warm_source = warm_source
    inwfun_deck = cnti%inwfun
    ibench_deck = ibench(1)
    ! The cheap guess is for a process that already holds c0. A new
    ! session's counter is 0; the latch says the basis still matches.
    ! Stored orbitals dropped for distance give way to the deck's guess:
    ! on the Si3N4 band the simple guess took 65 to 205 SCF steps where
    ! the deck's atomic guess took 37 to 52.
    IF (image%cfg_warm_steps > 0 .OR. embed_basis_latched) THEN
      IF (warm_source /= -1) cnti%inwfun = 3
      ibench(1) = 1
    END IF
    CALL embed_set_charges_snapshot(embed_charges_wanted())
    CALL embed_set_dipole_snapshot(.TRUE.)
    CALL embed_set_warm_orbitals(.TRUE.)
    CALL wfopts
    CALL note_scf_steps(iteropt%nfi)
    first_minimizer = 0_c_int
    first_steps = 0_c_int
    scf_minimizer = scf_minimizer_code()
    scf_steps = INT(iteropt%nfi, KIND=c_int)
    scf_gemax_final = REAL(scf_gemax, KIND=c_double)
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
      first_minimizer = scf_minimizer
      first_steps = scf_steps
      scf_minimizer = scf_minimizer_code()
      scf_steps = INT(iteropt%nfi, KIND=c_int)
      scf_gemax_final = REAL(scf_gemax, KIND=c_double)
      cntl%diis = was_diis
      cntl%pcg = was_pcg
      cntl%pcgmin = was_pcgmin
      cntl%prec = was_prec
    END IF
    cnti%inwfun = inwfun_deck
    ibench(1) = ibench_deck
    external_ion_energy = 0.0_real64
    CALL add_external_ion_terms(external_ion_energy)
    ! eext holds the deck's own EXTERNAL POTENTIAL energy unless session
    ! charges are active, so only that case owns the sum. A deck with
    ! EXTERNAL POTENTIAL ADD and no session charges is left as OpenCPMD
    ! reports it.
    IF (embed_external_count > 0) THEN
      ener_com%etot = ener_com%etot + ener_com%eext + external_ion_energy
      ener_com%eext = ener_com%eext + external_ion_energy
    END IF
    DO i = 1, 3 * n_external
      external_force(i) = REAL(embed_external_force(1 + MOD(i - 1, 3), &
           1 + (i - 1) / 3), KIND=c_double)
    END DO
    energy_h = REAL(ener_com%etot, KIND=c_double)
    image%energy%etot = energy_h
    image%energy%ekin = REAL(ener_com%ekin, KIND=c_double)
    image%energy%epseu = REAL(ener_com%epseu, KIND=c_double)
    image%energy%enl = REAL(ener_com%enl, KIND=c_double)
    image%energy%eht = REAL(ener_com%eht, KIND=c_double)
    image%energy%exc = REAL(ener_com%exc, KIND=c_double)
    image%energy%eext = REAL(ener_com%eext, KIND=c_double)
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
        CALL embed_project_gradient(n_atoms, recentre, grad, n_external, &
             external_force)
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
    CALL snapshot_prop_from_modules(image, n_atoms, grad, origin)
    ! rwfopt computes the ionic forces only for converged orbitals: an SCF
    ! that ran out of MAXITER leaves fion zero, which a caller would read as
    ! a stationary point. Report it as a failure instead, with the last
    ! pass's minimizer, step count, MAXITER and orbital gradient. A stopgm
    ! has its own report.
    IF (ropt_mod%convwf .AND. ABS(energy_h) > 1.0e-8_c_double) ok = 1_c_int
    IF (.NOT. ropt_mod%convwf) THEN
      IF (cpmdc_stop_code() == 0_c_int) &
          CALL cpmdc_note_scf_unconverged(scf_minimizer, scf_steps, &
          INT(cnti%nomore_iter, KIND=c_int), scf_gemax_final, &
          REAL(cntr%tolog, KIND=c_double), first_minimizer, first_steps)
    END IF
  END SUBROUTINE

  ! The minimizer the cntl flags select, in the codes of
  ! src/cpmdc_scf_report.h.
  INTEGER(c_int) FUNCTION scf_minimizer_code() RESULT(code)
    USE system, ONLY: cntl
    IF (cntl%diis) THEN
      code = 1_c_int
    ELSE IF (cntl%pcg .AND. cntl%pcgmin) THEN
      code = 3_c_int
    ELSE IF (cntl%pcg) THEN
      code = 2_c_int
    ELSE IF (cntl%tsde) THEN
      code = 4_c_int
    ELSE
      code = 0_c_int
    END IF
  END FUNCTION


  ! CPMDC_SCF_STEPS=1 reports the iteration count of the SCF that just
  ! returned. The error unit stays attached when unit 6 is a file, and
  ! libgfortran buffers it when it is a file too, so each report flushes.
  ! An SCF that started from EXTRAPOLATE WFN's prediction over two or
  ! more earlier calls is reported first, with that count.
  SUBROUTINE note_scf_steps(nfi)
    USE ISO_FORTRAN_ENV, ONLY: error_unit
    USE parac, ONLY: paral
    USE rwfopt_utils, ONLY: embed_extrapolation_used
    INTEGER, INTENT(IN) :: nfi
    CHARACTER(LEN=8) :: v
    INTEGER :: st
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_SCF_STEPS', v, STATUS=st)
    IF (st /= 0) RETURN
    IF (TRIM(v) /= '1') RETURN
    IF (.NOT. paral%io_parent) RETURN
    IF (embed_extrapolation_used() > 1) WRITE(error_unit, '(A,I0)') &
        'cpmdc_wfn_extrapolated ', embed_extrapolation_used()
    WRITE(error_unit, '(A,I0)') 'cpmdc_scf_steps ', nfi
    FLUSH(error_unit)
  END SUBROUTINE

  ! The largest root-mean-square distance per atom, in Angstrom, between
  ! the geometry stored orbitals converged at and the geometry of a call
  ! that may start from them. CPMDC_WARM_RMSD overrides it; a negative
  ! value starts from the stored orbitals at any distance.
  REAL(real64) FUNCTION warm_geometry_distance(a, b, cell, periodic) RESULT(d)
    REAL(real64), INTENT(IN) :: a(:), b(:), cell(9)
    LOGICAL, INTENT(IN) :: periodic
    INTERFACE
      FUNCTION cpmdc_geometry_rmsd(n, a, b, cell, periodic) BIND(C) RESULT(rmsd)
        IMPORT :: c_double, c_int
        INTEGER(c_int), VALUE :: n, periodic
        REAL(c_double), INTENT(IN) :: a(*), b(*), cell(*)
        REAL(c_double) :: rmsd
      END FUNCTION
    END INTERFACE
    d = HUGE(1.0_real64)
    IF (SIZE(a) /= SIZE(b) .OR. MOD(SIZE(a), 3) /= 0) RETURN
    d = cpmdc_geometry_rmsd(INT(SIZE(a)/3, c_int), a, b, cell, MERGE(1_c_int, 0_c_int, periodic))
  END FUNCTION warm_geometry_distance

  REAL(real64) FUNCTION warm_rmsd_limit()
    ! Default limit, Angstrom. On the Si3N4 BLYP 70 Ry band (7 atoms,
    ! 4 ranks), against the deck's guess in the same process, stored
    ! orbitals saved 4.5 s per call over five images at 0.12 Angstrom,
    ! 2.1 s at 0.15 and 1.2 s at 0.17, and cost 0.4 s at 0.19 and 2.4 s
    ! at 0.23. A band step capped at 0.1 Angstrom per atom stays inside.
    REAL(real64), PARAMETER :: warm_rmsd_default = 0.15_real64
    CHARACTER(LEN=32) :: v
    INTEGER :: st, ios
    REAL(real64) :: x
    warm_rmsd_limit = warm_rmsd_default
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_WARM_RMSD', v, STATUS=st)
    IF (st /= 0 .OR. LEN_TRIM(v) == 0) RETURN
    READ(v, *, IOSTAT=ios) x
    IF (ios == 0) warm_rmsd_limit = x
  END FUNCTION

  ! CPMDC_SCF_STEPS=1 also reports where each SCF started: source 0 the
  ! selected key's store, 1 another key's parked orbitals, -1 the guess
  ! because every stored geometry was beyond warm_rmsd_limit(), -2 the
  ! guess because nothing was stored; and the nearest distance, Angstrom.
  SUBROUTINE note_warm_start(source, dist, donor, periodic)
    USE parac, ONLY: paral
    INTEGER, INTENT(IN) :: source
    INTEGER(SELECTED_INT_KIND(18)), INTENT(IN) :: donor
    LOGICAL, INTENT(IN) :: periodic
    REAL(real64), INTENT(IN) :: dist
    CHARACTER(LEN=8) :: v
    INTEGER :: st
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_SCF_STEPS', v, STATUS=st)
    IF (st /= 0) RETURN
    IF (TRIM(v) /= '1') RETURN
    IF (.NOT. paral%io_parent) RETURN
    IF (source == -2) THEN
      WRITE(error_unit, '(A,I0,A,L1)') 'cpmdc_warm_start ', source, ' periodic=', periodic
    ELSE
      WRITE(error_unit, '(A,I0,1X,ES24.16E3,A,I0,A,L1)') &
          'cpmdc_warm_start ', source, dist, ' donor_key=', donor, ' periodic=', periodic
    END IF
    FLUSH(error_unit)
  END SUBROUTINE

  ! CPMDC_CHARGE_SNAPSHOT=0 skips the Hirshfeld snapshot. Any other
  ! value, or an unset variable, requests it. cpmd.x never calls the setter.
  LOGICAL FUNCTION embed_charges_wanted()
    CHARACTER(LEN=32) :: v
    INTEGER :: st
    embed_charges_wanted = .TRUE.
    CALL GET_ENVIRONMENT_VARIABLE('CPMDC_CHARGE_SNAPSHOT', v, STATUS=st)
    IF (st == 0 .AND. TRIM(v) == '0') embed_charges_wanted = .FALSE.
  END FUNCTION

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

  SUBROUTINE drop_stale_basis(image)
    USE rwfopt_utils, ONLY: embed_reset_warm_orbitals
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    CALL embed_reset_warm_orbitals()
    CALL clear_embed_basis_latch()
    image%cfg_warm_steps = 0
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
    embed_basis_changed = .FALSE.
    IF (.NOT. embed_basis_latched) RETURN
    ! A detached image has no cell of its own. The cell stored with the
    ! orbitals is the one a new session has to match.
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
    IF (embed_method_changed(n_atoms, z, knobs)) embed_basis_changed = .TRUE.
  END FUNCTION

  ! Every part of the basis but the cell. c0 is plane-wave coefficients,
  ! so the count of each element is the composition.
  LOGICAL FUNCTION embed_method_changed(n_atoms, z, knobs)
    INTEGER, INTENT(IN) :: n_atoms
    INTEGER(c_int), INTENT(IN) :: z(*)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    REAL(c_double) :: cutoff
    embed_method_changed = .FALSE.
    cutoff = REAL(knobs%cutoff_ry, KIND=c_double)
    IF (ABS(cutoff - embed_saved_cutoff) > 1.0e-8_c_double) &
        embed_method_changed = .TRUE.
    IF (knobs%charge /= embed_saved_charge) embed_method_changed = .TRUE.
    IF (knobs%mult /= embed_saved_mult) embed_method_changed = .TRUE.
    IF (TRIM(knobs%functional) /= TRIM(embed_saved_functional)) &
        embed_method_changed = .TRUE.
    IF (.NOT. ALLOCATED(embed_saved_deck)) THEN
      IF (ALLOCATED(knobs%input_deck)) THEN
        IF (LEN_TRIM(knobs%input_deck) > 0) embed_method_changed = .TRUE.
      END IF
    ELSE IF (.NOT. ALLOCATED(knobs%input_deck)) THEN
      IF (LEN_TRIM(embed_saved_deck) > 0) embed_method_changed = .TRUE.
    ELSE IF (knobs%input_deck /= embed_saved_deck) THEN
      embed_method_changed = .TRUE.
    END IF
    ! c0 is plane-wave coefficients, so the count of each element is the
    ! composition.
    IF (n_atoms /= embed_saved_natoms) embed_method_changed = .TRUE.
    IF (.NOT. same_element_counts(n_atoms, z)) embed_method_changed = .TRUE.
  END FUNCTION

  ! CPMDC_CELL_IN_PLACE is the largest principal strain of a cell from the
  ! setup cell, as a fraction of a length, that keeps the setup's plane waves
  ! on a constant-cutoff deck. The C entry gate has already read it on root.
  REAL(real64) FUNCTION embed_cell_strain_limit()
    embed_cell_strain_limit = REAL(embed_cell_in_place_limit, KIND=real64)
  END FUNCTION

  ! The largest eigenvalue of a symmetric 3x3 matrix (Smith 1961).
  REAL(real64) FUNCTION largest_symmetric_eigenvalue(s) RESULT(lmax)
    REAL(real64), INTENT(IN) :: s(3, 3)
    REAL(real64) :: p1, p2, p, q, r, b(3, 3), phi
    INTEGER :: i
    p1 = s(1,2)**2 + s(1,3)**2 + s(2,3)**2
    IF (p1 <= 0.0_real64) THEN
      lmax = MAX(s(1,1), s(2,2), s(3,3))
      RETURN
    END IF
    q = (s(1,1) + s(2,2) + s(3,3)) / 3.0_real64
    p2 = (s(1,1) - q)**2 + (s(2,2) - q)**2 + (s(3,3) - q)**2 + 2.0_real64 * p1
    p = SQRT(p2 / 6.0_real64)
    b = s
    DO i = 1, 3
      b(i,i) = b(i,i) - q
    END DO
    b = b / p
    r = 0.5_real64 * (b(1,1) * (b(2,2) * b(3,3) - b(2,3) * b(3,2)) - &
         b(1,2) * (b(2,1) * b(3,3) - b(2,3) * b(3,1)) + &
         b(1,3) * (b(2,1) * b(3,2) - b(2,2) * b(3,1)))
    r = MAX(-1.0_real64, MIN(1.0_real64, r))
    phi = ACOS(r) / 3.0_real64
    lmax = q + 2.0_real64 * p * COS(phi)
  END FUNCTION

  ! A constant-cutoff deck can keep its plane waves through a change of the
  ! periodic cell. The penalty of CONSTANT CUTOFF leaves the components
  ! above E0 + 2 sigma empty, so the G vectors of the setup cell still
  ! hold every occupied component of the new cell while they contain
  ! that sphere of it, and the density's G vectors the sphere of twice
  ! the radius. The plane waves between that sphere and the cutoff are
  ! then the setup's rather than the new cell's, which moves the stress
  ! away from that of a fresh setup in proportion to the strain, so the
  ! strain from the setup cell stays within CPMDC_CELL_IN_PLACE. Only the
  ! cell may differ from the latched basis: Gamma point, norm-conserving
  ! pseudopotentials, the setup's lattice class, and a minimum-image rule
  ! that embed_move_cell can follow.
  LOGICAL FUNCTION embed_cell_moves_in_place(n_atoms, z, cell, has_cell, knobs)
    USE prcp, ONLY: prcp_com
    USE gvec, ONLY: gvec_com
    USE kpts, ONLY: tkpts
    USE isos, ONLY: isos1
    USE system, ONLY: parm
    USE pslo, ONLY: pslo_com
    USE cnst, ONLY: fbohr
    INTEGER, INTENT(IN) :: n_atoms, has_cell
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(IN) :: cell(*)
    TYPE(embed_knobs), INTENT(IN) :: knobs
    REAL(real64) :: h_new(3, 3), t(3, 3), s(3, 3), stretch, shrink, sphere
    REAL(real64) :: limit
    REAL(c_double) :: setup_cell(9)
    INTEGER :: i, setup_code
    embed_cell_moves_in_place = .FALSE.
    limit = embed_cell_strain_limit()
    IF (limit <= 0.0_real64) RETURN
    IF (.NOT. cpmd_set_up .OR. .NOT. embed_setup_ht_set) RETURN
    IF (has_cell == 0 .OR. embed_saved_has_cell == 0) RETURN
    IF (isos1%tclust .OR. tkpts%tkpnt .OR. pslo_com%tivan) RETURN
    IF (prcp_com%akin <= 0.0_real64 .OR. prcp_com%skin <= 0.0_real64) RETURN
    IF (embed_method_changed(n_atoms, z, knobs)) RETURN
    DO i = 1, 3
      setup_cell(3*i-2:3*i) = embed_setup_ht(i,1:3) / fbohr
      h_new(i,1:3) = REAL(cell(3*i-2:3*i), KIND=real64) * fbohr
    END DO
    setup_code = periodic_symmetry_code(setup_cell, 1)
    IF (setup_code == 0) RETURN
    IF (periodic_symmetry_code(cell, has_cell) /= setup_code) RETURN
    SELECT CASE (parm%ibrav)
    CASE (1)
      IF (ABS(h_new(2,2) - h_new(1,1)) > 1.0e-10_real64 * h_new(1,1) .OR. &
          ABS(h_new(3,3) - h_new(1,1)) > 1.0e-10_real64 * h_new(1,1)) RETURN
    CASE (8)
    CASE (14)
      IF (ABS(parm%apbc(2)) > 0.0_real64) RETURN
    CASE DEFAULT
      RETURN
    END SELECT
    ! |G(m)|^2 = 4 pi^2 m' (H H')^-1 m for lattice rows H, so the ratio of
    ! the setup's to the new over m lies between the extreme eigenvalues
    ! of T T' with T = H_setup^-1 H_new, the squares of the principal
    ! stretches.
    t = MATMUL(inverse3(embed_setup_ht), h_new)
    s = MATMUL(t, TRANSPOSE(t))
    stretch = largest_symmetric_eigenvalue(s)
    shrink = largest_symmetric_eigenvalue(inverse3(s))
    IF (SQRT(stretch) - 1.0_real64 > limit) RETURN
    IF (1.0_real64 - 1.0_real64 / SQRT(shrink) > limit) RETURN
    sphere = prcp_com%gckin + 2.0_real64 * prcp_com%gskin
    IF (sphere <= 0.0_real64) RETURN
    IF (stretch * sphere > gvec_com%gcutw) RETURN
    IF (stretch * 4.0_real64 * sphere > gvec_com%gcut) RETURN
    embed_cell_moves_in_place = .TRUE.
  END FUNCTION

  FUNCTION inverse3(a) RESULT(b)
    REAL(real64), INTENT(IN) :: a(3, 3)
    REAL(real64) :: b(3, 3), det
    b(1,1) = a(2,2) * a(3,3) - a(2,3) * a(3,2)
    b(1,2) = a(1,3) * a(3,2) - a(1,2) * a(3,3)
    b(1,3) = a(1,2) * a(2,3) - a(1,3) * a(2,2)
    b(2,1) = a(2,3) * a(3,1) - a(2,1) * a(3,3)
    b(2,2) = a(1,1) * a(3,3) - a(1,3) * a(3,1)
    b(2,3) = a(1,3) * a(2,1) - a(1,1) * a(2,3)
    b(3,1) = a(2,1) * a(3,2) - a(2,2) * a(3,1)
    b(3,2) = a(1,2) * a(3,1) - a(1,1) * a(3,2)
    b(3,3) = a(1,1) * a(2,2) - a(1,2) * a(2,1)
    det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)
    b = b / det
  END FUNCTION

  ! Puts the caller's cell under the setup's G vectors. OpenCPMD's initrun
  ! does the same on each SCF with cntl%tpres, but after rinitwf, which
  ! builds a guess on the cell it finds. The minimum-image half lengths
  ! of rpiiint's pbc come from setsc, once per setup; a cell given by
  ! CELL VECTORS folds through htm1 instead (parm%apbc(2) = 0).
  SUBROUTINE embed_move_cell(cell)
    USE metr, ONLY: metr_com
    USE cnst, ONLY: fbohr
    USE system, ONLY: parm
    USE setsc_utils, ONLY: ihmat
    USE newcell_utils, ONLY: newcell
    REAL(c_double), INTENT(IN) :: cell(*)
    INTEGER :: i
    DO i = 1, 3
      metr_com%ht(i,1:3) = REAL(cell(3*i-2:3*i), KIND=real64) * fbohr
    END DO
    parm%a1 = metr_com%ht(1,:)
    parm%a2 = metr_com%ht(2,:)
    parm%a3 = metr_com%ht(3,:)
    CALL ihmat(metr_com%ht, metr_com%htm1, parm%omega)
    SELECT CASE (parm%ibrav)
    CASE (1)
      parm%apbc(1) = 0.5_real64 * metr_com%ht(1,1)
    CASE (8)
      DO i = 1, 3
        parm%apbc(i) = 0.5_real64 * metr_com%ht(i,i)
      END DO
    END SELECT
    CALL newcell
  END SUBROUTINE

  SUBROUTINE run_embed_scf(image, n_atoms, pos, z, cell, has_cell, energy_h, &
      grad, n_external, external_pos, external_q, external_width, &
      external_force, ok)
    USE cpmdc_embed_host_iface, ONLY: cpmdc_memfd_write, cpmdc_memfd_close, &
        cpmdc_pseudopotential_directory, cpmdc_prepare_pp_cwd, &
        cpmdc_restore_host_cwd, cpmdc_note_embed_failure, &
        cpmdc_embed_failure_pending, &
        cpmdc_note_deck_external_potential, &
        cpmdc_deck_defines_external_potential
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
    USE metr, ONLY: metr_com
    USE ropt, ONLY: init_pinf_pointers
    USE bicanonicalCpmd, ONLY: bicanonicalCpmdConfig, bicanonicalCpmdInputConfig, New
    USE bicanonicalConfig, ONLY: New
    TYPE(cpmdc_embed_image), INTENT(INOUT) :: image
    INTEGER, INTENT(IN) :: n_atoms, has_cell, n_external
    REAL(c_double), INTENT(IN) :: pos(*), cell(*), external_pos(*), &
        external_q(*), external_width(*)
    INTEGER(c_int), INTENT(IN) :: z(*)
    REAL(c_double), INTENT(OUT) :: energy_h
    REAL(c_double), INTENT(OUT) :: grad(*), external_force(*)
    INTEGER(c_int), INTENT(OUT) :: ok
    INTEGER :: ierr, idx, nmax, nlen, mfd
    TYPE(embed_knobs) :: knobs
    LOGICAL :: tinfo, warm, in_place
    CHARACTER(LEN=:), ALLOCATABLE :: deck
    CHARACTER(LEN=64) :: mempath
    CHARACTER(KIND=c_char) :: pp_c(1024)
    ok = 0_c_int
    energy_h = 0.0_c_double
    ! Close ignores a descriptor that was never opened.
    mfd = -1
    ! The count is checked before any collective. A rank that rejects it
    ! still agrees, or the others enter the warm evaluation or startpa.
    IF (n_atoms <= 0 .OR. n_atoms > HUGE(nmax) / 3) THEN
      CALL cpmdc_note_embed_failure('invalid arguments'//c_null_char)
    END IF
    IF (ranks_agree_failure()) RETURN
    nmax = n_atoms * 3
    DO idx = 1, nmax
      grad(idx) = 0.0_c_double
    END DO
    ! Saved orbitals match one method, cell, and elemental composition.
    ! A shape mismatch in embed_restore_orbitals calls stopgm, so drop
    ! the copy before rwfopt when any of those change. A new session and
    ! a reordering of the same atomic numbers are not a new basis, and
    ! neither is a cell that embed_cell_moves_in_place accepts.
    knobs = knobs_of(image)
    warm = .FALSE.
    in_place = .FALSE.
    IF (embed_basis_latched) THEN
      warm = .NOT. embed_basis_changed(image, n_atoms, z, cell, has_cell, knobs)
      IF (.NOT. warm) THEN
        in_place = embed_cell_moves_in_place(n_atoms, z, cell, has_cell, knobs)
        warm = in_place
      END IF
    END IF
    IF (warm) THEN
      ! embed_use_output_dir notes the directory failure on this rank.
      ! Agree before embed_eval_energy_grad, then return to the host
      ! directory. The setup stays; a warm failure does not tear it down.
      IF (.NOT. embed_use_output_dir(image)) ierr = 1
      IF (ranks_agree_failure()) RETURN
      ! The latched cell is the one metr_com%ht holds, even when the SCF
      ! that follows fails.
      IF (in_place) THEN
        CALL embed_move_cell(cell)
        CALL latch_warm_cell(image, cell, has_cell)
      END IF
      CALL embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, &
           n_external, external_pos, external_q, external_width, &
           external_force, ok)
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
    IF (embed_basis_latched) CALL drop_stale_basis(image)
    ! Cold: honor the rendered method deck.
    ! 1) &ATOMS that already has coordinates is kept.
    ! 2) Method sections keep their text. &ATOMS is rebuilt from the step.
    !    A '!SPECIES' entry supplies that element's file, LMAX, LOC, and
    !    KLEINMAN-BYLANDER. Other elements use the built-in table.
    ! 3) Else a minimal deck with the applied functional, cutoff, charge, and
    !    multiplicity.
    ! Geometry for forces always comes from the C arrays into TAU0 after parse.
    ! A basis change or a call after a CPMD stop: free the first setup
    ! before this one allocates the same module arrays again.
    IF (cpmd_set_up .OR. cpmd_setup_started) THEN
      IF (cpmdc_embed_teardown() == 0_c_int) THEN
        CALL cpmdc_note_embed_failure( &
             'cannot tear down the previous setup'//c_null_char)
      END IF
    END IF
    ! Compose, the memfd, and the pseudopotential directory are local.
    ! A rank that fails one of them skips the later local steps and
    ! meets the others at the agreement below. None of those steps is
    ! a collective, so the ranks still enter that agreement together.
    ! The deck is closed and the host directory restored only then.
    IF (cpmdc_embed_failure_pending() == 0_c_int) THEN
      CALL embed_compose_cold_deck(n_atoms, pos, z, cell, has_cell, deck, &
           nlen, ierr, knobs)
      IF (ierr /= 0 .OR. nlen < 1 .OR. .NOT. ALLOCATED(deck) .OR. &
          nlen > HUGE(0_c_int)) THEN
        IF (cpmdc_embed_failure_pending() == 0_c_int) THEN
          CALL cpmdc_note_embed_failure( &
               'cold deck composition failed'//c_null_char)
        END IF
      ELSE IF (n_external > 0 .AND. &
          cpmdc_deck_defines_external_potential(deck, &
          INT(nlen, KIND=c_size_t)) /= 0_c_int) THEN
        CALL cpmdc_note_deck_external_potential()
      ELSE
        mfd = INT(cpmdc_memfd_write(deck, INT(nlen, KIND=c_int), mempath, &
             INT(LEN(mempath), KIND=c_int)))
        IF (mfd < 0) THEN
          CALL cpmdc_note_embed_failure( &
               'cannot create the cold deck'//c_null_char)
        ELSE
          paral%io_parent = .TRUE.
          ! (compose path above honors method-only Cap'n decks.)
          ! C wrote a NUL-terminated /proc/self/fd/N path into mempath.
          DO idx = 1, LEN(mempath)
            IF (IACHAR(mempath(idx:idx)) == 0) THEN
              cnts%inputfile = mempath(1:idx-1)
              EXIT
            END IF
          END DO
          ! OpenCPMD get_pplib uses argv[2] as the PP library whenever
          ! argc>1, ignoring CPMD_PP_LIBRARY_PATH. Hosts like Catch2/eOn
          ! always pass filters so argc>1; relative *PP basenames then
          ! only resolve via recpnew's second-chance CWD lookup.
          ! cpmdc_prepare_pp_cwd chdirs to the library and also exports
          ! CPMD_PP_LIBRARY_PATH (trailing slash) for argc==1 hosts.
          ! The process leaves that directory again once the files are
          ! read. A failed getcwd and a failed chdir both fail that call.
          IF (cpmdc_pseudopotential_directory(pp_c, &
              INT(1024, KIND=c_size_t)) /= 0_c_int) THEN
            CALL cpmdc_note_embed_failure(pp_c)
          ELSE IF (cpmdc_prepare_pp_cwd(pp_c) /= 0_c_int) THEN
            CALL cpmdc_note_embed_failure( &
                 'cannot enter the pseudopotential directory'//c_null_char)
          END IF
        END IF
      END IF
    END IF
    IF (cold_setup_stopped()) RETURN
    ! A stopgm returns to the caller the hook is armed for. The next
    ! init in this list would run on undefined module state and can
    ! fault, and that fault never reaches the C handler. Check after
    ! every call. cpmd_set_up stays false until the list finishes.
    cpmd_setup_started = .TRUE.
    IF (.NOT. timer_started) THEN
      CALL tistart(tcpu0, twall0)
      timer_started = .TRUE.
    END IF
    IF (cold_setup_stopped()) RETURN
    CALL init_fileopen
    IF (cold_setup_stopped()) RETURN
    CALL startpa
    IF (cold_setup_stopped()) RETURN
    CALL New(bicanonicalCpmdInputConfig)
    IF (cold_setup_stopped()) RETURN
    tinfo = .TRUE.
    CALL init_pinf_pointers()
    IF (cold_setup_stopped()) RETURN
    CALL envir
    IF (cold_setup_stopped()) RETURN
    CALL setcnst
    IF (cold_setup_stopped()) RETURN
    CALL control
    IF (cold_setup_stopped()) RETURN
    CALL dftin
    IF (cold_setup_stopped()) RETURN
    CALL sysin
    IF (cold_setup_stopped()) RETURN
    ! PEF stress needs cntl%tpres before dqgalloc. Isolated/Hockney (tclust)
    ! must not set tpres: rinitwf→newcell→gf_periodic needs scg, which
    ! initclust only allocates for periodic cells.
    IF (.NOT. isos1%tclust .AND. embed_stress_wanted()) cntl%tpres = .TRUE.
    CALL setsc
    IF (cold_setup_stopped()) RETURN
    CALL detsp
    IF (cold_setup_stopped()) RETURN
    CALL mm_init
    IF (cold_setup_stopped()) RETURN
    CALL ratom
    IF (cold_setup_stopped()) RETURN
    ! Pseudopotential files are in memory. Leave that directory before the
    ! SCF so a stray write cannot land on the library.
    ! A directory failure is already noted. Agree it before this rank
    ! frees the group and the others enter the next setup call.
    IF (.NOT. embed_use_output_dir(image)) THEN
      ierr = 1
    END IF
    IF (cold_setup_stopped()) RETURN
    CALL vdwin
    IF (cold_setup_stopped()) RETURN
    CALL propin(tinfo)
    IF (cold_setup_stopped()) RETURN
    CALL setsys
    IF (cold_setup_stopped()) RETURN
    CALL New(bicanonicalCpmdConfig, bicanonicalCpmdInputConfig)
    IF (cold_setup_stopped()) RETURN
    CALL genxc
    IF (cold_setup_stopped()) RETURN
    CALL numpw
    IF (cold_setup_stopped()) RETURN
    CALL rinit
    IF (cold_setup_stopped()) RETURN
    CALL rinforce
    IF (cold_setup_stopped()) RETURN
    CALL fft_init()
    IF (cold_setup_stopped()) RETURN
    CALL ortho_init()
    IF (cold_setup_stopped()) RETURN
    CALL initclust
    IF (cold_setup_stopped()) RETURN
    CALL dg_init
    IF (cold_setup_stopped()) RETURN
    CALL nosalloc
    IF (cold_setup_stopped()) RETURN
    CALL exterp
    IF (cold_setup_stopped()) RETURN
    CALL setbasis
    IF (cold_setup_stopped()) RETURN
    CALL dqgalloc
    IF (cold_setup_stopped()) RETURN
    CALL prnginit
    IF (cold_setup_stopped()) RETURN
    CALL gle_alloc
    IF (cold_setup_stopped()) RETURN
    CALL vdw_wf_alloc
    IF (cold_setup_stopped()) RETURN
    ! The parsers have read the deck. Drop the anonymous fd before the SCF.
    ! A failed close is already noted. Agree it before any rank treats
    ! the setup as finished and enters the SCF.
    IF (close_cold_deck() /= 0) THEN
      ierr = 1
    END IF
    IF (cold_setup_stopped()) RETURN
    cpmd_set_up = .TRUE.
    cpmd_setup_started = .FALSE.
    embed_setup_ht = metr_com%ht
    embed_setup_ht_set = .TRUE.
    ! First and later forces: positions only from C arrays (nwchemc geom pattern).
    CALL embed_eval_energy_grad(image, n_atoms, pos, z, energy_h, grad, &
         n_external, external_pos, external_q, external_width, &
         external_force, ok)
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
  CONTAINS
    INTEGER FUNCTION close_cold_deck()
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int
      INTEGER(c_int) :: rc
      close_cold_deck = 0
      IF (mfd < 0) RETURN
      rc = cpmdc_memfd_close(INT(mfd, KIND=c_int))
      ! 0 closed. -2 released after an error: drop the number so a
      ! later close does not hit a reused descriptor. -1 still open:
      ! keep the number. mfd is local and the next run_embed_scf
      ! resets it, so only the same-call retry in cold_setup_stopped
      ! can close it.
      IF (rc == 0_c_int .OR. rc == -2_c_int) mfd = -1
      IF (rc == 0_c_int) RETURN
      close_cold_deck = 1
      ! Keep a failure that was already noted. This one is the reason
      ! only when the setup had not failed yet.
      IF (cpmdc_embed_failure_pending() /= 0_c_int) RETURN
      CALL cpmdc_note_embed_failure( &
           'cannot close the cold deck'//c_null_char)
    END FUNCTION close_cold_deck

    ! A noted failure with no partial setup to free. Every rank calls
    ! this, the text is copied, and only then does this rank return to
    ! the host directory. The cold setup list is not torn down.
    LOGICAL FUNCTION ranks_agree_failure()
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int
      INTERFACE
        FUNCTION cpmdc_agree_cold_stop() BIND(C, NAME='cpmdc_agree_cold_stop') &
            RESULT(flag)
          IMPORT :: c_int
          INTEGER(c_int) :: flag
        END FUNCTION cpmdc_agree_cold_stop
      END INTERFACE
      INTEGER(c_int) :: agreed, rc
      agreed = cpmdc_agree_cold_stop()
      ranks_agree_failure = agreed /= 0_c_int
      IF (.NOT. ranks_agree_failure) RETURN
      rc = cpmdc_restore_host_cwd()
    END FUNCTION ranks_agree_failure

    LOGICAL FUNCTION cold_setup_stopped()
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int
      INTERFACE
        FUNCTION cpmdc_agree_cold_stop() BIND(C, NAME='cpmdc_agree_cold_stop') &
            RESULT(flag)
          IMPORT :: c_int
          INTEGER(c_int) :: flag
        END FUNCTION cpmdc_agree_cold_stop
      END INTERFACE
      INTEGER(c_int) :: agreed, rc
      ! stopgm may return on one rank only. Agree before this rank
      ! leaves the list or frees the group, or the others enter the
      ! next setup collective and the job waits.
      agreed = cpmdc_agree_cold_stop()
      cold_setup_stopped = agreed /= 0_c_int
      IF (.NOT. cold_setup_stopped) RETURN
      ! Partial module state is undefined. Free it and leave cpmd_set_up
      ! false so the C handler reports the stop. The deck fd goes too.
      ! Every rank takes this path once the stop has been agreed.
      rc = cpmdc_restore_host_cwd()
      ierr = cpmdc_embed_teardown()
      ierr = close_cold_deck()
    END FUNCTION cold_setup_stopped
  END SUBROUTINE

#endif
END MODULE cpmd_embed_c_api
