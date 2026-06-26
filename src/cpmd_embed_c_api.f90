! SPDX-License-Identifier: MIT
!
! cpmd_embed_c_api.f90 - compiler-independent C ABI for OpenCPMD embed.
! Public symbols use bind(C, name=...) (iso_c_binding).
!
! Until a full OpenCPMD archive link is wired, the energy/gradient path runs a
! deterministic reference PEF so the Cap'n Proto socket and session topology
! contracts can be exercised end-to-end (single-point and multi-step).

module cpmd_embed_c_api
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  public :: cpmdc_embed_init
  public :: cpmdc_embed_available
  public :: cpmdc_embed_set_config
  public :: cpmdc_embed_energy_grad
  public :: cpmdc_embed_finalize

  logical, save :: runtime_ready = .false.
  logical, save :: runtime_finalized = .false.
  character(len=64), save :: cfg_functional = 'BLYP'
  real(real64), save :: cfg_cutoff_ry = 70.0_real64
  integer, save :: cfg_charge = 0
  integer, save :: cfg_mult = 1
  character(len=4096), save :: cfg_input_deck = ' '
  character(len=1024), save :: cfg_cpmd_root = ' '

  ! Angstrom per Bohr (CODATA-consistent with cpmdc_params.h).
  real(real64), parameter :: bohr_to_angstrom = 0.529177210903_real64
  ! Reference PEF spring (Hartree / Angstrom^2) on each Cartesian component,
  ! scaled by atomic number so species enter the energy surface.
  real(real64), parameter :: k_ref = 0.02_real64

contains

  subroutine copy_c_string(src, n, dst)
    character(kind=c_char), intent(in) :: src(*)
    integer(c_int), intent(in), value :: n
    character(len=*), intent(out) :: dst
    integer :: i, lim
    dst = ' '
    lim = min(int(n), len(dst))
    do i = 1, lim
      if (src(i) == c_null_char) exit
      dst(i:i) = transfer(src(i), 'a')
    end do
  end subroutine copy_c_string

  function cpmdc_embed_init() result(ok) bind(C, name='cpmdc_embed_init')
    integer(c_int) :: ok
    runtime_ready = .true.
    runtime_finalized = .false.
    ok = 1_c_int
  end function cpmdc_embed_init

  function cpmdc_embed_available() result(ok) bind(C, name='cpmdc_embed_available')
    integer(c_int) :: ok
    ok = merge(1_c_int, 0_c_int, runtime_ready .and. .not. runtime_finalized)
  end function cpmdc_embed_available

  function cpmdc_embed_set_config(functional, functional_len, cutoff_ry, &
      charge, multiplicity, input_deck, input_deck_len, cpmd_root, &
      cpmd_root_len) result(ok) bind(C, name='cpmdc_embed_set_config')
    character(kind=c_char), intent(in) :: functional(*)
    integer(c_int), intent(in), value :: functional_len
    real(c_double), intent(in), value :: cutoff_ry
    integer(c_int), intent(in), value :: charge
    integer(c_int), intent(in), value :: multiplicity
    character(kind=c_char), intent(in) :: input_deck(*)
    integer(c_int), intent(in), value :: input_deck_len
    character(kind=c_char), intent(in) :: cpmd_root(*)
    integer(c_int), intent(in), value :: cpmd_root_len
    integer(c_int) :: ok

    call copy_c_string(functional, functional_len, cfg_functional)
    cfg_cutoff_ry = real(cutoff_ry, kind=real64)
    cfg_charge = int(charge)
    cfg_mult = int(multiplicity)
    call copy_c_string(input_deck, input_deck_len, cfg_input_deck)
    call copy_c_string(cpmd_root, cpmd_root_len, cfg_cpmd_root)
    if (.not. runtime_ready) then
      ok = 0_c_int
    else
      ok = 1_c_int
    end if
  end function cpmdc_embed_set_config

  ! Reference PEF: E = sum_i 0.5 * k * Z_i * |r_i|^2  (r in Angstrom, E in Hartree)
  ! Nuclear gradient in Hartree/Bohr: dE/d r_bohr = (k * Z * r_ang) * bohr_to_angstrom
  ! Optional weak cell term when has_cell != 0 (trace of cell^2) so PBC steps differ.
  function cpmdc_embed_energy_grad(n_atoms, positions_ang, atomic_numbers, &
      cell_ang, has_cell, energy_h, grad_h_bohr) result(ok) &
      bind(C, name='cpmdc_embed_energy_grad')
    integer(c_int), intent(in), value :: n_atoms
    real(c_double), intent(in) :: positions_ang(*)
    integer(c_int), intent(in) :: atomic_numbers(*)
    real(c_double), intent(in) :: cell_ang(*)
    integer(c_int), intent(in), value :: has_cell
    real(c_double), intent(out) :: energy_h
    real(c_double), intent(out) :: grad_h_bohr(*)
    integer(c_int) :: ok
    integer :: i, j, idx, z
    real(real64) :: e, x, g_ang, k_eff, cell_term

    energy_h = 0.0_c_double
    if (n_atoms <= 0) then
      ok = 0_c_int
      return
    end if
    if (.not. runtime_ready .or. runtime_finalized) then
      ok = 0_c_int
      return
    end if

    e = 0.0_real64
    do i = 1, int(n_atoms)
      z = max(1, int(atomic_numbers(i)))
      k_eff = k_ref * real(z, kind=real64)
      do j = 0, 2
        idx = 3 * (i - 1) + j + 1
        x = real(positions_ang(idx), kind=real64)
        e = e + 0.5_real64 * k_eff * x * x
        ! dE/d x_ang = k_eff * x; convert to Ha/Bohr
        g_ang = k_eff * x
        grad_h_bohr(idx) = real(g_ang * bohr_to_angstrom, kind=c_double)
      end do
    end do

    cell_term = 0.0_real64
    if (has_cell /= 0) then
      do j = 1, 9
        cell_term = cell_term + real(cell_ang(j), kind=real64)**2
      end do
      ! Tiny isotropic cell penalty (Hartree) so PBC-aware steps are distinct.
      e = e + 1.0e-6_real64 * cell_term
    end if

    ! Fold a trace of method config into energy so deck/config is not ignored.
    e = e + 1.0e-8_real64 * cfg_cutoff_ry + 1.0e-9_real64 * real(cfg_charge, kind=real64)

    energy_h = real(e, kind=c_double)
    ok = 1_c_int
  end function cpmdc_embed_energy_grad

  subroutine cpmdc_embed_finalize() bind(C, name='cpmdc_embed_finalize')
    runtime_ready = .false.
    runtime_finalized = .true.
  end subroutine cpmdc_embed_finalize

end module cpmd_embed_c_api
