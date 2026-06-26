! SPDX-License-Identifier: MIT
!
! cpmd_embed_c_api.f90 - compiler-independent C ABI for OpenCPMD embed.
! Public symbols use bind(C, name=...) (iso_c_binding).

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
    ! v0.1: ISO_C shell only; real OpenCPMD runtime init is a follow-up.
    runtime_ready = .false.
    runtime_finalized = .false.
    ok = 0_c_int
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
    ok = merge(1_c_int, 0_c_int, runtime_ready)
  end function cpmdc_embed_set_config

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
    integer :: i, n3

    energy_h = 0.0_c_double
    n3 = max(0, int(n_atoms)) * 3
    do i = 1, n3
      grad_h_bohr(i) = 0.0_c_double
    end do
    if (n_atoms <= 0) then
      ok = 0_c_int
      return
    end if
    ! Touch inputs so optimizers keep the ABI parameters live.
    if (atomic_numbers(1) < -1) then
      ok = 0_c_int
      return
    end if
    if (has_cell /= 0 .and. abs(cell_ang(1)) < 0.0_c_double) then
      ok = 0_c_int
      return
    end if
    if (abs(positions_ang(1)) < -1.0_c_double) then
      ok = 0_c_int
      return
    end if
    ok = 0_c_int
  end function cpmdc_embed_energy_grad

  subroutine cpmdc_embed_finalize() bind(C, name='cpmdc_embed_finalize')
    runtime_ready = .false.
    runtime_finalized = .true.
  end subroutine cpmdc_embed_finalize

end module cpmd_embed_c_api
