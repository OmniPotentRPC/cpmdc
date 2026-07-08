! SPDX-License-Identifier: MIT
!
! Decode CPMDParams with public capnp-fortran and apply config knobs into
! embed state (functional / cutOffRy / charge / multiplicity / deck / root).

module cpmdc_embed_apply_params_mod
  use, intrinsic :: iso_c_binding, only: &
      c_char, c_double, c_int, c_ptr, c_size_t, c_null_char, &
      c_associated, c_f_pointer
  use, intrinsic :: iso_fortran_env, only: int8, int64, real64
  use capnp
  use potentials_capnp
  implicit none
  private

  public :: cpmdc_embed_apply_params
  public :: cpmdc_embed_get_config
  public :: applied_functional, applied_cutoff_ry
  public :: applied_charge, applied_mult
  public :: applied_input_deck, applied_cpmd_root

  character(len=64), save :: applied_functional = 'BLYP'
  real(real64), save :: applied_cutoff_ry = 70.0_real64
  integer, save :: applied_charge = 0
  integer, save :: applied_mult = 1
  character(len=4096), save :: applied_input_deck = ' '
  character(len=1024), save :: applied_cpmd_root = ' '

contains

  subroutine reset_applied()
    applied_functional = 'BLYP'
    applied_cutoff_ry = 70.0_real64
    applied_charge = 0
    applied_mult = 1
    applied_input_deck = ' '
    applied_cpmd_root = ' '
  end subroutine reset_applied

  subroutine copy_alloc_text(src, dst, maxlen)
    character(len=:), allocatable, intent(in) :: src
    character(len=*), intent(out) :: dst
    integer, intent(in) :: maxlen
    integer :: n
    dst = ' '
    if (.not. allocated(src)) return
    n = min(len_trim(src), maxlen)
    if (n > 0) dst(1:n) = src(1:n)
  end subroutine copy_alloc_text

  subroutine cchars_to_f(cbuf, clen, fstr)
    character(kind=c_char), intent(in) :: cbuf(*)
    integer(c_int), intent(in) :: clen
    character(len=*), intent(out) :: fstr
    integer :: i, n
    fstr = ' '
    if (clen <= 0) return
    n = min(int(clen), len(fstr))
    do i = 1, n
      if (cbuf(i) == c_null_char) exit
      fstr(i:i) = cbuf(i)
    end do
  end subroutine cchars_to_f

  subroutine f_to_c_chars(fstr, cbuf, clen)
    character(len=*), intent(in) :: fstr
    character(kind=c_char), intent(out) :: cbuf(*)
    integer(c_int), intent(in) :: clen
    integer :: i, n
    n = min(len_trim(fstr), max(0, int(clen) - 1))
    do i = 1, n
      cbuf(i) = fstr(i:i)
    end do
    if (clen > 0) cbuf(n + 1) = c_null_char
  end subroutine f_to_c_chars

  !> Decode CPMDParams wire bytes into embed knobs; accept C-rendered deck.
  !> Optional common-overlay scalars: empty functional / cutoff<=0 / has_* =0
  !> leave the wire effective config unchanged for that field.
  function cpmdc_embed_apply_params(params_capnp, params_capnp_size, &
      input_deck, input_deck_len, functional_ov, functional_ov_len, &
      cutoff_ov, has_charge_ov, charge_ov, has_mult_ov, mult_ov) result(rc) &
      bind(C, name='cpmdc_embed_apply_params')
    type(c_ptr), intent(in), value :: params_capnp
    integer(c_size_t), intent(in), value :: params_capnp_size
    character(kind=c_char), intent(in) :: input_deck(*)
    integer(c_int), intent(in), value :: input_deck_len
    character(kind=c_char), intent(in) :: functional_ov(*)
    integer(c_int), intent(in), value :: functional_ov_len
    real(c_double), intent(in), value :: cutoff_ov
    integer(c_int), intent(in), value :: has_charge_ov, charge_ov
    integer(c_int), intent(in), value :: has_mult_ov, mult_ov
    integer(c_int) :: rc
    integer(int8), pointer :: raw(:)
    integer(int8), allocatable :: bytes(:)
    type(capnp_message_t), target :: msg
    type(c_p_m_d_params_t) :: params
    type(c_p_m_d_input_section_t) :: sec
    type(c_p_m_d_system_section_t) :: sys
    type(c_p_m_d_dft_section_t) :: dft
    type(capnp_ptr_t) :: sections
    character(len=:), allocatable :: s
    character(len=64) :: fov
    integer :: err, i, n, tag, ib
    integer(int64) :: n64
    real(real64) :: cut

    rc = -1_c_int
    call reset_applied()
    if (input_deck_len > 0) then
      n = min(int(input_deck_len), len(applied_input_deck))
      do ib = 1, n
        if (input_deck(ib) == c_null_char) exit
        applied_input_deck(ib:ib) = input_deck(ib)
      end do
    end if
    if (.not. c_associated(params_capnp) .or. params_capnp_size <= 0) return
    n = int(params_capnp_size)
    call c_f_pointer(params_capnp, raw, [n])
    allocate (bytes(0:n - 1))
    bytes(0:n - 1) = raw(1:n)
    call capnp_deserialize_bytes(bytes, msg, err)
    if (err /= CAPNP_OK) then
      deallocate (bytes)
      return
    end if
    params = c_p_m_d_params_read_root(msg, err)
    if (err /= CAPNP_OK) then
      call capnp_message_free(msg)
      deallocate (bytes)
      return
    end if

    ! Top-level effective config (matches C cpmdc_params_effective_config).
    call c_p_m_d_params_functional_get(params, s, err)
    if (err == CAPNP_OK) call copy_alloc_text(s, applied_functional, 64)
    if (len_trim(applied_functional) == 0) applied_functional = 'BLYP'
    cut = c_p_m_d_params_cut_off_ry_get(params)
    applied_cutoff_ry = cut
    if (applied_cutoff_ry <= 0.0_real64) applied_cutoff_ry = 70.0_real64
    applied_charge = int(c_p_m_d_params_charge_get(params))
    applied_mult = max(1, int(c_p_m_d_params_multiplicity_get(params)))
    call c_p_m_d_params_cpmd_root_get(params, s, err)
    if (err == CAPNP_OK) call copy_alloc_text(s, applied_cpmd_root, 1024)

    ! Section overrides: system (cutoff/charge/mult) and dft (functional).
    sections = c_p_m_d_params_input_sections_get(params, err)
    n64 = 0_int64
    if (err == CAPNP_OK) n64 = capnp_list_len(sections)
    do i = 0, int(n64) - 1
      sec = c_p_m_d_params_input_sections_get_elem(params, i, err)
      if (err /= CAPNP_OK) cycle
      tag = c_p_m_d_input_section_which(sec)
      if (tag == C_P_M_D_INPUT_SECTION_SYSTEM_TAG) then
        sys = c_p_m_d_input_section_system_get(sec, err)
        if (err /= CAPNP_OK .or. sys%p%kind == CAPNP_PK_NULL) cycle
        cut = c_p_m_d_system_section_cut_off_ry_get(sys)
        if (cut > 0.0_real64) applied_cutoff_ry = cut
        applied_charge = int(c_p_m_d_system_section_charge_get(sys))
        n = int(c_p_m_d_system_section_multiplicity_get(sys))
        if (n > 0) applied_mult = n
      else if (tag == C_P_M_D_INPUT_SECTION_DFT_TAG) then
        dft = c_p_m_d_input_section_dft_get(sec, err)
        if (err /= CAPNP_OK .or. dft%p%kind == CAPNP_PK_NULL) cycle
        call c_p_m_d_dft_section_functional_get(dft, s, err)
        if (err == CAPNP_OK .and. allocated(s)) then
          if (len_trim(s) > 0) call copy_alloc_text(s, applied_functional, 64)
        end if
      end if
    end do

    ! CommonMethodSpec overlay scalars (after section walk).
    if (functional_ov_len > 0) then
      call cchars_to_f(functional_ov, functional_ov_len, fov)
      if (len_trim(fov) > 0) applied_functional = fov
    end if
    if (cutoff_ov > 0.0_c_double) applied_cutoff_ry = real(cutoff_ov, real64)
    if (has_charge_ov /= 0) applied_charge = int(charge_ov)
    if (has_mult_ov /= 0) applied_mult = max(1, int(mult_ov))

    call capnp_message_free(msg)
    deallocate (bytes)
    rc = 0_c_int
  end function cpmdc_embed_apply_params

  function cpmdc_embed_get_config(functional, functional_len, cutoff_ry, &
      charge, mult, input_deck, input_deck_len, cpmd_root, cpmd_root_len) &
      result(rc) bind(C, name='cpmdc_embed_get_config')
    character(kind=c_char), intent(out) :: functional(*)
    integer(c_int), intent(in), value :: functional_len
    real(c_double), intent(out) :: cutoff_ry
    integer(c_int), intent(out) :: charge, mult
    character(kind=c_char), intent(out) :: input_deck(*)
    integer(c_int), intent(in), value :: input_deck_len
    character(kind=c_char), intent(out) :: cpmd_root(*)
    integer(c_int), intent(in), value :: cpmd_root_len
    integer(c_int) :: rc
    call f_to_c_chars(applied_functional, functional, functional_len)
    cutoff_ry = real(applied_cutoff_ry, c_double)
    charge = int(applied_charge, c_int)
    mult = int(applied_mult, c_int)
    call f_to_c_chars(applied_input_deck, input_deck, input_deck_len)
    call f_to_c_chars(applied_cpmd_root, cpmd_root, cpmd_root_len)
    rc = 0_c_int
  end function cpmdc_embed_get_config

end module cpmdc_embed_apply_params_mod
