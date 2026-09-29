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

contains

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
    n = min(len_trim(fstr), max(0, int(clen)))
    do i = 1, n
      cbuf(i) = fstr(i:i)
    end do
    if (n < int(clen)) cbuf(n + 1) = c_null_char
  end subroutine f_to_c_chars

  !> Decode CPMDParams wire bytes into embed knobs; accept C-rendered deck.
  !> Optional common-overlay scalars: empty functional / cutoff<=0 / has_* =0
  !> leave the wire effective config unchanged for that field.
  function cpmdc_embed_apply_params(params_capnp, params_capnp_size, &
      input_deck, input_deck_len, functional_ov, functional_ov_len, &
      cutoff_ov, has_charge_ov, charge_ov, has_mult_ov, mult_ov, &
      functional_out, functional_cap, cutoff_out, charge_out, mult_out, &
      deck_out, deck_cap, root_out, root_cap) result(rc) &
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
    character(kind=c_char), intent(out) :: functional_out(*)
    integer(c_int), intent(in), value :: functional_cap
    real(c_double), intent(out) :: cutoff_out
    integer(c_int), intent(out) :: charge_out, mult_out
    character(kind=c_char), intent(out) :: deck_out(*)
    integer(c_int), intent(in), value :: deck_cap
    character(kind=c_char), intent(out) :: root_out(*)
    integer(c_int), intent(in), value :: root_cap
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
    character(len=64) :: functional_l
    real(real64) :: cutoff_l
    integer :: charge_l, mult_l
    character(len=4096) :: deck_l
    character(len=1024) :: root_l

    rc = -1_c_int
    functional_l = 'BLYP'
    cutoff_l = 70.0_real64
    charge_l = 0
    mult_l = 1
    deck_l = ' '
    root_l = ' '
    if (input_deck_len > 0) then
      n = min(int(input_deck_len), len(deck_l))
      do ib = 1, n
        if (input_deck(ib) == c_null_char) exit
        deck_l(ib:ib) = input_deck(ib)
      end do
    end if
    if (.not. c_associated(params_capnp) .or. params_capnp_size <= 0) return
    if (params_capnp_size > int(huge(n), kind=c_size_t)) return
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
    if (err == CAPNP_OK) call copy_alloc_text(s, functional_l, 64)
    if (len_trim(functional_l) == 0) functional_l = 'BLYP'
    cut = c_p_m_d_params_cut_off_ry_get(params)
    cutoff_l = cut
    if (cutoff_l <= 0.0_real64) cutoff_l = 70.0_real64
    charge_l = int(c_p_m_d_params_charge_get(params))
    mult_l = max(1, int(c_p_m_d_params_multiplicity_get(params)))
    call c_p_m_d_params_cpmd_root_get(params, s, err)
    if (err == CAPNP_OK) call copy_alloc_text(s, root_l, 1024)

    ! Section overrides: system (cutoff/charge/mult) and dft (functional).
    sections = c_p_m_d_params_input_sections_get(params, err)
    n64 = 0_int64
    if (err == CAPNP_OK) n64 = capnp_list_len(sections)
    ! INT of a value that does not fit the default integer is undefined.
    if (n64 < 0_int64 .or. n64 > int(huge(i), kind=int64)) then
      call capnp_message_free(msg)
      deallocate (bytes)
      return
    end if
    do i = 0, int(n64) - 1
      sec = c_p_m_d_params_input_sections_get_elem(params, i, err)
      if (err /= CAPNP_OK) cycle
      tag = c_p_m_d_input_section_which(sec)
      if (tag == C_P_M_D_INPUT_SECTION_SYSTEM_TAG) then
        sys = c_p_m_d_input_section_system_get(sec, err)
        if (err /= CAPNP_OK .or. sys%p%kind == CAPNP_PK_NULL) cycle
        cut = c_p_m_d_system_section_cut_off_ry_get(sys)
        if (cut > 0.0_real64) cutoff_l = cut
        charge_l = int(c_p_m_d_system_section_charge_get(sys))
        n = int(c_p_m_d_system_section_multiplicity_get(sys))
        if (n > 0) mult_l = n
      else if (tag == C_P_M_D_INPUT_SECTION_DFT_TAG) then
        dft = c_p_m_d_input_section_dft_get(sec, err)
        if (err /= CAPNP_OK .or. dft%p%kind == CAPNP_PK_NULL) cycle
        call c_p_m_d_dft_section_functional_get(dft, s, err)
        if (err == CAPNP_OK .and. allocated(s)) then
          if (len_trim(s) > 0) call copy_alloc_text(s, functional_l, 64)
        end if
      end if
    end do

    ! CommonMethodSpec overlay scalars (after section walk).
    if (functional_ov_len > 0) then
      call cchars_to_f(functional_ov, functional_ov_len, fov)
      if (len_trim(fov) > 0) functional_l = fov
    end if
    if (cutoff_ov > 0.0_c_double) cutoff_l = real(cutoff_ov, real64)
    if (has_charge_ov /= 0) charge_l = int(charge_ov)
    if (has_mult_ov /= 0) mult_l = max(1, int(mult_ov))

    call capnp_message_free(msg)
    deallocate (bytes)
    call f_to_c_chars(functional_l, functional_out, functional_cap)
    cutoff_out = real(cutoff_l, c_double)
    charge_out = int(charge_l, c_int)
    mult_out = int(mult_l, c_int)
    call f_to_c_chars(deck_l, deck_out, deck_cap)
    call f_to_c_chars(root_l, root_out, root_cap)
    rc = 0_c_int
  end function cpmdc_embed_apply_params

end module cpmdc_embed_apply_params_mod
