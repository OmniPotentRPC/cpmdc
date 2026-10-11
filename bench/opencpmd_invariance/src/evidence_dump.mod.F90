! Evidence only, not for upstream as is.
! Writes raw f64 results when CPMDC_EVIDENCE_BIN is set.
! fion is the nuclear gradient (printed as GRADIENTS (-FORCES)).
! The file stores the physical force, which is -fion.
! Stress is paiu, column major.

module evidence_dump
  use ener, only: ener_com
  use ions, only: ions0, ions1
  use kinds, only: real_8
  use norm, only: cnorm, gemax
  use parac, only: parai, paral
  use strs, only: paiu
  use system, only: ncpw, parm, parap
  implicit none
  public :: evidence_write_results
  public :: evidence_log_step
contains

  subroutine evidence_write_results(fion)
    real(real_8), intent(in) :: fion(:, :, :)
    character(len=1024) :: path
    character(len=16) :: nm
    character(len=16), parameter :: names(24) = (/ &
         'etot            ', 'ekin            ', 'epseu           ', &
         'enl             ', 'eht             ', 'ehep            ', &
         'ehee            ', 'ehii            ', 'exc             ', &
         'vxc             ', 'egc             ', 'esr             ', &
         'eself           ', 'eeig            ', 'eband           ', &
         'entropy         ', 'eext            ', 'ehsic           ', &
         'ecnstr          ', 'erestr          ', 'eefield         ', &
         'ebogo           ', 'amu             ', 'etddft          ' /)
    real(real_8) :: values(24)
    integer :: iu, is, ia, k, r, length, stat, nforce
    integer(4) :: i4
    logical :: opened

    call get_environment_variable('CPMDC_EVIDENCE_BIN', path, length, stat)
    if (stat /= 0 .or. length <= 0) return
    if (.not. paral%io_parent) return

    values(1) = ener_com%etot
    values(2) = ener_com%ekin
    values(3) = ener_com%epseu
    values(4) = ener_com%enl
    values(5) = ener_com%eht
    values(6) = ener_com%ehep
    values(7) = ener_com%ehee
    values(8) = ener_com%ehii
    values(9) = ener_com%exc
    values(10) = ener_com%vxc
    values(11) = ener_com%egc
    values(12) = ener_com%esr
    values(13) = ener_com%eself
    values(14) = ener_com%eeig
    values(15) = ener_com%eband
    values(16) = ener_com%entropy
    values(17) = ener_com%eext
    values(18) = ener_com%ehsic
    values(19) = ener_com%ecnstr
    values(20) = ener_com%erestr
    values(21) = ener_com%eefield
    values(22) = ener_com%ebogo
    values(23) = ener_com%amu
    values(24) = ener_com%etddft

    nforce = 0
    do is = 1, ions1%nsp
      nforce = nforce + ions0%na(is) * 3
    end do

    open(newunit=iu, file=trim(path(1:length)), access='stream', &
         form='unformatted', status='replace', iostat=stat)
    if (stat /= 0) return
    opened = .true.
    i4 = int(z'45564431', 4)
    write(iu) i4
    i4 = 1_4
    write(iu) i4
    i4 = int(ions1%nat, 4)
    write(iu) i4
    i4 = 24_4
    write(iu) i4
    do k = 1, 24
      nm = names(k)
      write(iu) nm
      write(iu) values(k)
    end do
    i4 = int(nforce, 4)
    write(iu) i4
    do is = 1, ions1%nsp
      do ia = 1, ions0%na(is)
        do k = 1, 3
          write(iu) -fion(k, ia, is)
        end do
      end do
    end do
    i4 = 9_4
    write(iu) i4
    do ia = 1, 3
      do is = 1, 3
        write(iu) paiu(is, ia)
      end do
    end do
    write(iu) parm%omega
    write(iu) gemax
    write(iu) cnorm
    i4 = int(parai%nproc, 4)
    write(iu) i4
    if (allocated(parap%sparm) .and. allocated(parap%nrxpl)) then
      do r = 0, parai%nproc - 1
        i4 = int(parap%sparm(3, r), 4)
        write(iu) i4
        i4 = int(parap%nrxpl(r, 1), 4)
        write(iu) i4
        i4 = int(parap%nrxpl(r, 2), 4)
        write(iu) i4
      end do
    else
      do r = 0, parai%nproc - 1
        i4 = int(ncpw%ngw, 4)
        write(iu) i4
        i4 = 0_4
        write(iu) i4
        write(iu) i4
      end do
    end if
    close(iu)
    if (.not. opened) return
  end subroutine evidence_write_results

  subroutine evidence_log_step(tcpu)
    real(real_8), intent(in) :: tcpu
    character(len=1024) :: path
    integer :: iu, length, stat

    call get_environment_variable('CPMDC_EVIDENCE_STEPS', path, length, stat)
    if (stat /= 0 .or. length <= 0) return
    if (.not. paral%io_parent) return
    open(newunit=iu, file=trim(path(1:length)), access='stream', &
         form='unformatted', status='old', position='append', iostat=stat)
    if (stat /= 0) then
      open(newunit=iu, file=trim(path(1:length)), access='stream', &
           form='unformatted', status='replace', iostat=stat)
      if (stat /= 0) return
    end if
    write(iu) tcpu
    close(iu)
  end subroutine evidence_log_step

end module evidence_dump
