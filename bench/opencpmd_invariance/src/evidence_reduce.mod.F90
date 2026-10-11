! Evidence only, not for upstream as is.
! Rank-ordered gather and an int128 fixed-point MPI reduction for mp_sum.
! The fixed point is x * 2^90. Bits that do not fit, including
! subnormals, are summed afterwards in rank order as an f64 residual.
! A one-rank reduction therefore returns the input bits.

#include "evidence_flags.h"

module evidence_reduce
#ifdef __PARALLEL
  use mpi
#endif
  implicit none
  public :: evidence_allreduce_r8
  public :: evidence_allreduce_r8_s
  public :: evidence_allreduce_r8_inplace
  public :: evidence_allreduce_c8
  public :: evidence_allreduce_c8_s
  public :: evidence_allreduce_c8_inplace
  public :: evidence_reduce_r8
  public :: evidence_reduce_r8_s
  public :: evidence_reduce_c8
  public :: evidence_reduce_c8_s
  public :: evidence_gorder_sum
  public :: evidence_pack_f64
  public :: evidence_unpack_i128
  integer, parameter :: evidence_i16 = selected_int_kind(38)
  integer, parameter :: evidence_scale_bias = 985
  integer, parameter :: evidence_tag = 44021
contains

  function evidence_pack_f64(x) result(q)
    real(8), intent(in) :: x
    integer(evidence_i16) :: q
    integer(8) :: bits, frac
    integer :: expb, sh, sgn
    bits = transfer(x, 0_8)
    sgn = int(shiftr(bits, 63))
    expb = int(ibits(bits, 52, 11))
    frac = ibits(bits, 0, 52)
    q = 0_evidence_i16
    if (expb == 0 .or. expb == 2047) return
    frac = ior(frac, shiftl(1_8, 52))
    sh = expb - evidence_scale_bias
    if (sh >= 0) then
      if (sh > 74) return
      q = shiftl(int(frac, evidence_i16), sh)
    else if (sh > -53) then
      q = shiftr(int(frac, evidence_i16), -sh)
    else
      return
    end if
    if (sgn == 1) q = -q
  end function evidence_pack_f64

  function evidence_unpack_i128(q) result(x)
    integer(evidence_i16), intent(in) :: q
    real(8) :: x
    integer(evidence_i16) :: aq, mant
    integer(8) :: bits, frac
    integer :: p, sh, expb, rb, sticky, sgn
    if (q == 0_evidence_i16) then
      x = 0.0d0
      return
    end if
    sgn = 0
    aq = q
    if (q < 0_evidence_i16) then
      sgn = 1
      aq = -q
    end if
    p = 127 - leadz(aq)
    sh = p - 52
    if (sh >= 0) then
      mant = shiftr(aq, sh)
    else
      mant = shiftl(aq, -sh)
    end if
    rb = 0
    sticky = 0
    if (sh > 0) then
      rb = int(ibits(aq, sh - 1, 1))
      if (sh > 1) then
        if (iand(aq, shiftl(1_evidence_i16, sh - 1) - 1_evidence_i16) /= 0_evidence_i16) sticky = 1
      end if
    end if
    if (rb == 1) then
      if (sticky == 1 .or. iand(mant, 1_evidence_i16) == 1_evidence_i16) mant = mant + 1_evidence_i16
    end if
    if (mant >= shiftl(1_evidence_i16, 53)) then
      mant = shiftr(mant, 1)
      p = p + 1
    end if
    frac = int(iand(mant, shiftl(1_evidence_i16, 52) - 1_evidence_i16), 8)
    expb = p - 90 + 1023
    if (expb <= 0) then
      x = 0.0d0
      return
    else if (expb >= 2047) then
      bits = ior(shiftl(int(sgn, 8), 63), shiftl(2047_8, 52))
      x = transfer(bits, x)
      return
    end if
    bits = ior(shiftl(int(sgn, 8), 63), ior(shiftl(int(expb, 8), 52), frac))
    x = transfer(bits, x)
  end function evidence_unpack_i128

#ifdef __PARALLEL
  subroutine evidence_ensure_exact(dtype, op, ierr)
    integer, intent(out) :: dtype, op, ierr
    integer, save :: saved_dtype = -1, saved_op = -1
    logical, save :: ready = .false.
    external evidence_sum_i128
    if (.not. ready) then
      call mpi_type_contiguous(2, MPI_INTEGER8, saved_dtype, ierr)
      if (ierr /= 0) return
      call mpi_type_commit(saved_dtype, ierr)
      if (ierr /= 0) return
      call mpi_op_create(evidence_sum_i128, .true., saved_op, ierr)
      if (ierr /= 0) return
      ready = .true.
    end if
    dtype = saved_dtype
    op = saved_op
    ierr = 0
  end subroutine evidence_ensure_exact

  subroutine evidence_exact_allreduce_r8(sendbuf, recvbuf, n, comm, ierr)
    real(8), intent(in) :: sendbuf(*)
    real(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
    integer(evidence_i16), allocatable :: qs(:), qr(:)
    real(8), allocatable :: residual(:), summed(:)
    integer :: dtype, op, i
    real(8) :: y
    if (n <= 0) then
      ierr = 0
      return
    end if
    call evidence_ensure_exact(dtype, op, ierr)
    if (ierr /= 0) return
    allocate(qs(n), qr(n), residual(n), summed(n))
    do i = 1, n
      qs(i) = evidence_pack_f64(sendbuf(i))
      y = evidence_unpack_i128(qs(i))
      residual(i) = sendbuf(i) - y
    end do
    call mpi_allreduce(qs, qr, n, dtype, op, comm, ierr)
    if (ierr /= 0) return
    call evidence_gather_allreduce_r8(residual, summed, n, comm, ierr)
    if (ierr /= 0) return
    do i = 1, n
      recvbuf(i) = evidence_unpack_i128(qr(i)) + summed(i)
    end do
  end subroutine evidence_exact_allreduce_r8

  subroutine evidence_gather_allreduce_r8(sendbuf, recvbuf, n, comm, ierr)
    real(8), intent(in) :: sendbuf(*)
    real(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
    real(8), allocatable :: tmp(:)
    integer :: me, np, r, i
    integer :: stat(MPI_STATUS_SIZE)
    call mpi_comm_rank(comm, me, ierr)
    if (ierr /= 0) return
    call mpi_comm_size(comm, np, ierr)
    if (ierr /= 0) return
    if (n <= 0) return
    allocate(tmp(n))
    if (me == 0) then
      do i = 1, n
        recvbuf(i) = sendbuf(i)
      end do
      do r = 1, np - 1
        call mpi_recv(tmp, n, MPI_DOUBLE_PRECISION, r, evidence_tag, comm, stat, ierr)
        if (ierr /= 0) return
        do i = 1, n
          recvbuf(i) = recvbuf(i) + tmp(i)
        end do
      end do
    else
      call mpi_send(sendbuf, n, MPI_DOUBLE_PRECISION, 0, evidence_tag, comm, ierr)
      if (ierr /= 0) return
    end if
    call mpi_bcast(recvbuf, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  end subroutine evidence_gather_allreduce_r8

  subroutine evidence_gather_reduce_r8(sendbuf, recvbuf, n, root, comm, ierr)
    real(8), intent(in) :: sendbuf(*)
    real(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, root, comm
    integer, intent(out) :: ierr
    real(8), allocatable :: tmp(:)
    integer :: me, np, r, i
    integer :: stat(MPI_STATUS_SIZE)
    call mpi_comm_rank(comm, me, ierr)
    if (ierr /= 0) return
    call mpi_comm_size(comm, np, ierr)
    if (ierr /= 0) return
    if (n <= 0) return
    allocate(tmp(n))
    if (me == root) then
      do i = 1, n
        recvbuf(i) = sendbuf(i)
      end do
      do r = 0, np - 1
        if (r == root) cycle
        call mpi_recv(tmp, n, MPI_DOUBLE_PRECISION, r, evidence_tag, comm, stat, ierr)
        if (ierr /= 0) return
        do i = 1, n
          recvbuf(i) = recvbuf(i) + tmp(i)
        end do
      end do
    else
      call mpi_send(sendbuf, n, MPI_DOUBLE_PRECISION, root, evidence_tag, comm, ierr)
    end if
  end subroutine evidence_gather_reduce_r8
#endif

  subroutine evidence_allreduce_r8(sendbuf, recvbuf, n, comm, ierr)
    real(8), intent(in) :: sendbuf(*)
    real(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
    integer :: i
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT)
    call evidence_exact_allreduce_r8(sendbuf, recvbuf, n, comm, ierr)
#elif defined(CPMDC_EVIDENCE_MP_GATHER)
    call evidence_gather_allreduce_r8(sendbuf, recvbuf, n, comm, ierr)
#else
    call mpi_allreduce(sendbuf, recvbuf, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
#else
    ierr = 0
    do i = 1, n
      recvbuf(i) = sendbuf(i)
    end do
#endif
  end subroutine evidence_allreduce_r8

  subroutine evidence_allreduce_r8_s(sendbuf, recvbuf, comm, ierr)
    real(8), intent(in) :: sendbuf
    real(8), intent(inout) :: recvbuf
    integer, intent(in) :: comm
    integer, intent(out) :: ierr
    real(8) :: send1(1), recv1(1)
    send1(1) = sendbuf
    call evidence_allreduce_r8(send1, recv1, 1, comm, ierr)
    recvbuf = recv1(1)
  end subroutine evidence_allreduce_r8_s

  subroutine evidence_allreduce_r8_inplace(buf, n, comm, ierr)
    real(8), intent(inout) :: buf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT)
    call evidence_exact_allreduce_r8(buf, buf, n, comm, ierr)
#elif defined(CPMDC_EVIDENCE_MP_GATHER)
    block
      real(8), allocatable :: local(:)
      integer :: i
      allocate(local(max(n, 1)))
      do i = 1, n
        local(i) = buf(i)
      end do
      call evidence_gather_allreduce_r8(local, buf, n, comm, ierr)
    end block
#else
    call mpi_allreduce(MPI_IN_PLACE, buf, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
#endif
#else
    ierr = 0
#endif
  end subroutine evidence_allreduce_r8_inplace

  subroutine evidence_copy_c8_to_r8(z, r, n)
    complex(8), intent(in) :: z(*)
    real(8), intent(out) :: r(*)
    integer, intent(in) :: n
    integer :: i
    do i = 1, n
      r(2 * i - 1) = real(z(i), 8)
      r(2 * i) = aimag(z(i))
    end do
  end subroutine evidence_copy_c8_to_r8

  subroutine evidence_copy_r8_to_c8(r, z, n)
    real(8), intent(in) :: r(*)
    complex(8), intent(out) :: z(*)
    integer, intent(in) :: n
    integer :: i
    do i = 1, n
      z(i) = cmplx(r(2 * i - 1), r(2 * i), 8)
    end do
  end subroutine evidence_copy_r8_to_c8

  subroutine evidence_allreduce_c8(sendbuf, recvbuf, n, comm, ierr)
    complex(8), intent(in) :: sendbuf(*)
    complex(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
    integer :: i
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT) || defined(CPMDC_EVIDENCE_MP_GATHER)
    block
      real(8), allocatable :: sr(:), rr(:)
      allocate(sr(max(2 * n, 1)), rr(max(2 * n, 1)))
      call evidence_copy_c8_to_r8(sendbuf, sr, n)
      call evidence_allreduce_r8(sr, rr, 2 * n, comm, ierr)
      if (ierr == 0) call evidence_copy_r8_to_c8(rr, recvbuf, n)
    end block
#else
    call mpi_allreduce(sendbuf, recvbuf, n, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
#else
    ierr = 0
    do i = 1, n
      recvbuf(i) = sendbuf(i)
    end do
#endif
  end subroutine evidence_allreduce_c8

  subroutine evidence_allreduce_c8_s(sendbuf, recvbuf, comm, ierr)
    complex(8), intent(in) :: sendbuf
    complex(8), intent(inout) :: recvbuf
    integer, intent(in) :: comm
    integer, intent(out) :: ierr
    complex(8) :: send1(1), recv1(1)
    send1(1) = sendbuf
    call evidence_allreduce_c8(send1, recv1, 1, comm, ierr)
    recvbuf = recv1(1)
  end subroutine evidence_allreduce_c8_s

  subroutine evidence_allreduce_c8_inplace(buf, n, comm, ierr)
    complex(8), intent(inout) :: buf(*)
    integer, intent(in) :: n, comm
    integer, intent(out) :: ierr
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT) || defined(CPMDC_EVIDENCE_MP_GATHER)
    block
      complex(8), allocatable :: local(:)
      integer :: i
      allocate(local(max(n, 1)))
      do i = 1, n
        local(i) = buf(i)
      end do
      call evidence_allreduce_c8(local, buf, n, comm, ierr)
    end block
#else
    call mpi_allreduce(MPI_IN_PLACE, buf, n, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
#endif
#else
    ierr = 0
#endif
  end subroutine evidence_allreduce_c8_inplace

  subroutine evidence_reduce_r8(sendbuf, recvbuf, n, root, comm, ierr)
    real(8), intent(in) :: sendbuf(*)
    real(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, root, comm
    integer, intent(out) :: ierr
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT)
    block
      real(8), allocatable :: tmp(:)
      integer :: me, i
      allocate(tmp(max(n, 1)))
      call evidence_exact_allreduce_r8(sendbuf, tmp, n, comm, ierr)
      call mpi_comm_rank(comm, me, ierr)
      if (me == root) then
        do i = 1, n
          recvbuf(i) = tmp(i)
        end do
      end if
    end block
#elif defined(CPMDC_EVIDENCE_MP_GATHER)
    call evidence_gather_reduce_r8(sendbuf, recvbuf, n, root, comm, ierr)
#else
    call mpi_reduce(sendbuf, recvbuf, n, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, ierr)
#endif
#else
    integer :: i
    ierr = 0
    do i = 1, n
      recvbuf(i) = sendbuf(i)
    end do
#endif
  end subroutine evidence_reduce_r8

  subroutine evidence_reduce_r8_s(sendbuf, recvbuf, root, comm, ierr)
    real(8), intent(in) :: sendbuf
    real(8), intent(inout) :: recvbuf
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    real(8) :: send1(1), recv1(1)
    send1(1) = sendbuf
    call evidence_reduce_r8(send1, recv1, 1, root, comm, ierr)
    recvbuf = recv1(1)
  end subroutine evidence_reduce_r8_s

  subroutine evidence_reduce_c8(sendbuf, recvbuf, n, root, comm, ierr)
    complex(8), intent(in) :: sendbuf(*)
    complex(8), intent(inout) :: recvbuf(*)
    integer, intent(in) :: n, root, comm
    integer, intent(out) :: ierr
#ifdef __PARALLEL
#if defined(CPMDC_EVIDENCE_MP_EXACT) || defined(CPMDC_EVIDENCE_MP_GATHER)
    block
      real(8), allocatable :: sr(:), rr(:)
      integer :: me, i
      allocate(sr(max(2 * n, 1)), rr(max(2 * n, 1)))
      call evidence_copy_c8_to_r8(sendbuf, sr, n)
      call evidence_reduce_r8(sr, rr, 2 * n, root, comm, ierr)
      call mpi_comm_rank(comm, me, ierr)
      if (me == root .and. ierr == 0) then
        call evidence_copy_r8_to_c8(rr, recvbuf, n)
      end if
    end block
#else
    call mpi_reduce(sendbuf, recvbuf, n, MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, ierr)
#endif
#else
    integer :: i
    ierr = 0
    do i = 1, n
      recvbuf(i) = sendbuf(i)
    end do
#endif
  end subroutine evidence_reduce_c8

  subroutine evidence_reduce_c8_s(sendbuf, recvbuf, root, comm, ierr)
    complex(8), intent(in) :: sendbuf
    complex(8), intent(inout) :: recvbuf
    integer, intent(in) :: root, comm
    integer, intent(out) :: ierr
    complex(8) :: send1(1), recv1(1)
    send1(1) = sendbuf
    call evidence_reduce_c8(send1, recv1, 1, root, comm, ierr)
    recvbuf = recv1(1)
  end subroutine evidence_reduce_c8_s

  subroutine evidence_gorder_sum(n, key, val, comm, total)
    integer, intent(in) :: n, comm
    integer(8), intent(in) :: key(*)
    real(8), intent(in) :: val(*)
    real(8), intent(out) :: total
#ifdef __PARALLEL
    integer :: ierr, np, i, nt
    integer, allocatable :: counts(:), displs(:)
    integer(8), allocatable :: allkey(:)
    real(8), allocatable :: allval(:)
    integer, allocatable :: ord(:)
    integer(8) :: ksend(1)
    real(8) :: vsend(1)
    call mpi_comm_size(comm, np, ierr)
    if (ierr /= 0) then
      total = 0.0d0
      return
    end if
    allocate(counts(0:np - 1), displs(0:np - 1))
    call mpi_allgather(n, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, comm, ierr)
    displs(0) = 0
    do i = 1, np - 1
      displs(i) = displs(i - 1) + counts(i - 1)
    end do
    nt = displs(np - 1) + counts(np - 1)
    allocate(allkey(max(nt, 1)), allval(max(nt, 1)), ord(max(nt, 1)))
    if (n > 0) then
      call mpi_allgatherv(key, n, MPI_INTEGER8, allkey, counts, displs, MPI_INTEGER8, comm, ierr)
      call mpi_allgatherv(val, n, MPI_DOUBLE_PRECISION, allval, counts, displs, MPI_DOUBLE_PRECISION, comm, ierr)
    else
      ksend(1) = 0_8
      vsend(1) = 0.0d0
      call mpi_allgatherv(ksend, 0, MPI_INTEGER8, allkey, counts, displs, MPI_INTEGER8, comm, ierr)
      call mpi_allgatherv(vsend, 0, MPI_DOUBLE_PRECISION, allval, counts, displs, MPI_DOUBLE_PRECISION, comm, ierr)
    end if
    call evidence_order_keys(nt, allkey, ord)
    total = 0.0d0
    do i = 1, nt
      total = total + allval(ord(i))
    end do
#else
    integer :: i
    integer, allocatable :: ord(:)
    allocate(ord(max(n, 1)))
    call evidence_order_keys(n, key, ord)
    total = 0.0d0
    do i = 1, n
      total = total + val(ord(i))
    end do
#endif
  end subroutine evidence_gorder_sum

  subroutine evidence_order_keys(n, key, ord)
    integer, intent(in) :: n
    integer(8), intent(in) :: key(*)
    integer, intent(out) :: ord(*)
    integer :: i, nn, tmp
    do i = 1, n
      ord(i) = i
    end do
    if (n <= 1) return
    do i = n / 2, 1, -1
      tmp = ord(i)
      call evidence_sift(key, ord, i, n, tmp)
    end do
    nn = n
    do while (nn > 1)
      tmp = ord(nn)
      ord(nn) = ord(1)
      nn = nn - 1
      call evidence_sift(key, ord, 1, nn, tmp)
    end do
  end subroutine evidence_order_keys

  subroutine evidence_sift(key, ord, lo, hi, tmp)
    integer(8), intent(in) :: key(*)
    integer, intent(inout) :: ord(*)
    integer, intent(in) :: lo, hi, tmp
    integer :: root, child
    root = lo
    do
      child = root * 2
      if (child > hi) exit
      if (child < hi) then
        if (evidence_key_before(key, ord(child + 1), ord(child))) child = child + 1
      end if
      if (.not. evidence_key_before(key, ord(child), tmp)) exit
      ord(root) = ord(child)
      root = child
    end do
    ord(root) = tmp
  end subroutine evidence_sift

  function evidence_key_before(key, ia, ib) result(yes)
    integer(8), intent(in) :: key(*)
    integer, intent(in) :: ia, ib
    logical :: yes
    if (key(ia) /= key(ib)) then
      yes = key(ia) > key(ib)
    else
      yes = ia > ib
    end if
  end function evidence_key_before

end module evidence_reduce

#ifdef __PARALLEL
subroutine evidence_sum_i128(invec, inoutvec, len, datatype)
  use evidence_reduce, only: evidence_i16
  implicit none
  integer :: len, datatype
  integer(8) :: invec(2, len), inoutvec(2, len)
  integer(evidence_i16) :: a, b
  integer :: i
  do i = 1, len
    a = transfer(invec(:, i), a)
    b = transfer(inoutvec(:, i), b)
    b = a + b
    inoutvec(:, i) = transfer(b, inoutvec(:, i))
  end do
end subroutine evidence_sum_i128
#endif
