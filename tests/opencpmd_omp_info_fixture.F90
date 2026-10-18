module mp_interface
  use mpi
  implicit none
  integer, parameter :: mp_max_processor_name=256
  interface mp_bcast
    module procedure bcast_integer, bcast_integers, bcast_text
  end interface
contains
  subroutine mp_get_node_env(comm,n,task)
    integer,intent(in)::comm
    integer,intent(out)::n,task
    integer::rank,total,split,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,total,ierr)
    split=max(1,total/3)
    if(rank<split)then
      n=split;task=rank
    else
      n=total-split;task=rank-split
    endif
  end subroutine
  subroutine mp_get_processor_name(name)
    character(*),intent(out)::name
    integer::rank,total,ierr
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,total,ierr)
    name='fixture-node-b'
    if(rank<max(1,total/3))name='fixture-node-a'
  end subroutine
  subroutine mp_sum(values,n,comm)
    integer,intent(in)::n,comm
    integer,intent(inout)::values(*)
    integer::low,high,ierr
    call MPI_Allreduce(n,low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(n,high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(low/=high)then
      write(*,*)'inconsistent native report counts:',low,high
      call MPI_Abort(comm,92,ierr)
    endif
    call MPI_Allreduce(MPI_IN_PLACE,values,n,MPI_INTEGER,MPI_SUM,comm,ierr)
  end subroutine
  subroutine bcast_integer(value,root,comm)
    integer,intent(inout)::value
    integer,intent(in)::root,comm
    integer::ierr
    call MPI_Bcast(value,1,MPI_INTEGER,root,comm,ierr)
  end subroutine
  subroutine bcast_integers(values,n,root,comm)
    integer,intent(in)::n,root,comm
    integer,intent(inout)::values(*)
    integer::ierr
    call MPI_Bcast(values,n,MPI_INTEGER,root,comm,ierr)
  end subroutine
  subroutine bcast_text(value,root,comm)
    character(*),intent(inout)::value
    integer,intent(in)::root,comm
    integer::ierr
    call MPI_Bcast(value,len(value),MPI_CHARACTER,root,comm,ierr)
  end subroutine
end module
module parac
  implicit none
  type integer_parameters
    integer::ncpus=1,cp_grp,io_source=0,cp_me
  end type
  type logical_parameters
    logical::io_parent
  end type
  type(integer_parameters)::parai
  type(logical_parameters)::paral
end module
module native_report
  use mp_interface
  use parac
  !$ use omp_lib
  implicit none
contains
#include "native_print_omp_info.inc"
end module
program check_report
  use mpi
  use native_report
  implicit none
  integer::ierr,total
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,total,ierr)
  if(total<3)call MPI_Abort(MPI_COMM_WORLD,93,ierr)
  parai%cp_grp=MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD,parai%cp_me,ierr)
  paral%io_parent=parai%cp_me==parai%io_source
  call print_omp_info()
  call MPI_Finalize(ierr)
end program
