module mpi_mod

  use :: mpi
  ! include 'mpif.h'

  implicit none

contains

  subroutine forced_exit(id)

    integer,intent(in) :: id
    integer :: ierr
  
    print '(i3," forced exit of rank: ")',id
    call MPI_Finalize(ierr)
    stop
  
  end subroutine forced_exit

end module mpi_mod
