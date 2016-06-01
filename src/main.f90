!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is free software: you can redistribute it and/or modify        |
!|  it under the terms of the GNU General Public License as published by  |
!|  the Free Software Foundation, either version 3 of the License, or     |
!|  (at your option) any later version.                                   |
!|                                                                        |
!|  BDpack is distributed in the hope that it will be useful,             |
!|  but WITHOUT ANY WARRANTY; without even the implied warranty of        |
!|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
!|  GNU General Public License for more details.                          |
!|                                                                        |
!|  You should have received a copy of the GNU General Public License     |
!|  along with BDpack.  If not, see <http://www.gnu.org/licenses/>.       |
!%------------------------------------------------------------------------%
program BDpack

  use :: mpi
  use :: inp_mod, only: read_inp,driver
  use :: dlt_mod, only: dlt_bs

  implicit none

  ! MPI variables
  integer :: ierr,p,id,narg
  character(len=20) :: inpFile

  ! Initialize MPI
  call MPI_Init(ierr)
  ! Get the number of processes
  call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)
  ! Get the individual process ID
  call MPI_Comm_rank(MPI_COMM_WORLD,id,ierr)

  
  narg=command_argument_count()
  if (narg /= 0) then
    call get_command_argument(1,inpFile)
  else
    inpFile='input.dat'
  end if

  call read_inp(id,inpFile)

  select case (driver)
    case ('dilute_bs')
      call dlt_bs(p,id)
    case ('semidilute_bs')
!      call smdlt_bs(p,id)
    case default
      print '(" Error in main: Incorrect driver.")'
      stop
  end select

  ! Finalizing MPI
  call MPI_Finalize(ierr)

  if (id == 0) print *,'*** NORMAL END OF THE PROGRAM ***'

end program BDpack
