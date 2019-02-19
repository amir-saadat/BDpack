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

  use :: dlt_mod, only: dlt_bs
  use :: smdlt_mod, only: smdlt_bs
  use :: mpi
#ifdef USE_GPU
  use :: dev_cumod, only: init_dev
#ifdef USE_MAGMA
  use :: magma_cumod, only: init_magma
#endif
#endif

  implicit none


  ! MPI variables
  integer :: ierr,p,id,narg
  character(len=20) :: inpFile
  character(len=20) :: driver


  ! Initialize MPI
  call MPI_Init(ierr)
  ! Get the number of processes
  call MPI_Comm_size(MPI_COMM_WORLD,p,ierr)
  ! Get the individual process ID
  call MPI_Comm_rank(MPI_COMM_WORLD,id,ierr)

#ifdef USE_GPU
  ! Initialize GPU devices
  call init_dev()
#ifdef USE_MAGMA
  call init_magma()
#endif
#endif

  narg=command_argument_count()
  if (narg /= 0) then
    call get_command_argument(1,inpFile)
  else
    inpFile='input.dat'
  end if

  call read_inp()

  select case (driver)
    case ('dilute_bs')
      call dlt_bs(p,id)
    case ('semidilute_bs')
      call smdlt_bs(p,id)
    case default
      print '(" Error in main: Incorrect driver.")'
      stop
  end select

  ! Finalizing MPI
  call MPI_Finalize(ierr)

  if (id == 0) then
    print *
    print '(" *** NORMAL END OF THE PROGRAM ***")'
    select case (driver)
      case ('dilute_bs')
        print *
        print '("%---------------------------------------------------------------%")'
        print '(" | ***Please cite the following article if BDpack washelpful*** |")'
        print '(" |                                                              |")'
        print '(" |  A. Saadat and B. Khomami, J. Chem. Phys., 2014, 140, 184903 |")'
        print '(" |                                                              |")'
        print '("%---------------------------------------------------------------%")'
      case ('semidilute_bs')
        print *
        print '("%---------------------------------------------------------------%")'
        print '(" | ***Please cite the following article if BDpack washelpful*** |")'
        print '(" |                                                              |")'
        print '(" |   A. Saadat and B. Khomami, Phys. Rev. E, 2015, 92, 033307   |")'
        print '(" |                                                              |")'
        print '("%---------------------------------------------------------------%")'
    end select
  end if

contains

  subroutine read_inp()

    use :: prcn_mod
    use,intrinsic :: iso_fortran_env
    use :: strg_mod, only: parse,value
    use :: inp_dlt, only: read_dlt
!    use :: smdlt_inp, only: smdlt_read

    integer :: il,j,ntokens,u1,stat,ios
    real(wp),allocatable :: u(:)
    character(len=1024) :: line
    character(len=100) :: tokens(10)

    ! Default values
    driver='dilute_bs'

    open (newunit=u1,action='read',file=inpFile,status='old')
    il=1
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" io_mod: Error reading line ", i0, " Process ID ", i0)', il,id
        stop
      elseif (line(1:1) == '#') then
        il=il+1
        cycle ef ! commented line
      else
        il=il+1
      end if
      call parse(line,': ',tokens,ntokens)
      if (ntokens > 0) then
        do j=1,ntokens
          if (trim(adjustl(tokens(j))) == 'driver') then
            driver=trim(adjustl(tokens(j+1)))
          end if
        end do ! j
      end if ! ntokens
    end do ef

    close(u1)

    select case (driver)
      case ('dilute_bs')
        call read_dlt(id,inpFile)
      case ('semidilute_bs')
!        call smdlt_read(id,inpFile)
    end select

  end subroutine read_inp

end program BDpack
