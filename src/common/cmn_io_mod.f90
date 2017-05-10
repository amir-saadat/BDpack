!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2017:                                            |
!|  Fluid Mechanics Laboratory (Shaqfeh's Group)                          |
!|  Stanford University                                                   |
!|  Material Research and Innovation Laboratory                           |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat        <asaadat@stanford.edu>                  |
!|  Advisor:   Eric S. G. Shaqfeh <esgs@stanford.edu>                     |
!|             Bamin Khomami      <bkhomami@utk.edu>                      |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is a free software: you can redistribute it and/or modify      |
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
!--------------------------------------------------------------------
!
! MODULE: 
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: 
!> Applies homogeneous flow to the entities inside the box
!--------------------------------------------------------------------

module cmn_io_mod

  use :: strg_mod
  use :: iso_fortran_env
  use :: prcn_mod

  public :: read_input
  interface read_input
    module procedure read_input_int
    module procedure read_input_wp
    module procedure read_input_ch
  end interface read_input

contains

  !> Reads the input for character
  !! \param id the MPI id
  !! \param inp_ph the input phrase used for the variable
  !! \param os offset from the starting argument, first argument for var
  !! \param nvar the number of variables
  !! \param var the actual variables
  subroutine read_input_ch(inp_ph,os,var,def)
 
    integer,intent(in) :: os
    character(len=*),intent(in) :: inp_ph
    character(len=*),intent(inout) :: var
    character(len=*),intent(in),optional :: def
    integer :: i,j,ntokens,u1,il,stat
    character(len=1024) :: line
    character(len=100) :: tokens(10)
    character(len=20) :: inpFile

    narg=command_argument_count()

    if (narg /= 0) then
      call get_command_argument(1,inpFile)
    else
      inpFile='input.dat'
    end if

    open (newunit=u1,action='read',file=inpFile,status='old')
    il=1

    if (present(def)) then
      var=def
    endif
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        if (.not.present(def)) then
          print '(" Error: variable ",a," was not found.")',inp_ph
          stop ! should be changed
        endif
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" cmn_io_mod: Error reading line ", i0)',il
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
          if(trim(adjustl(tokens(j))) == inp_ph) then
            var=trim(adjustl(tokens(j+os+1)))
            print '(" ",a," = ",a)',inp_ph,var
            exit ef
          end if
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

  end subroutine read_input_ch

  !> Reads the input for character
  !! \param id the MPI id
  !! \param inp_ph the input phrase used for the variable
  !! \param os offset from the starting argument, first argument for var
  !! \param var the value of variable
  !! \param def the default value of var
  subroutine read_input_int(inp_ph,os,var,def)
    
    integer,intent(in) :: os
    character(len=*),intent(in) :: inp_ph
    integer,intent(inout) :: var
    integer,intent(in),optional :: def
    integer :: i,j,ntokens,u1,il,stat
    character(len=1024) :: line
    character(len=100) :: tokens(10)
    character(len=20) :: inpFile

    narg=command_argument_count()

    if (narg /= 0) then
      call get_command_argument(1,inpFile)
    else
      inpFile='input.dat'
    end if

    open (newunit=u1,action='read',file=inpFile,status='old')
    il=1

    if (present(def)) then
      var=def
    endif
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        if (.not.present(def)) then
          print '(" Error: variable ",a," was not found.")',inp_ph
          stop ! should be changed
        endif
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" cmn_io_mod: Error reading line ", i0)',il
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
          if(trim(adjustl(tokens(j))) == inp_ph) then
            call value(tokens(j+os+1),var,ios)
            print '(" ",a," = ",i)',inp_ph,var
            exit ef
          end if
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

  end subroutine read_input_int


  !> Reads the input for character
  !! \param id the MPI id
  !! \param inp_ph the input phrase used for the variable
  !! \param os offset from the starting argument, first argument for var
  !! \param nvar the number of variables for var
  !! \param var the actual variables
  subroutine read_input_wp(inp_ph,os,var,def)

    use :: strg_mod
    use :: iso_fortran_env
    
    integer,intent(in) :: os
    character(len=*),intent(in) :: inp_ph
    real(wp),intent(inout) :: var
    real(wp),intent(in),optional :: def
    integer :: i,j,ntokens,u1,il,stat
    character(len=1024) :: line
    character(len=100) :: tokens(10)
    character(len=20) :: inpFile

    narg=command_argument_count()

    if (narg /= 0) then
      call get_command_argument(1,inpFile)
    else
      inpFile='input.dat'
    end if

    open (newunit=u1,action='read',file=inpFile,status='old')
    il=1

    if (present(def)) then
      var=def
    endif
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        if (.not.present(def)) then
          print '(" Error: variable ",a," was not found.")',inp_ph
          stop ! should be changed
        endif
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" cmn_io_mod: Error reading line ", i0)',il
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
          if(trim(adjustl(tokens(j))) == inp_ph) then
            call value(tokens(j+os+1),var,ios)
            print '(" ",a," = ",f20.5)',inp_ph,var
            exit ef
          end if
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

  end subroutine read_input_wp

end module cmn_io_mod