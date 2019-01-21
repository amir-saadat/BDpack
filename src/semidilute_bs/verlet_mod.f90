!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2018:                                            |
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
! MODULE: verlet
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: construction of verlet list
!> 
!
!--------------------------------------------------------------------

module verlet_mod
  
  use :: prcn_mod

  implicit none

  !> A public type for constructing verlet list
  type verlet
   
    private

  contains

  end type verlet

  ! Protected module variables:
  protected :: nnc,shifts

  !> Number of neighbering cells
  integer,save :: nnc
  !> The neighboring cells offset
  integer,allocatable,save :: shifts(:,:)
  !> The coordinates for neighboring cells
  integer,allocatable :: j_clx(:),j_cly(:),j_clz(:),j_cll(:)

contains 

  !> Initializes the verlet module
  !! \param id The rank of the process
  subroutine init_verlet(id)

    use :: arry_mod, only: print_vector
    use :: flow_mod, only: FlowType
    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id


    select case (FlowType)

      case ('Equil','PSF')
        ! nnc=13
        nnc=27
        allocate(shifts(nnc,3))
        shifts(1,:) =[ 0, 0,-1]
        shifts(2,:) =[ 1, 0,-1]
        shifts(3,:) =[ 1, 0, 0]
        shifts(4,:) =[ 1, 0, 1]
        shifts(5,:) =[-1, 1,-1]
        shifts(6,:) =[ 0, 1,-1]
        shifts(7,:) =[ 1, 1,-1]
        shifts(8,:) =[-1, 1, 0]
        shifts(9,:) =[ 0, 1, 0]
        shifts(10,:)=[ 1, 1, 0]
        shifts(11,:)=[-1, 1, 1]
        shifts(12,:)=[ 0, 1, 1]
        shifts(13,:)=[ 1, 1, 1]

        shifts(14,:)=[ 0, 0, 0]

        shifts(15:27,:) = - shifts(1:13,:)

      case ('PEF')
        ! nnc=31
        nnc=63
        allocate(shifts(nnc,3))
        shifts(1,:) =[ 0, 0,-1]
        shifts(2,:) =[ 1, 0,-1]
        shifts(3,:) =[ 2, 0,-1]
        shifts(4,:) =[ 3, 0,-1]
        shifts(5,:) =[ 1, 0, 0]
        shifts(6,:) =[ 2, 0, 0]
        shifts(7,:) =[ 3, 0, 0]
        shifts(8,:) =[ 1, 0, 1]
        shifts(9,:) =[ 2, 0, 1]
        shifts(10,:)=[ 3, 0, 1]
        shifts(11,:)=[-3, 1,-1]
        shifts(12,:)=[-2, 1,-1]
        shifts(13,:)=[-1, 1,-1]
        shifts(14,:)=[ 0, 1,-1]
        shifts(15,:)=[ 1, 1,-1]
        shifts(16,:)=[ 2, 1,-1]
        shifts(17,:)=[ 3, 1,-1]
        shifts(18,:)=[-3, 1, 0]
        shifts(19,:)=[-2, 1, 0]
        shifts(20,:)=[-1, 1, 0]
        shifts(21,:)=[ 0, 1, 0]
        shifts(22,:)=[ 1, 1, 0]
        shifts(23,:)=[ 2, 1, 0]
        shifts(24,:)=[ 3, 1, 0]
        shifts(25,:)=[-3, 1, 1]
        shifts(26,:)=[-2, 1, 1]
        shifts(27,:)=[-1, 1, 1]
        shifts(28,:)=[ 0, 1, 1]
        shifts(29,:)=[ 1, 1, 1]
        shifts(30,:)=[ 2, 1, 1]
        shifts(31,:)=[ 3, 1, 1]

        shifts(32,:)=[ 0, 0, 0]

        shifts(33:63,:) = - shifts(1:31,:)
!        this%ncps(1:2)=bs(1:2)/(sqrt(10._wp)*rc_F)
!        this%ncps(3)=bs(3)/rc_F
    end select

    allocate(j_clx(nnc))
    allocate(j_cly(nnc))
    allocate(j_clz(nnc))
    allocate(j_cll(nnc))

  end subroutine init_verlet

  subroutine del_verlet()

    deallocate(shifts)
    deallocate(j_clx,j_cly,j_clz,j_cll)

  end subroutine del_verlet

end module verlet_mod
