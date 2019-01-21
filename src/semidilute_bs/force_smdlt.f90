!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
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
! MODULE: force
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Apr 2014
!
! DESCRIPTION: 
!> Calculating the conservative forces
!--------------------------------------------------------------------
module force_smdlt
  
  use :: prcn_mod

  implicit none

!  private :: 

  type,abstract :: force

    private

  contains
    
    procedure(updateforce),deferred :: update

  end type force

  abstract interface
 
    subroutine updateforce(this,Rbx,Rby,Rbz,bs,invbs,itime,nchain,nseg,nbead,&
                          ntotseg,ntotsegx3,ntotbead,ntotbeadx3,Qt)
      import :: force
      import :: wp
      implicit none
      class(force),intent(inout) :: this
      real(wp),intent(in) :: Rbx(:)
      real(wp),intent(in) :: Rby(:)
      real(wp),intent(in) :: Rbz(:)
      real(wp),intent(in) :: bs(3),invbs(3)
!      real(wp),intent(inout) :: F(:)
      integer,intent(in) :: itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbead,ntotbeadx3
      real(wp),intent(in) :: Qt(:)
      
    end subroutine updateforce

  end interface

  ! Protected module variables:
!  protected :: Fphi,Fx,Fy,Fz,rFphi

  !> Total force acting on the particles
  real(wp),allocatable,target,save :: Fphi(:)
  !> x coordinate of total force
  real(wp),pointer,save :: Fx(:)
  !> y coordinate of total force
  real(wp),pointer,save :: Fy(:)
  !> z coordinate of total force
  real(wp),pointer,save :: Fz(:)
  !> Interparticle distance times total force
  real(wp),save :: rFphi(4)

contains

  !> Initializes force_mod module variables
  subroutine init_force(ntotbeadx3)

    integer,intent(in) :: ntotbeadx3

    allocate(Fphi(ntotbeadx3))
    Fx => Fphi(1:ntotbeadx3-2:3)
    Fy => Fphi(2:ntotbeadx3-1:3)
    Fz => Fphi(3:ntotbeadx3:3)

  end subroutine init_force

  !> Deletion of the force_mod module variables
  subroutine del_force()

    use :: verlet_mod, only: del_verlet

    nullify(Fx,Fy,Fz)
    deallocate(Fphi)

    call del_verlet()

  end subroutine del_force

end module force_smdlt