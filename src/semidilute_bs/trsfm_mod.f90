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
! MODULE: Transformation
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: 
!> (1) Interchanging the configurational state from Q to R
!! (2) Converting {R,rc} to position vector Rb
!! (3) Converting Q to segmental force Fseg
!--------------------------------------------------------------------

module trsfm_mod

  use :: prcn_mod

  implicit none

  ! Private module procedures:
  private :: init_trsfm_t  ,&
             applypbc_glob ,&
             applypbc_rec  ,&
             map           ,&
             remap         ,&
             unwrap_box
             
  !> A public type for configurational transformation
  type trsfm

    private

    !> The x-component of transformed position of the beads
    real(wp),pointer,public :: Rbtrx(:) => null()
    !> The y-component of transformed position of the beads
    real(wp),pointer,public :: Rbtry(:) => null()
    !> The x-component of transformed position of the center of mass
    real(wp),pointer,public :: rcmtrx(:) => null()
    !> The y-component of transformed position of the center of mass
    real(wp),pointer,public :: rcmtry(:) => null()

  contains

    procedure,pass(this) :: init => init_trsfm_t
    procedure,pass(this) :: applypbc => applypbc_glob
    procedure,pass(this) :: unwrap => unwrap_box
    final :: del_trsfm_t

  end type trsfm

  ! Private module variables:
  private :: ieps,eps_mx,L1_0,L2_0,L1p,L2p,sinth0,costh0,eps_p
  ! Protected module variables:
  protected :: theta_0,theta,bsx,bsy,invbsx,invbsy,sinth,costh,&
               tanb,L1,L2,reArng
  !> The number of periodic strains applied to the box
  integer,save :: ieps
  !> if the box should retain to its original state
  logical,save :: reArng
  ! Parameters used in shear flow:
  !{
  !> The maximum strain corresponding to maximum deflection angle
  real(wp),save :: eps_mx
  !> modulo(strain,max(strain))
  real(wp),save :: eps_m
  !> Deflection of the box in the top layer
  real(wp),save :: delrx_L
  !> Deflection of the box in the middle height of the box
  real(wp),save :: delrx_m
  !}
  ! Parameters used in elongational flow:
  !{
  !> The initial angle of the box with respect to x-axis in PEF
  real(wp),save :: theta_0
  !> The current angle of the box with respect to x-axis in PEF
  real(wp),save :: theta
  !> The periodic strain in planar elongational flow
  real(wp),save :: eps_p
  !> The dimension of the deformed box in x-direction
  real(wp),save :: bsx
  !> The dimension of the deformed box in y-direction
  real(wp),save :: bsy
  !> The inversed dimension of the deformed box in x-direction
  real(wp),save :: invbsx
  !> The inversed dimension of the deformed box in y-direction
  real(wp),save :: invbsy
  !> The initial box basis vector, L1,0
  real(wp),save :: L1_0(2)
  !> The initial box basis vector, L2,0
  real(wp),save :: L2_0(2)
  !> The current box basis vector, L1
  real(wp),save :: L1(2)
  !> The current box basis vector, L2
  real(wp),save :: L2(2)
  !> The current box rotated basis vector, L1'
  real(wp),save :: L1p(2)
  !> The current box rotated basis vector, L2'
  real(wp),save :: L2p(2)
  !> The sin(theta)
  real(wp),save :: sinth
  !> The cos(theta)
  real(wp),save :: costh
  !> The sin(theta0)
  real(wp),save :: sinth0
  !> The cos(theta0)
  real(wp),save :: costh0
  !> The tan(beta), where beta is the angle of L2' with y-axis
  real(wp),save :: tanb
  !}

contains

  !> Initializes the trsfm module
  subroutine init_trsfm()

    use :: flow_mod, only: FlowType

    real(wp),parameter :: PI=3.1415926535897958648_wp!4*atan(1.0_wp)
    ! The maximum deflection angle in planar shear flow
    real(wp) :: theta_mx
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:init_trsfm"
#endif
    select case (FlowType)
      case ('PSF')
        theta_mx=PI/4
        eps_mx=tan(theta_mx)
      case ('PEF')
        theta_0=atan(-0.5_wp*(1-sqrt(5._wp)))
        sinth0=sin(theta_0)
        costh0=cos(theta_0)
        eps_p=log(0.5_wp*(3+sqrt(5._wp)))
    end select

  end subroutine init_trsfm

  !> Initializes trsfm module variables when itime=1
  !! \param bs the initial dimension of the box
  subroutine init_trsfm_tm(bs)

    use :: flow_mod, only: FlowType

    real(wp),intent(in) :: bs(3)
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:init_trsfm_tm"
#endif
    ieps=1
    reArng=.false.
    select case (FlowType)
      case ('PSF')
        eps_m=0._wp
        delrx_L=0._wp
        delrx_m=0._wp
      case ('PEF')
        theta=theta_0
        sinth=sinth0
        costh=costh0
        tanb=0._wp
        ! Note that the initial box should be square in x-y plane
        L1_0=bs(1)*[ costh0,sinth0]
        L2_0=bs(1)*[-sinth0,costh0]
    end select

  end subroutine init_trsfm_tm

  !> Updates rearrangement logical variable at each time step
  !! \param eps applied strain
  !! \param bs the dimension of the box
  subroutine update_arng(eps)

    use :: flow_mod, only: FlowType

    real(wp) :: M(2,2),eps,eps_r
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:update_arng"
#endif
    select case (FlowType)
      case ('PSF')
        eps_m=mod(eps,eps_mx)
        eps_r=eps/eps_mx
      case ('PEF')
        eps_m=mod(eps,eps_p)
        eps_r=eps/eps_p
    end select

    reArng=.false.
    if (FlowType == 'PSF' .or. FlowType == 'PEF') then
      if (floor(eps_r) == ieps) then
        ieps=ieps+1
        reArng=.true.
      else
        reArng=.false.
      end if
    endif

  end subroutine update_arng

  !> Updates trsfm module variables at each time step
  !! \param eps applied strain
  !! \param bs the dimension of the box
  subroutine update_trsfm(bs)

    use :: flow_mod, only: FlowType

    real(wp),intent(in) :: bs(3)
    real(wp) :: M(2,2)
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:update_trsfm"
#endif
    select case (FlowType)
      case ('PSF')
        delrx_L=eps_m*bs(2)
        delrx_m=delrx_L/2
      case ('PEF')
        L1=[L1_0(1)*exp(eps_m),L1_0(2)*exp(-eps_m)]
        L2=[L2_0(1)*exp(eps_m),L2_0(2)*exp(-eps_m)]
        theta=atan(L1(2)/L1(1))
        sinth=sin(theta)
        costh=cos(theta)
        M(1,:)=[ costh,sinth]
        M(2,:)=[-sinth,costh]
        L1p=matmul(M,L1)
        L2p=matmul(M,L2)
        bsx=sqrt(L1(1)**2+L1(2)**2)
        bsy=L2p(2)
        invbsx=1/bsx
        invbsy=1/bsy
        tanb=L2p(1)/L2p(2)
    end select

  end subroutine update_trsfm

  !> Constructor for trsfm type
  !! \param ntotbead total number of beads inside the box
  !! \param nchain the number of chain inside the box
  subroutine init_trsfm_t(this,Rbtr,rcmtr)

!    use :: inp_smdlt, only: ntotbead,nchain
    use :: flow_mod, only: FlowType

    class(trsfm),intent(inout) :: this
    ! integer,intent(in) :: ntotbead,nchain
    real(wp),intent(in),target,contiguous :: Rbtr(:,:)
    real(wp),intent(in),target,contiguous :: rcmtr(:,:)
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:init_trsfm_t"
#endif
    select case (FlowType)
      case ('Equil')
        ! Rbtr and rcmtr in this case are zero-sized
        this%Rbtrx(1:size(Rbtr)) => Rbtr
        this%rcmtrx(1:size(rcmtr)) => rcmtr
      case ('PSF')
        this%Rbtrx => Rbtr(:,1)
        this%rcmtrx => rcmtr(:,1)
      case ('PEF')
        this%Rbtrx => Rbtr(:,1)
        this%Rbtry => Rbtr(:,2)
        this%rcmtrx => rcmtr(:,1)
        this%rcmtry => rcmtr(:,2)
    end select
        
  end subroutine init_trsfm_t

  !> Applying periodic boundary condition, global call
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param str the modulo(strain,max(strain))
  !! \param b_img the image of the beads inside the primary box
  !! \param cm_img the image of the center of mass inside the primary box
  subroutine applypbc_glob(this,bs,invbs,Rbx,Rby,Rbz,rcm,b_img,cm_img,itime)

    use :: flow_mod, only: FlowType

    class(trsfm),intent(inout) :: this
    integer,intent(in) :: itime
    real(wp),intent(in) :: bs(3),invbs(3)
    real(wp),intent(inout) :: Rbx(:)
    real(wp),intent(inout) :: Rby(:)
    real(wp),intent(inout) :: Rbz(:)
    real(wp),intent(inout) :: rcm(:,:)
    integer,intent(inout) :: b_img(:,:)
    integer,intent(inout) :: cm_img(:,:)
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:applypbc_glob"
#endif
    if (FlowType /= 'Equil') call map(this,Rbx,Rby,rcm,itime)
    call applypbc_rec(this,bs,invbs,Rbx,Rby,Rbz,rcm,b_img,cm_img,itime)
    if (FlowType /= 'Equil') call remap(this,bs,invbs,Rbx,Rby,rcm,b_img,cm_img,itime)

  end subroutine applypbc_glob

  !> Applying periodic boundary condition on a rectangular box
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param b_img the image of the beads inside the primary box
  !! \param cm_img the image of the center of mass inside the primary box
  subroutine applypbc_rec(this,bs,invbs,Rbx,Rby,Rbz,rcm,b_img,cm_img,itime)

    use :: flow_mod, only: FlowType
    use :: arry_mod, only: print_vector,print_matrix

    class(trsfm),intent(inout) :: this
    real(wp),intent(in) :: bs(3),invbs(3)
    integer,intent(in) :: itime
    real(wp),intent(inout) :: Rbx(:)
    real(wp),intent(inout) :: Rby(:)
    real(wp),intent(inout) :: Rbz(:)
    real(wp),intent(inout) :: rcm(:,:)
    integer,intent(inout) :: b_img(:,:)
    integer,intent(inout) :: cm_img(:,:)
    integer :: igb,ich
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:applypbc_rec"
#endif
    select case (FlowType)
      case ('Equil')
!$omp parallel default(private) shared(Rbx,Rby,Rbz,b_img,bs,invbs,cm_img,rcm)
!$omp do simd
        do igb=1, size(Rbx,1)
          ! calculating the image of the beads
          b_img(igb,1)=b_img(igb,1)-nint(Rbx(igb)*invbs(1)-0.5_wp)
          b_img(igb,2)=b_img(igb,2)-nint(Rby(igb)*invbs(2)-0.5_wp)
          b_img(igb,3)=b_img(igb,3)-nint(Rbz(igb)*invbs(3)-0.5_wp)
          Rbx(igb)=Rbx(igb)-nint(Rbx(igb)*invbs(1)-0.5_wp)*bs(1)
          Rby(igb)=Rby(igb)-nint(Rby(igb)*invbs(2)-0.5_wp)*bs(2)
          Rbz(igb)=Rbz(igb)-nint(Rbz(igb)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd

! call print_matrix(rcm,'rcmtrsfm')
!$omp do simd
        do ich=1, size(rcm,1)
          cm_img(ich,1)=-nint(rcm(ich,1)*invbs(1)-0.5_wp)
          cm_img(ich,2)=-nint(rcm(ich,2)*invbs(2)-0.5_wp)
          cm_img(ich,3)=-nint(rcm(ich,3)*invbs(3)-0.5_wp)
          rcm(ich,1)=rcm(ich,1)-nint(rcm(ich,1)*invbs(1)-0.5_wp)*bs(1)
          rcm(ich,2)=rcm(ich,2)-nint(rcm(ich,2)*invbs(2)-0.5_wp)*bs(2)
          rcm(ich,3)=rcm(ich,3)-nint(rcm(ich,3)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd
!$omp end parallel
! call print_matrix(cm_img,'cm_imgtrsfm')
! call print_matrix(rcm,'rcmtrsfmafter')
      case ('PSF')
!$omp parallel default(private) shared(this,Rby,Rbz,b_img,bs,invbs,cm_img,rcm)
!$omp do simd
        do igb=1, size(Rby,1)
          ! calculating the image of the beads
          b_img(igb,1)=b_img(igb,1)-nint(this%Rbtrx(igb)*invbs(1)-0.5_wp)
          b_img(igb,2)=b_img(igb,2)-nint(Rby(igb)*invbs(2)-0.5_wp)
          b_img(igb,3)=b_img(igb,3)-nint(Rbz(igb)*invbs(3)-0.5_wp)
          this%Rbtrx(igb)=this%Rbtrx(igb)-nint(this%Rbtrx(igb)*invbs(1)-0.5_wp)*bs(1)
          Rby(igb)=Rby(igb)-nint(Rby(igb)*invbs(2)-0.5_wp)*bs(2)
          Rbz(igb)=Rbz(igb)-nint(Rbz(igb)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          cm_img(ich,1)=-nint(this%rcmtrx(ich)*invbs(1)-0.5_wp)
          cm_img(ich,2)=-nint(rcm(ich,2)*invbs(2)-0.5_wp)
          cm_img(ich,3)=-nint(rcm(ich,3)*invbs(3)-0.5_wp)
          this%rcmtrx(ich)=this%rcmtrx(ich)-nint(this%rcmtrx(ich)*invbs(1)-0.5_wp)*bs(1)
          rcm(ich,2)=rcm(ich,2)-nint(rcm(ich,2)*invbs(2)-0.5_wp)*bs(2)
          rcm(ich,3)=rcm(ich,3)-nint(rcm(ich,3)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd
!$omp end parallel
      case ('PEF')
!$omp parallel default(private) shared(this,Rbz,b_img,bs,invbs,bsx,invbsx,bsy,invbsy,reArng,cm_img,rcm)
!$omp do simd
        do igb=1, size(Rbz,1)
          b_img(igb,1)=b_img(igb,1)-nint(this%Rbtrx(igb)*invbsx-0.5_wp)
          b_img(igb,2)=b_img(igb,2)-nint(this%Rbtry(igb)*invbsy-0.5_wp)
          b_img(igb,3)=b_img(igb,3)-nint(Rbz(igb)*invbs(3)-0.5_wp)
          this%Rbtrx(igb)=this%Rbtrx(igb)-nint(this%Rbtrx(igb)*invbsx-0.5_wp)*bsx
          this%Rbtry(igb)=this%Rbtry(igb)-nint(this%Rbtry(igb)*invbsy-0.5_wp)*bsy
          Rbz(igb)=Rbz(igb)-nint(Rbz(igb)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          cm_img(ich,1)=-nint(this%rcmtrx(ich)*invbsx-0.5_wp)
          cm_img(ich,2)=-nint(this%rcmtry(ich)*invbsy-0.5_wp)
          cm_img(ich,3)=-nint(rcm(ich,3)*invbs(3)-0.5_wp)
          this%rcmtrx(ich)=this%rcmtrx(ich)-nint(this%rcmtrx(ich)*invbsx-0.5_wp)*bsx
          this%rcmtry(ich)=this%rcmtry(ich)-nint(this%rcmtry(ich)*invbsy-0.5_wp)*bsy
          rcm(ich,3)=rcm(ich,3)-nint(rcm(ich,3)*invbs(3)-0.5_wp)*bs(3)
        end do
!!$omp end do simd
!$omp end parallel

    end select

  end subroutine applypbc_rec

  !> Mapping the particles to a rectangular box
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbtr the beads position transdormed to a rectangular box
  subroutine map(this,Rbx,Rby,rcm,itime)
  
    use :: flow_mod, only: FlowType

    class(trsfm),intent(inout) :: this
    integer,intent(in) :: itime
    real(wp),intent(in) :: Rbx(:)
    real(wp),intent(in) :: Rby(:)
    real(wp),intent(in) :: rcm(:,:)
    integer :: igb,ich
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:map"
#endif
    select case (FlowType)
      case ('PSF')
!$omp parallel default(private) shared(this,Rbx,Rby,eps_m,delrx_m,rcm)
!$omp do simd
        do igb=1, size(Rbx,1)
          this%Rbtrx(igb)=Rbx(igb)-eps_m*Rby(igb)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          this%rcmtrx(ich)=rcm(ich,1)-eps_m*rcm(ich,2)
        end do
!!$omp end do simd
!$omp end parallel
      case ('PEF') 
!$omp parallel default(private) shared(this,Rbx,Rby,sinth,costh,tanb,rcm)
!$omp do simd
        do igb=1, size(Rbx,1)
          this%Rbtry(igb)=-sinth*Rbx(igb)+costh*Rby(igb)
          this%Rbtrx(igb)= costh*Rbx(igb)+sinth*Rby(igb)-tanb*this%Rbtry(igb)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          this%rcmtry(ich)=-sinth*rcm(ich,1)+costh*rcm(ich,2)
          this%rcmtrx(ich)= costh*rcm(ich,1)+sinth*rcm(ich,2)-tanb*this%rcmtry(ich)
        end do
!!$omp end do simd
!$omp end parallel
    end select

  end subroutine map

  !> Mapping back the particles from the rectangular box to deformed box
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param b_img the image of the beads inside the primary box
  !! \param cm_img the image of the center of mass inside the primary box
  subroutine remap(this,bs,invbs,Rbx,Rby,rcm,b_img,cm_img,itime)

    use :: flow_mod, only: FlowType

    class(trsfm),intent(inout) :: this
    real(wp),intent(in) :: bs(3),invbs(3)
    real(wp),intent(inout),contiguous :: Rbx(:)
    real(wp),intent(inout),contiguous :: Rby(:)
    real(wp),intent(inout) :: rcm(:,:)
    integer,intent(inout) :: b_img(:,:)
    integer,intent(inout) :: cm_img(:,:)
    integer :: igb,ich
    real(wp) :: Rbxtmp,Rbytmp,rcmxtmp,rcmytmp
    integer :: itime
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:remap"
#endif
    select case (FlowType)
      case ('PSF')
!$omp parallel default(private) shared(this,Rbx,Rby,eps_m,delrx_m,invbs,bs,reArng,rcm)
!$omp do simd
        do igb=1, size(Rbx,1)
          Rbx(igb)=this%Rbtrx(igb)+eps_m*Rby(igb)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          rcm(ich,1)=this%rcmtrx(ich)+eps_m*rcm(ich,2)
        end do
!!$omp end do simd
!$omp end parallel
      case ('PEF')
!$omp parallel default(private) &
!$omp shared(this,Rbx,Rby,tanb,sinth,costh,sinth0,costh0,invbs,bs,reArng,b_img,rcm)
!$omp do simd
        do igb=1, size(Rbx,1)
          Rbx(igb)=this%Rbtrx(igb)+tanb*this%Rbtry(igb)
          Rby(igb)=sinth*Rbx(igb)+costh*this%Rbtry(igb)
          Rbx(igb)=costh*Rbx(igb)-sinth*this%Rbtry(igb)
        end do
!!$omp end do simd
!$omp do simd
        do ich=1, size(rcm,1)
          rcm(ich,1)=this%rcmtrx(ich)+tanb*this%rcmtry(ich)
          rcm(ich,2)=sinth*rcm(ich,1)+costh*this%rcmtry(ich)
          rcm(ich,1)=costh*rcm(ich,1)-sinth*this%rcmtry(ich)
        end do
!!$omp end do simd
!$omp end parallel
    end select

  end subroutine remap

  subroutine unwrap_box(this,bs,Rbx,Rby,b_img,itime)

    use :: flow_mod, only: FlowType

    class(trsfm),intent(inout) :: this
    real(wp),intent(in) :: bs(3)
    integer,intent(in) :: itime
    real(wp),intent(inout),contiguous :: Rbx(:)
    real(wp),intent(inout),contiguous :: Rby(:)
    integer,intent(inout) :: b_img(:,:)
    integer :: igb
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:unwrap_box"
#endif
    select case (FlowType)
      case ('PSF')
!$omp parallel default(private) &
!$omp shared(this,Rbx,eps_m,b_img,bs)
!$omp do simd
        do igb=1, size(Rbx,1)
          ! map
          this%Rbtrx(igb)=Rbx(igb)-eps_m*Rby(igb)
          ! rectangualr unwrap
          this%Rbtrx(igb)=this%Rbtrx(igb)-b_img(igb,1)*bs(1)
          b_img(igb,1)=0
          ! remap
          Rbx(igb)=this%Rbtrx(igb)+eps_m*Rby(igb)
        end do
!!$omp end do simd
!$omp end parallel
      case ('PEF')
!$omp parallel default(private) &
!$omp shared(this,Rbx,Rby,tanb,sinth,costh,b_img,bsx,bsy)
!$omp do simd
        do igb=1, size(Rbx,1)
          ! map
          this%Rbtry(igb)=-sinth*Rbx(igb)+costh*Rby(igb)
          this%Rbtrx(igb)= costh*Rbx(igb)+sinth*Rby(igb)-tanb*this%Rbtry(igb)
          ! rectangualr unwrap
          this%Rbtrx(igb)=this%Rbtrx(igb)-b_img(igb,1)*bsx
          this%Rbtry(igb)=this%Rbtry(igb)-b_img(igb,2)*bsy
          b_img(igb,1:2)=0
          ! remap
          Rbx(igb)=this%Rbtrx(igb)+tanb*this%Rbtry(igb)
          Rby(igb)=sinth*Rbx(igb)+costh*this%Rbtry(igb)
          Rbx(igb)=costh*Rbx(igb)-sinth*this%Rbtry(igb)
        end do
!!$omp end do simd
!$omp end parallel
    end select

  end subroutine unwrap_box

  !> Destructor for trsfm type
  subroutine del_trsfm_t(this)

    use :: flow_mod, only: FlowType

    type(trsfm) :: this
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:del_trsfm_t"
#endif
    select case (FlowType)
      case ('PSF')
        nullify(this%Rbtrx,this%rcmtrx)
      case ('PEF')
        nullify(this%Rbtrx,this%Rbtry)
        nullify(this%rcmtrx,this%rcmtry)
    end select

  end subroutine del_trsfm_t

  !> Rotates a vector around z-axis
  !! \param V the vector in three dimension
  !! \param th the angle of rotation
  subroutine zrotate(V,th)

    real(wp),intent(inout) :: V(:)
    real(wp),intent(in) :: th
    integer :: i,os
    real(wp) :: Vx,Vy
#ifdef Debuge_sequence
	write(*,*) "module:trsfm_mod:zrotate"
#endif
!$omp parallel default(private) shared(V,th)
!$omp do schedule(auto)
    do i=1, size(V)/3
      os=(i-1)*3
      Vx=V(os+1)
      Vy=V(os+2)
      V(os+1)= cos(th)*Vx+sin(th)*Vy
      V(os+2)=-sin(th)*Vx+cos(th)*Vy
    end do
!!$omp end do
!$omp end parallel

  end subroutine zrotate

end module trsfm_mod
