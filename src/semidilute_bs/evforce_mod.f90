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
! MODULE: EV force
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Apr 2014
!
! DESCRIPTION: 
!> Calculating the excluded volume force
!--------------------------------------------------------------------
module evforce_mod

  use :: prcn_mod
  use :: force_smdlt, only: force
  use :: evverlet_mod, only: evverlet

  implicit none

  ! Private module procedures:
  private :: init_evforce_t ,&
             update_vlt_lst ,&
             update_force

  !> A public type for EV force calcualtion
  type, extends(force) :: evforce

    private
    
    !> An object for construcing neighbor list for ev force
    type(evverlet) :: evvlt
    !> The position vector in the previous list update iteration
    real(wp),allocatable :: Rb0(:)
    !> The neighbor list
    integer,allocatable :: nlst(:,:)
    !> Gaussian ev force parameters
    !> @{
    real(wp) :: fctr,efctr
    !> @}

  contains
     
    procedure,pass(this) :: init => init_evforce_t
    procedure,pass(this) :: update_vlt => update_vlt_lst
    procedure,pass(this) :: update => update_force
    final :: del_evforce

  end type evforce
  
  ! Private module variables:
  ! private :: 
  ! Protected module variables:
  protected :: EVForceLaw,rc_F,zstar,dstar
  
  !> The type of ev force
  character(len=10),save :: EVForceLaw
  !> The cutoff radius for ev force
  real(wp),save :: rc_F
  !> The search radius for verlet list
  real(wp),save :: rs_F
  !> The ev strength parameter z*
  real(wp),save :: zstar
  !> The ev broadness parameter d*
  real(wp),save :: dstar

contains

  !> Initialization of the evforce module
  !! \param id The rank of the process
  subroutine init_evforce(id)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id
    integer :: il,j,ntokens,u1,stat,ios
    character(len=1024) :: line
    character(len=100) :: tokens(10)
    character(len=10) :: dstarCalc
    real(wp) :: s_F
    
    ! default settings:
    EVForceLaw='NoEV' 
    dstar=1._wp;dstarCalc='Kumar'
    rc_F=7._wp;s_F=0.2_wp
    open (newunit=u1,action='read',file='input.dat',status='old')
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
          select case (trim(adjustl(tokens(j))))
            case ('EVForceLaw')
              EVForceLaw=trim(adjustl(tokens(j+1)))
            case ('zstar')
              call value(tokens(j+1),zstar,ios)
            case ('dstar')
              call value(tokens(j+1),dstar,ios)
              dstarCalc=trim(adjustl(tokens(j+2)))
            case ('rc_F')
              call value(tokens(j+1),rc_F,ios)
              call value(tokens(j+2),s_F,ios)
          end select 
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

    rs_F=rc_F*(1+s_F)

    ! The parameters used for EV potentials
    select case (EVForceLaw)
      case ('Gauss')
        ! zstar and dstar (EV potential parameters) setting
        select case (dstarCalc)
          case ('Kumar')
            if (zstar == 0) then
              if (id == 0) print '(" zstar should be nonzero if Gauss is used.")'
              stop
            end if
            dstar=dstar*zstar**(1._wp/5)
          case default
            ! No change to dstar is necessary.
        end select

        if (id == 0) then
          print *
          print '(" Polymers are in good solvent with zstar=",f12.3,f12.3)',zstar
          print '(" Cutoff radius for EV force: ",f12.3)',rc_F
          print '(" Search radius for EV force: ",f12.3)',rs_F
        end if

      case ('LJ')
      case ('NoEV')
        if (id == 0) then
          print '(" Polymers are in theta solvent.")'
        end if
    end select
  
  end subroutine init_evforce

  !> Constructor for evforce type
  subroutine init_evforce_t(this,id,bs,ntotbead,ntotbeadx3)

    class(evforce),intent(inout) :: this
    integer,intent(in) :: id,ntotbead,ntotbeadx3
    real(wp),intent(in) :: bs(3)

    ! The parameters used for EV potentials
    select case (EVForceLaw)
      case ('Gauss')
        this%fctr=zstar/dstar**5
        this%efctr=1/(2*dstar**2)
        allocate(this%Rb0(ntotbeadx3))
      case ('LJ')
      case ('NoEV')
    end select

    call this%evvlt%init(rs_F,bs,ntotbead)
    if (id == 0) call this%evvlt%get_ncps()
  
  end subroutine init_evforce_t

  !> Updates the verlet list
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
!  subroutine update_vlt_lst(this,Rbx,Rby,Rbz,bs,invbs,itime,itrst,ntotbead,ntotbeadx3)
  subroutine update_vlt_lst(this,Rbx,Rby,Rbz,Rb,bs,invbs,itime,itrst,ntotbead,ntotbeadx3)

    use :: flow_mod, only: FlowType

    class(evforce),intent(inout) :: this
    real(wp),intent(in) :: Rbx(:),Rby(:),Rbz(:)
    real(wp),intent(in) :: Rb(:)
    real(wp),intent(in) :: bs(3),invbs(3)
    integer,intent(in) :: itime,itrst,ntotbead,ntotbeadx3
    real(wp) :: dispmax
    logical :: update

    if (itime == itrst+1) then
      update=.true.
    else
      ! Calculate maximum displacement since last update:
      dispmax=maxval(abs(Rb-this%Rb0))
      ! A conservative testing of the list skin crossing:
      dispmax=2*sqrt(3*dispmax*dispmax)
      update=dispmax > (rs_F-rc_F)
    end if
    if (update) then
      ! Save positions for next evaluations:
      this%Rb0=Rb
      if ((FlowType == 'PEF').and.(itime /= itrst+1)) then
        call this%evvlt%init_cll(rs_F,bs,ntotbead)
      end if
      call this%evvlt%cnstr_cll(Rbx,Rby,Rbz,itime,ntotbead,ntotbeadx3)
      call this%evvlt%cnstr_nab(Rbx,Rby,Rbz,rs_F,bs,invbs,this%nlst,itime,ntotbead,ntotbeadx3)
    end if

  end subroutine update_vlt_lst

  !> Updates the force by adding ev force contribution
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  subroutine update_force(this,Rbx,Rby,Rbz,bs,invbs,itime,nchain,nseg,nbead,&
                          ntotseg,ntotsegx3,ntotbead,ntotbeadx3,Qt)

    use :: arry_mod, only: print_vector
    use :: trsfm_mod, only: eps_m,tanb,sinth,costh
    use :: flow_mod, only: FlowType
    use :: force_smdlt, only: Fx,Fy,Fz,rFphi,Fphi

    class(evforce),intent(inout) :: this    
    real(wp),intent(in) :: Rbx(:)
    real(wp),intent(in) :: Rby(:)
    real(wp),intent(in) :: Rbz(:)
    real(wp),intent(in) :: bs(3),invbs(3) 
!    real(wp),intent(inout) :: F(:)
    integer,intent(in) :: itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbead,ntotbeadx3
    real(wp),intent(in) :: Qt(:)
    integer :: iint,i,j
    real(wp) :: rijx,rijy,rijz,rijsq,rijytmp,Fevij(3)

!    this%Fev=0._wp
    do iint=1, size(this%nlst,1)
      i=this%nlst(iint,1)
      j=this%nlst(iint,2)
      rijx=Rbx(i)-Rbx(j)
      rijy=Rby(i)-Rby(j)
      rijz=Rbz(i)-Rbz(j)
      rijx=rijx-nint(rijx*invbs(1))*bs(1)
      rijy=rijy-nint(rijy*invbs(2))*bs(2)
      rijz=rijz-nint(rijz*invbs(3))*bs(3)
      select case (FlowType)
        case ('PSF')
          rijx=rijx+eps_m*rijy
        case ('PEF')
          rijytmp=rijy
          rijx=rijx+tanb*rijytmp
          rijy=sinth*rijx+costh*rijytmp
          rijx=costh*rijx-sinth*rijytmp
      end select
      rijsq=rijx*rijx+rijy*rijy+rijz*rijz
      if (rijsq <= rc_F**2) then
        Fevij=this%fctr*[rijx,rijy,rijz]*exp(-rijsq*this%efctr)
        Fx(i)=Fx(i)+Fevij(1)
        Fy(i)=Fy(i)+Fevij(2)
        Fz(i)=Fz(i)+Fevij(3)
        Fx(j)=Fx(j)-Fevij(1)
        Fy(j)=Fy(j)-Fevij(2)
        Fz(j)=Fz(j)-Fevij(3)
!        this%Fevx(i)=this%Fevx(i)+this%fctr*rijx*exp(-rijsq*this%efctr)
!        this%Fevy(i)=this%Fevy(i)+this%fctr*rijy*exp(-rijsq*this%efctr)
!        this%Fevz(i)=this%Fevz(i)+this%fctr*rijz*exp(-rijsq*this%efctr)
!        this%Fevx(j)=this%Fevx(j)-this%fctr*rijx*exp(-rijsq*this%efctr)
!        this%Fevy(j)=this%Fevy(j)-this%fctr*rijy*exp(-rijsq*this%efctr)
!        this%Fevz(j)=this%Fevz(j)-this%fctr*rijz*exp(-rijsq*this%efctr)
        rFphi(1)=rFphi(1)+rijx*Fevij(1)
        rFphi(2)=rFphi(2)+rijx*Fevij(2)
        rFphi(3)=rFphi(3)+rijy*Fevij(2)
        rFphi(4)=rFphi(4)+rijz*Fevij(3)
      
      end if
    end do

!    F=F+this%Fev

  end subroutine update_force

  !> Destructor for evforce type
  subroutine del_evforce(this)

    type(evforce),intent(inout) :: this
   
  end subroutine del_evforce

end module evforce_mod
