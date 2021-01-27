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
! MODULE: Spring force
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Apr 2014
!
! DESCRIPTION: 
!> Calculating the net spring forces on the beads
!--------------------------------------------------------------------
module sprforce_mod

  use :: prcn_mod
  use :: force_smdlt, only: force

  implicit none

  ! Private module procedures:
  private :: init_sprforce_t ,&
             update_force

  !> A public type for EV force calcualtion
  type, extends(force) :: sprforce

    private
    !> The spring force
    real(wp),allocatable :: Fs(:)

  contains
     
    procedure,pass(this) :: init => init_sprforce_t
    procedure,pass(this) :: update => update_force
    final :: del_sprforce

  end type sprforce
  
  ! Private module variables:
!  private ::
  ! Protected module variables:
  protected :: ForceLaw,b,qmx,WLC_v,WLC_A,WLC_B
  
  !> The type of ev force
  character(len=10),save :: ForceLaw
  integer,save :: ForceLaw_i
  !> Different types of force law
  integer,parameter :: Hookean=1     ,&
                       FENE   =2     ,&
                       ILCCP  =3     ,&
                       WLC_MS =4     ,&
                       WLC_UD =5
  !> The maximum dimensionless squared length of a spring
  real(wp),save :: b
  !> The maximum dimensionless length of a spring
  real(wp),save :: qmx
  !> The number of Kuhn steps per spring for WLC model
  real(wp),save :: WLC_v
  !> Force parameter for WLC-UD model
  real(wp),save :: WLC_A
  !> Force parameter for WLC-UD model
  real(wp),save :: WLC_B

contains

  !> Initialization of the sprforce module
  !! \param id The rank of the process
  subroutine init_sprforce(id)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id
    integer :: il,j,ntokens,u1,stat,ios
    character(len=1024) :: line
    character(len=100) :: tokens(10)

    ! default values:
    ForceLaw='Hookean'

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
            case ('SPR-Force')
              ForceLaw=trim(adjustl(tokens(j+1)))
            case ('b')
              call value(tokens(j+1),b,ios)
            case ('N_Ks')
              call value(tokens(j+1),WLC_v,ios)
          end select
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

    qmx=sqrt(b)

    select case (ForceLaw)
    case('Hookean')
      ForceLaw_i=Hookean
    case('FENE')
      ForceLaw_i=FENE
    case('ILCCP')
      ForceLaw_i=ILCCP
    case('WLC_MS')
      ForceLaw_i=WLC_MS
    case('WLC_UD')
      ForceLaw_i=WLC_UD
    case default
      print('(" Force law not properly chosen.")')
    end select
  
    select case (ForceLaw)
      case ('WLC_UD')
        WLC_A=3._wp/32-3/(8*WLC_v)-3/(2*WLC_v**2)
        WLC_B=(13._wp/32+0.4086_wp/WLC_v-14.79_wp/(4*WLC_v**2))/ &
                    (1-4.225_wp/(2*WLC_v)+4.87_wp/(4*WLC_v**2))
    end select
  end subroutine init_sprforce

  !> Constructor for sprforce type
  !! \param id The rank of the process
  subroutine init_sprforce_t(this,id,ntotsegx3)

    class(sprforce),intent(inout) :: this
    integer,intent(in) :: id,ntotsegx3

    allocate(this%Fs(ntotsegx3))
  
  end subroutine init_sprforce_t

  !> Updates the force by adding spring force contribution
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param F totoal force on particles
  subroutine update_force(this,Rbx,Rby,Rbz,bs,invbs,itime,nchain,nseg,nbead,&
                          ntotseg,ntotsegx3,ntotbead,ntotbeadx3,Qt)

    use :: arry_mod, only: print_vector
    use :: conv_mod, only: Bbar_vals,Bbar_cols,Bbar_rowInd
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m,tanb,sinth,costh
    use :: force_smdlt, only: Fphi,rFphi
    
    class(sprforce),intent(inout) :: this
    real(wp),intent(in) :: Rbx(:)
    real(wp),intent(in) :: Rby(:)
    real(wp),intent(in) :: Rbz(:)
    real(wp),intent(in) :: bs(3),invbs(3)
!    real(wp),intent(inout) :: F(:)
    integer,intent(in) :: itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbead,ntotbeadx3
    real(wp),intent(in) :: Qt(:)
    integer :: its,ich,osb,oss,is
    real(wp) :: qx,qy,qz,qsq,q,Ftmp,qytmp


!!$omp parallel default(private) &
!!$omp shared(this,ntotseg,nchain,nbead,nseg,Rbx,Rby,Rbz,bs,invbs)  &
!!$omp shared(ForceLaw,b,qmx,FlowType,eps_m,tanb,sinth,costh,itime) &
!!$omp shared(WLC_v,WLC_A,WLC_B) reduction(-:rFphi) 
!!$omp do schedule(auto)
    do its=1, ntotseg
      ! ich=(its-1)/nseg+1
      ! oss=(ich-1)*nseg
      ! osb=(ich-1)*nbead
      ! is=its-oss
      ! qx=Rbx(osb+is+1)-Rbx(osb+is)
      ! qy=Rby(osb+is+1)-Rby(osb+is)
      ! qz=Rbz(osb+is+1)-Rbz(osb+is)
      qx=Qt((its-1)*3+1)
      qy=Qt((its-1)*3+2)
      qz=Qt((its-1)*3+3)

      qx=qx-nint(qx*invbs(1))*bs(1)
      qy=qy-nint(qy*invbs(2))*bs(2)
      qz=qz-nint(qz*invbs(3))*bs(3)
      select case (FlowType)
        case ('PSF')
          qx=qx+eps_m*qy
        case ('PEF')
          qytmp=qy
          qx=qx+tanb*qytmp
          qy=sinth*qx+costh*qytmp
          qx=costh*qx-sinth*qytmp
      end select
      qsq=qx*qx+qy*qy+qz*qz
!if ((itime==1).or.(itime==100)) then
!if (its<=10) then
!!print *,'eps_m',eps_m
!print *,'its',its
!print *,'qv',qx,qy,qz
!print *,'qsq',qsq
!end if
!end if
      select case (ForceLaw)
        case ('FENE')
          Ftmp = 1/(1-qsq/b)
        case ('WLC_MS')
          q=sqrt(qsq)
          Ftmp = 2*qmx/(3*q)*(0.25/((1-q/qmx)**2)-0.25+q/qmx)
        case ('WLC_UD')
          Ftmp=2._wp/3*(1/(1-qsq/b)**2-7/(2*WLC_v*(1-qsq/b))+WLC_A+WLC_B*(1-qsq/b))
        case ('ILCCP')
          Ftmp = (1-qsq/b/3)/(1-qsq/b)
        case ('Hookean')
          Ftmp = 1._wp
      end select
      ! this%Fs((oss+is-1)*3+1)=Ftmp*qx
      ! this%Fs((oss+is-1)*3+2)=Ftmp*qy
      ! this%Fs((oss+is-1)*3+3)=Ftmp*qz
      this%Fs((its-1)*3+1)=Ftmp*qx
      this%Fs((its-1)*3+2)=Ftmp*qy
      this%Fs((its-1)*3+3)=Ftmp*qz

      rFphi(1)=rFphi(1)-qx*Ftmp*qx
      rFphi(2)=rFphi(2)-qx*Ftmp*qy
      rFphi(3)=rFphi(3)-qy*Ftmp*qy
      rFphi(4)=rFphi(4)-qz*Ftmp*qz
    end do
!!$omp end do
!!$omp end parallel
#ifdef USE_DP
      call mkl_dcsrmv('T',ntotsegx3,ntotbeadx3,-1._wp,'GIIF',Bbar_vals,&
                      Bbar_cols,Bbar_rowInd,Bbar_rowInd(2),this%Fs,1._wp,Fphi)
#elif USE_SP
      call mkl_scsrmv('T',ntotsegx3,ntotbeadx3,-1._wp,'GIIF',Bbar_vals,&
                      Bbar_cols,Bbar_rowInd,Bbar_rowInd(2),this%Fs,1._wp,Fphi)
#endif

  end subroutine update_force

  !> Destructor for spring force type
  subroutine del_sprforce(this)

    type(sprforce),intent(inout) :: this

  end subroutine del_sprforce

end module sprforce_mod
