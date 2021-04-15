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
             update_force    ,&
             update_bendforce

  !> A public type for EV force calcualtion
  type, extends(force) :: sprforce

    private
    !> The spring force
    real(wp),allocatable :: Fs(:)
    !> Total Bending force acting on the particles
    real(wp),allocatable :: Fbnd(:)

  contains

    procedure,pass(this) :: init => init_sprforce_t
    procedure,pass(this) :: update => update_force
    procedure,pass(this) :: updatebend => update_bendforce
    final :: del_sprforce

  end type sprforce

  ! Private module variables:
!  private ::
  ! Protected module variables:

  protected :: ForceLaw_i,ForceLaw,b,qmx,WLC_v,WLC_A,WLC_B, WLC_C, RWS_v,RWS_C,RWS_D


  !> The type of ev force
  character(len=10),save :: ForceLaw
  integer,save :: ForceLaw_i
  !> Different types of force law
  integer,parameter :: Hookean=1     ,&
                       FENE   =2     ,&
                       ILCCP  =3     ,&
                       WLC_MS =4     ,&
                       WLC_UD =5     ,&
                       WLC_GEN=6    ,&
                       WLC_SK =6     ,&
                       RWS    =7

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
  !> Force parameter for WLC-UD model
  real(wp),save :: WLC_C
  !> Force parameter for WLC-UD model
  real(wp),save :: RWS_v
  !> Force parameter for WLC-UD model
  real(wp),save :: RWS_C
  !> Force parameter for WLC-UD model
  real(wp),save :: RWS_D

contains

  !> Initialization of the sprforce module
  !! \param id The rank of the process
  subroutine init_sprforce(id)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id
    integer :: il,j,ntokens,u1,stat,ios
    character(len=1024) :: line
    character(len=100) :: tokens(50)
#ifdef Debuge_sequence
    write(*,*) "module:sprforce_mod:init_sprforce"
#endif
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
              call value(tokens(j+1),RWS_v,ios)
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
    case('RWS')
      ForceLaw_i=RWS
    case('ILCCP')
      ForceLaw_i=ILCCP
    case('WLC_MS')
      ForceLaw_i=WLC_MS
    case('WLC_UD')
      ForceLaw_i=WLC_UD
    case('WLC_SK')
      ForceLaw_i=WLC_SK
    case('WLC_GEN')
      ForceLaw_i=WLC_GEN
    case default
      print('(" The selected force law is not available.")')
      stop
    end select

    select case (ForceLaw)
      case ('WLC_UD','WLC_SK','WLC_GEN')
        WLC_A=3._wp/32-3/(8*WLC_v)-3/(2*WLC_v**2)
        WLC_B=(13._wp/32+0.4086_wp/WLC_v-14.79_wp/(4*WLC_v**2))/ &
                    (1-4.225_wp/(2*WLC_v)+4.87_wp/(4*WLC_v**2))
        WLC_C=(1-1.2370_wp*(2*WLC_v)+0.8105_wp*(4*WLC_v**2))/ &
                    (2*WLC_v- 1.0243_wp*(4*WLC_v**2)+0.4595_wp*(8*WLC_v**3))
      case ('RWS')
        RWS_C=3-10._wp/(3*RWS_v)+10._wp/(27*RWS_v*RWS_v)
        RWS_D=1+2._wp/(3*RWS_v)+10._wp/(27*RWS_v*RWS_v)
    end select
  end subroutine init_sprforce


  !> Constructor for sprforce type
  !! \param id The rank of the process
  subroutine init_sprforce_t(this,id,ntotsegx3)

    class(sprforce),intent(inout) :: this
    integer,intent(in) :: id,ntotsegx3
#ifdef Debuge_sequence
    write(*,*) "module:sprforce_mod:init_sprforce_t"
#endif
    allocate(this%Fs(ntotsegx3))

  end subroutine init_sprforce_t



  !> Updates the force by adding spring force contribution
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param F totoal force on particles
  !!> Called by CalcForce
  subroutine update_force(this,Rbx,Rby,Rbz,bs,invbs,itime,nchain,nseg,nbead,&
                          ntotseg,ntotsegx3,ntotbead,ntotbeadx3,Qt)
    !MB
    !  subroutine update_force(this,Rbx,Rby,Rbz,bs,invbs,itime,ntotseg,ntotsegx3,ntotbead,ntotbeadx3,Qt)

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
    integer,intent(in) :: itime,ntotseg,ntotsegx3,ntotbead,ntotbeadx3
    integer,intent(in) :: nchain,nseg,nbead !not needed here
    real(wp),intent(in) :: Qt(:)
    integer :: its,ich,osb,oss,is
    real(wp) :: qx,qy,qz,qsq,q,Ftmp,qytmp
    logical :: Qbcomptest
    Qbcomptest=.FALSE.
#ifdef Debuge_sequence
    write(*,*) "module:sprforce_mod:update_force"
#endif
!!$omp parallel default(private) &
!!$omp shared(this,ntotseg,nchain,nbead,nseg,Rbx,Rby,Rbz,bs,invbs)  &
!!$omp shared(ForceLaw,b,qmx,FlowType,eps_m,tanb,sinth,costh,itime) &
!!$omp shared(WLC_v,WLC_A,WLC_B,RWS_D,RWS_C) reduction(-:rFphi)
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

      !Write(*,*) 'Qtspr',its, Qt((its-1)*3+1:(its-1)*3+3)

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

     if (qsq >= b) then
        Qbcomptest=.TRUE.
     end if
    !if ((itime==1).or.(itime==100)) then
    !if (its<=10) then
    !!print *,'eps_m',eps_m
    !print *,'its',its
    !print *,'qv',qx,qy,qz
    !print *,'qsq',qsq
    !end if
    !end if
    !   qmx=sqrt(b)
      select case (ForceLaw)
        case ('FENE')
            Ftmp = 1/(1-qsq/b)
        case ('WLC_MS')
            q=sqrt(qsq)
            Ftmp = 2*qmx/(3*q)*(0.25/((1-q/qmx)**2)-0.25+q/qmx)
        case ('WLC_UD','WLC_SK','WLC_GEN')
            Ftmp=2._wp/3*(1/(1-qsq/b)**2-7/(2*WLC_v*(1-qsq/b))+WLC_A+WLC_B*(1-qsq/b))
        case ('ILCCP')
            Ftmp = (1-qsq/b/3)/(1-qsq/b)
        case ('Hookean')
            Ftmp = 1._wp
        case ('RWS')
       ! Eq52, JOR49, 2005 Underhill and Doyle
            Ftmp =(RWS_C/qmx)*(1-(RWS_D/RWS_C)*(qsq/b))/(1-qsq/b)
      end select
      ! this%Fs((oss+is-1)*3+1)=Ftmp*qx
      ! this%Fs((oss+is-1)*3+2)=Ftmp*qy
      ! this%Fs((oss+is-1)*3+3)=Ftmp*qz
!      Write(*,*) its,'Fseg', Ftmp
      !! \Spring Force(X,Y,Z)= Ftmp* q(X,Y,Z)
      this%Fs((its-1)*3+1)=Ftmp*qx
      this%Fs((its-1)*3+2)=Ftmp*qy
      this%Fs((its-1)*3+3)=Ftmp*qz
      !! \R_v.Fs_v=(r_v-r_(v-1))*(Fc_v-Fc_(v-1))----> =-Q_v*Fc_v
      rFphi(1)=rFphi(1)-qx*Ftmp*qx
      rFphi(2)=rFphi(2)-qx*Ftmp*qy
      rFphi(3)=rFphi(3)-qy*Ftmp*qy
      rFphi(4)=rFphi(4)-qz*Ftmp*qz
    end do
!!$omp end do
!!$omp end parallel
    if (Qbcomptest) then
        write(*,*) " ==> Warning: [Q_spr > b] has happened <=="
    end if
#ifdef USE_DP
      call mkl_dcsrmv('T',ntotsegx3,ntotbeadx3,-1._wp,'GIIF',Bbar_vals,&
                      Bbar_cols,Bbar_rowInd,Bbar_rowInd(2),this%Fs,1._wp,Fphi)
#elif USE_SP
      call mkl_scsrmv('T',ntotsegx3,ntotbeadx3,-1._wp,'GIIF',Bbar_vals,&
                      Bbar_cols,Bbar_rowInd,Bbar_rowInd(2),this%Fs,1._wp,Fphi)
#endif
!    do its=1,ntotseg+1
!        Write(*,*) its,'Fphi', Fphi(its*3-2:its*3)
!    end do



  end subroutine update_force


  !> Updates the bending force For WLC_SK/GEN
  !! \param R coordinate of the position vector
  !! \param  this%Fbnd totoal Bending force on particles
  !!> Called by Move_box
  subroutine update_bendforce(this,Qt,R,nchain,nseg,nbead,nchain_cmb,nseg_cmb,ntotbead,nseg_cmbbb,nseg_cmbar,add_cmb,Na,Ia,bs,invbs)

    use :: arry_mod, only: print_vector
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m,tanb,sinth,costh
    use :: force_smdlt, only: rFphi,Fphi

    class(sprforce),intent(inout) :: this

    !real(wp),ALLOCATABLE ::  this%Fbnd(:)
    real(wp),intent(in) :: R(:),Qt(:)
	integer,intent(in) :: Ia(:)
    real(wp),intent(in) :: bs(3),invbs(3)
    integer,intent(in) :: nchain,nseg,nbead,nseg_cmb,nseg_cmbar,nchain_cmb,nSeg_cmbbb,ntotbead,Na
    real(wp) :: thta(-2:0),cost(-2:0),thtal,thtar,costl,costr
    real(wp) :: thta_s,cost_s
    real(wp) :: qtmp(3,-2:0),qmg(-2:0),ehat(3,-2:0)
    real(wp) :: qtmpl(3),qtmpr(3),qmgl,qmgr,ehatl(3),ehatr(3)
    integer  :: os,nbead_cmb,osl,iarm,osSbb,ntotbeadx3,nbead_cmbbb
    integer  :: osb,oss,Osb1,OsS1,oslbbb,ibead,ichain,ibead_arm
    logical  :: add_cmb

    nbead_cmb=nseg_cmb+1
    nbead_cmbbb=nseg_cmbbb+1
	!write(*,*) "Bending called",ForceLaw
    select case (ForceLaw)
      case ('WLC_SK','WLC_GEN')
#ifdef Debuge_sequence
        write(*,*) "Bending !",ForceLaw
	    write(*,*) "module:sprforce_mod:update_bendforce_SK"
#endif
        ntotbeadx3=ntotbead*3
        allocate(this%Fbnd(ntotbeadx3))
        this%Fbnd=0._wp
        !do ibead=1,(nchain*nseg+nchain_cmb*nseg_cmb)
		! Write(*,*) 'Qt',ibead, Qt((ibead-1)*3+1:(ibead-1)*3+3)
		!end do
!!$omp parallel default(private) &
!!$omp shared(this,nchain,nbead,nseg,nseg_cmbar,nchain_cmb,nbead_cmbbb,nseg_cmbbb,nseg_cmb,R,Q,Na,Ia,sinth,tanb,costh,eps_m,invbs,bs,FlowType,add_cmb)  &
!!$omp shared(WLC_C) reduction(-:rFphi)
!!$omp do schedule(auto)

        if (nbead >= 3 .and. nchain /= 0) then
         do ichain=1, nchain

          Osb1=(ichain-1)*nbead
          OsS1=(ichain-1)*nSeg

          do ibead=3, nbead

            osS=OsS1+(ibead-1)
            osb=Osb1+ ibead

            qtmp(:,-1)=Qt((osS-1)*3+1:(osS-1)*3+3)
            qtmp(:,-2)=Qt((osS-2)*3+1:(osS-2)*3+3)

            qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
			qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
            select case (FlowType)
              case ('PSF')
                qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
				qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
              case ('PEF')
                qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
				qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
            end select

            qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2)))
            qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1)))
            ehat(:,-2)=qtmp(:,-2)/qmg(-2)
            ehat(:,-1)=qtmp(:,-1)/qmg(-1)
            thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
            cost(-1)=cos(thta(-1))
            !!Fi
            this%Fbnd(osb*3-2:osb*3)= this%Fbnd(osb*3-2:osb*3)+&
                         WLC_C*(1/qmg(-1))*(ehat(1:3,-2)-cost(-1)*ehat(1:3,-1))
            !!Fi-1
             this%Fbnd((osb-1)*3-2:(osb-1)*3)= this%Fbnd((osb-1)*3-2:(osb-1)*3)+&
                        WLC_C*( ehat(1:3,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))-&
                                ehat(1:3,-2)*(1/qmg(-1)+cost(-1)/qmg(-2)) )
            !!Fi-2
             this%Fbnd((osb-2)*3-2:(osb-2)*3)= this%Fbnd((osb-2)*3-2:(osb-2)*3)+&
                        WLC_C*(1/qmg(-2))*( cost(-1)*ehat(1:3,-2)-ehat(1:3,-1))
            end do !ibead
         end do !ichain
        end if


        if (add_cmb) then
         !write(*,*) "nbead_cmbbb",nbead_cmbbb,"nchain_cmb",nchain_cmb
         do ichain=1, nchain_cmb

          Osb1=nchain*nbead+(ichain-1)*(nSeg_cmb+1)
          OsS1=nchain*nSeg+(ichain-1)*nSeg_cmb

          ! Loop over backbone
          do ibead=3, nbead_cmbbb

            osS=OsS1+(ibead-1)
            osb=Osb1+ ibead
            qtmp(:,-1)=Qt((osS-1)*3+1:(osS-1)*3+3)
            qtmp(:,-2)=Qt((osS-2)*3+1:(osS-2)*3+3)

            !write(*,*) 'before',qtmp(:,-1)
            qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
			!write(*,*) 'after',qtmp(:,-1)
			qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
            select case (FlowType)
              case ('PSF')
                qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
				qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
              case ('PEF')
                qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
				qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
            end select

            qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2)))
            qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1)))
            ehat(:,-2)=qtmp(:,-2)/qmg(-2)
            ehat(:,-1)=qtmp(:,-1)/qmg(-1)
            thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
            cost(-1)=cos(thta(-1))

            !!Fi
             this%Fbnd(osb*3-2:osb*3)= this%Fbnd(osb*3-2:osb*3)+&
                        WLC_C*(1/qmg(-1))*(ehat(:,-2)-cost(-1)*ehat(:,-1))

            !!Fi-1
            this%Fbnd((osb-1)*3-2:(osb-1)*3)= this%Fbnd((osb-1)*3-2:(osb-1)*3)+&
                        WLC_C*( ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))-&
                                ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2)) )
            !!Fi-2
            this%Fbnd((osb-2)*3-2:(osb-2)*3+0)= this%Fbnd((osb-2)*3-2:(osb-2)*3)+&
                        WLC_C*(1/qmg(-2))*(cost(-1)*ehat(:,-2)-ehat(:,-1))

            !Write(*,*) ichain,ibead,osb-2,this%Fbnd((osb-2)*3-2:(osb-2)*3+0)
            !Write(*,*) ichain,ibead,osb-1,this%Fbnd((osb-1)*3-2:(osb-1)*3+0)
			!Write(*,*) ichain,ibead,osb , this%Fbnd(osb*3-2:osb*3)



            end do !ibead_cmbbb
          !go to 12
          ! Loop over Arms
          ! Bead[end of BB]
          Osb1 = nchain*nbead + (ichain-1)*nbead_cmb+nbead_cmbbb
          !Segment Qt[end of BB]
          OsS1= nchain*nSeg + (ichain-1)*nSeg_cmb+nSeg_cmbbb

          if (nseg_cmbar == 1) then
            !write(*,*) "nseg_cmbar == 1"
            do iarm=1, Na
              !Backbone _Arm
              !Segment of the backbone connected to the arm
              !Ia(:) Backbone bead place of arm
              osSbb=nchain*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)
              !bead id for iarm start  !Ia: Backbone bead place of arm
              oslbbb=nchain*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))
              !Seg-Arm
              osS=OsS1+(iarm-1)*nseg_cmbar
              ! osb+1 =Bead# start arm iarm
              osb=osb1+(iarm-1)*(nseg_cmbar)   ! nbead_arm = seg_cmbarm
              !-------------------------------------
              ! Effect of first bead of arm on the Backbone.
              ! Left and right segments
              ! ...O-L-O-R-O...
              !        | arm Seg 1
              !        O arm bead 1
              ! -------------------------------
              qtmpl(:)= -Qt(osSbb*3-2:osSbb*3)
              qtmpr(:)= Qt(osSbb*3+1:osSbb*3+3)
              qtmp(:,-1)= Qt(osS*3+1:osS*3+3)  ! 1st segment of the arm

              qtmpl=QtEq(qtmpl,invbs,bs)
              qtmpr=QtEq(qtmpr,invbs,bs)
              qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
              !qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
              select case (FlowType)
                case ('PSF')
                  qtmpl=QtPSF(qtmpl,eps_m)
				  qtmpr=QtPSF(qtmpr,eps_m)
				  qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
                  !qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
                case ('PEF')
				  qtmpl=QtPEF(qtmpl,sinth,tanb,costh)
				  qtmpr=QtPEF(qtmpr,sinth,tanb,costh)
                  qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
                  !qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
              end select

              qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
              qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
              ehatl(:)=qtmpl(:)/qmgl
              ehatr(:)=qtmpr(:)/qmgr

              qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) !1st_seg arm
              ehat(:,-1)=qtmp(:,-1)/qmg(-1)  !1st_seg arm
              thtal=acos(dot(qtmpl,qtmp(:,-1))/(qmgl*qmg(-1)))
              thtar=acos(dot(qtmpr,qtmp(:,-1))/(qmgr*qmg(-1)))
              costl=cos(thtal)
              costr=cos(thtar)
              ! force on BB Bead Left  n-1
               this%Fbnd(oslbbb*3-5:oslbbb*3-3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+&
                                          WLC_C*(1/qmgl)*(costl*ehatl(:)-ehat(:,-1))
              ! force on BB Bead Right n+1
               this%Fbnd(oslbbb*3+1:oslbbb*3+3)= this%Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                          WLC_C*(1/qmgr)*(costr*ehatr(:)-ehat(:,-1))
              ! force on BB Bead Right n
               this%Fbnd(oslbbb*3-2:oslbbb*3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C*&
                            (  (ehat(:,-1)*(1/qmgl  +costl/qmg(-1)) -ehatl(:)*(1/qmg(-1)+costl/qmgl))-&
                               (ehat(:,-1)*(1/qmgr  +costr/qmg(-1)) -ehatr(:)*(1/qmg(-1)+costr/qmgr)) )
              ! 1st bead of the arm F(v,v-1)+F(v,v+1)
               this%Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( (1/qmg(-1))*(ehatr(:)-costr*ehat(:,-1))-&
                                                      (1/qmg(-1))*(ehatl(:)-costl*ehat(:,-1)) )
              end do

          else if (nseg_cmbar == 2) then
            !write(*,*) "nseg_cmbar == 2"
            do iarm=1, Na
              !Backbone _Arm
              !Segment of the backbone connected to the arm
              !Ia: Backbone bead place of arm
              osSbb=nchain*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)
              !bead id for iarm start  !Ia: Backbone bead place of arm
              oslbbb=nchain*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))
              !Seg-Arm
              osS=OsS1+(iarm-1)*nseg_cmbar
              ! osb/osS +1 = First Beadtot#/Segtot# start arm iarm
              osb=osb1+(iarm-1)*(nseg_cmbar)   ! nbead_arm = nseg_cmbarm

              qtmpl(:)= Qt(osSbb*3-2:osSbb*3)
              qtmpr(:)= Qt(osSbb*3+1:osSbb*3+3)

              qtmp(:,-1)=Qt(osS*3+1:osS*3+3)          ! First segment of the arm
              qtmp(:,-2)=Qt((osS+1)*3+1:(osS+1)*3+3)

              qtmpl=QtEq(qtmpl,invbs,bs)
              qtmpr=QtEq(qtmpr,invbs,bs)
              qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
              qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
              select case (FlowType)
                case ('PSF')
                  qtmpl=QtPSF(qtmpl,eps_m)
				  qtmpr=QtPSF(qtmpr,eps_m)
				  qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
                  qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
                case ('PEF')
				  qtmpl=QtPEF(qtmpl,sinth,tanb,costh)
				  qtmpr=QtPEF(qtmpr,sinth,tanb,costh)
                  qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
                  qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
              end select

			  qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
              qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
              ehatl(:)=qtmpl(:)/qmgl
              ehatr(:)=qtmpr(:)/qmgr
              qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) ! 1St_seg arm
              qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2)))
              ehat(:,-1)=qtmp(:,-1)/qmg(-1)  !1st_seg arm
              ehat(:,-2)=qtmp(:,-2)/qmg(-2)
              thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2))) !Angle Seg1&Seg2
              cost(-1)=cos(thta(-1))
              ! force on BB Bead Left  n-1
               this%Fbnd(oslbbb*3-5:oslbbb*3-3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+&
                                      WLC_C*(1/qmgl)*(costl*ehatl(:)-ehat(:,-1))
              ! force on BB Bead Right n+1
               this%Fbnd(oslbbb*3+1:oslbbb*3+3)= this%Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                        WLC_C*(1/qmgr)*(costr*ehatr(:)-ehat(:,-1))
              ! force on BB Bead Right n
               this%Fbnd(oslbbb*3-2:oslbbb*3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C* &
                            ( (ehat(:,-1)*(1/qmgl  +costl/qmg(-1)) -ehatl(:)*(1/qmg(-1)+costl/qmgl))-&
                              (ehat(:,-1)*(1/qmgr  +costr/qmg(-1)) -ehatr(:)*(1/qmg(-1)+costr/qmgr))+&
                              (1/qmg(-1))*(cost(-1)*ehat(:,-1)- ehat(:,-2)) )
              ! 1st bead of the arm 2@F(v,v-1)+F(v,v)
               this%Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( (1/qmg(-1))*(ehatr(:)-costr*ehat(:,-1))-&
                                                      (1/qmg(-1))*(ehatl(:)-costl*ehat(:,-1))+&
                                                      (ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2))-&
                                                       ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1)) ) )
              ! 2ed Bead F(v,v-1)
               this%Fbnd((osb+1)*3+1:(osb+1)*3+3)=WLC_C*(1/qmg(-2)*(ehat(:,-1)-cost(-1)*ehat(:,-2)))
              end do
          else if (nseg_cmbar >= 3) then
            !write(*,*) "nseg_cmbar >= 3", Na
            do iarm=1, Na
              !Backbone _Arm
              !Segment of the backbone connected to the arm
              !Ia: Backbone bead place of arm
              osSbb=nchain*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)
              !bead id for iarm start  !Ia: Backbone bead place of arm
              oslbbb=nchain*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))
              !Seg-Arm
              osS=OsS1+(iarm-1)*nseg_cmbar    ! +1 =1st Seg/bead iarm
              osb=osb1+(iarm-1)*(nseg_cmbar)   ! nbead_arm = nseg_cmbarm
              qtmpl(:)= Qt(osSbb*3-2:osSbb*3)
              qtmpr(:)= Qt(osSbb*3+1:osSbb*3+3)

              qtmp(:,-1)=Qt(osS*3+1:osS*3+3)  ! 1st segment of the arm
              qtmp(:,-2)=Qt(osS*3+4:osS*3+6)  ! 2ed segment of the arm

              qtmpl=QtEq(qtmpl,invbs,bs)
              qtmpr=QtEq(qtmpr,invbs,bs)
              qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
              qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
              select case (FlowType)
                case ('PSF')
                  qtmpl=QtPSF(qtmpl,eps_m)
				  qtmpr=QtPSF(qtmpr,eps_m)
				  qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
                  qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
                case ('PEF')
				  qtmpl=QtPEF(qtmpl,sinth,tanb,costh)
				  qtmpr=QtPEF(qtmpr,sinth,tanb,costh)
                  qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
                  qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
              end select
              qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
              qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
              ehatl(:)=qtmpl(:)/qmgl
              ehatr(:)=qtmpr(:)/qmgr
              qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) !1st_seg arm
              ehat(:,-1)=qtmp(:,-1)/qmg(-1)  !1st_seg arm
              qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2)))
              ehat(:,-2)=qtmp(:,-2)/qmg(-2)
              thtal=acos(dot(qtmpl,qtmp(:,-1))/(qmgl*qmg(-1)))
              thtar=acos(dot(qtmpr,qtmp(:,-1))/(qmgr*qmg(-1)))
              costl=cos(thtal)
              costr=cos(thtar)
              thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
              cost(-1)=cos(thta(-1))

              ! force on BB Bead Left  n-1
               this%Fbnd(oslbbb*3-5:oslbbb*3-3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+&
                                          WLC_C*(1/qmgl)*(costl*ehatl(:)-ehat(:,-1))
              ! force on BB Bead Right n+1
               this%Fbnd(oslbbb*3+1:oslbbb*3+3)= this%Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                          WLC_C*(1/qmgr)*(costr*ehatr(:)-ehat(:,-1))
              ! force on BB Bead Right n
               this%Fbnd(oslbbb*3-2:oslbbb*3)= this%Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C* &
                            ( (ehat(:,-1)*(1/qmgl  +costl/qmg(-1)) -ehatl(:)*(1/qmg(-1)+costl/qmgl))-&
                              (ehat(:,-1)*(1/qmgr  +costr/qmg(-1)) -ehatr(:)*(1/qmg(-1)+costr/qmgr))+&
                              (1/qmg(-1))*(cost(-1)*ehat(:,-1)-ehat(:,-2)) )

              !1st bead of the arm 2@F(v,v-1)+F(v,v) !!                     +  F(v,v+1) on the loop
               this%Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( (1/qmg(-1))*(ehatr(:)-costr*ehat(:,-1))-&
                                                      (1/qmg(-1))*(ehatl(:)-costl*ehat(:,-1))+&
                                                      (ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2))-&
                                                       ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))) )
              ! 2ed Bead F(v,v-1) !!!            +F(v,v)+ F(v,v+1) on the loop
               this%Fbnd((osb+1)*3+1:(osb+1)*3+3)=WLC_C*(1/qmg(-2))*(ehat(:,-1)-cost(-1)*ehat(:,-2))


              do ibead_arm=3,nseg_cmbar
                osS=OsS1+(iarm-1)*nseg_cmbar +ibead_arm
                osb=osb1+(iarm-1)*(nseg_cmbar)+ibead_arm ! seg_cmbarm=nbead_arm
				! Index not like before!
                qtmp(:, 0)=Qt((osS-1)*3+1:(osS-1)*3+3)
                qtmp(:,-1)=Qt((osS-2)*3+1:(osS-2)*3+3)
                qtmp(:,-2)=Qt((osS-3)*3+1:(osS-3)*3+3)

				qtmp(:, 0)=QtEq(qtmp(:, 0),invbs,bs)
                qtmp(:,-1)=QtEq(qtmp(:,-1),invbs,bs)
			    qtmp(:,-2)=QtEq(qtmp(:,-2),invbs,bs)
                select case (FlowType)
                  case ('PSF')
                    qtmp(:, 0)=QtPSF(qtmp(:, 0),eps_m)
				    qtmp(:,-1)=QtPSF(qtmp(:,-1),eps_m)
                    qtmp(:,-2)=QtPSF(qtmp(:,-2),eps_m)
                  case ('PEF')
				    qtmp(:,0) =QtPEF(qtmp(:,0),sinth,tanb,costh)
                    qtmp(:,-1)=QtPEF(qtmp(:,-1),sinth,tanb,costh)
                    qtmp(:,-2)=QtPEF(qtmp(:,-2),sinth,tanb,costh)
                end select

                qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2))) !qmg(-1)
                qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) !qmg( 0)
                qmg( 0)=sqrt(dot(qtmp(:,0),qtmp(:,0)))
                ehat(:,-2)=qtmp(:,-2)/qmg(-2)
                ehat(:,-1)=qtmp(:,-1)/qmg(-1)
                ehat(:, 0)=qtmp(:,0) /qmg(0)
				thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
                cost(-1)=cos(thta(-1))

                thta(-2)=acos(dot(qtmp(:,-1),qtmp(:,0))/(qmg(-1)*qmg(0)))
                cost(-2)=cos(thta(-2))
                !Fi = +Fi,i-1
                 this%Fbnd(osb*3-2:osb*3)= this%Fbnd(osb*3-2:osb*3)+&
                                      WLC_C*(1/qmg(0))*(ehat(:,-1)-cost(-1)*ehat(:,0))
                !Fi-1 = +F(i-1,i-1)
                 this%Fbnd((osb-1)*3-2:(osb-1)*3+0)= this%Fbnd((osb-1)*3-2:(osb-1)*3)+&
                                      WLC_C*( ehat(:,0)*(1/qmg(-1)+cost(-1)/qmg(0))-&
                                              ehat(:,-1)*(1/qmg(0)+cost(-1)/qmg(-1)) )
                !Fi-2 = +F(i-2,i-2+1)
                 this%Fbnd((osb-2)*3-2:(osb-2)*3)= this%Fbnd((osb-2)*3-2:(osb-2)*3)+&
                                      WLC_C*(1/qmg(-1))*(cost(-1)*ehat(:,-1)-ehat(:,0))

                end do
            end do
12        end if !seg_cmbar
         end do !ichain_comb
        end if !add comb


!!$omp end do
!!$omp do schedule(auto)


        !rFphi, Fphi
        do ibead=1, ntotbead

          !Write(*,*) ibead,'Fphi',Fphi(ibead*3-2:ibead*3)
		  !Write(*,*) ibead,'Fbnd',this%Fbnd(ibead*3-2:ibead*3)

          Fphi(ibead*3-2)=Fphi(ibead*3-2)+ this%Fbnd(ibead*3-2)
          Fphi(ibead*3-1)=Fphi(ibead*3-1)+ this%Fbnd(ibead*3-1)
          Fphi(ibead*3) = Fphi(ibead*3)+ this%Fbnd(ibead*3)

          rFphi(1)=rFphi(1)+ R(ibead*3-2)* this%Fbnd(ibead*3-2) !xx
          rFphi(2)=rFphi(2)+ R(ibead*3-1)* this%Fbnd(ibead*3-2) !yx
          rFphi(3)=rFphi(3)+ R(ibead*3-1)* this%Fbnd(ibead*3-1) !yy
          rFphi(4)=rFphi(4)+ R(ibead*3) *  this%Fbnd(ibead*3)   !zz
         end do

!!$omp end do
!!$omp end parallel

      deallocate(this%Fbnd)
    case default
    !    write(*,*) "Bending Not Working!"
    !Do nothing
    end select

 end subroutine update_bendforce


  function QtEq(q,invbs,bs) result(q1)

    real(wp) :: q1(3)
    real(wp),INTENT(in) :: invbs(:),bs(:),q(:)

    q1(1)=q(1)-nint(q(1)*invbs(1))*bs(1)
    q1(2)=q(2)-nint(q(2)*invbs(2))*bs(2)
    q1(3)=q(3)-nint(q(3)*invbs(3))*bs(3)

  end function QtEq

  function QtPSF(q,eps_m) result(q1)

    real(wp) :: q1(3)
    real(wp),INTENT(in) :: q(:)
    real(wp),INTENT(in) :: eps_m

    q1(1)=q(1)+eps_m*q(2)
	q1(2)=q(2)
	q1(3)=q(3)

  end function QtPSF

  function QtPEF(q,sinth,tanb,costh) result(q1)

    real(wp) :: q1(3)
    real(wp),INTENT(in) :: q(:)
    real(wp),INTENT(in) :: sinth,tanb,costh
	real(wp):: qytmp

    qytmp=q(2)
    q1(1)=q(1)+tanb*qytmp
    q1(2)=sinth*q1(1)+costh*qytmp
    q1(1)=costh*q1(1)-sinth*qytmp
    q1(3)=q(3)

  end function QtPEF


  !> Destructor for spring force type
  subroutine del_sprforce(this)

    type(sprforce),intent(inout) :: this

#ifdef Debuge_sequence
    write(*,*) "module:sprforce_mod:del_sprforce"
#endif


  end subroutine del_sprforce




end module sprforce_mod
