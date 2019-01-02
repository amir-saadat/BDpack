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
module box_mod

  use,intrinsic :: iso_c_binding
  use :: mkl_dfti
  use :: types
  use :: prcn_mod
  use :: flow_mod, only: flow
  use :: chain_mod
  use :: io_mod
  use :: verlet_mod, only: verlet
  use :: sprforce_mod, only: sprforce
  use :: evforce_mod, only: evforce
  use :: trsfm_mod, only: trsfm
#ifdef USE_GPU
  use :: cudafor
  use :: trsfm_cumod, only: trsfm_cu_t
  use :: conv_cumod, only: conv_cu_t
  use :: chain_cumod, only: chain_cu_t
  use :: rndm_cumod, only: rndm_cu_t
  use :: hi_cumod, only: hi_cu_t
  use :: sprforce_cumod, only: sprforce_cu_t
#endif
!$ use :: omp_lib
  
  implicit none

  ! Private module procedures:
  private :: init_box_t   ,&
             move_box     ,&
             read_dmp_cnf ,&
             read_init_cnf,&
             write_cnf    ,&
             calc_mat_fcn ,&
             del_box_t    ,&
             calcHI       ,&
             calcForce

  type :: box

    !> The simulation box ID (process ID)
    integer :: ID
    !> The box initial dimensions
    real(wp) :: size(3)
    !> The inverse of box initial dimensions
    real(wp) :: invsize(3)
    !> The origin of the box
    real(wp) :: origin(3)
    !> An object for applying flow to the entities inside the box
    type(flow) :: Boxflow
    type(decomp) :: decompRes
    !> An object for memory transaction of configurational info
    type(conf_io) :: BoxIO
    !> An object for calculating net spring force
    type(sprforce) :: Boxsprf
    !> An object for calculating ev force
    type(evforce) :: Boxevf
    !> An object for performing transformation  
    type(trsfm) :: Boxtrsfm
    !> The position vector of all beads
    real(wp),allocatable :: Rb_tilde(:)
    !> The x-component of position vector of all beads
    real(wp),allocatable :: Rbx(:)
    !> The y-component of position vector of all beads
    real(wp),allocatable :: Rby(:)
    !> The z-component of position vector of all beads
    real(wp),allocatable :: Rbz(:)
    !> The connectivity vector of all polymer segments
    real(wp),allocatable :: Qdagger_tilde(:)
    !> The center of mass of all chains
    real(wp),pointer :: rcm_tilde(:,:)
    !> The image for the center of mass of all chains
    integer,pointer :: cmif(:,:)
    !> The array for storing the beads transformed position
    real(wp),allocatable :: Rbtr(:,:)
    !> The array for storing the center of mass transformed position
    real(wp),allocatable :: rcmtr(:,:)
    !> The bead to center of mass distance of all polymer beads
    real(wp),allocatable :: R_tilde(:)
    !> The image of the beads inside the primary box
    integer,allocatable :: b_img(:,:)
    !> The image of the center of mass inside the primary box
    integer,allocatable :: cm_img(:,:)
    !> The pointer to the chain section  of configurational arrays
    type(chain),allocatable :: BoxChains(:)
    #ifdef USE_GPU
      !> The position vector of all beads
      real(wp),device,allocatable :: Rb_d(:)
      !> The x-component of position vector of all beads
      real(wp),device,pointer :: Rbx_d(:)
      !> The y-component of position vector of all beads
      real(wp),device,pointer :: Rby_d(:)
      !> The z-component of position vector of all beads
      real(wp),device,pointer :: Rbz_d(:)
      !> The bead to center of mass distance of all polymer beads
      real(wp),device,allocatable :: R_d(:)
      !> The connectivity vector of all polymer segments
      real(wp),device,allocatable :: Q_d(:)
      !> The center of mass of all chains
      real(wp),device,pointer :: rcm_d(:,:)
      !> The image of the beads inside the primary box
      integer,device,pointer :: b_img_d(:,:)
      !> The image of the center of mass inside the primary box
      integer,device,pointer :: cm_img_d(:,:)
      !> The pointer to the chain section of configurational arrays
      type(chain_cu_t),allocatable :: BoxChains_d(:)
      !> An object for performing transformation on device
      type(trsfm_cu_t) :: Boxtrsfm_d
      !> An object for performing conversions on device
      type(conv_cu_t) :: Boxconv_d
      !> An object for random number generation on device
      type(rndm_cu_t) :: Boxrndm_d
      !> An object for hi calculation on device
      type(hi_cu_t) :: Boxhi_d
      !> An object for calculating net spring force on device
      type(sprforce_cu_t) :: Boxsprf_d
    #endif
    
    contains

      procedure,pass(this) :: init => init_box_t
      procedure,pass(this) :: move => move_box
      procedure,pass(this) :: read_init => read_init_cnf
      procedure,pass(this) :: read_dmp => read_dmp_cnf
      procedure,pass(this) :: write => write_cnf 
      procedure,pass(this) :: calc_mf => calc_mat_fcn
      final :: del_box_t

  end type box

  ! real(wp),device,pointer :: Rbx_dp(:)
  ! real(wp),device,pointer :: Rby_dp(:)
  ! real(wp),device,pointer :: Rbz_dp(:)
  ! real(wp),device,pointer :: rcm_dp(:,:)
  ! integer,device,pointer :: b_img_dp(:,:)
  ! integer,device,pointer :: cm_img_dp(:,:)

  
  ! Protected module variables:
  protected :: nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3,lambda

  !> The number of chains in the box
  integer,save :: nchain
  !> The number of segments in a chain
  integer,save :: nseg
  !> The number of beads in a chain
  integer,save :: nbead
  !> The number of segment components in a chain
  integer,save :: nsegx3
  !> The number of bead components in a chain
  integer,save :: nbeadx3
  !> The number of segments in the box
  integer,save :: ntotseg
  !> The number of beads in the box
  integer,save :: ntotbead
  !> The number of segment components in the box
  integer,save :: ntotsegx3
  !> The number of bead components in the box
  integer,save :: ntotbeadx3
  !> The precision of real parameters
  integer,save :: MPI_REAL_WP
  !> The relaxation time of the polymer
  real(wp),save :: lambda
  !> The initial dimension of the box
  real(wp),save :: Lbox(3)

contains

  !> Initializes the box module
  subroutine init_box(myrank,nprun)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env
    use :: conv_mod, only: init_conv
    use :: flow_mod, only: init_flow
    use :: trsfm_mod, only: init_trsfm
    use :: hi_mod, only: init_hi
    use :: verlet_mod, only: init_verlet
    use :: force_smdlt, only: init_force
    use :: sprforce_mod, only: init_sprforce
    use :: evforce_mod, only: init_evforce
    use :: mpi
    #ifdef USE_GPU
      use :: rndm_cumod, only: init_rndm_d
      use :: flow_cumod, only: init_flow_d
      use :: conv_cumod, only: init_conv_d
      use :: hiverlet_cumod, only: init_hiverlet_d
      use :: force_cumod, only: init_force_d
    #endif
    !include 'mpif.h'
    
    integer,intent(in) :: myrank,nprun
    integer :: j,ntokens,u1,il,stat,size_sp,size_dp,ierr,ios
    character(len=1024) :: line 
    character(len=100) :: tokens(10)
    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp) :: b,hstar
    character(len=10) :: LambdaMethod

    ! default setting:
    LambdaMethod='Rouse'
    hstar=0._wp

    open (newunit=u1,action='read',file='input.dat',status='old')
    il=1
    ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" box_mod: Error reading line ",i0, " Process ID ",i0)',il,myrank
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
            case ('nchain')
              call value(tokens(j+1),nchain,ios)
            case ('nseg')
              call value(tokens(j+1),nseg,ios)
            case ('Rel-Model')
              LambdaMethod=trim(adjustl(tokens(j+1)))
              if (LambdaMethod == 'Self') then
                call value(tokens(j+2),lambda,ios)
              end if
            case ('hstar')
              call value(tokens(j+1),hstar,ios)
            case ('b')
              call value(tokens(j+1),b,ios)
            case ('Box-dim')
              call value(tokens(j+1),Lbox(1),ios)
              call value(tokens(j+2),Lbox(2),ios)
              call value(tokens(j+3),Lbox(3),ios)
          end select
        enddo
      end if
    end do ef
    close(u1)

    nbead=nseg+1
    nsegx3=nseg*3
    nbeadx3=nbead*3
    ntotseg=nseg*nchain
    ntotbead=nbead*nchain
    ntotsegx3=ntotseg*3
    ntotbeadx3=ntotbead*3
    ! Longest dimensionless relaxation time setting
    if (LambdaMethod /= 'Self') lambda=0._wp
    select case (LambdaMethod)
      case ('Tanner')
        lambda=(b*nseg+7)/(15*b*nseg)*b/(b+5)*(2*(nseg+1)**2+7-&
               12*((nseg+1)**2+1)/real(nseg+1,kind=wp)/(b+7))
      case ('Zimm')
        if (hstar /= 0._wp) then
          ! lambda_1:
          lambda=nbead**1.5*1.22_wp/(hstar*PI**2)
          ! lambda_eta:
         ! lambda=2.39_wp*(nbead**1.5*1.22_wp/(hstar*PI**2))
        else          
          if (myrank == 0) print '(" Zimm relaxation time is only for non-zero hstar.")'
          stop
        end if
      case ('Rouse')
        ! lambda_1:
        lambda=0.5*1.0/(sin(PI/(2*nbead)))**2
        ! lambda_eta:
       ! lambda=1.64_wp*0.5*1.0/(sin(PI/(2*nbead)))**2
      case ('Prabhakar')
        if (hstar /= 0.0_wp) then
          lambda=0.42_wp*nseg**(1.5)*sqrt(12._wp)/(hstar*PI**1.5)
        end if
      case ('Self')
      case default
        if (myrank == 0) print '(" Error: Incorrect Rel-Model.")'
        stop
    end select
    if (myrank == 0) then
      print *
      print '(" Longest Relaxation Time:",f10.4)',lambda
    end if

    if (lambda == 0._wp) then
      if (myrank == 0) print '(" Error: Incorrect LambdaMethod.")'
      stop
    end if

    ! Set MPI working precision WP:
    call MPI_Type_size(MPI_REAL,size_sp,ierr)
    call MPI_Type_size(MPI_DOUBLE_PRECISION,size_dp,ierr)
    if (wp == size_sp) then
      MPI_REAL_WP=MPI_REAL
    else if (wp == size_dp) then
      MPI_REAL_WP=MPI_DOUBLE_PRECISION
    end if

    !----------------------------------------------
    !>>> Initialization of the modules:
    !----------------------------------------------
    call init_flow(myrank)
    call init_trsfm()
    call init_conv(nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3)
    ! if (hstar /= 0._wp) then
      call init_hi(myrank,ntotbead,ntotbeadx3,Lbox)
    ! end if
    call init_io(myrank,nchain,nsegx3,nbeadx3,ntotbeadx3,nprun)
    call init_verlet(myrank)
    call init_force(ntotbeadx3)
    call init_sprforce(myrank,nseg)
    call init_evforce(myrank)

    #ifdef USE_GPU
      call init_rndm_d(ntotbeadx3)
      call init_flow_d()
      call init_conv_d(nseg,ntotsegx3,ntotbeadx3)
      call init_hiverlet_d()
      call init_force_d(ntotbeadx3)
    #endif

  end subroutine init_box

  subroutine init_box_t(this,myrank,nproc,nprun,runrst)

    use :: arry_mod, only: print_matrix,print_vector
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m
    use :: evforce_mod, only: rc_F
    use :: sprforce_mod, only: qmx
    use :: mpi
    !include 'mpif.h'

    class(box),intent(inout) :: this
    integer,intent(in) :: myrank,nproc,nprun,runrst
    integer :: ichain,offseti,offsetj,info
    integer :: ndiag,isegx3,ibead,jseg,maxnz_K,ierr,offsetch1
    integer :: job_Bbar(8),job_B(8),job_K(8)

    this%ID=myrank
    this%size=Lbox
    this%invsize(:)=1/(this%size)
    this%origin=0._wp

    ! Instantiation of Boxflow:
    call this%Boxflow%init(nchain,nbead)

    select case (FlowType)
      case ('Equil')
        if (CoMDiff) then
          allocate(this%cmif(nchain,3))
        else
          allocate(this%cmif(0,0))
        endif
        allocate(this%Rbtr(0,0))
        allocate(this%rcmtr(0,0))
      case ('PSF')
        allocate(this%Rbtr(ntotbead,1))
        allocate(this%rcmtr(nchain,1))
        allocate(this%cmif(0,0))
      case ('PEF')
        allocate(this%Rbtr(ntotbead,2))
        allocate(this%rcmtr(nchain,2))
        allocate(this%cmif(0,0))
    end select
    ! Instantiation of Boxtrsfm:
    call this%Boxtrsfm%init(nchain,ntotbead,this%Rbtr,this%rcmtr)

    ! Allocation of general "tilde" variables:
    allocate(this%R_tilde(ntotbeadx3))
    allocate(this%rcm_tilde(nchain,3))
    allocate(this%Qdagger_tilde(ntotsegx3))

    ! Array for bead positions:
    allocate(this%Rb_tilde(ntotbeadx3))
    allocate(this%Rbx(ntotbead),&
             this%Rby(ntotbead),&
             this%Rbz(ntotbead))

    allocate(this%b_img(ntotbead,3))
    allocate(this%cm_img(nchain,3))

    ! Allocation of box chains:
    allocate(this%BoxChains(nchain))

    do ichain=1, nchain
      call this%BoxChains(ichain)%init(ichain,nsegx3,nbead,nbeadx3,this%Rb_tilde,this%Rbx,&
                         this%Rby,this%Rbz,this%Qdagger_tilde,this%R_tilde,this%rcm_tilde,&
                         this%cmif,this%Rbtr,this%rcmtr,this%b_img,this%cm_img)
    end do

    ! Instantiation of BoxIO:
    call this%BoxIO%init(myrank,nproc,nchain,nsegx3,nbeadx3,MPI_REAL_WP)

    call this%BoxIO%read(myrank,nproc,this%size,nchain,nseg,nbead,nsegx3,&
              ntotbead,ntotsegx3,ntotbeadx3,nprun,runrst,qmx,MPI_REAL_WP)

    ! Instantiation of Boxevf:
    call this%Boxevf%init(myrank,this%size,ntotbead,ntotbeadx3)

    ! Instantiation of Boxsprf:
    call this%Boxsprf%init(myrank,ntotsegx3)

    #ifdef USE_GPU
      ! Allocation of configurational variables:
      allocate(this%R_d(ntotbeadx3))
      allocate(this%rcm_d(nchain,3))
      allocate(this%Q_d(ntotsegx3))
      ! Array for bead positions:
      allocate(this%Rb_d(ntotbeadx3))
      allocate(this%Rbx_d(ntotbead),this%Rby_d(ntotbead),this%Rbz_d(ntotbead))
      ! beads and com image:
      allocate(this%b_img_d(ntotbead,3))
      allocate(this%cm_img_d(nchain,3))
      ! Allocation of box chains:
      allocate(this%BoxChains_d(nchain))
      do ichain=1, nchain
        call this%BoxChains_d(ichain)%init(ichain,nsegx3,nbead,nbeadx3,&
          this%Rb_d,this%Rbx_d,this%Rby_d,this%Rbz_d,this%Q_d,this%R_d,&
          this%rcm_d,this%b_img_d,this%cm_img_d)
      end do
      ! Instantiation of Boxtrsfm_d
      call this%Boxtrsfm_d%init()
      ! Instantiation of Boxtrsfm_d
      call this%Boxconv_d%init()
      ! Instantiation of Boxrndm_d
      call this%Boxrndm_d%init(myrank)
      ! Instantiation of Boxhi_d
      call this%Boxhi_d%init(myrank,ntotbead,this%size)
      ! Instantiation of Boxsprf_d:
      call this%Boxsprf_d%init(ntotsegx3)
    #endif
    
    if (myrank == 0) then
      print *
      if ((this%size(1) == this%size(2)) .and. &
          (this%size(2) == this%size(3))) then
        print '(" Periodic Box with Equal dimensions:")'
      else
        print '(" Periodic Box with Non-equal dimensions:")'
      end if
      print '(3(f14.5,1x))',this%size(:)
      if ((FlowType == 'PEF') .and. &
          (this%size(1) /= this%size(2))) then
        print'(" In PEF, box should be square in x-y plane.")'
        stop
      end if
    end if

  end subroutine init_box_t

  !> Moving the particles ...
  !!
  !! \param eps current value of applied strain
  subroutine move_box(this,itime,ntime,irun,Pe,dt,col,myrank,eps,itrst)

    use :: diffcalc_mod
    use :: hi_mod, only: Rb0,ncols,dw_bltmp,dispmax,rlist_D,HIcalc_mode,hstar,&
                         rc_D,Diff_tens,DF_tot,dw_bltmpP,dw_bl,update_lst,Coeff_tens
    use :: arry_mod, only: print_vector,print_matrix
    use :: conv_mod, only: RtoRbc,RbctoRb,RbtoRbc,QtoR,RbctoQ
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: zrotate,theta_0,bsx,bsy,invbsx,invbsy,init_trsfm_tm,&
                     tanb,update_trsfm,eps_m,delrx_L,delrx_m,reArng,sinth,costh,&
                     L1,L2,update_arng
    use :: tmng_mod, only: et_CF,et_HI,et_DT,et_DEC,et_PR,et_QR,et_PME,et_FFT,&
           et_IFFT,et_SPR,et_INT,et_INF,et_R,et_K,et_EIKR,et_EW,et_DCR,et_DCK,&
           HIcount,PMEcount,et_CM,tick,tock,doTiming,et_vlt,et_mklfft,et_mklsparse,&
           et_cufft,et_cusparse,et_cuinfl
!    use :: conv_mod, only: Bbar_vals,Bbar_cols,Bbar_rowInd
    use :: io_mod, only: Qst,rcmst,Rbst,CoMDiff
    use :: force_smdlt, only: Fphi,rFphi
    use :: evforce_mod, only: EVForceLaw
    #ifdef USE_GPU
      use :: cublas
      use :: trsfm_cumod, only: init_trsfm_tm_d,update_trsfm_d,update_arng_d
      use :: rndm_cumod, only: dw_bl_d
      use :: force_cumod, only: Fphi_d
      use :: chain_cumod
      use :: diffcalc_cumod, only: PME_d
    #endif
    
    class(box),intent(inout) :: this
    real(wp),intent(in) :: Pe,dt,eps
    integer,intent(in) :: itime,ntime,irun,myrank,itrst
    real(wp) :: rbead(3),Delrx_y,rbeadx_tr,rijtmp(3),ri(3),rj(3),rijx_tr,rcmx,rcmy
    integer :: ichain,ibead,offsetch,offsetb,offset,col,iglobbead,ich
    integer(long) :: count0,count1
    logical :: updateList
    real(wp),pointer :: rcmP(:) => null()
    real(wp),parameter :: coeff=1/sqrt(2.0_wp)
    #ifdef USE_GPU
      type(dim3) :: dimGrid, dimBlock
      integer :: istat,ierr
      real(wp) :: et
    #endif

real(wp),allocatable :: rbtest(:)



    allocate(rbtest(ntotbeadx3))

    !----------------------
    !>>> Initial time step:
    !----------------------

    if (itime == itrst+1) then

      do ichain=1, nchain
        this%BoxChains(ichain)%chain_Q(:)=Qst(:,ichain,irun)
        this%BoxChains(ichain)%chain_rcm(:)=rcmst(ichain,:,irun)
        if ((FlowType == 'Equil').and.CoMDiff) then
          ! this%BoxChains(ichain)%rcm_ImFlag(:)=cmifst(ichain,:,irun)
          this%BoxChains(ichain)%chain_cmif(:)=cmifst(ichain,:,irun)
        end if
      end do
      this%Rb_tilde=Rbst(:,irun)
      this%b_img=0
      
      if (itrst /= 0) then
        call update_arng(eps)
        call update_trsfm(this%size)
      else
        call init_trsfm_tm(this%size)
      end if

      if ((FlowType == 'PEF') .and. (this%BoxIO%initmode == 'st')) then
        call zrotate(this%Qdagger_tilde,-theta_0)
        do ich=1, nchain
          rcmx=this%rcm_tilde(ich,1)
          rcmy=this%rcm_tilde(ich,2)
          this%rcm_tilde(ich,1)= cos(-theta_0)*rcmx+sin(-theta_0)*rcmy
          this%rcm_tilde(ich,2)=-sin(-theta_0)*rcmx+cos(-theta_0)*rcmy
        end do
      end if

      call QtoR(this%Qdagger_tilde,this%R_tilde,ntotsegx3,ntotbeadx3)
      call RtoRbc(this%R_tilde,this%rcm_tilde,this%Rbx,this%Rby,this%Rbz,&
                  nchain,nbead,ntotbead,ntotbeadx3)

      #ifdef USE_GPU
        this%Q_d=this%Qdagger_tilde
        this%rcm_d=this%rcm_tilde
        this%Rb_d=this%Rb_tilde
        this%R_d=this%R_tilde
        this%Rbx_d=this%Rbx
        this%Rby_d=this%Rby
        this%Rbz_d=this%Rbz

        this%b_img_d=0

        if (itrst /= 0) then
          call update_arng_d()
          call update_trsfm_d()
        else
          call init_trsfm_tm_d()
        end if

      #endif

      this%decompRes%Success=.true.

      ! Parameters for timing:
      HIcount=0;PMEcount=0
      et_CF=0._wp;et_HI=0._wp;et_DT=0._wp;et_DEC=0._wp
      et_PR=0._wp;et_R=0._wp;et_K=0._wp;et_CM=0._wp
      et_vlt=0._wp
      select case (HIcalc_mode)
        case ('PME')
          et_PME=0._wp;et_FFT=0._wp;et_IFFT=0._wp;et_SPR=0._wp
          et_INT=0._wp;et_INF=0._wp;et_DCR=0._wp;et_DCK=0._wp

          et_cufft=0._wp
          et_mklfft=0._wp
          et_cusparse=0._wp
          et_mklsparse=0._wp
          et_cuinfl=0._wp

        case ('Ewald')
          et_EW=0._wp;et_EIKR=0._wp
      end select

    end if ! itime == itrst+1

    !-----------------------------
    !>>> Random number generation:
    !-----------------------------
! print*,'itime',itime
! call print_vector(dw_bl(:,1),'dworigin')
    #ifdef USE_GPU
      call this%Boxrndm_d%gen(nchain,ntotbeadx3,dt)
      dw_bl_d=dw_bl
    #endif

    !-----------------------------------------
    !>>> Applying Periodic Boundary Condition:
    !-----------------------------------------

    ! Note that PBC is applied in the beginning of the time step to stay
    ! consistent with calculation of forces.


    #ifdef USE_GPU

      call this%Boxtrsfm_d%applypbc(this%size,this%invsize,this%Rbx_d,this%Rby_d,&
        this%Rbz_d,this%rcm_d,this%b_img_d,this%cm_img_d,nbead,itime)
      call this%Boxconv_d%RbctoRb(this%Rbx_d,this%Rby_d,this%Rbz_d,this%Rb_d,ntotbead)

      this%Qdagger_tilde=this%Q_d
      this%rcm_tilde=this%rcm_d
      this%R_tilde=this%R_d
      this%Rb_tilde=this%Rb_d
      this%Rbx=this%Rbx_d
      this%Rby=this%Rby_d
      this%Rbz=this%Rbz_d


    #else

      call this%Boxtrsfm%applypbc(this%size,this%invsize,this%Rbx,this%Rby,this%Rbz,&
        this%rcm_tilde,this%b_img,this%cm_img,nbead,itime)
      call RbctoRb(this%Rbx,this%Rby,this%Rbz,this%Rb_tilde,ntotbead)

    #endif

    ! rbtest=this%Rb_d
    ! print*,'rb1--0',rbtest
    ! print*,'rb2--0',this%Rb_tilde

    !---------------------------------------------------------
    !>>> Constructing Verlet neighbor list for EV calculation:
    !---------------------------------------------------------
    if (doTiming) call tick(count0)
    if (EVForceLaw /= 'NoEV') then
      #ifdef USE_GPU
        select case (FlowType)
          case ('Equil')
            ! call this%Boxevf_d%update_vlt(this%Rbx,this%Rby,this%Rbz,this%Rb_tilde,&
            !                this%size,this%invsize,itime,itrst,ntotbead,ntotbeadx3)
          case ('PSF')
          case ('PEF')
        end select
      #else
        select case (FlowType)
          case ('Equil')
            call this%Boxevf%update_vlt(this%Rbx,this%Rby,this%Rbz,this%Rb_tilde,&
                           this%size,this%invsize,itime,itrst,ntotbead,ntotbeadx3)
          case ('PSF')
            call this%Boxevf%update_vlt(this%Boxtrsfm%Rbtrx,this%Rby,this%Rbz,&
                                   this%Rb_tilde,this%size,this%invsize,itime,&
                                   itrst,ntotbead,ntotbeadx3)
          case ('PEF')
            call this%Boxevf%update_vlt(this%Boxtrsfm%Rbtrx,this%Boxtrsfm%Rbtry,this%Rbz,&
                    this%Rb_tilde,[bsx,bsy,this%size(3)],[invbsx,invbsy,this%invsize(3)],&
                    itime,itrst,ntotbead,ntotbeadx3)
        end select
      #endif
    end if
    if (doTiming) et_vlt=et_vlt+tock(count0)

    !----------------------------------------------------------
    !>>> Constructing Verlet neighbour-list for HI calculation:
    !----------------------------------------------------------

    if (HIcalc_mode == 'PME') then
      #ifdef USE_GPU

    call update_lst(this%Rb_tilde,this%Boxtrsfm%Rbtrx,itime,itrst,nchain,nbead,&
      nbeadx3,ntotbead,this%size,this%origin)
    call this%Boxhi_d%updatelst(this%Rb_tilde,this%Boxtrsfm%Rbtrx,itime,&
        itrst,nchain,nbead,nbeadx3,ntotbead,this%size,this%origin)


        call this%Boxhi_d%updatelst_(this%Rbx_d,this%Rby_d,this%Rbz_d,this%Rb_d,itime,&
          itrst,nchain,nbead,nbeadx3,ntotbead,ntotbeadx3,this%size,this%origin)


      #else
        call update_lst(this%Rb_tilde,this%Boxtrsfm%Rbtrx,itime,itrst,nchain,nbead,&
          nbeadx3,ntotbead,this%size,this%origin)
      #endif
    end if

    !---------------------------------------------------
    !>>> Calculating conservative forces; spring and EV:
    !---------------------------------------------------

    if (doTiming) call tick(count0)
    call calcForce(this,itime)
    if (doTiming) et_CF=et_CF+tock(count0)

!if ((itime==24).or.(itime==1000)) then
!print*,'itime',itime
!call print_vector(Fphi,'fphi')
!!call print_vector(this%Rb_tilde,'rb')
!!call print_matrix(this%rcm_tilde,'rcm')
!end if
    !-----------------------------------------
    !>>> Calculating hydrodynamic interaction:
    !-----------------------------------------

! print*,'box1','size',size(this%Boxhi_d%P_vals)

    if (doTiming) call tick(count0)
!    if (hstar /= 0._wp) then
      if ( (mod(itime,ncols) == 1) .or. (ncols == 1) ) then
        call calcHI(this,itime,dt)
      end if
! print*,'box2'

!    end if
    if (doTiming) et_HI=et_HI+tock(count0)

    !--------------------------------
    !>>> Integration in Euler scheme:
    !--------------------------------

    ! rbtest=dw_bl_d(:,1)
    ! print*,'dw cuda',rbtest
    ! rbtest=Fphi_d
    ! print*,'fphi cuda',rbtest
    ! rbtest=this%Rb_d
    ! print*,'rb1',rbtest

    if (doTiming) call tick(count0)

    #ifdef USE_GPU


      if (hstar == 0._wp) then
        call cublasDaxpy(ntotbeadx3,0.25_wp*dt,Fphi_d,1,this%Rb_d,1)
        call cublasDaxpy(ntotbeadx3,coeff,dw_bl_d(:,col),1,this%Rb_d,1)
      else
        call PME_d(this%Boxhi_d,Fphi_d,ntotbead,ntotbeadx3,this%size)

        call cublasDaxpy(ntotbeadx3,0.25_wp*dt,Fphi_d,1,this%Rb_d,1)
        ! dw_bl_d(:,col)=dw_bltmp(:,col)
        call cublasDaxpy(ntotbeadx3,coeff,dw_bl_d(:,col),1,this%Rb_d,1)
      endif

      call this%Boxconv_d%RbtoRbc(this%Rb_d,this%Rbx_d,this%Rby_d,this%Rbz_d,ntotbead)

      this%Qdagger_tilde=this%Q_d
      this%rcm_tilde=this%rcm_d
      this%R_tilde=this%R_d
      this%Rb_tilde=this%Rb_d
      this%Rbx=this%Rbx_d
      this%Rby=this%Rby_d
      this%Rbz=this%Rbz_d

    #else

      if (FlowType /= 'Equil') call this%Boxflow%apply(Pe,dt,this%Rb_tilde,ntotbeadx3)

      if (hstar == 0._wp) then
        this%Rb_tilde=this%Rb_tilde+0.25_wp*dt*Fphi+coeff*dw_bl(:,col)
      else
        select case (HIcalc_mode)
          case ('Ewald')
            call symv(Diff_tens,Fphi,this%Rb_tilde,alpha=0.25_wp*dt,beta=1._wp)
          case ('PME')
            ! #ifdef USE_GPU
            !   call PME_dev(this%Boxhi_d,Fphi,ntotbead,this%size,DF_tot)
            ! #else
              call PME_cpu(Fphi,ntotbead,this%size,DF_tot)
            ! #endif
            this%Rb_tilde=this%Rb_tilde+0.25_wp*dt*DF_tot
        end select

        ! #ifdef USE_GPU
        !   this%Rb_tilde=this%Rb_tilde+coeff*dw_bltmp(:,col)
        ! #else
        #ifdef USE_DP
          this%Rb_tilde=this%Rb_tilde+coeff*dw_bltmp(:,col)
        #elif USE_SP
          this%Rb_tilde=this%Rb_tilde+coeff*real(dw_bltmp(:,col),kind=wp)
        #endif
        ! #endif

      end if

      call RbtoRbc(this%Rb_tilde,this%Rbx,this%Rby,this%Rbz,ntotbead)


    #endif

        ! this%b_img=this%b_img_d
        ! this%cm_img=this%cm_img_d
        ! print*,'b_img',this%b_img
        ! print*,'cm_img',this%cm_img
        ! call print_vector(this%Rb_tilde,'rb')
        ! call print_matrix(this%rcm_tilde,'rcm')

        if (itime == 1000 .or. itime==2000) then

        ! #ifdef USE_GPU
        ! dw_bltmp=dw_bl_d
        ! #endif
        ! call print_vector(dw_bltmp(:,col),'dw_bl')
        call print_vector(this%Rb_tilde,'rb')
        ! call print_matrix(this%rcm_tilde,'rcm')

        endif

    if (doTiming) et_PR=et_PR+tock(count0)

    ! print*,'dw_bl',dw_bltmp(:,col)
    ! print*,'df',DF_tot
    ! print*,'rb2',this%Rb_tilde
    ! stop

    call tick(count0)

    #ifdef USE_GPU

    !   call update_chain_krnl <<<(ntotbead+255)/256,256>>>(this%b_img_d,this%cm_img_d,this%rcm_d,&
    !     this%Rbx_d,this%Rby_d,this%Rbz_d,nbead,ntotbead,nchain,this%size(1),this%size(2),this%size(3))
    !   ierr = cudaGetLastError()
    !   if (ierr /= cudaSuccess) then
    !     print '(" update chain error: ",a)', cudaGetErrorString(ierr)
    !     stop
    !   endif

      do ichain=1, nchain
        call this%BoxChains_d(ichain)%update(this%size)
      end do

      this%rcm_tilde=this%rcm_d
      ! call print_matrix(this%rcm_tilde,'rcm')

    #else

      !$omp parallel default(private) shared(nchain,nbead,this)
      !$omp do schedule(auto)
        do ichain=1, nchain
          call this%BoxChains(ichain)%update(nbead,this%size,this%invsize)
        end do
      !$omp end do
      !$omp end parallel

      ! call print_matrix(this%rcm_tilde,'rcm')

    #endif

    if (doTiming) et_CM=et_CM+tock(count0)

    !------------------------------------------------
    !>>> Updating the dimension of the deformed box :
    !------------------------------------------------

    #ifdef USE_GPU

      call update_arng_d()
      call update_trsfm_d()

    #else

      call update_arng(eps)
      if (reArng) then
        call this%Boxtrsfm%unwrap(this%size,this%Rbx,this%Rby,this%b_img,itime)
      end if
      call update_trsfm(this%size)

    #endif

  end subroutine move_box

  subroutine read_init_cnf(this,p,irun)

    use :: conv_mod, only: QtoR

    class(box),intent(inout) :: this
    integer,intent(in) :: p,irun

    call this%BoxIO%read_init(p,irun,nchain,MPI_REAL_WP)

  end subroutine read_init_cnf

  subroutine read_dmp_cnf(this,p,irun,idmp,ndmp)

    use :: conv_mod, only: QtoR

    class(box),intent(inout) :: this
    integer,intent(in) :: p,irun,idmp,ndmp

    call this%BoxIO%read_dmp(p,irun,idmp,this%Qdagger_tilde,this%rcm_tilde,this%cmif,&
                             nchain,ntotsegx3,ndmp,MPI_REAL_WP)
    call QtoR(this%Qdagger_tilde,this%R_tilde,ntotsegx3,ntotbeadx3)

  end subroutine read_dmp_cnf

  subroutine write_cnf(this,id,p,itime,ntime,irun,idmp,time,Wi,dt,ndmp)

    use :: flow_mod, only: FlowType
    use :: conv_mod, only: RbctoQ,QtoR
    use :: trsfm_mod, only: bsx,bsy,invbsx,invbsy,eps_m,sinth,costh,tanb

    class(box),intent(inout) :: this
    integer,intent(in) :: id,p,itime,ntime,irun,idmp,ndmp
    real(wp),intent(in) :: time,Wi,dt
    integer :: igb

    
    select case (FlowType)
      case ('Equil')
        call RbctoQ(this%Rbx,this%Rby,this%Rbz,this%Qdagger_tilde,this%size,&
                                 this%invsize,nseg,nbead,ntotseg)
      case ('PSF')
        !$omp parallel default(private) shared(this,ntotbead,eps_m)
        !$omp do simd
        do igb=1, ntotbead
          this%Boxtrsfm%Rbtrx(igb)=this%Rbx(igb)-eps_m*this%Rby(igb)
        end do
        ! !$omp end do simd
        !$omp end parallel
        call RbctoQ(this%Boxtrsfm%Rbtrx,this%Rby,this%Rbz,this%Qdagger_tilde,&
                             this%size,this%invsize,nseg,nbead,ntotseg)
      case ('PEF')
        !$omp parallel default(private) shared(this,ntotbead,sinth,costh,tanb)
        !$omp do simd
        do igb=1, ntotbead
          this%Boxtrsfm%Rbtry(igb)=-sinth*this%Rbx(igb)+costh*this%Rby(igb)
          this%Boxtrsfm%Rbtrx(igb)= costh*this%Rbx(igb)+sinth*this%Rby(igb)-tanb*this%Boxtrsfm%Rbtry(igb)
        end do
        ! !$omp end do simd
        !$omp end parallel
        call RbctoQ(this%Boxtrsfm%Rbtrx,this%Boxtrsfm%Rbtry,this%Rbz,this%Qdagger_tilde,&
                            [bsx,bsy,this%size(3)],[invbsx,invbsy,this%invsize(3)],nseg,&
                            nbead,ntotseg)
    end select
    call QtoR(this%Qdagger_tilde,this%R_tilde,ntotsegx3,ntotbeadx3)

    call this%BoxIO%write(id,p,itime,ntime,irun,idmp,time,Wi,dt,nchain,nbead,nsegx3,nbeadx3,ntotsegx3,&
                       ntotbeadx3,ndmp,lambda,MPI_REAL_WP,this%Qdagger_tilde,this%rcm_tilde,this%cmif,&
                       this%Rb_tilde,this%R_tilde)

  end subroutine write_cnf

  subroutine calc_mat_fcn(this,id,irun,itime,time,Wi,Pe,dt,tgap,ntime,tss,trst,nrun,nprun,tend)

    use :: pp_smdlt, only: material_func

    class(box),intent(inout) :: this
    integer,intent(in) :: id,irun,itime,tgap,ntime,nrun,nprun
    real(wp),intent(in) :: time,Wi,Pe,dt,tss,trst,tend

    call material_func(id,irun,itime,time,Wi,Pe,dt,tgap,ntime,nchain,nseg,nbead,nsegx3,nbeadx3,&
              lambda,tss,trst,MPI_REAL_WP,nrun,nprun,tend,this%size,this%BoxChains,this%R_tilde)

  end subroutine calc_mat_fcn

  subroutine del_box(myrank)

    use :: conv_mod, only: del_conv
    use :: hi_mod, only: del_hi
    use :: force_smdlt, only: del_force

    integer,intent(in) :: myrank

    call del_conv()

    call del_force()

    call del_hi(myrank,Lbox)

  end subroutine del_box

  !> Box destructor
  subroutine del_box_t(this)

    use :: arry_mod, only: print_matrix,print_vector
    use :: flow_mod, only: FlowType
    use :: force_smdlt, only: del_force

    type(box) :: this
    integer :: ichain

    select case (FlowType)
      case ('Equil')
        if (CoMDiff) deallocate(this%cmif)
      case ('PSF')
        deallocate(this%Rbtr)
        deallocate(this%rcmtr)
      case ('PEF')
        deallocate(this%Rbtr)
        deallocate(this%rcmtr)
    end select
    deallocate(this%BoxChains)
    deallocate(this%Rb_tilde)
    deallocate(this%Rbx,this%Rby,this%Rbz)
    deallocate(this%R_tilde)
    deallocate(this%rcm_tilde)!,this%Fseg_tilde)
    deallocate(this%Qdagger_tilde)!,Fbead_tilde)
    deallocate(this%b_img,this%cm_img)

    #ifdef USE_GPU
      deallocate(this%Rb_d)
      deallocate(this%Rbx_d,this%Rby_d,this%Rbz_d)
      deallocate(this%Q_d)
      deallocate(this%rcm_d)
    #endif
    
  end subroutine del_box_t

  !> Calculating HI 
  subroutine calcHI(this,itime,dt)

    use :: hi_mod, only: HI_M,Mstart,restartHIpar,updateHIpar,hstar
    use :: tmng_mod, only: et_DT,et_DEC,HIcount,doTiming,tick,tock
    use :: arry_mod, only: print_vector,print_matrix
    #ifdef USE_GPU
      use :: diffcalc_mod, only: calcDiffTens_cpu,calcDiffTens_dev
      use :: diffdcmp_mod, only: calcBrownNoise_cpu,calcBrownNoise_dev
    #else
      use :: diffcalc_mod, only: calcDiffTens_cpu
      use :: diffdcmp_mod, only: calcBrownNoise_cpu
    #endif


    #ifdef USE_GPU
      use :: diffcalc_cumod, only: calcDiffTens_d
      use :: diffdcmp_cumod, only: calcBrownNoise_d
    #endif
    use :: mpi
    !include 'mpif.h'


    #ifdef USE_GPU
    use :: rndm_cumod, only: dw_bl_d
    #endif
    use :: hi_mod, only: dw_bltmp

    class(box),intent(inout) :: this
    integer,intent(in) :: itime
    real(wp),intent(in) :: dt
    integer :: offsetch,offsetbead,offsetseg,ibead,iseg,ichain
    integer(long) :: count0
    real(wp) :: t1,t2,rcm(3),R(3),rv(3)
    ! After 1/10th of lambda:
    if ((mod(itime,ceiling(lambda/(dt*100))) == 0) .and. &
        (HI_M /= Mstart) .and. (hstar /= 0._wp))   then
      call restartHIpar(this%ID,this%size,ntotbead)
    end if
    HIlp: do
      HIcount=HIcount+1
      



      if (doTiming) call tick(count0)


      ! #ifdef USE_GPU
      !   call calcDiffTens_d(this%Boxhi_d,this%Rb_d,ntotbead,this%size,this%origin)
      ! #endif

      #ifdef USE_GPU

        if (hstar /= 0._wp) then
          ! call calcDiffTens_dev(this%Boxhi_d,this%Rb_tilde,ntotbead,this%size,this%origin,itime)
          call calcDiffTens_d(this%Boxhi_d,this%Rb_d,ntotbead,this%size,this%origin)
        endif
        if (doTiming) then
          et_DT=et_DT+tock(count0)
          call tick(count0)
        end if
        dw_bltmp=dw_bl_d
        ! remember that the following routine is where PME_d is called for the first time
        call calcBrownNoise_d(this%Boxhi_d,this%decompRes,itime,ntotbeadx3,this%size)
      #else

        if (hstar /= 0._wp) then
          ! #ifdef USE_GPU
          !   call calcDiffTens_dev(this%Boxhi_d,this%Rb_tilde,ntotbead,this%size,this%origin,itime)
          ! #else
            call calcDiffTens_cpu(this%Rb_tilde,ntotbead,this%size,this%origin,itime)
          ! #endif
        endif
        if (doTiming) then
          et_DT=et_DT+tock(count0)
          call tick(count0)
        end if
        ! #ifdef USE_GPU
          ! call calcBrownNoise_dev(this%Boxhi_d,this%decompRes,itime,ntotbeadx3,this%size)
        ! #else

          call calcBrownNoise_cpu(this%decompRes,itime,ntotbeadx3,this%size)
        ! #endif


      #endif

        ! print*,'nnzb_h',this%Boxhi_d%nnzb
        ! print*,'col_h',this%Boxhi_d%P_ColInd_h
        ! print*,'val_h',this%Boxhi_d%P_Val_h
        ! print*,'row',this%Boxhi_d%P_RowPtr_h


      
      if (doTiming) et_DEC=et_DEC+tock(count0)
      if (this%decompRes%Success) then
        exit HIlp
      else
        call updateHIpar(this%ID,this%size,ntotbead)
        this%decompRes%Success=.true.
        cycle HIlp
      end if
    end do HIlp

  end subroutine calcHI

  !> Calculating the force on the beads 
  subroutine calcForce(this,itime)

    use :: arry_mod, only: print_vector
    use :: force_smdlt, only: Fphi,rFphi
    use :: evforce_mod, only: EVForceLaw
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: bsx,bsy,invbsx,invbsy,eps_m,reArng
    #ifdef USE_GPU
      use :: force_cumod, only: Fphi_d
    #endif

    class(box),intent(inout) :: this
    integer,intent(in) :: itime

    Fphi=0._wp
    rFphi=0._wp

    select case (FlowType)
      case ('Equil')
        call this%Boxsprf%update(this%Rbx,this%Rby,this%Rbz,this%size,&
                                 this%invsize,itime,nchain,nseg,nbead,&
                                 ntotseg,ntotsegx3,ntotbeadx3)
        if (EVForceLaw /= 'NoEV') then
          call this%Boxevf%update(this%Rbx,this%Rby,this%Rbz,this%size,&
                                  this%invsize,itime,nchain,nseg,nbead,&
                                  ntotseg,ntotsegx3,ntotbeadx3)
        end if
      case ('PSF')
        call this%Boxsprf%update(this%Boxtrsfm%Rbtrx,this%Rby,this%Rbz,&
                             this%size,this%invsize,itime,nchain,nseg,&
                             nbead,ntotseg,ntotsegx3,ntotbeadx3)
        if (EVForceLaw /= 'NoEV') then
          call this%Boxevf%update(this%Boxtrsfm%Rbtrx,this%Rby,this%Rbz,&
                               this%size,this%invsize,itime,nchain,nseg,&
                               nbead,ntotseg,ntotsegx3,ntotbeadx3)
        end if
      case ('PEF')
        call this%Boxsprf%update(this%Boxtrsfm%Rbtrx,this%Boxtrsfm%Rbtry,this%Rbz,&
                           [bsx,bsy,this%size(3)],[invbsx,invbsy,this%invsize(3)],&
                           itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbeadx3)
        if (EVForceLaw /= 'NoEV') then
          call this%Boxevf%update(this%Boxtrsfm%Rbtrx,this%Boxtrsfm%Rbtry,this%Rbz,&
                            [bsx,bsy,this%size(3)],[invbsx,invbsy,this%invsize(3)],&
                            itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbeadx3)
        end if
    end select

    #ifdef USE_GPU

      Fphi_d=0._wp

      call this%Boxsprf_d%update(this%Rbx_d,this%Rby_d,this%Rbz_d,this%size,&
        this%invsize,itime,nchain,nseg,nbead,ntotseg,ntotsegx3,ntotbeadx3)
    #endif

  end subroutine calcForce

end module box_mod
