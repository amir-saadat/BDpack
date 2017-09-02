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
module dlt_mod

  use :: mpi
  use :: iso_fortran_env
  use :: inp_dlt
  use :: prcn_mod
  use :: cmn_io_mod, only: read_input
  use :: arry_mod, only: logspace,linspace,print_vector,print_matrix,sort
  use :: force_mod, only: sprforce,sprupdate,bndforce,bndupdate,tetforce,tetupdate
  use :: dcmp_mod, only: Lanczos,BlockLanczos,MKLsyevr,BlockChebyshev
  use :: pp_mod, only: pp_init,pp_init_tm,data_prcs,conf_sort,del_pp
  use :: intrn_mod, only: intrn_t,wall_rflc,print_wcll

  implicit none

  contains

  subroutine dlt_bs(p,id)

    integer,intent(in) :: p,id
    ! MPI variables
    integer :: ierr,SizeSP,SizeDP,MPI_REAL_WP
    real(wp) :: realvar
    integer  :: intvar
    ! Variables used for making output
    integer :: q_recvsubarray,q_resizedrecvsubarray
    integer :: R_recvsubarray,R_resizedrecvsubarray
    integer :: rc_recvsubarray,rc_resizedrecvsubarray
    integer :: ia_recvsubarray,ia_resizedrecvsubarray
    integer :: cnf_tp_recvsubarray,cnf_tp_resizedrecvsubarray
    integer,dimension(3) :: q_starts,q_sizes,q_subsizes
    integer,dimension(3) :: R_starts,R_sizes,R_subsizes
    integer,dimension(3) :: rc_starts,rc_sizes,rc_subsizes
    integer,dimension(2) :: ia_starts,ia_sizes,ia_subsizes
    integer,dimension(2) :: cnf_tp_starts,cnf_tp_sizes,cnf_tp_subsizes
    integer,allocatable  :: q_counts(:),q_disps(:)
    integer,allocatable  :: R_counts(:),R_disps(:)
    integer,allocatable  :: rc_counts(:),rc_disps(:)
    integer,allocatable  :: ia_counts(:),ia_disps(:)
    integer,allocatable  :: cnf_tp_counts(:),cnf_tp_disps(:)
    integer(kind=MPI_ADDRESS_KIND) :: q_start,q_extent
    integer(kind=MPI_ADDRESS_KIND) :: R_start,R_extent
    integer(kind=MPI_ADDRESS_KIND) :: rc_start,rc_extent
    integer(kind=MPI_ADDRESS_KIND) :: ia_start,ia_extent
    integer(kind=MPI_ADDRESS_KIND) :: cnf_tp_start,cnf_tp_extent
    real(wp),parameter :: PI=4*atan(1._wp)
    real(wp),parameter :: sqrtPI=sqrt(PI)
    ! Parameters used for calculating Brownian displacement
    real(wp),parameter :: c1=14.14858378_wp,c2=1.21569221_wp
    ! Constant paramter used in SDE
    real(wp),parameter :: coeff=1/sqrt(2._wp)
    ! Characters:
    character(len=1024) :: file1,format_str,fstat,fmt1
    ! Integers:
    integer :: offset,offseti,offsetj,idx,is,os
    integer :: jcheck,kcheck,iseed,jp,iseg,jbead,i,j,k,nseg_bb,iarm,ku,kl
    integer :: nu,mu,ibead,jseg,ichain,ip,iread,iPe,idt,iAdjSeq,itime
    integer :: icol,jcol,kcol,lcol,jchain,istr,info,mrestart,Lrestart
    integer :: icount,iglob,jglob,nbead_bb,tgap_dmp
    ! Reals:
    real(wp) :: xxkappa,xykappa,yykappa,zzkappa,strain,fctr,lambdaminFixman
    real(wp) :: lambdamaxFixman,time,time_check1,time_check2,time_check3
    real(wp) :: time_check4,sqrtdt,rtpassed,wx,wy,wz,eps
    !    real(double) :: tA0,tA1
    real(single) :: time_begin,time_end
    ! Logicals:
    logical :: newSeq
    !    logical :: newSeq,EVcalcd
    ! Non-allocatable flow arrays
    integer,dimension(3) :: ipiv
    real(wp),dimension(2) :: lambdaBE
    real(wp),dimension(3) :: SumBdotw,SumDdotF,LdotBdotw,Ftet,Fbartet,rf_in
    real(wp),dimension(3,3) :: kappareg,totMobilTens,invtotMobilTens
    ! Allocatable arrays:
    integer,allocatable,dimension(:) :: mch,Lch
    real(wp),allocatable,dimension(:) :: qstar,qbar,Fbead,Fbarseg
    real(wp),allocatable,dimension(:) :: Fbarbead,Fev,w,FBr,Kdotq,ADFev,RHScnt
    real(wp),allocatable,dimension(:) :: Fbarev,root_f,Ybar,qctemp,uminus,uplus
    real(wp),allocatable,dimension(:) :: Ddotuminus,Ddotuplus
    real(wp),allocatable,dimension(:) :: Fbnd,Fbarbnd,Fbar,divD
    real(wp),allocatable,dimension(:),target :: RHS,RHSbase,qc,Fseg,DdotF
    real(wp),dimension(:),pointer  :: RHSP,RHSbaseP,qcP,FsegP,FBrblP,rcmP,rcmPy
    real(wp),dimension(:),pointer :: rchrP,FphiP,DdotFPx,DdotFPy,DdotFPz,qP
    real(wp),dimension(:),pointer :: rvmrcP,rvmrcPx,rvmrcPy,rvmrcPz,qPx,qPy,qPz
    real(double),dimension(:),pointer :: wbltempP2,BdotwP,BdotwPx,BdotwPy,BdotwPz
    real(wp),allocatable,dimension(:,:) :: Kappa,Amat,Bmat,wbl,aBlLan,WBlLan
    real(wp),allocatable,dimension(:,:) :: VBlLan,VcntBlLan,Eye,Dsh
    real(wp),allocatable,dimension(:,:) :: KappaBF,AmatBF,WeightTenstmp,rf0
    real(wp),allocatable,dimension(:,:),target :: q,rcm,rchr,Fphi,MobilTens,rvmrc
    real(wp),allocatable,dimension(:,:),target :: WeightTens,qstart,rcmstart
    real(wp),allocatable,dimension(:,:),target :: rchrstart
    real(wp),dimension(:,:),pointer :: AdotDP2,MobilTensP1,MobilTensP2
    real(wp),dimension(:,:),pointer :: DiffTensP,AdotDP1
    real(double),dimension(:,:),pointer :: wbltempP1
    real(double),dimension(:,:),pointer :: CoeffTensP
    real(wp),allocatable,dimension(:,:,:) :: rdnt,rdn
    ! For gathering data by rank 0
    real(wp),allocatable,dimension(:,:,:),target :: qxT,qyT,qzT,qTot,RTot,rcmT,FphiTot
    real(wp),allocatable,dimension(:,:,:),target :: rchrT,FBrbl,DiffTens,AdotD
    real(double),allocatable,dimension(:,:,:),target :: CoeffTens,wbltemp
    integer,allocatable,dimension(:) :: cnf_tp
    integer,allocatable,dimension(:,:) :: iaTot,cnf_tpTot
    ! File units
    integer :: u2,u3,u4,u21,u22,u23,u24,u25,u26,u27,u34,u39,u40,u41,u42
    ! objects
    type(intrn_t) :: myintrn


    call prcs_inp(id,p)

    ! Foramats used:
    1   format(3(f20.8,2x))
    2   format(i4,1x,f8.2,1x,e11.3,1x,f14.7,1x,3(i4,1x),f14.5)

    call pp_init(id)

    ! Controlling index (0 or -1):
    jcheck=0;kcheck=0

    ! Set MPI working precision - WP
    call MPI_Type_size(MPI_REAL,SizeSP,ierr)
    call MPI_Type_size(MPI_DOUBLE_PRECISION,SizeDP,ierr)
    if (wp == SizeSP) then
      MPI_REAL_WP=MPI_REAL
    else if (wp == SizeDP) then
      MPI_REAL_WP=MPI_DOUBLE_PRECISION
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    !--------------------------------------------------------------
    !>>>>> Initializations and allocations in master process:
    !--------------------------------------------------------------

    ! Specifying parameters by rank 0
    if (id == 0) then
      ! Initialization of random number generator in rank 0
      ! Choice of iseed: 0 <= iseed <= 2000000000 (2E+9);
      iseed=657483726
      call ranils(iseed)
      write (*,*)
      write (*,*) "%------------------------------------------------------------%"
      write (*,*) "| ***Start of BDpack program to perform Brownian dynamics*** |"
      write (*,*) "|         simulation for infinitely dilute solution          |"
      write (*,*) "%------------------------------------------------------------%"
      ! For making output
      ! q
      q_sizes(1)=nsegx3;q_sizes(2)=npchain;q_sizes(3)=p
      q_subsizes(1)=nsegx3;q_subsizes(2)=npchain;q_subsizes(3)=1
      q_starts(1)=0;q_starts(2)=0;q_starts(3)=0
      call MPI_Type_create_subarray(3,q_sizes,q_subsizes,q_starts,MPI_ORDER_&
        &FORTRAN,MPI_REAL_WP,q_recvsubarray,ierr)
      call MPI_Type_commit(q_recvsubarray,ierr)
      q_extent=sizeof(realvar)
      q_start=0
      call MPI_Type_create_resized(q_recvsubarray,q_start,q_extent,q_resized&
       &recvsubarray,ierr)
      call MPI_Type_commit(q_resizedrecvsubarray,ierr)
      ! rvmrc or R
      R_sizes=[nbeadx3,npchain,p]
      R_subsizes=[nbeadx3,npchain,1]
      R_starts=0
      call MPI_Type_create_subarray(3,R_sizes,R_subsizes,R_starts,MPI_ORDER_&
        &FORTRAN,MPI_REAL_WP,R_recvsubarray,ierr)
      call MPI_Type_commit(R_recvsubarray,ierr)
      R_extent=sizeof(realvar)
      R_start=0
      call MPI_Type_create_resized(R_recvsubarray,R_start,R_extent,R_resized&
       &recvsubarray,ierr)
      call MPI_Type_commit(R_resizedrecvsubarray,ierr)
      ! rcm and rchr
      if ((CoM) .or. (CoHR)) then
        rc_sizes(1)=3;rc_sizes(2)=npchain;rc_sizes(3)=p
        rc_subsizes(1)=3;rc_subsizes(2)=npchain;rc_subsizes(3)=1
        rc_starts(1)=0;rc_starts(2)=0;rc_starts(3)=0
        call MPI_Type_create_subarray(3,rc_sizes,rc_subsizes,rc_starts,MPI_OR&
          &DER_FORTRAN,MPI_REAL_WP,rc_recvsubarra&
          &y,ierr)
        call MPI_Type_commit(rc_recvsubarray,ierr)
        rc_extent=sizeof(realvar)
        rc_start=0
        call MPI_Type_create_resized(rc_recvsubarray,rc_start,rc_extent,rc_re&
         &sizedrecvsubarray,ierr)
        call MPI_Type_commit(rc_resizedrecvsubarray,ierr)
      end if
      ! Ia
      if ((tplgy == 'Comb').and.(arm_plc /= 'Fixed')) then
        ia_sizes=[Na,p]
        ia_subsizes=[Na,1]
        ia_starts=0
        call MPI_Type_create_subarray(2,ia_sizes,ia_subsizes,ia_starts,MPI_ORDER_&
          &FORTRAN,MPI_INTEGER,ia_recvsubarray,ierr)
        call MPI_Type_commit(ia_recvsubarray,ierr)
        ia_extent=sizeof(intvar)
        ia_start=0
        call MPI_Type_create_resized(ia_recvsubarray,ia_start,ia_extent,ia_resized&
         &recvsubarray,ierr)
        call MPI_Type_commit(ia_resizedrecvsubarray,ierr)
        allocate(iaTot(Na,p))
      end if
      ! cnf_tp
      if (cnf_srt) then
        cnf_tp_sizes=[npchain,p]
        cnf_tp_subsizes=[npchain,1]
        cnf_tp_starts=0
        call MPI_Type_create_subarray(2,cnf_tp_sizes,cnf_tp_subsizes,cnf_tp_starts,&
          MPI_ORDER_FORTRAN,MPI_INTEGER,cnf_tp_recvsubarray,ierr)
        call MPI_Type_commit(cnf_tp_recvsubarray,ierr)
        cnf_tp_extent=sizeof(intvar)
        cnf_tp_start=0
        call MPI_Type_create_resized(cnf_tp_recvsubarray,cnf_tp_start,cnf_tp_extent,&
         cnf_tp_resizedrecvsubarray,ierr)
        call MPI_Type_commit(cnf_tp_resizedrecvsubarray,ierr)
        allocate(cnf_tpTot(npchain,p))
        if (initmode == 'rst') then
          open(newunit=u41,file='data/cnf_tp.dat',status='unknown',position='append')
        else
          open(newunit=u41,file='data/cnf_tp.dat',status='replace',position='append')
        end if
      end if
      ! ia_time, inter-arrival time
      ! if (EV_bw == 'Rflc_bc') then
      !   ia_tm_sizes=[nseg,ntime,nchain,p]
      !   ia_tm_subsizes=[nseg,ntime,nchain,1]
      !   ia_tm_starts=0
      !   call MPI_Type_create_subarray(2,ia_tm_sizes,cnf_subsizes,cnf_tp_starts,&
      !                                 MPI_ORDER_FORTRAN,MPI_INTEGER,cnf_tp_recvsuba&
      !                                 &rray,ierr)
      !   call MPI_Type_commit(cnf_tp_recvsubarray,ierr)
      !   cnf_tp_extent=sizeof(intvar)
      !   cnf_tp_start=0
      !   call MPI_Type_create_resized(cnf_tp_recvsubarray,cnf_tp_start,cnf_tp_extent,&
      !                                cnf_tp_resizedrecvsubarray,ierr)
      !   call MPI_Type_commit(cnf_tp_resizedrecvsubarray,ierr)
      !   allocate(cnf_tpTot(npchain,p))
      !   if (initmode == 'rst') then
      !     open(newunit=u41,file='data/cnf_tp.dat',status='unknown',position='append')
      !   else
      !     open(newunit=u41,file='data/cnf_tp.dat',status='replace',position='append')
      !   end if
      ! end if
      ! allocation of total random number
      allocate(rdnt(nbeadx3,ncols,nchain))
      allocate(qTot(nsegx3,npchain,p))
      allocate(qxT(npchain,nseg,p))
      allocate(qyT(npchain,nseg,p))
      allocate(qzT(npchain,nseg,p))
      allocate(RTot(nbeadx3,npchain,p))
      allocate(FphiTot(nbeadx3,npchain,p))
      if (CoM) allocate(rcmT(3,npchain,p))
      if (CoHR) allocate(rchrT(3,npchain,p))
    end if ! id.eq.0

    !----------------------------------------------------------
    !>>>>> Allocating arrays:
    !----------------------------------------------------------

    ! Allocating the local random arrays
    allocate(rdn(nbeadx3,ncols,npchain))
    ! Note: the order nseg,nchain is selected because Fortran is column major.
    allocate (DiffTens(nbeadx3,nbeadx3,npchain),CoeffTens(nbeadx3,nbeadx3,npchain))
    allocate (qc(nsegx3),qstar(nsegx3),Fseg(nsegx3),wbl(nbeadx3,ncols))
    allocate (wbltemp(nbeadx3,ncols,npchain),w(nbeadx3),Kappa(nsegx3,nsegx3))
    allocate (Amat(nsegx3,nbeadx3),Kdotq(nsegx3),FBr(nsegx3))
    allocate (FBrbl(nsegx3,ncols,npchain),AdotD(nsegx3,nbeadx3,npchain))
    allocate (ADFev(nsegx3),Fev(nbeadx3),Fbead(nbeadx3),RHS(nsegx3),RHScnt(nsegx3))
    allocate (RHSbase(nsegx3),Fbarseg(nsegx3),qbar(nsegx3),Fbarbead(nbeadx3))
    allocate (q(nsegx3,npchain),qstart(nsegx3,npchain),Fphi(nbeadx3,npchain))
    allocate (Bmat(nbeadx3,nsegx3),Fbarev(nbeadx3),KappaBF(2,nsegx3),Fbar(nbeadx3))
    allocate (qctemp(nsegx3),DdotF(nbeadx3),Fbnd(nbeadx3),Fbarbnd(nbeadx3))
    select case (tplgy)
    case ('Linear')
      allocate(AmatBF(4,nbeadx3))
    case ('Comb')
    end select
    allocate(rvmrc(nbeadx3,npchain))
    if ((hstar.eq.0.0_wp).or.(DecompMeth.eq.'Chebyshev')) allocate(Eye(nbeadx3,nbeadx3))
    if (DecompMeth == 'Chebyshev') then
      allocate(Dsh(nbeadx3,nbeadx3),Ybar(nbeadx3),uminus(nbeadx3),uplus(nbeadx3))
      allocate(Ddotuminus(nbeadx3),Ddotuplus(nbeadx3),Lch(npchain))
      ! For calculating the average of iteration number
      Lch(:)=LCheb
    elseif (DecompMeth == 'Lanczos') then
      allocate(aBlLan(nbeadx3,ncols),WBlLan(nbeadx3,ncols),VBlLan(nbeadx3,mBlLan*ncols))
      allocate(Ybar(nbeadx3),VcntBlLan(nbeadx3,ncols),mch(npchain))
      ! For calculating the average of iteration number
      mch(:)=mBlLan
    end if
    if ((hstar /= 0._WP) .and. (HITens == 'Blake')) allocate(divD(nbead))
    ! allocate(divD(nbead))
    if (CoM) allocate(rcm(3,npchain),rcmstart(3,npchain))
    if (CoHR) allocate(rchr(3,npchain),rchrstart(3,npchain),&
     MobilTens(nbeadx3,nbeadx3),WeightTens&
     &(3,nbeadx3),weightTenstmp(3,nbeadx3))
    if (srf_tet) then
      allocate(rf0(3,npchain))
      call read_input('rf-in',0,rf_in(1),def=0._wp)
      call read_input('rf-in',1,rf_in(2),def=0._wp)
      call read_input('rf-in',2,rf_in(3),def=0._wp)
    endif
    ! For making the output
    allocate (q_counts(p),q_disps(p))
    allocate (R_counts(p),R_disps(p))
    if (CoM .or. CoHR) then
      allocate(rc_counts(p),rc_disps(p))
    end if
    if ((tplgy == 'Comb').and.(arm_plc /= 'Fixed')) then
      allocate(ia_counts(p),ia_disps(p))
    end if
    if (cnf_srt) then
      allocate(cnf_tp(npchain),cnf_tp_counts(p),cnf_tp_disps(p))
    end if
    do jp=1, p
      q_counts(jp)=1
      q_disps(jp)=(jp-1)*nsegx3*npchain
      R_counts(jp)=1
      R_disps(jp)=(jp-1)*nbeadx3*npchain
      if (CoM .or. CoHR) then
        rc_counts(jp)=1
        rc_disps(jp)=(jp-1)*3*npchain
      end if
      if ((tplgy == 'Comb').and.(arm_plc /= 'Fixed')) then
        ia_counts(jp)=1
        ia_disps(jp)=(jp-1)*Na
      end if
      if (cnf_srt) then
        cnf_tp_counts(jp)=1
        cnf_tp_disps(jp)=(jp-1)*npchain
      end if
    end do
    if ((tplgy == 'Comb').and.(arm_plc /= 'Fixed')) then
      write(fmt1,'("(",i3.3,"(i4,1x))")') Na
      if ((iflow == 1).and.(initmode == 'st').and.(arm_plc /= 'Stored')) then
        call MPI_Gatherv(Ia(2:Na+1),Na,MPI_INTEGER,iaTot,ia_counts,ia_disps,&
         ia_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
        if (id == 0) then
          open(newunit=u39,file='data/ia.dat',status='replace')
          do ip=1, p
            write(u39,fmt1) iaTot(:,ip)
          end do
        end if ! id==0
      else
        Ia(1)=1
        do ip=1, p
          if (ip == 1) then
            open(newunit=u39,file='data/ia.dat',status='old',position='rewind')
            do iread=1, id
              read(u39,*)
            end do
          end if
          if (id == ip-1) then
            read(u39,*) Ia(2:Na+1)
          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
        end do
      end if
      close(u39)
    end if
    call sort(Ia)

    !----------------------------------------------------------------
    !>>>>> Constant tensors in SDE:
    !----------------------------------------------------------------

    ! Specifying kappa based on type of flow
    if (iflow == 1) then ! Finding Equilibrium
      xxkappa=0._wp
      xykappa=0._wp
      yykappa=0._wp
      zzkappa=0._wp
    else if (iflow == 2) then ! Shear Flow
      xxkappa=0._wp
      xykappa=1._wp
      yykappa=0._wp
      zzkappa=0._wp
    else if (iflow == 3) then ! Uniaxial Extension
      xxkappa=1._wp
      xykappa=0._wp
      yykappa=-0.5_wp
      zzkappa=-0.5_wp
    else if (iflow == 4) then ! biaxial Extension
      xxkappa=1._wp
      xykappa=0._wp
      yykappa=1._wp
      zzkappa=-2.0_wp
    else if (iflow == 5) then ! Planar Extension
      xxkappa=1._wp
      xykappa=0._wp
      yykappa=-1._wp
      zzkappa=0._wp
    end if

    ! To be used in Predictor-Corrector step
    Kappa=0._wp
    forall (iseg=1:3*(nseg-1)+1:3)
      Kappa(iseg,iseg)=xxkappa
      Kappa(iseg,iseg+1)=xykappa
      Kappa(iseg+1,iseg+1)=yykappa
      Kappa(iseg+2,iseg+2)=zzkappa
    end forall
    Kappareg(1,1:3)=(/xxkappa,xykappa,0.0_wp/)
    Kappareg(2,1:3)=(/0.0_wp ,yykappa,0.0_wp/)
    Kappareg(3,1:3)=(/0.0_wp ,0.0_wp ,zzkappa/)
    ! Amat is Bbar and Bmat is B in "DPL" Bird et al.
    Amat=0.0_wp
    select case (tplgy)
    case ('Linear')
      nseg_bb=nseg
      nbead_bb=nbead
      do iseg=1, nseg
        offseti=3*(iseg-1)
        do jbead=1, nbead
          offsetj=3*(jbead-1)
          if (iseg == jbead) then
            forall (i=1:3) Amat(offseti+i,offsetj+i)=-1._wp
          elseif (iseg == jbead-1) then
            forall (i=1:3) Amat(offseti+i,offsetj+i)= 1._wp
          end if
        end do
      end do
    case ('Comb')
      nseg_bb=nseg-Na*nseg_ar
      nbead_bb=nseg_bb+1
      iarm=1
      do iseg=1, nseg
        offseti=3*(iseg-1)
        do jbead=1, nbead
          offsetj=3*(jbead-1)
          if (iseg <= nseg_bb) then
            if (jbead == iseg) then
              forall (i=1:3) Amat(offseti+i,offsetj+i)=-1._wp
            elseif (jbead == iseg+1) then
              forall (i=1:3) Amat(offseti+i,offsetj+i)= 1._wp
            end if
          else ! iseg > nseg_bb
            if (iseg-nseg_bb-(iarm-1)*nseg_ar == 1) then
              if (jbead == Ia(iarm+1)) then
                forall (i=1:3) Amat(offseti+i,offsetj+i)=-1._wp
              elseif (jbead == iseg+1) then
                forall (i=1:3) Amat(offseti+i,offsetj+i)= 1._wp
              elseif (jbead == nbead) then
                iarm=iarm+1
              end if
            else
              if (jbead == iseg) then
                forall (i=1:3) Amat(offseti+i,offsetj+i)=-1._wp
              elseif (jbead == iseg+1) then
                forall (i=1:3) Amat(offseti+i,offsetj+i)= 1._wp
              end if
            end if
          end if
        end do
      end do
    end select
    ! Constructing banded form of Kappa
    KappaBF=0._wp
    ku=1;kl=0
    do j=1, nsegx3
      k=ku+1-j
      do i=max(1,j-ku),min(nsegx3,j+kl)
        KappaBF(k+i,j)=Kappa(i,j)
      end do
    end do
    if (tplgy == 'Linear') then
      AmatBF=0._wp
      ! Constructing banded form of Amat
      ku=3;kl=0
      do j=1, nbeadx3
        k=ku+1-j
        do i=max(1,j-ku),min(nsegx3,j+kl)
          AmatBF(k+i,j)=Amat(i,j)
        end do
      end do
    end if
    Bmat=0.0_wp
    select case (tplgy)
    case ('Linear')
      do ibead=1, nbead
        offseti=3*(ibead-1)
        do jseg=1, nseg
          offsetj=3*(jseg-1)
          if (ibead > jseg) then
            forall (i=1:3)
              Bmat(offseti+i,offsetj+i)=jseg/real(nbead,kind=wp)
            end forall
          else
            forall (i=1:3)
              Bmat(offseti+i,offsetj+i)=-(1-jseg/real(nbead,kind=wp))
            end forall
          end if
        end do
      end do
    case ('Comb')
     ! Constructing the elements of the first row of B
     do k=1, nseg_bb
       forall (i=1:3) Bmat(i,3*(k-1)+i)=-(nseg_bb-k+1)/real(nbead,kind=wp)
     end do
     do iarm=1, Na
       fctr=(Na-iarm+1)*nseg_ar/real(nbead,kind=wp)
       do k=Ia(iarm), Ia(iarm+1)-1
         forall (i=1:3)
           Bmat(i,3*(k-1)+i)=Bmat(i,3*(k-1)+i)-fctr
         end forall
       end do ! k
       do k=1, nseg_ar
         idx=nseg_bb+(iarm-1)*nseg_ar+k
         forall (i=1:3)
           Bmat(i,3*(idx-1)+i)=Bmat(i,3*(idx-1)+i)-&
           (nseg_ar-k+1)/real(nbead,kind=wp)
         end forall
       end do ! k
     end do ! iarm
     ! Constructing the rest of the rows in backbone
     do nu=2, nseg_bb+1
       forall (i=1:3) Bmat(3*(nu-1)+i,:)=Bmat(i,:)
       do k=1, nu-1
         forall (i=1:3)
           Bmat(3*(nu-1)+i,3*(k-1)+i)=Bmat(3*(nu-1)+i,3*(k-1)+i)+1
         end forall
       end do ! k
     end do ! nu
     ! Constructing the rows for the arms
     do iarm=1, Na
       do mu=1, nseg_ar
         nu=nseg_bb+1+(iarm-1)*nseg_ar+mu
         forall (i=1:3) Bmat(3*(nu-1)+i,:)=Bmat(i,:)
         do k=1, Ia(iarm+1)-1
           forall (i=1:3)
             Bmat(3*(nu-1)+i,3*(k-1)+i)=Bmat(3*(nu-1)+i,3*(k-1)+i)+1
           end forall
         end do ! k
         do k=1, mu
           idx=nseg_bb+(iarm-1)*nseg_ar+k
           forall (i=1:3)
             Bmat(3*(nu-1)+i,3*(idx-1)+i)=Bmat(3*(nu-1)+i,3*(idx-1)+i)+1
           end forall
         end do ! k
       end do ! mu
     end do ! iarm
   end select
   if ((hstar == 0._wp) .or. (DecompMeth == 'Chebyshev')) then
    Eye=0._wp
    forall (i=1:nbeadx3) Eye(i,i)=1._wp
  end if
  if (hstar == 0._wp) then
    !      Eye=0._wp
    do ichain=1, npchain
      DiffTens(:,:,ichain)=Eye
      CoeffTens(:,:,ichain)=Eye
    end do
  end if
  if (EV_bb == 'NoEV') Fev=0._wp;Fbarev=0._wp
  if (ForceLaw /= 'WLC_GEN') Fbnd=0._wp;Fbarbnd=0._wp
  if (DecompMeth == 'Chebyshev') then
    !      Eye=0._wp
    forall (i=1:nbeadx3) uminus(i)=real((-1)**i,kind=wp)
    forall (i=1:nbeadx3) uplus(i)=1.0_wp
  end if


  !----------------------------------------------------------------
  !>>>>> Configuration initialization:
  !----------------------------------------------------------------

  if (iflow == 1) then
    if (initmode == 'st') then
      do ichain=1, npchain
        do iseg=1, nseg
          offset=3*(iseg-1)
          select case (tplgy)
          case ('Linear')
            qstart(offset+1:offset+3,ichain)=(/infrx*qmax,infry*qmax,infrz*qmax/)
            !                qstart(offset+1:offset+3,ichain)=[0.4_wp*qmax,0.1_wp*qmax,0.2_wp*qmax]
          case ('Comb')
            if (iseg <= nseg_bb) then
              qstart(offset+1:offset+3,ichain)=(/0.9_wp*qmax,0._wp,0._wp/)
            else
              if (mod((iseg-nseg_bb-1)/nseg_ar+1,2) == 0) then
                qstart(offset+1:offset+3,ichain)=(/0._wp,-0.9_wp*qmax,0._wp/)
              else
                qstart(offset+1:offset+3,ichain)=(/0._wp,0.9_wp*qmax,0._wp/)
              end if
            end if
          end select
        end do
      end do
      if (CoM) rcmstart=0._wp
      !        if (CoM) rcmstart(1:3,1)=[-4.517_wp,2.308_wp,-3.395_wp]
      !        if (CoM) rcmstart(1:3,2)=[250-0.034_wp,250-3.173_wp,250-1.454_wp]
      if (CoHR) rchrstart=0._wp
    elseif (initmode == 'rst') then
      ! In order to prevent probable race condition
      do ip=1, p
        if (ip == 1) then
          open(newunit=u2,file='data/q.rst.dat',status='old',position='rewind')
          if (CoM) then
            open(newunit=u3,file='data/CoM.rst.dat',status='old',position='rewind')
          end if
          if (CoHR) then
            open(newunit=u4,file='data/CoHR.rst.dat',status='old',position='rewind')
          end if
          do iread=1, id*nseg*npchain
            read(u2,*)
          end do
          do iread=1, id*npchain
            if (CoM) read(u3,*);if (CoHR) read(u4,*)
          end do
        end if
        if (id == ip-1) then
          do ichain=1, npchain
            do iseg=1, nseg
              offset=3*(iseg-1)
              read(u2,*) qstart(offset+1:offset+3,ichain)
            end do
            if (CoM) read(u3,*) rcmstart(1:3,ichain)
            if (CoHR) read(u4,*) rchrstart(1:3,ichain)
          end do
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end do
    end if ! initmode
  else ! iflow /= 1
    do ip=1, p
      if (ip == 1) then
        if (initmode == 'st') then
          open(newunit=u2,file='data/q.st.dat',status='old',position='rewind')
          if (CoM) then
            open(newunit=u3,file='data/CoM.st.dat',status='old',position='rewind')
          end if
          if (CoHR) then
            open(newunit=u4,file='data/CoHR.st.dat',status='old',position='rewind')
          end if
        elseif (initmode == 'rst') then
          open(newunit=u2,file='data/q.rst.dat',status='old',position='rewind')
          if (CoM) then
            open(newunit=u3,file='data/CoM.rst.dat',status='old',position='rewind')
          end if
          if (CoHR) then
            open(newunit=u4,file='data/CoHR.rst.dat',status='old',position='rewind')
          end if
        end if
        do iread=1, id*nseg*npchain
          read(u2,*)
        end do
        do iread=1, id*npchain
          if (CoM) read(u3,*);if (CoHR) read(u4,*)
        end do
      end if
      if (id == ip-1) then
        do ichain=1, npchain
          do iseg=1, nseg
            offset=3*(iseg-1)
            read(u2,*) qstart(offset+1:offset+3,ichain)
          end do
          if (CoM) read(u3,*) rcmstart(1:3,ichain)
          if (CoHR) read(u4,*) rchrstart(1:3,ichain)
        end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end do
  end if
  close (u2);if (CoM) close (u3);if (CoHR) close(u4)

  allocate(root_f(PrScale*nroots))

  !----------------------------------------------------------------
  !>>>>> HI and EV-bb initialization:
  !----------------------------------------------------------------

  call myintrn%init(id)

  ! call hi_init()
  ! call ev_init()
  ! call evbw_init(id)

  !----------------------------------------------------------------
  !>>>>> Time integration of SDE:
  !----------------------------------------------------------------

  ! Loop over Pe number
  do iPe=1, nWi
    ! Loop over dt
    do idt=1, ndt

      if (.not.Adjust_dt) call lookup_tab(dt(iPe,idt))

      if (id == 0) then
        write (*,*)
        write (*,*) "%--------------------------------------------------%"
        write (*,*) "| ***Start of time integration in all processes*** |"
        write (*,*) "%--------------------------------------------------%"
        write (*,*)
        write(*,'(7x,a)') 'Wi            dt            ntime'
        write(*,'(7x,a)') '---------------------------------'
        write(*,'(f14.7,1x,e10.2,1x,i10)') Wi(iPe),dt(iPe,idt),ntime(iPe,idt)
        write (*,*)
      end if

      if (TimerA) then
        !          if (id == 0) tA0=MPI_Wtime()
        if (id.eq.0) call cpu_time(time_begin)

      end if

      newSeq=.true.
      iAdjSeq=1
      time=0._wp
      time_check1=frm_rt_rep*lambda ! For reporting the time passed.
      time_check2=tss*lambda+frm_rt_dmp*lambda ! For making dump files.

      tgap_dmp=ceiling(frm_rt_dmp*lambda/dt(iPe,idt))


      time_check3=frm_rt_pp*lambda ! For post processing.
      time_check4=frm_rt_rst*lambda ! For making restart files.

      if (id == 0) then
        print *
        print '(" Number of iterations for dumping:",i10)',tgap_dmp
        print *
      end if

      do itime=1, ntime(iPe,idt)
        if (Adjust_dt) then
          if ( (itime == itime_AdjSeq(iPe,idt,iAdjSeq)) .and. &
           (iAdjSeq <= (nAdjSeq-1)) ) newSeq=.true.
          if (newSeq) then
            ! increment the sequence after the first adjustment itime.
            if (itime /= 1) iAdjSeq=iAdjSeq+1
            dt(iPe,idt)=dt_tmp(iPe,idt)*AdjFact(iAdjSeq)
            newSeq=.false.
            sqrtdt=sqrt(dt(iPe,idt))
            if (id == 0) then
              print *
              print '(1x,a,1x,i2)', 'Adjusting Sequence:',iAdjSeq
              print *, 'Note!!: Time step adjustment has changed as follows:'
              print *
              print '(7x,a)', 'Wi            Pe            dt          ntime'
              print '(7x,a)', '---------------------------------------------'
              print '(f14.7,1x,f14.7,1x,e10.2,1x,i10)', Wi(iPe),Pe(iPe),dt(iPe,idt),&
                  ntime(iPe,idt)
              print *
            end if
            call lookup_tab(dt(iPe,idt))
          end if ! newSeq
        end if ! Adjust_dt
        ! Calculating time passed based on time step:
        time=time+dt(iPe,idt)
        if (id == 0) then
          ! Constructing a block of random numbers for ncols time steps by rank 0
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
            do ichain=1, nchain
              do icol=1, ncols
                do ibead=1, 3*nbead
                  rdnt(ibead,icol,ichain)=ranuls()-0.5
                end do
              end do
            end do
          end if
        end if
        ! Scattering the generated random numbers from rank 0 the owner.
        if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
          call MPI_Scatter(rdnt,3*nbead*ncols*npchain,MPI_REAL_WP,rdn,&
           3*nbead*ncols*npchain,MPI_REAL_WP,0,MPI_COMM_WORLD,ierr)
          jcol=1
        end if
        if ((time >= time_check1) .and. (id == 0)) then
          if (initmode == 'rst') then
            rtpassed=(time+trst*lambda)/lambda
          else
            rtpassed=time/lambda
          end if
          print '(f7.3," Chain-Relaxation-Time(s) Passed")',rtpassed
          time_check1=time_check1+frm_rt_rep*lambda
        end if

        ! Things to do in the first time step:
        if (itime == 1) then
          if (.not.Adjust_dt) sqrtdt=sqrt(dt(iPe,idt))
          q=qstart
          do jchain=1, npchain
            qP => q(:,jchain)
            rvmrcP => rvmrc(:,jchain)
            call gemv(Bmat,qP,rvmrcP)
            if (srf_tet) then
              rf0(:,jchain)=rf_in(1:3)!rvmrc(1:3,jchain)+rcmstart(:,jchain)
            end if
          end do
          if (CoM) rcm=rcmstart
          if (CoHR) rchr=rchrstart
          call pp_init_tm(id)
          istr=1 ! the index for dumpstr()
        end if

        do ichain=1, npchain
          ! Important Note: each chain should construct its own block,
          ! if the vectorial components (i.e. a(:,ichain)) are going
          ! to be used in the rest of the program. These tensors are
          ! used for "ncols" time steps for EACH CHAIN:
          DiffTensP => DiffTens(:,:,ichain)
          CoeffTensP => CoeffTens(:,:,ichain)
          wbltempP1 => wbltemp(:,:,ichain)
          AdotDP1 => AdotD(:,:,ichain)
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
            ! Setting the random numbers for the brownian force
            do kcol=1, ncols
              ! In case of having overlap for the #cols and adjusting sequence:
              if ( Adjust_dt .and. (ncols > 1) ) then
                if ((itime+kcol) >= itime_AdjSeq(iPe,idt,iAdjSeq)) then
                  if ((iAdjSeq+1) > nAdjSeq) then
                    sqrtdt=sqrt(AdjFact(iAdjSeq)*dt_tmp(iPe,idt))
                  else
                    sqrtdt=sqrt(AdjFact(iAdjSeq+1)*dt_tmp(iPe,idt))
                  end if
                end if
              end if
              do ibead=1, nbead
                offset=3*(ibead-1)
                wx=rdn(offset+1,kcol,ichain);wx=sqrtdt*wx*(c1*wx**2+c2)
                wy=rdn(offset+2,kcol,ichain);wy=sqrtdt*wy*(c1*wy**2+c2)
                wz=rdn(offset+3,kcol,ichain);wz=sqrtdt*wz*(c1*wz**2+c2)
                wbl(offset+1,kcol)=wx;wbl(offset+2,kcol)=wy;wbl(offset+3,kcol)=wz
              end do
            end do
          end if
          w(:)=wbl(:,jcol)


          ! Constructing an array for each individual chain
          qc(:)=q(:,ichain)
          qP => q(:,ichain)
          rvmrcP => rvmrc(:,ichain)
          qPx => qP(1:nsegx3-2:3)
          qPy => qP(2:nsegx3-1:3)
          qPz => qP(3:nsegx3:3)
          rvmrcPx => rvmrcP(1:nbeadx3-2:3)
          rvmrcPy => rvmrcP(2:nbeadx3-1:3)
          rvmrcPz => rvmrcP(3:nbeadx3:3)
          if (CoM) then
            rcmP => rcm(:,ichain)
          end if


          ! if (EV_bw == 'Rflc_bc') then
          !   call wall_rflc(rvmrcPy,rcmP(2))
          ! end if


          call sprforce(id,qc,nseg,ForceLaw,TruncMethod,Fseg)


          if (ForceLaw == 'WLC_GEN') call bndforce(nbead_bb,qc,Fbnd,itime)
          ! Calculation of Diffusion Tensor and Excluded Volume Force
          !            EVcalcd=.false.
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
            if ((hstar /= 0._wp).or.(DecompMeth /= 'Cholesky')) then
              !                if (EV_bb /= 'NoEV') EVcalcd=.true.
              !                if (EV_bw /= 'NoEV') EVbwcalcd=.true.
              !call hicalc2(rvmrcP,nseg,DiffTensP,Fev)
              !call print_matrix(DiffTensP,'d1')
              !                  call HICalc(rvmrcP,nseg,HITens,DiffTensP,EV_bb,Fev)
              !call print_matrix(DiffTensP,'d2')

              call myintrn%calc(rvmrcP,rcmP,nseg,DiffTensP,divD,Fev,Fbarev,&
                calchi=.true.,calcdiv=.true.,calcevbb=.true.,calcevbw=.true.)

            end if
            if (DecompMeth == 'Cholesky') then
              if (hstar /= 0._wp) then
                CoeffTensP=real(DiffTensP,kind=double)
                wbltempP1=real(wbl,kind=double)

                !! I need to change the storage scheme for symmetric tensors
                !! for dilute solutions in the absence of the wall

                if (HITens == 'Blake') then
                  call potrf(CoeffTensP,info=info)
                else
                  call potrf(CoeffTensP,info=info)
                endif

                if (info /= 0) then
                  print '(" Unsuccessful Cholesky fact. of D in main.")'
                  print '(" info: ",i3)',info
                  stop
                end if


                if (HITens == 'Blake') then
                  call trmm(CoeffTensP,wbltempP1,transa='T')
                else
                  call trmm(CoeffTensP,wbltempP1,transa='T')
                endif

              else
                wbltempP1=real(wbl,kind=double)
              end if
            elseif (DecompMeth == 'Lanczos') then
              if ((mod(itime,upfactr*ncols) == 1) .or. (upfactr == 1)) then
                mrestart=mBlLan
                if (ncols == 1) then
                  call Lanczos(real(DiffTensP,kind=double),real(wbl,kind=double),real(WBlLan,kind=double),&
                    real(Ybar,kind=double),nbeadx3,real(errormin,kind=double),mubBlLan,mrestart,&
                    wbltempP1,msetinp=mset)
                else
                  call BlockLanczos(real(DiffTensP,kind=double),real(wbl,kind=double),real(aBlLan,kind=double),&
                    real(WBlLan,kind=double),real(Ybar,kind=double),nbeadx3,ncols,real(errormin,kind=double),&
                    mubBlLan,mrestart,wbltempP1,msetinp=mset)
                endif
                mch(ichain)=mrestart
              else
                if (ncols == 1) then
                  call Lanczos(real(DiffTensP,kind=double),real(wbl,kind=double),real(WBlLan,kind=double),&
                   real(Ybar,kind=double),nbeadx3,real(errormin,kind=double),mubBlLan,mch(ichain),&
                   wbltempP1,msetinp=mset)
                else
                  call BlockLanczos(real(DiffTensP,kind=double),real(wbl,kind=double),real(aBlLan,kind=double),&
                    real(WBlLan,kind=double),real(Ybar,kind=double),nbeadx3,ncols,real(errormin,kind=double),&
                    mubBlLan,mch(ichain),wbltempP1,msetinp=mset)
                end if
              end if
            elseif (DecompMeth == 'Chebyshev') then
              if ((mod(itime,upfactr*ncols) == 1) .or. (upfactr == 1)) then
                Lrestart=LCheb
                ! Calculation of dmin and dmax passed with lambdain to
                ! BlockChebyshev routine
                call symv(DiffTensP,uminus,Ddotuminus)
                call symv(DiffTensP,uplus,Ddotuplus)
                lambdaminFixman=dot(uminus,Ddotuminus)/nbeadx3
                lambdamaxFixman=dot(uplus,Ddotuplus)/nbeadx3
                lambdaBE=(/lambdaminFixman/2,2*lambdamaxFixman/)
                call BlockChebyshev(real(DiffTensP,kind=double),real(Eye,kind=double),real(Dsh,kind=double),&
                  real(wbl,kind=double),nbeadx3,ncols,Lub,Lrestart,wbltempP1,real(Ybar,kind=double),&
                  MKLsyevr,real(errormin,kind=double),lambdainp=real(lambdaBE,kind=double),Lsetinp=Lset)
                Lch(ichain)=Lrestart
              else
                ! Calculation of dmin and dmax passed with lambdain to
                ! BlockChebyshev routine
                call symv(DiffTensP,uminus,Ddotuminus)
                call symv(DiffTensP,uplus,Ddotuplus)
                lambdaminFixman=dot(uminus,Ddotuminus)/nbeadx3
                lambdamaxFixman=dot(uplus,Ddotuplus)/nbeadx3
                lambdaBE=(/lambdaminFixman/2,2*lambdamaxFixman/)
                call BlockChebyshev(real(DiffTensP,kind=double),real(Eye,kind=double),real(Dsh,kind=double),&
                 real(wbl,kind=double),nbeadx3,ncols,Lub,Lch(ichain),wbltempP1,real(Ybar,kind=double),&
                 MKLsyevr,real(errormin,kind=double),lambdainp=real(lambdaBE,kind=double),Lsetinp=Lset)
              end if
            end if
            do lcol=1, ncols
              wbltempP2 => wbltemp(:,lcol,ichain)
              FBrblP => FBrbl(:,lcol,ichain)
              if (tplgy == 'Linear') then
                call gbmv(AmatBF,real(wbltempP2,kind=wp),FBrblP,kl=0,m=nsegx3,alpha=coeff)

                !TYL: HI for tethered bead -------------------------------------
                !print *, 'Brownian force calculaton------------'
                !call print_matrix(Amat,'Amat')
                !call print_vector(real(wbltempP2,kind=wp),'W')
                !call print_vector(FBrblP, 'Br force BEFORE')
                if (srf_tet) then
                  call gemv(Amat(:,1:3),real(wbltempP2(1:3),kind=wp),FBrblP,alpha=-coeff,&
                    beta=1._wp)
                  !call print_vector(FBrblP, 'Br force AFTER')
                end if
                !TYL: HI for tethered bead -------------------------------------

              else
                call gemv(Amat,real(wbltempP2,kind=wp),FBrblP,alpha=coeff)
              end if
            end do
            ! Calculation of AdotD=Amat.D, to be used in Predictor-Corrector
            if (hstar /= 0._wp) then
              call symm(DiffTensP,Amat,AdotDP1,side='R')
            else
              AdotDP1=Amat
            end if
          end if
          FBr=FBrbl(:,jcol,ichain)


          ! rflc doesn't need this
          !            if (hstar == 0._wp .and. EV_bw /= 'Rflc_bc') then
          if (hstar == 0._wp) then
            !                call EVCalc(rvmrcP,nseg,EV_bb,Fev)
            call myintrn%calc(rvmrcP,rcmP,nseg,DiffTensP,divD,Fev,Fbarev,&
              calcevbb=.true.,calcevbw=.true.)
            !call evcalc2(rvmrcP,nseg,Fev)
          end if


          !============ Predictor-Corrector =============!
          !--------Predictor Algorithm----------!
          ! Kdotq=dt*Pe*(Kappa.q)               !
          ! Fbead=-A'.Fseg                      !
          ! qstar=q                             !
          ! qstar:=qstar+Kdotq                  !
          ! Fphi=Fbead+Fev+Fbnd                 !
          ! qstar:=qstar+(1/4)*dt*(AdotD.Fbead) !
          ! qstar:=qstar+(1/4)*dt*(AdotD.Fev)   !
          ! qstar:=qstar+(1/4)*dt*(AdotD.Fbnd)  !
          ! qstar:=qstar+FBr                    !
          !-------------------------------------!
          call gbmv(KappaBF,qc,Kdotq,kl=0,alpha=Pe(iPe)*dt(iPe,idt))
          if (tplgy == 'Linear') then
            call gbmv(AmatBF,Fseg,Fbead,kl=0,m=nsegx3,alpha=-1.0_wp,trans='T')
          else
            call gemv(Amat,Fseg,Fbead,alpha=-1._wp,trans='T')
          end if
          call copy(qc,qstar)
          call axpy(Kdotq,qstar)
          Fphi(:,ichain)=Fbead+Fev+Fbnd
          if (applFext) then
            Fphi(1,ichain)=Fphi(1,ichain)-Fext0
            Fphi(nbeadx3-2,ichain)=Fphi(nbeadx3-2,ichain)+Fext0
          end if
          if (srf_tet) then
            call tetforce(rvmrcP,rcm(:,ichain),DiffTensP,dt(iPe,idt),Ftet,&
              rf0(:,ichain),itime)
            Fphi(1:3,ichain)=Fphi(1:3,ichain)+Ftet(1:3)
          end if
          call gemv(AdotDP1,Fphi(:,ichain),qstar,alpha=0.25*dt(iPe,idt),&
            beta=1._wp)

          !TYL: HI for tethered bead -------------------------------------------
          !print *, 'Predictor calculaton------------'
          !call print_matrix(AdotDP1,'AdotDP1')
          !call print_vector(Fphi(:,ichain),'Fphi(:,ichain)')
          !call print_vector(qstar,'qstar BEFORE')
          if (srf_tet) then
            call gemv(AdotDP1(:,1:3),Fphi(1:3,ichain),qstar,&
              alpha=-0.25*dt(iPe,idt),beta=1._wp)
            !call print_vector(qstar,'qstar AFTER')
          end if
          !TYL: HI for tethered bead -------------------------------------------

          !! Blake's part
          if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
            do is=1, nseg
              os=(is-1)*3
              qstar(os+2)=qstar(os+2)+(divD(is+1)-divD(is))*0.25*dt(iPe,idt)
            enddo
          endif
          !!-------------


          call axpy(FBr,qstar)
          !-------First Corrector Algorithm-------!
          ! RHS=q                                 !
          ! RHS:=RHS+1/2*Kdotq (from Predictor)   !
          ! Fbarev:=Fev+Fstarev                   !
          ! RHS:=RHS+(1/4)dt*(AdotD.Fev)          !
          ! RHS:=RHS+FBr (from Predictor)         !
          ! RHScnt=RHS(part of it for 2ndCorr.)   !
          ! RHS:=RHS+1/2*dt*(Pe*Kappa.qstar)      !
          ! RHS:=RHS+(1/2)*dt*Fseg                !
          ! Fbarseg=Fseg; Fbarbead=Fbead          !
          ! Inside the loop:                      !
          ! RHSP:=RHSP+(1/4)*dt*(AdotDP.Fbarbead) !
          ! Fbarbead=-A'.Fbarseg                  !
          !---------------------------------------!
          call gemv(Bmat,qstar,rvmrcP)
          call copy(qc,RHS)
          call axpy(Kdotq,RHS,a=0.5_wp)

          ! rflc doesn't need this
          !            if ((EV_bb/='NoEV').or.(EV_bw/='NoEV') .and. EV_bw /= 'Rflc_bc') then
          if ((EV_bb/='NoEV').or.(EV_bw/='NoEV')) then
            !              call EVUpdate(Fev,rvmrcP,Fbarev)
            call myintrn%calc(rvmrcP,rcmP,nseg,DiffTensP,divD,Fev,Fbarev,&
              updtevbb=.true.,updtevbw=.true.)
            !call print_vector(Fev,'fev3')
            !call evupdate2(Fev,rvmrcP,nseg,Fbarev)
            !call print_vector(Fev,'fev4')
            !stop
          end if


          if (ForceLaw == 'WLC_GEN') then
            call bndupdate(nbead_bb,Fbnd,qstar,Fbarbnd,itime)
          end if
          Fbar=Fbarev+Fbarbnd
          if (applFext) then
            Fbar(1)=Fbar(1)-Fext0
            Fbar(nbeadx3-2)=Fbar(nbeadx3-2)+Fext0
          end if
          if (srf_tet) then
            call tetupdate(Ftet,rvmrcP,rcm(:,ichain),DiffTensP,dt(iPe,idt),&
             Fbartet,rf0(:,ichain),itime)
            Fbar(1:3)=Fbar(1:3)+Fbartet(1:3)
            !call print_vector(Fbartet(1:3),'Tether force')
          end if
          call gemv(AdotDP1,Fbar,RHS,alpha=0.25*dt(iPe,idt),beta=1._wp)

          !TYL: HI for tethered bead -------------------------------------------
          !print *, 'First corrector calculaton------------'
          !call print_matrix(AdotDP1,'AdotDP1')
          !call print_vector(Fbar,'Fbar')
          !call print_vector(RHS,'RHS BEFORE')
          if (srf_tet) then
            call gemv(AdotDP1(:,1:3),Fbar(1:3),RHS,alpha=-0.25*dt(iPe,idt),&
              beta=1._wp)
            !call print_vector(RHS,'RHS AFTER')
          end if
          !TYL: HI for tethered bead -------------------------------------------


          !! Blake's part
          if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
            do is=1, nseg
              os=(is-1)*3
              RHS(os+2)=RHS(os+2)+(divD(is+1)-divD(is))*0.25*dt(iPe,idt)
            enddo
          endif
          !!-------------


          call axpy(FBr,RHS)
          call copy(RHS,RHScnt)
          call gemv(Kappa,qstar,RHS,alpha=0.5*Pe(iPe)*dt(iPe,idt),beta=1._wp)
          call axpy(Fseg,RHS,a=0.5*dt(iPe,idt))
          call copy(Fseg,Fbarseg)
          call copy(Fbead,Fbarbead)
          do iseg=1, nseg
            offset=3*(iseg-1)
            RHSP => RHS(offset+1:offset+3)
            AdotDP2 => AdotD(offset+1:offset+3,:,ichain)
            call gemv(AdotDP2,Fbarbead,RHSP,alpha=0.25*dt(iPe,idt),beta=1._wp)

            !TYL: HI for tethered bead -----------------------------------------
            !print *, 'First corrector calculaton per segment------------'
            !call print_matrix(AdotDP2,'AdotDP2')
            !call print_vector(Fbarbead,'Fbarbead')
            !call print_vector(RHSP,'RHSP BEFORE')
            if (srf_tet) then
              call gemv(AdotDP2(:,1:3),Fbarbead(1:3),RHSP,&
              alpha=-0.25*dt(iPe,idt),beta=1._wp)
              !call print_vector(RHSP,'RHSP AFTER')
            end if
            !TYL: HI for tethered bead -----------------------------------------

            call sprupdate(id,root_f,PrScale,nroots,dt(iPe,idt),RHSP,qstar,iseg,&
              nseg,ForceLaw,TruncMethod,qbar,Fbarseg,Fbarbead,tplgy,Amat,nseg_bb,&
              nseg_ar,Ia,Na,itime)
          end do
          !----------Second Corrector Algorithm----------!
          ! q=qbar;Fseg=Fbarseg;Fbead=Fbarbead           !
          ! RHSbase=RHScnt(from 1stCorr.)for while loop. !
          ! While Loop,do loop:                          !
          ! RHSP=RHSbaseP                                !
          ! RHSP:=RHSP+(1/2)*dt*(Pe*Kappareg.qP)         !
          ! RHSP:=RHSP+(1/2)*dt*FsegP                    !
          ! RHSP:=RHSP+(1/4)dt*(AdotDP.Fbead)            !
          ! Updating q based on Seg. Cubic Eq.           !
          ! Fbead=-A'.Fseg                               !
          !----------------------------------------------!
          call copy(qbar,qc)
          call copy(Fbarseg,Fseg)
          call copy(Fbarbead,Fbead)
          call copy(RHScnt,RHSbase)
          icount=0;eps=1.0_wp
          do while (eps >= tol)
            eps=0.0_wp
            qctemp=qc
            do iseg=1, nseg
              offset=3*(iseg-1)
              RHSP => RHS(offset+1:offset+3);RHSbaseP => RHSbase(offset+1:offset+3)
              call copy(RHSbaseP,RHSP)
              qcP => qc(offset+1:offset+3);FsegP => Fseg(offset+1:offset+3)
              AdotDP2 => AdotD(offset+1:offset+3,:,ichain)
              call gemv(Kappareg,qcP,RHSP,alpha=0.5*Pe(iPe)*dt(iPe,idt),beta=1.0_wp)
              call axpy(FsegP,RHSP,a=0.5*dt(iPe,idt))
              call gemv(AdotDP2,Fbead,RHSP,alpha=0.25*dt(iPe,idt),beta=1.0_wp)

              !TYL: HI for tethered bead ---------------------------------------
              !print *, 'Second corrector calculaton perbead------------'
              !call print_matrix(AdotDP2,'AdotDP2')
              !call print_vector(Fbead,'Fbead')
              !call print_vector(RHSP,'RHSP BEFORE')
              if (srf_tet) then
                call gemv(AdotDP2(:,1:3),Fbead(1:3),RHSP,&
                alpha=-0.25*dt(iPe,idt),beta=1._wp)
                !call print_vector(RHSP,'RHSP AFTER')
              end if
              !TYL: HI for tethered bead ---------------------------------------

              call sprupdate(id,root_f,PrScale,nroots,dt(iPe,idt),RHSP,qbar,iseg,&
                nseg,ForceLaw,TruncMethod,qc,Fseg,Fbead,tplgy,Amat,nseg_bb,nseg_ar,&
                Ia,Na,itime)

            end do
            eps=nrm2(qc-qctemp)/nrm2(qctemp)
            icount=icount+1
            if (icount > 5000) then
              print *
              print '(" Convergance Problem in 2nd Corrector.")'
              print '(" time index: ",i10)',itime
              print '(" Total iterations: ",i10," Residual: ",f14.7)',icount,eps
              if (hstar /= 0._wp) then
                if (DecompMeth == 'Lanczos') then
                  print '(" No. iterations in (block) Lanczos algorithm: ",i4)',&
                  mch(ichain)
                elseif (DecompMeth == 'Chebyshev') then
                  print '(" Eigen value range for  diffusion tensor: ",2(f14.7))',&
                  lambdaBE(:)
                  print '(" No. iterations in Chebyshev algorithm: ",i4)',Lch(ichain)
                end if
              end if
              stop
            end if
          end do ! while loop
          !==================================================!

          ! Inserting back the final result to original arrays
          q(:,ichain)=qc(:)
          Fphi(:,ichain)=Fbead(:)+Fbar(:)
          call gemv(Bmat,qc,rvmrcP)
          ! Calculating center of mass and/or center of hydrodynamic resistance movement
          if (CoM) then
            FphiP => Fphi(:,ichain)
            rcmP => rcm(:,ichain) ! might be redundant
            BdotwPx => wbltemp(1:nbeadx3-2:3,jcol,ichain)
            BdotwPy => wbltemp(2:nbeadx3-1:3,jcol,ichain)
            BdotwPz => wbltemp(3:nbeadx3:3,jcol,ichain)
            SumBdotw=[sum(real(BdotwPx,kind=wp)),&
            sum(real(BdotwPy,kind=wp)),&
            sum(real(BdotwPz,kind=wp))]

            !TYL: HI for tethered bead ---------------------------------------
            if (srf_tet) then
              SumBdotw = SumBdotw - [real(BdotwPx(1),kind=wp),&
              real(BdotwPy(1),kind=wp),&
              real(BdotwPz(1),kind=wp)]
            end if
            !TYL: HI for tethered bead ---------------------------------------

            if (hstar.ne.0._wp) then
              call symv(DiffTensP,FphiP,DdotF)
            else
              DdotF=FphiP
            end if

            !TYL: HI for tethered bead ---------------------------------------
            if (srf_tet) then
              call gemv(DiffTensP(:,1:3),FphiP(1:3),DdotF,alpha=-1._wp,beta=1._wp)
            end if
            !TYL: HI for tethered bead ---------------------------------------

            DdotFPx => DdotF(1:nbeadx3-2:3)
            DdotFPy => DdotF(2:nbeadx3-1:3)
            DdotFPz => DdotF(3:nbeadx3:3)
            SumDdotF=(/sum(DdotFPx),sum(DdotFPy),sum(DdotFPz)/)
            rcmP=rcmP+(Pe(iPe)*matmul(Kappareg,rcmP)+1._wp/(4*nbead)*SumDdotF)*&
            dt(iPe,idt)+coeff/nbead*SumBdotw



            !! Blake's part
            if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
              rcmP(2)=rcmP(2)+1._wp/(4*nbead)*sum(divD)*dt(iPe,idt)
            endif
            !!-------------


          end if
          if (CoHR) then
            if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
              if (DecompMeth == 'Cholesky') then
                MobilTens=CoeffTensP
              else
                MobilTens=DiffTensP
                call potrf(MobilTens,info=info)
                if (info /= 0) then
                  print '(" Unsuccessful Cholesky fact. of Mobility Tensor in main.")'
                  print '(" info: ",i3)',info
                  stop
                end if
              end if
              call potri(MobilTens,info=info)
              if (info /= 0) then
                print '(" Unsuccessful diffusion-matrix inversion in HIEVCalc.")'
                print '(" info: ",i3)',info
                stop
              end if
              do i=1, nbeadx3
                do j=i+1, nbeadx3
                  MobilTens(j,i)=MobilTens(i,j)
                end do
              end do
              WeightTenstmp=0._wp;totMobilTens=0._wp
              do ibead=1, nbead
                ! Calculating the global location
                iglob=3*(ibead-1)
                MobilTensP1 => MobilTens(iglob+1:iglob+3,:)
                WeightTenstmp=WeightTenstmp+MobilTensP1
                do jbead=ibead, nbead
                  jglob=3*(jbead-1)
                  MobilTensP2 => MobilTens(iglob+1:iglob+3,jglob+1:jglob+3)
                  if (jbead == ibead) then
                    totMobilTens=totMobilTens+MobilTensP2
                  else
                    totMobilTens=totMobilTens+MobilTensP2+transpose(MobilTensP2)
                  end if
                end do
              end do
              invtotMobilTens=totMobilTens
              call getrf(invtotMobilTens,ipiv=ipiv,info=info)
              if (info /= 0) then
                print '(" Unsuccessful LU factorization in HIEVCalc.")'
                print '(" info: ",i3)',info
                stop
              end if
              call getri(invtotMobilTens,ipiv,info=info)
              if (info /= 0) then
                print '(" Unsuccessful Total-Mobility-Tensor inversion in HIEVCalc.")'
                print '(" info: ",i3)',info
                stop
              end if
              call gemm(invtotMobilTens,WeightTenstmp,WeightTens)
            end if ! mod(itime,ncol)==1 or ncol==1
            rchrP => rchr(:,ichain)
            BdotwP => wbltemp(:,jcol,ichain)
            call gemv(WeightTens,real(BdotwP,kind=wp),LdotBdotw)
            rchrP=rchrP+Pe(iPe)*matmul(kappareg,rchrP)*dt(iPe,idt)+coeff*LdotBdotw
          end if


          ! rflc part
          if (EV_bw == 'Rflc_bc') then
            !call wall_rflc(dt(iPe,idt),itime,time,id,ichain,qPy,rvmrcPy,rcmP(2),rf_in)
            call wall_rflc(myintrn%evbw,dt(iPe,idt),itime,time,id,ichain,&
              qPx,qPy,qPz,rvmrcPx,rvmrcPy,rvmrcPz,rcmP(1),rcmP(2),rcmP(3),&
              rf_in)
          end if
          !-----


        end do ! ichain loop

        !----------------------------------------------------------------
        !>>>>> data processing and outputs:
        !----------------------------------------------------------------

        if ( (time >= time_check3) .or. (itime == ntime(iPe,idt)) ) then
          time_check3=time_check3+frm_rt_pp*lambda


          if (EV_bw == 'Rflc_bc') then
            if (itime == ntime(iPe,idt)) call print_wcll(myintrn%evbw,id,p,MPI_REAL_WP,time)
          endif


          call data_prcs(id,itime,time,idt,iPe,q,rvmrc,Fphi,rcm,rcmstart,rchr,&
           nseg_bb,mch,Lch,MPI_REAL_WP)
          if (cnf_srt) then
            call conf_sort(q,rvmrc,nseg_bb,nbead_bb,cnf_tp)
            call MPI_Gatherv(cnf_tp,npchain,MPI_INTEGER,cnf_tpTot,cnf_tp_counts,&
             cnf_tp_disps,cnf_tp_resizedrecvsubarray,0,MPI_COMM_&
             &WORLD,ierr)
            if (id == 0) then
              do ip = 1, p
                do ichain=1, npchain
                  write(u41,'(" ",i1)') cnf_tpTot(ichain,ip)
                end do
              end do
            end if
          end if
        end if
        ! To be done at specified strains:
        if (initmode == 'rst') then
          strain=Pe(iPe)*(time+trst*lambda)
        else
          strain=Pe(iPe)*time
        end if
        if (DumpstrCalc) then
          if (istr <= nstr) then
            if (strain >= dumpstr(istr)) then
              ! For writing data to the output file:
              call MPI_Gatherv(q,nsegx3*npchain,MPI_REAL_WP,qTot,q_counts,q_disps,&
               q_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
              if (id == 0) then
                format_str="(A,f6.2,'.dat')"
                if (initmode.eq.'rst') then
                  write(file1,format_str)'data/qdump_str',Pe(iPe)*(time+trst*lambda)
                else
                  write(file1,format_str)'data/qdump_str',Pe(iPe)*time
                end if
                open(newunit=u40,file=trim(adjustl(file1)),status='replace')
                do ip = 1, p
                  do ichain=1, npchain
                    do iseg= 1, nseg
                      offset=3*(iseg-1)
                      write(u40,1) qTot(offset+1:offset+3,ichain,ip)
                    end do
                  end do
                end do
              end if ! id == 0
              istr=istr+1
            end if ! strain ...
          end if ! istr ...
        end if ! DumpstrCalc
        if (time >= time_check4) then
          time_check4=time_check4+frm_rt_rst*lambda
          ! For writing restart data to the output file
          call MPI_Gatherv(q,nsegx3*npchain,MPI_REAL_WP,qTot,q_counts,q_disps,&
           q_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          if (CoM) then
            call MPI_Gatherv(rcm,3*npchain,MPI_REAL_WP,rcmT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (CoHR) then
            call MPI_Gatherv(rchr,3*npchain,MPI_REAL_WP,rchrT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (id == 0) then
            open(newunit=u21,file='data/q.rst.dat',status='replace')
            if (CoM) open(newunit=u22,file='data/CoM.rst.dat',status='replace')
            if (CoHR) open(newunit=u23,file='data/CoHR.rst.dat',status='replace')
            do ip = 1, p
              do ichain =1, npchain
                do iseg = 1, nseg
                  offset=3*(iseg-1)
                  write(u21,1) qTot(offset+1:offset+3,ichain,ip)
                end do
                if (CoM) write(u22,1) rcmT(1:3,ichain,ip)
                if (CoHR) write(u23,1) rchrT(1:3,ichain,ip)
              end do
            end do
            if (initmode == 'rst') then
              rtpassed=(time+trst*lambda)/lambda
            else
              rtpassed=time/lambda
            end if
            write(u21,'(f7.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed";close(u21)
            if (CoM) then
              write(u22,'(f7.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed";close(u22)
            end if
            if (CoHR) then
              write(u23,'(f7.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed";close(u23)
            end if
          end if ! id == 0
        end if ! mod(itime,lambda/dt)==0



        ! if (time >= time_check2) then
        if ((time >= tss*lambda) .and. (mod(itime,tgap_dmp) == 0)) then

          time_check2=time_check2+frm_rt_dmp*lambda
          ! For writing equilibrium or final data to the output file
          call MPI_Gatherv(q,nsegx3*npchain,MPI_REAL_WP,qTot,q_counts,q_disps,&
           q_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          call MPI_Gatherv(rvmrc,nbeadx3*npchain,MPI_REAL_WP,RTot,R_counts,R_disps,&
           R_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          call MPI_Gatherv(Fphi,nbeadx3*npchain,MPI_REAL_WP,FphiTot,R_counts,R_disps,&
           R_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          if (CoM) then
            call MPI_Gatherv(rcm,3*npchain,MPI_REAL_WP,rcmT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (CoHR) then
            call MPI_Gatherv(rchr,3*npchain,MPI_REAL_WP,rchrT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (id == 0) then
            if (TimerA) then
              open(newunit=u24,file='data/TimerA.dat',status='unknown',position='append')
              write(u24,*) "iflow, Wi, dt, After ~ Chain RelTime, Nbeads, Nchains, Nproc, ExecutionTime"
              write(u24,*) "---------------------------------------------------------------------------"
              !                tA1=MPI_Wtime()
              call cpu_time(time_end)
              !                write(24,2) iflow,Wi(iPe),dt(iPe,idt),tend,nbead,nchain,p,(tA1-tA0)
              write(u24,2) iflow,Wi(iPe),dt(iPe,idt),tend,nbead,nchain,p,(time_end-time_begin)
            end if
            ! Writing equilibrium or final data to the output file
            if (jcheck == 0) then
              jcheck=-1
              if (initmode == 'st') then
                fstat='replace'
              else
                fstat='unknown'
              end if
              if (iflow == 1) then
                open(newunit=u34,file='data/q.equil.dat',status=fstat,position='append')
                open(newunit=u25,file='data/R.equil.dat',status=fstat,position='append')
                open(newunit=u42,file='data/fphi.equil.dat',status=fstat,position='append')
                if (CoM) open(newunit=u26,file='data/CoM.equil.dat',status=fstat,position='append')
                if (CoHR) open(newunit=u27,file='data/CoHR.equil.dat',status=fstat,position='append')
              else
                open(newunit=u34,file='data/q.flow.dat',status=fstat,position='append')
                open(newunit=u25,file='data/R.flow.dat',status=fstat,position='append')
                open(newunit=u42,file='data/fphi.flow.dat',status=fstat,position='append')
                if (CoM) open(newunit=u26,file='data/CoM.flow.dat',status=fstat,position='append')
                if (CoHR) open(newunit=u27,file='data/CoHR.flow.dat',status=fstat,position='append')

              end if ! iflow == 1
            end if ! jcheck
            do ip = 1, p
              do ichain =1, npchain
                do iseg = 1, nseg
                  offset=3*(iseg-1)
                  write(u34,1) qTot(offset+1:offset+3,ichain,ip)
                end do
                do ibead = 1, nbead
                  offset=3*(ibead-1)
                  write(u25,1) RTot(offset+1:offset+3,ichain,ip)
                  write(u42,1) FphiTot(offset+1:offset+3,ichain,ip)
                end do
                if (CoM) write(u26,1) rcmT(1:3,ichain,ip)
                if (CoHR) write(u27,1) rchrT(1:3,ichain,ip)
              end do
            end do
          end if ! id==0
        end if ! time >= ...

        jcol=jcol+1 ! col in the block of random number columns
      end do ! time loop
      ! resetting restart time
      if (initmode == 'rst') trst=0._wp
    end do ! dt loop
  end do ! Pe loop

  !----------------------------------------------------------------
  !>>>>> Deallocation of arrays and closing files:
  !----------------------------------------------------------------
  if (id == 0) then
    close(u25);close(u34)
    if (TimerA) close (u24)
    if (CoM) close(u26)
    if (CoHR) close(u27)
    if (DumpstrCalc) close(u40)
    if (cnf_srt) close(u41)
  end if

  if (id == 0) then
    deallocate (rdnt)
    deallocate (qTot,qxT,qyT,qzT,RTot,FphiTot)
    if (CoM) deallocate(rcmT)
    if (CoHR) deallocate(rchrT)
  end if

  deallocate(rdn)
  deallocate(qc,qstar,Fseg,w,wbl,wbltemp,Kappa,Amat,FBr,FBrbl,Kdotq,AdotD)
  deallocate(Fbead,RHS,RHScnt,RHSbase,Fbarseg,ADFev,qbar,Fbarbead,Fev,Fphi)
  deallocate(Bmat,q,qstart,Fbarev,KappaBF,qctemp,DdotF,rvmrc,DiffTens)
  deallocate(CoeffTens,Fbnd,Fbarbnd)
  if (tplgy == 'Linear') deallocate(AmatBF)
  if (DecompMeth == 'Chebyshev') then
    deallocate(Dsh,Eye,Ybar,uminus,uplus,Ddotuminus,Ddotuplus,Lch)
  elseif (DecompMeth == 'Lanczos') then
    deallocate(aBlLan,WBlLan,VBlLan,Ybar,VcntBlLan,mch)
  end if
  if (CoM) deallocate(rcm,rcmstart)
  if (CoHR) deallocate(rchr,rchrstart,MobilTens,WeightTens,WeightTenstmp)
  deallocate(q_counts,q_disps)
  if ((CoM) .or. (CoHR)) deallocate(rc_counts,rc_disps)
  if ((tplgy == 'Comb').and.(arm_plc /= 'Fixed')) then
    deallocate(ia_counts,ia_disps)
  end if
  if (cnf_srt) then
    deallocate(cnf_tp,cnf_tp_counts,cnf_tp_disps)
  end if
  call del_inp()
  call del_pp(id)

  contains

  !----------------------------------------------------------------
  !>>>>> Inline subroutines:
  !----------------------------------------------------------------

  ! Setting up the lookup table For using in Corrector:
  subroutine lookup_tab(dttmp)

    use :: root_mod, only: CubeRoot,root_fndr

    real(wp) :: dttmp,a1,a2,a3,coeffs(8),rhsmag,denom
    integer :: nr

    if (ForceLaw /= 'Hookean') then
      do nr=1, (PrScale*nroots)-1
        ! Note!!: the root of rhsmag=(0.01/PrScale) is root_f(2).
        rhsmag=nr*(0.01_wp/PrScale)
        select case (ForceLaw)
        case ('WLC_MS')
          ! WLC, Worm Like Chain proposed by Marko and Siggia
          denom=1+dttmp/3
          a1=-(((2+0.75*dttmp)*qmax)+rhsmag)/denom
          a2=(((1+(0.5*dttmp))*b)+(2*qmax*rhsmag))/denom
          a3=-rhsmag*b/denom
        case ('WLC_UD','WLC_GEN')
          ! WLC model by Underhill and Doyle and the generalized version
          coeffs(1)=rhsmag
          coeffs(2)=-1-dttmp/3+7*dttmp/(6*WLC_v)-dttmp/3*(WLC_A+WLC_B)
          coeffs(3)=-2*rhsmag/b
          coeffs(4)=2/b-7*dttmp/(6*WLC_v*b)+dttmp/b*(2*WLC_A/3+WLC_B)
          coeffs(5)=rhsmag/b**2
          coeffs(6)=-1/b**2-dttmp/b**2*(WLC_A/3+WLC_B)
          coeffs(7)=0._wp
          coeffs(8)=WLC_B*dttmp/(3*b**3)
        case ('ILCCP')
          ! ILCCP, Inverse Langevin Chain (Cohen-Pade approximation)
          denom=1+dttmp/6
          a1=-rhsmag/denom
          a2=-(1+0.5*dttmp)*b/denom
          a3=rhsmag*b/denom
        case ('FENE')
          ! FENE
          a1=-rhsmag
          a2=-b*(1+dttmp/2)
          a3=b*rhsmag
        case ('RWS')
          ! RWS, From Underhill and Doyle
          denom=1+RWS_D*dttmp/6
          a1=-rhsmag/denom
          a2=-(1+RWS_C/6*dttmp)*b/denom
          a3=rhsmag*b/denom
        end select
        if ((ForceLaw == 'WLC_UD').or.(ForceLaw == 'WLC_GEN')) then
          call root_fndr(real(coeffs,kind=double),real(qmax,kind=double),root_f(nr+1))
        else
          call CubeRoot(real(a1,kind=double),real(a2,kind=double)  ,&
            real(a3,kind=double),real(qmax,kind=double),&
            root_f(nr+1))
        end if
      end do
      root_f(1)=0._wp
    end if

  end subroutine lookup_tab

  ! Random numeber seeding (from H. C. Ottinger):
  subroutine ranils(iseed)

    integer,intent(in) :: iseed
    integer,parameter :: in=2147483563,ik=40014,iq=53668,ir=12211,ntab=32
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Initial seeds for two random number generators
    idum=iseed+123456789
    idum2=idum

    ! Load the shuffle table (after 8 warm-ups)
    do 10 j=ntab+8,1,-1
     k=idum/iq
     idum=ik*(idum-k*iq)-k*ir
     if(idum < 0) idum=idum+in
     if(j <= ntab) iv(j)=idum
     10 continue
     iy=iv(1)
     return

   end subroutine ranils

   ! Uniform random number generator (from H. C. Ottinger):
   real(wp) function ranuls()

   integer,parameter :: in1=2147483563,ik1=40014,iq1=53668,ir1=12211,&
   in2=2147483399,ik2=40692,iq2=52774,ir2=3791 ,&
   ntab=32,inm1=in1-1,ndiv=1+inm1/ntab
   real(wp),parameter :: an=1./in1
   integer :: iv(ntab),idum,idum2,iy
   integer :: k,j

   common /ranbls/ idum,idum2,iy,iv

   ! Linear congruential generator 1
   k=idum/iq1
   idum=ik1*(idum-k*iq1)-k*ir1
   if(idum < 0._wp) idum=idum+in1

   ! Linear congruential generator 2
   k=idum2/iq2
   idum2=ik2*(idum2-k*iq2)-k*ir2
   if(idum2 < 0._wp) idum2=idum2+in2

   !Shuffling and subtracting
   j=1+iy/ndiv
   iy=iv(j)-idum2
   iv(j)=idum
   if(iy < 1) iy=iy+inm1
   ranuls=an*iy
   return

 end function ranuls

 ! Gaussian random number generator (from H. C. Ottinger):
 real(wp) function rangls()

 integer :: iflag
 real(wp) :: gauss2,x1,x2,xsq,aux

 save iflag,gauss2
 data iflag/0/

 if(iflag == 0) then
  10 continue

  ! pair of uniform random numbers in [-1,1]x[-1,1]
  x1=2*ranuls()-1
  x2=2*ranuls()-1

  ! if not in the unit circle, try again
  xsq=x1*x1+x2*x2
  if(xsq >= 1._wp .or. xsq == 0._wp) goto 10
  ! pair of gaussian random numbers; return one and
  ! save the other for next time
  aux=sqrt(-2*log(xsq)/xsq)
  rangls=x1*aux
  gauss2=x2*aux
  iflag=1
else
  rangls=gauss2
  iflag=0
endif
return

end function rangls

end subroutine dlt_bs

end module dlt_mod
