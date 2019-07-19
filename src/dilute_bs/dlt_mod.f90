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
  use :: intrn_mod, only: intrn_t,wall_rflc_sph,wall_rflc,print_wcll
  use :: rand_mod, only: ranils,ranuls,rangls
  use :: sde_mod, only: sde_t
  use :: sph_sde_mod, only: sph_sde_t

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
    integer :: rsph_recvsubarray,rsph_resizedrecvsubarray
    integer :: psph_recvsubarray,psph_resizedrecvsubarray
    integer :: qsph_recvsubarray,qsph_resizedrecvsubarray
    integer :: rcm_recvsubarray,rcm_resizedrecvsubarray
    integer :: ia_recvsubarray,ia_resizedrecvsubarray
    integer :: cnf_tp_recvsubarray,cnf_tp_resizedrecvsubarray
    integer,dimension(3) :: q_starts,q_sizes,q_subsizes
    integer,dimension(3) :: R_starts,R_sizes,R_subsizes
    integer,dimension(3) :: rc_starts,rc_sizes,rc_subsizes
    integer,dimension(3) :: rsph_starts,rsph_sizes,rsph_subsizes
    integer,dimension(3) :: psph_starts,psph_sizes,psph_subsizes
    integer,dimension(3) :: qsph_starts,qsph_sizes,qsph_subsizes
    integer,dimension(4) :: rcm_starts,rcm_sizes,rcm_subsizes
    integer,dimension(2) :: ia_starts,ia_sizes,ia_subsizes
    integer,dimension(2) :: cnf_tp_starts,cnf_tp_sizes,cnf_tp_subsizes
    integer,allocatable  :: q_counts(:),q_disps(:)
    integer,allocatable  :: R_counts(:),R_disps(:)
    integer,allocatable  :: rc_counts(:),rc_disps(:)
    integer,allocatable  :: rsph_counts(:),rsph_disps(:)
    integer,allocatable  :: psph_counts(:),psph_disps(:)
    integer,allocatable  :: qsph_counts(:),qsph_disps(:)
    integer,allocatable  :: rcm_counts(:),rcm_disps(:)
    integer,allocatable  :: ia_counts(:),ia_disps(:)
    integer,allocatable  :: cnf_tp_counts(:),cnf_tp_disps(:)
    integer(kind=MPI_ADDRESS_KIND) :: q_start,q_extent
    integer(kind=MPI_ADDRESS_KIND) :: R_start,R_extent
    integer(kind=MPI_ADDRESS_KIND) :: rc_start,rc_extent
    integer(kind=MPI_ADDRESS_KIND) :: rsph_start,rsph_extent
    integer(kind=MPI_ADDRESS_KIND) :: psph_start,psph_extent
    integer(kind=MPI_ADDRESS_KIND) :: qsph_start,qsph_extent
    integer(kind=MPI_ADDRESS_KIND) :: rcm_start,rcm_extent
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
    integer :: nu,mu,ibead,jseg,ichain,ip,iread,iPe,idt,iAdjSeq,itime,ichain_pp
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
    real(wp),dimension(3) :: SumBdotw,SumDdotF,LdotBdotw!,Ftet,Fbartet!,rf_in
    real(wp),dimension(3,3) :: kappareg,totMobilTens,invtotMobilTens
    real(wp),dimension(3) :: U_unif
    ! Allocatable arrays:
    integer,allocatable,dimension(:) :: mch,Lch
    real(wp),allocatable,dimension(:) :: qstar,qbar,Fbead,Fbarseg
    real(wp),allocatable,dimension(:) :: Fbarbead,Fev,w,FBr,Kdotq,ADFev,RHScnt,Fev_sph
    real(wp),allocatable,dimension(:) :: Fbarev,root_f,Ybar,qctemp,uminus,uplus
    real(wp),allocatable,dimension(:) :: Ddotuminus,Ddotuplus
    real(wp),allocatable,dimension(:) :: Fbnd,Fbarbnd,Fbar,divD
    real(wp),allocatable,dimension(:),target :: RHS,RHSbase,qc,Fseg,DdotF
    real(wp),dimension(:),pointer  :: RHSP,RHSbaseP,qcP,FsegP,FBrblP!,rcmP,rcmPy
    real(wp),dimension(:),pointer :: rchrP,FphiP,DdotFPx,DdotFPy,DdotFPz,qP
    real(wp),dimension(:),pointer :: rvmrcP,rvmrcPx,rvmrcPy,rvmrcPz,qPx,qPy,qPz
    real(double),dimension(:),pointer :: wbltempP2,BdotwP,BdotwPx,BdotwPy,BdotwPz
    real(wp),allocatable,dimension(:,:) :: Kappa,Amat,Bmat,wbl,aBlLan,WBlLan,Amat_sph
    real(wp),allocatable,dimension(:,:) :: VBlLan,VcntBlLan,Eye,Dsh
    real(wp),allocatable,dimension(:,:) :: KappaBF,AmatBF,WeightTenstmp!,rf0
    real(wp),allocatable,dimension(:,:),target :: q,rchr,Fphi,MobilTens,rvmrc!,rcm
    real(wp),allocatable,dimension(:,:),target :: r_sph
    real(wp),allocatable,dimension(:),target :: F_sph,Fphi_all,Fphi_all_temp
    real(wp),allocatable,dimension(:,:),target :: p_sph, q_sph
    real(wp),allocatable,dimension(:,:),target :: WeightTens,qstart!,rcmstart
    real(wp),allocatable,dimension(:,:),target :: rchrstart
    real(wp),dimension(:,:),pointer :: AdotDP2,MobilTensP1,MobilTensP2
    real(wp),dimension(:,:),pointer :: DiffTensP,AdotDP1
    real(double),dimension(:,:),pointer :: wbltempP1
    real(double),dimension(:,:),pointer :: CoeffTensP
    real(wp),allocatable,dimension(:,:,:) :: rdnt,rdn
    real(wp),allocatable,dimension(:,:,:) :: rdnt_sph,rdn_sph
    real(wp),allocatable,dimension(:,:) :: wbl_sph
    real(wp),allocatable,dimension(:,:,:) :: rdnt_sph_or,rdn_sph_or
    real(wp),allocatable,dimension(:,:) :: wbl_sph_or
    ! For gathering data by rank 0
    real(wp),allocatable,dimension(:,:,:),target :: qxT,qyT,qzT,qTot,RTot,FphiTot!,rcmT
    real(wp),allocatable,dimension(:,:,:),target :: rchrT,FBrbl,DiffTens,AdotD
    real(wp),allocatable,dimension(:,:,:),target :: r_sphT
    real(wp),allocatable,dimension(:,:,:),target :: p_sphT, q_sphT
    real(double),allocatable,dimension(:,:,:),target :: CoeffTens,wbltemp
    integer,allocatable,dimension(:) :: cnf_tp
    integer,allocatable,dimension(:,:) :: iaTot,cnf_tpTot
    ! File units
    integer :: u2,u3,u4,u21,u22,u23,u24,u25,u26,u27,u34,u39,u40,u41,u42,u_rf_in
    integer :: u_rsph_rst,u_pqsph_rst,u_rsph_flow,u_pqsph_flow
    ! objects
    type(intrn_t) :: myintrn
    type(sde_t) :: mysde
    type(sph_sde_t) :: mysphsde

    !need to organize later
    !real(wp),dimension(nsegx3) :: U_seg
    !real(wp),dimension(nbeadx3) :: U_bead
    real(wp),allocatable,dimension(:) :: U_seg
    real(wp),allocatable,dimension(:) :: U_bead

    real(wp),allocatable,dimension(:,:,:) :: rf_in
    real(wp),allocatable,dimension(:,:,:) :: rf0
    real(wp),dimension(3) :: rf0_rel_unit
    real(wp),dimension(3) :: dr_sph_rflc
    real(wp),allocatable,dimension(:,:,:),target :: rcm,rcmstart
    real(wp),dimension(:,:),pointer :: rcmP,rcmPy
    real(wp),dimension(:),pointer :: r_sphP
    real(wp),dimension(:),pointer :: p_sphP,q_sphP
    real(wp),allocatable,dimension(:,:,:,:),target ::rcmT
    integer :: idx_pp_seg,idx_pp_bead
    real(wp),allocatable,dimension(:) :: Ftet,Fbartet
    logical :: debug_TYL

    debug_TYL = .false.

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
      ! rchr
      if (CoHR) then
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! r_sph
      !rsph_sizes, rsph_subsizes, rsph_starts, rsph_recvsubarray,rsph_extent
      !rsph_start,rsph_resizedrecvsubarray
      if (sph_move) then
        rsph_sizes(1)=3;rsph_sizes(2)=npchain;rsph_sizes(3)=p
        rsph_subsizes(1)=3;rsph_subsizes(2)=npchain;rsph_subsizes(3)=1
        rsph_starts(1)=0;rsph_starts(2)=0;rsph_starts(3)=0
        call MPI_Type_create_subarray(3,rsph_sizes,rsph_subsizes,rsph_starts,MPI_OR&
          &DER_FORTRAN,MPI_REAL_WP,rsph_recvsubarra&
          &y,ierr)
        call MPI_Type_commit(rsph_recvsubarray,ierr)
        rsph_extent=sizeof(realvar)
        rsph_start=0
        call MPI_Type_create_resized(rsph_recvsubarray,rsph_start,rsph_extent,rsph_re&
         &sizedrecvsubarray,ierr)
        call MPI_Type_commit(rsph_resizedrecvsubarray,ierr)

        psph_sizes(1)=3;psph_sizes(2)=npchain;psph_sizes(3)=p
        psph_subsizes(1)=3;psph_subsizes(2)=npchain;psph_subsizes(3)=1
        psph_starts(1)=0;psph_starts(2)=0;psph_starts(3)=0
        call MPI_Type_create_subarray(3,psph_sizes,psph_subsizes,psph_starts,MPI_OR&
          &DER_FORTRAN,MPI_REAL_WP,psph_recvsubarra&
          &y,ierr)
        call MPI_Type_commit(psph_recvsubarray,ierr)
        psph_extent=sizeof(realvar)
        psph_start=0
        call MPI_Type_create_resized(psph_recvsubarray,psph_start,psph_extent,psph_re&
         &sizedrecvsubarray,ierr)
        call MPI_Type_commit(psph_resizedrecvsubarray,ierr)

        qsph_sizes(1)=3;qsph_sizes(2)=npchain;qsph_sizes(3)=p
        qsph_subsizes(1)=3;qsph_subsizes(2)=npchain;qsph_subsizes(3)=1
        qsph_starts(1)=0;qsph_starts(2)=0;qsph_starts(3)=0
        call MPI_Type_create_subarray(3,qsph_sizes,qsph_subsizes,qsph_starts,MPI_OR&
          &DER_FORTRAN,MPI_REAL_WP,qsph_recvsubarra&
          &y,ierr)
        call MPI_Type_commit(qsph_recvsubarray,ierr)
        qsph_extent=sizeof(realvar)
        qsph_start=0
        call MPI_Type_create_resized(qsph_recvsubarray,qsph_start,qsph_extent,qsph_re&
         &sizedrecvsubarray,ierr)
        call MPI_Type_commit(qsph_resizedrecvsubarray,ierr)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! rcm
      if (CoM) then
        rcm_sizes(1)=3;rcm_sizes(2)=npchain;rcm_sizes(3)=nchain_pp;rcm_sizes(4)=p
        rcm_subsizes(1)=3;rcm_subsizes(2)=npchain;rcm_subsizes(3)=nchain_pp;rcm_subsizes(4)=1
        rcm_starts(1)=0;rcm_starts(2)=0;rcm_starts(3)=0;rcm_starts(4)=0
        call MPI_Type_create_subarray(4,rcm_sizes,rcm_subsizes,rcm_starts,MPI_OR&
          &DER_FORTRAN,MPI_REAL_WP,rcm_recvsubarra&
          &y,ierr)
        call MPI_Type_commit(rcm_recvsubarray,ierr)
        rcm_extent=sizeof(realvar)
        rcm_start=0
        call MPI_Type_create_resized(rcm_recvsubarray,rcm_start,rcm_extent,rcm_re&
         &sizedrecvsubarray,ierr)
        call MPI_Type_commit(rcm_resizedrecvsubarray,ierr)
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
      allocate(rdnt(nbeadx3+3,ncols,nchain))
      if (sph_move) then
        allocate(rdnt_sph(3,ncols,nchain))
        allocate(rdnt_sph_or(3,ncols,nchain))
      endif
      allocate(qTot(nsegx3,npchain,p))
      allocate(qxT(npchain,nseg,p))
      allocate(qyT(npchain,nseg,p))
      allocate(qzT(npchain,nseg,p))
      allocate(RTot(nbeadx3,npchain,p))
      allocate(FphiTot(nbeadx3,npchain,p))
      if (CoM) allocate(rcmT(3,npchain,nchain_pp,p))
      if (CoHR) allocate(rchrT(3,npchain,p))
      if (sph_move) then
        allocate(r_sphT(3,npchain,p))
        allocate(p_sphT(3,npchain,p))
        allocate(q_sphT(3,npchain,p))
      end if

    end if ! id.eq.0

    !----------------------------------------------------------
    !>>>>> Allocating arrays:
    !----------------------------------------------------------

    ! Allocating the local random arrays
    allocate(rdn(nbeadx3+3,ncols,npchain))
    if (sph_move) then
      allocate(rdn_sph(3,ncols,npchain))
      allocate(rdn_sph_or(3,ncols,npchain))
    endif
    allocate (r_sph(3,npchain))
    allocate (p_sph(3,npchain))
    allocate (q_sph(3,npchain))
    allocate (F_sph(3))
    ! Note: the order nseg,nchain is selected because Fortran is column major.
    allocate (DiffTens(nbeadx3+3,nbeadx3+3,npchain),CoeffTens(nbeadx3+3,nbeadx3+3,npchain))
    allocate (qc(nsegx3),qstar(nsegx3),Fseg(nsegx3),wbl(nbeadx3+3,ncols),wbl_sph(3,ncols),wbl_sph_or(3,ncols))
    allocate (wbltemp(nbeadx3+3,ncols,npchain),w(nbeadx3+3),Kappa(nsegx3,nsegx3))
    allocate (Amat(nsegx3,nbeadx3),Kdotq(nsegx3),FBr(nsegx3),Amat_sph(nsegx3,nbeadx3+3))
    allocate (FBrbl(nsegx3,ncols,npchain),AdotD(nsegx3,nbeadx3+3,npchain))
    allocate (ADFev(nsegx3),Fev(nbeadx3),Fbead(nbeadx3),RHS(nsegx3),RHScnt(nsegx3),Fev_sph(3))
    allocate (RHSbase(nsegx3),Fbarseg(nsegx3),qbar(nsegx3),Fbarbead(nbeadx3))
    allocate (q(nsegx3,npchain),qstart(nsegx3,npchain),Fphi(nbeadx3,npchain),Fphi_all(nbeadx3+3),Fphi_all_temp(nbeadx3+3))
    allocate (Bmat(nbeadx3,nsegx3),Fbarev(nbeadx3),KappaBF(2,nsegx3),Fbar(nbeadx3))
    allocate (qctemp(nsegx3),DdotF(nbeadx3+3),Fbnd(nbeadx3),Fbarbnd(nbeadx3))
    select case (tplgy)
    case ('Linear')
      allocate(AmatBF(4,nbeadx3))
    case ('Comb')
    end select
    allocate(rvmrc(nbeadx3,npchain))

    if (sph_flow) allocate(U_seg(nsegx3))
    if (sph_flow) allocate(U_bead(nbeadx3))


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
    if (CoM) allocate(rcm(3,npchain,nchain_pp),rcmstart(3,npchain,nchain_pp))
    if (CoHR) allocate(rchr(3,npchain),rchrstart(3,npchain),&
     MobilTens(nbeadx3,nbeadx3),WeightTens&
     &(3,nbeadx3),weightTenstmp(3,nbeadx3))
    if (srf_tet) then
      allocate(rf0(3,npchain,nchain_pp))
      allocate(rf_in(3,npchain,nchain_pp))
      allocate(Ftet(3*nchain_pp))
      allocate(Fbartet(3*nchain_pp))
      ! call read_input('rf-in',0,rf_in(1),def=0._wp)
      ! call read_input('rf-in',1,rf_in(2),def=0._wp)
      ! call read_input('rf-in',2,rf_in(3),def=0._wp)
    endif
    if (unif_flow) then
      call read_input('U-Unif',0,U_unif(1),def=0._wp)
      call read_input('U-Unif',1,U_unif(2),def=0._wp)
      call read_input('U-Unif',2,U_unif(3),def=0._wp)
    endif
    ! For making the output
    allocate (q_counts(p),q_disps(p))
    allocate (R_counts(p),R_disps(p))
    if (CoHR) then
      allocate(rc_counts(p),rc_disps(p))
    end if
    if (sph_move) then
      allocate(rsph_counts(p),rsph_disps(p))
      allocate(psph_counts(p),psph_disps(p))
      allocate(qsph_counts(p),qsph_disps(p))
    end if
    if (CoM) then
      allocate(rcm_counts(p),rcm_disps(p))
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
      if (CoHR) then
        rc_counts(jp)=1
        rc_disps(jp)=(jp-1)*3*npchain
      end if
      if (sph_move) then
        rsph_counts(jp)=1
        rsph_disps(jp)=(jp-1)*3*npchain
        psph_counts(jp)=1
        psph_disps(jp)=(jp-1)*3*npchain
        qsph_counts(jp)=1
        qsph_disps(jp)=(jp-1)*3*npchain
      end if
      if (CoM) then
        rcm_counts(jp)=1
        rcm_disps(jp)=(jp-1)*3*npchain*nchain_pp
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

    if (sph_move) then
      call mysphsde%init()
    end if
    call mysde%init(Kappareg,Kappa,Amat,Amat_sph,Bmat,KappaBF,AmatBF,nbead_bb,nseg_bb)

    call print_matrix(Amat_sph,'Amat_sph from dlt')

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
            if (CoHR) read(u4,*)
          end do
          do iread=1, id*npchain*nchain_pp
            if (CoM) read(u3,*)
          end do
        end if
        if (id == ip-1) then
          do ichain=1, npchain
            do iseg=1, nseg
              offset=3*(iseg-1)
              read(u2,*) qstart(offset+1:offset+3,ichain)
            end do
            if (CoM) then
              do ichain_pp=1,nchain_pp
                read(u3,*) rcmstart(1:3,ichain,ichain_pp)
              end do
            end if
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
          if (srf_tet) then
            open(newunit=u_rf_in,file='data/rfin.st.dat',status='old',position='rewind')
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
        do iread=1, id*npchain*nchain_pp
          if (srf_tet) read(u_rf_in,*)
        end do
        do iread=1, id*npchain
          if (CoHR) read(u4,*)
        end do
        do iread=1, id*npchain*nchain_pp
          if (CoM) read(u3,*)
        end do
      end if
      if (id == ip-1) then
        do ichain=1, npchain
          do iseg=1, nseg
            offset=3*(iseg-1)
            read(u2,*) qstart(offset+1:offset+3,ichain)
          end do
          if (srf_tet) then
            do ichain_pp=1,nchain_pp
              read(u_rf_in,*) rf_in(1:3,ichain,ichain_pp)
            end do
          end if
          if (CoM) then
            do ichain_pp=1,nchain_pp
              read(u3,*) rcmstart(1:3,ichain,ichain_pp)
            end do
          end if
          if (CoHR) read(u4,*) rchrstart(1:3,ichain)
        end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end do
    ! if (srf_tet) then
    !   do ichain_pp=1, nchain_pp
    !     read(u_rf_in,*) rf_in(1:3,ichain_pp)
    !   end do
    ! end if
  end if
  close (u2);if (CoM) close (u3);if (CoHR) close(u4);if (srf_tet) close(u_rf_in)

  ! do ip=1,p
  !   if (id == ip-1) then
  !     print *, '---------------'
  !     print *, 'id = ',id
  !     do ichain=1,npchain
  !       call print_matrix(rf_in(1:3,ichain,:),'rf_in(1:3,ichain,:)')
  !       call print_matrix(rcmstart(1:3,ichain,:),'rcmstart(1:3,ichain,:)')
  !       call print_vector(qstart(:,ichain),'qstart(:,ichain)')
  !     end do
  !
  !     print *, '---------------'
  !   end if
  !   call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ! end do
  !call print_matrix(rf_in,'rf_in')
  !call print_matrix(rcmstart(:,1,:),'rcmstart')
  !call print_vector(qstart(:,1),'qstart')
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
        if (debug_TYL) then
          print *, '----------------------------------'
          print *, 'itime = ', itime
          print *, 'time = ', time
          print *, '----------------------------------'
        end if
        if (id == 0) then
          ! Constructing a block of random numbers for ncols time steps by rank 0
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
            do ichain=1, nchain
              do icol=1, ncols
                if (sph_move) then
                  rdnt_sph(1,icol,ichain)=rangls()
                  rdnt_sph(2,icol,ichain)=rangls()
                  rdnt_sph(3,icol,ichain)=rangls()
                  rdnt_sph_or(1,icol,ichain)=rangls()
                  rdnt_sph_or(2,icol,ichain)=rangls()
                  rdnt_sph_or(3,icol,ichain)=rangls()
                end if

                do ibead=1, 3*nbead+3
                  rdnt(ibead,icol,ichain)=ranuls()-0.5
                end do
              end do
            end do
          end if
        end if
        ! Scattering the generated random numbers from rank 0 the owner.
        if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
          call MPI_Scatter(rdnt,3*(nbead+1)*ncols*npchain,MPI_REAL_WP,rdn,&
           3*(nbead+1)*ncols*npchain,MPI_REAL_WP,0,MPI_COMM_WORLD,ierr)
          if (sph_move) then
            call MPI_Scatter(rdnt_sph,3*ncols*npchain,MPI_REAL_WP,rdn_sph,&
             3*ncols*npchain,MPI_REAL_WP,0,MPI_COMM_WORLD,ierr)
            call MPI_Scatter(rdnt_sph_or,3*ncols*npchain,MPI_REAL_WP,rdn_sph_or,&
             3*ncols*npchain,MPI_REAL_WP,0,MPI_COMM_WORLD,ierr)
          end if
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
              do ichain_pp=1,nchain_pp
                !rf0(:,jchain)=rf_in(1:3)
                rf0(:,jchain,ichain_pp) = rf_in(1:3,jchain,ichain_pp)
              end do
            end if
          end do

          ! do ip=1,p
          !   if (id == ip-1) then
          !     print *, '---------------'
          !     print *, 'id = ',id
          !     do ichain=1,npchain
          !       call print_matrix(rf0(1:3,ichain,:),'rf0(1:3,ichain,:)')
          !     end do
          !     print *, '---------------'
          !   end if
          !   call MPI_Barrier(MPI_COMM_WORLD,ierr)
          ! end do


          !call print_matrix(rf0(:,1,:),'rf0')
          r_sph=0
          p_sph=0
          p_sph(1,:) = 1._wp
          q_sph=0
          q_sph(2,:) = 1._wp
          F_sph=0

          if (sph_move) then
            Ftet=0
          end if
          if (CoM) rcm=rcmstart
          if (CoHR) rchr=rchrstart
          call pp_init_tm(id)
          istr=1 ! the index for dumpstr()
        end if

        !call print_vector(rvmrc(:,1),'rvmrc')

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
              do ibead=1, nbead+1
                offset=3*(ibead-1)
                wx=rdn(offset+1,kcol,ichain);wx=sqrtdt*wx*(c1*wx**2+c2)
                wy=rdn(offset+2,kcol,ichain);wy=sqrtdt*wy*(c1*wy**2+c2)
                wz=rdn(offset+3,kcol,ichain);wz=sqrtdt*wz*(c1*wz**2+c2)
                if ((MOD(ibead-1,nbead_ind)==0) .and. (ibead < nbead+1)) then !tethered beads have no brownian motion
                  wbl(offset+1,kcol)=0._wp;wbl(offset+2,kcol)=0._wp;wbl(offset+3,kcol)=0._wp
                else
                  wbl(offset+1,kcol)=wx;wbl(offset+2,kcol)=wy;wbl(offset+3,kcol)=wz
                end if
              end do
              if (sph_move) then
                wbl_sph(1:3,kcol) = sqrtdt*rdn_sph(1:3,kcol,ichain)
                wbl_sph_or(1:3,kcol) = sqrtdt*rdn_sph_or(1:3,kcol,ichain)
              end if
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
            rcmP => rcm(:,ichain,:)
            !call print_matrix(rcmP(:,:),'rcmP')
          end if
          r_sphP => r_sph(:,ichain)
          p_sphP => p_sph(:,ichain)
          q_sphP => q_sph(:,ichain)

          ! if (EV_bw == 'Rflc_bc') then
          !   call wall_rflc(rvmrcPy,rcmP(2))
          ! end if


          call sprforce(id,qc,nseg,ForceLaw,TruncMethod,Fseg)

          if (debug_TYL) then
            call print_vector(Fseg(:),'Fseg(922): ')
            call print_vector(qc(:),'qc(923): ')
          end if
          ! !TYL - correcting the force in the tethered springs---
          ! do ichain_pp=1,nchain_pp
          !   offset = nseg_indx3*(ichain_pp-1)
          !
          !   rf0_rel_unit(:) = rf0(:,ichain,ichain_pp)  - r_sph(:,ichain)
          !   rf0_rel_unit(:) = rf0_rel_unit(:) / (sqrt(rf0_rel_unit(1)**2+ &
          !     rf0_rel_unit(2)**2 + rf0_rel_unit(3)**2))
          !
          !   Fseg(offset+1:offset+3) = Fseg(offset+1:offset+3) - &
          !     rf0_rel_unit(1:3)/dot(rf0_rel_unit(1:3),qc(offset+1:offset+3))
          ! end do
          ! !------------------------------------------------------

          !call print_vector(Fseg(:),'Fseg after: ')

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

  ! if (ichain==1) then
  !   if (itime>4250 .and. mod(itime,5)==0) then
  !     print*,'itime',itime
  !     call print_vector(rvmrcP,'rvmrc')
  !     read(*,*)
  !   endif
  ! endif
              !print *, 'before first calc--------------------------------'

              call print_matrix(DiffTensP,'DiffTensP before calc:')

              call myintrn%calc(id,itime,rvmrcP,rcmP,r_sphP,nseg,DiffTensP,divD,Fev,Fbarev,Fev_sph,&
                calchi=.true.,calcdiv=.true.,calcevbb=.true.,calcevbw=.true.)

              call print_matrix(DiffTensP,'DiffTensP after calc:')

              if (debug_TYL) then
                call print_vector(r_sphP,'r_sphP to calc D')
                call print_vector(rvmrcP(1:9),'rvmrcP 1 to calc D')
                call print_vector(rcmP(:,1),'rcmP 1 to calc D')
                call print_vector(rvmrcP(10:18),'rvmrcP 2 to calc D')
                call print_vector(rcmP(:,2),'rcmP 2 to calc D')
                call print_matrix(DiffTensP,'DiffTensP (966):')
              end if
              !print *, 'after first calc--------------------------------'

            end if
            if (DecompMeth == 'Cholesky') then
              if (hstar /= 0._wp) then
                CoeffTensP=real(DiffTensP,kind=double)
                wbltempP1=real(wbl,kind=double)

                !! I need to change the storage scheme for symmetric tensors
                !! for dilute solutions in the absence of the wall
                !print *, 'before Cholesky'
                if (HITens == 'Blake') then
                  call potrf(CoeffTensP,info=info)
                else



                  call print_matrix(CoeffTensP,'CoeffTensP before potrf (DiffTensP)')
                  call potrf(CoeffTensP,info=info)
                  call print_matrix(CoeffTensP,'CoeffTensP after potrf')

                endif
                !print *, 'after Cholesky'
                if (info /= 0) then
                  print '(" Unsuccessful Cholesky fact. of D in main.")'
                  print '(" info: ",i3)',info
                  stop
                end if
                !print *, 'after Cholesky check'


                if (HITens == 'Blake') then
                  call trmm(CoeffTensP,wbltempP1,transa='T')
                else
                  call print_matrix(wbltempP1,'wbltempP1 (dW)')
                  call trmm(CoeffTensP,wbltempP1,transa='T')
                  call print_matrix(wbltempP1,'wbltempP1 (C*dW)')
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
                !AmatBF is not right
                !call gbmv(AmatBF,real(wbltempP2,kind=wp),FBrblP,kl=0,m=nsegx3,alpha=coeff)
                call print_vector(wbltempP2,'wbltempP2 (C*dW)')
                call gemv(Amat_sph,real(wbltempP2,kind=wp),FBrblP,alpha=coeff)
                call print_vector(FBrblP,'FBrblP ((1/sqrt(2))*A*C*dW)')

                ! !TYL: HI for tethered bead -------------------------------------
                ! !print *, 'Brownian force calculaton------------'
                ! if (srf_tet) then
                !   do ichain_pp=1,nchain_pp !(ichain_pp-1)*nbead_indx3+1
                !     !nseg_indx3*(ichain-1)+1:nseg_indx3*ichain,nbead_indx3*(ichain-1)+1:nbead_indx3*ichain
                !     call gemv(Amat(:,(ichain_pp-1)*nbead_indx3+1:(ichain_pp-1)*nbead_indx3+3),&
                !       real(wbltempP2((ichain_pp-1)*nbead_indx3+1:(ichain_pp-1)*nbead_indx3+3),&
                !       kind=wp),FBrblP,alpha=-coeff,beta=1._wp)
                !   end do
                ! end if
                ! !TYL: HI for tethered bead -------------------------------------

              else
                call gemv(Amat,real(wbltempP2,kind=wp),FBrblP,alpha=coeff)
              end if
            end do
            ! Calculation of AdotD=Amat.D, to be used in Predictor-Corrector
            if (hstar /= 0._wp) then
              call symm(DiffTensP,Amat_sph,AdotDP1,side='R')

              call print_matrix(AdotDP1,'AdotDP1 (A*D):')

            else
              AdotDP1=Amat
            end if
          end if
          FBr=FBrbl(:,jcol,ichain)


          ! rflc doesn't need this
          !            if (hstar == 0._wp .and. EV_bw /= 'Rflc_bc') then
          if (hstar == 0._wp) then
            !                call EVCalc(rvmrcP,nseg,EV_bb,Fev)
            call myintrn%calc(id,itime,rvmrcP,rcmP,r_sphP,nseg,DiffTensP,divD,Fev,Fbarev,Fev_sph,&
              calcevbb=.true.,calcevbw=.true.)
            !call evcalc2(rvmrcP,nseg,Fev)
          end if

          !------------ advancing the configuration -------------------------!
          if (sph_move) then
            call mysphsde%advance(r_sph,p_sph,q_sph,rf0,iPe,idt,ichain,Fseg,Fbead,wbl_sph,wbl_sph_or,jcol,&
            Fev_sph,F_sph,Fev,Fbnd,DiffTensP,Amat,Fphi,Fphi_all,Fphi_all_temp,wbltempP1,coeff)

            if (debug_TYL) then
              call print_vector(r_sph(:,1),'r_sph (1100) = ')
              call print_vector(rf0(:,1,1),'rf0(:,1,1) (1101) = ')
              call print_vector(rf0(:,1,2),'rf0(:,1,2) (1101) = ')
            end if
            !print *, '------------------------'
            !call print_vector(Ftet(:),'Ftet = ')
            !call print_vector(qc, 'Spring vector')
            !call mysphsde%advance(r_sph,rf0,iPe,idt,ichain,Ftet,wbl_sph,jcol)
          end if

          call mysde%advance(myintrn,id,iPe,idt,ichain,itime,Kdotq,qc,Fseg,Fbead,Fev,Fbnd,qstar,Fphi,&
           rvmrcP,rcm,DiffTensP,Ftet,rf0,AdotDP1,divD,FBr,RHS,rcmP,r_sphP,Fbarev,nbead_bb,Fbarbnd,Fbar,Fbartet,&
           RHScnt,Fbarseg,Fbarbead,root_f,qbar,nseg_bb,AdotD,RHSbase,qctemp,mch,Lch,lambdaBE,r_sph,F_sph,Fphi_all)
          !------------ advancing the configuration -------------------------!


          ! Inserting back the final result to original arrays
          q(:,ichain)=qc(:)

          ! commented below line to test Euler scheme 3/11/19
          !Fphi(:,ichain)=Fbead(:)+Fbar(:)
          !note: Fbar includes Fbarev + Fbarbnd + Fbartet
          !      Fbead includes the spring forces

          !call print_vector(Fphi(:,ichain),'Fphi(after sde advance)')

          call gemv(Bmat,qc,rvmrcP)
          ! Calculating center of mass and/or center of hydrodynamic resistance movement
          if (CoM) then

            print*, '----------------'
            print*, 'calculation of CoM'
            print*, '----------------'

            call print_vector(Fphi(:,ichain),'Fphi(:,ichain)')
            call print_vector(Fphi_all(:),'Fphi_all(:)')
            FphiP => Fphi(:,ichain)
            rcmP => rcm(:,ichain,:) ! might be redundant
            ! BdotwPx => wbltemp(1:nbeadx3-2:3,jcol,ichain)
            ! BdotwPy => wbltemp(2:nbeadx3-1:3,jcol,ichain)
            ! BdotwPz => wbltemp(3:nbeadx3:3,jcol,ichain)
            ! SumBdotw=[sum(real(BdotwPx,kind=wp)),&
            ! sum(real(BdotwPy,kind=wp)),&
            ! sum(real(BdotwPz,kind=wp))]

            if (hstar.ne.0._wp) then
              call symv(DiffTensP,Fphi_all(:),DdotF)


              call print_vector(DdotF,'DdotF')

            else
              DdotF=Fphi_all(:)
            end if

            do ichain_pp=1,nchain_pp



              !TYL: HI for tethered bead ---------------------------------------
              if (srf_tet) then
                !call gemv(DiffTensP(:,1:3),FphiP(1:3),DdotF,alpha=-1._wp,beta=1._wp)
                call gemv(DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,:),FphiP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),DdotF,alpha=-1._wp,&
                  beta=1._wp,trans='T')

                if (debug_TYL) then
                  call print_matrix(DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,:),'DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,:)')
                  call print_vector(FphiP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),'FphiP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3)')
                end if


                !TYL: adding back self-mobility of first bead-------------------------
                call gemv(DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,&
                                    nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),&
                          FphiP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),&
                          DdotF(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),&
                          alpha=1._wp,beta=1._wp)

                if (debug_TYL) then
                  call print_matrix(DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,&
                                    nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3),'DiffTensP(nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3,nbead_indx3*(ichain_pp-1) + 1:nbead_indx3*(ichain_pp-1) + 3)')
                end if
                !TYL: adding back self-mobility of first bead-------------------------


                call print_vector(DdotF,'fixing DdotF')

              end if
              !TYL: HI for tethered bead ---------------------------------------


              if (unif_flow) then
                rcmP(:,ichain_pp)=rcmP(:,ichain_pp)+U_unif*dt(iPe,idt)
              endif

            end do

            do ichain_pp=1,nchain_pp
              BdotwPx => wbltemp(nbead_indx3*(ichain_pp-1) + 1 : nbead_indx3*ichain_pp - 2 : 3,jcol,ichain)
              BdotwPy => wbltemp(nbead_indx3*(ichain_pp-1) + 2 : nbead_indx3*ichain_pp - 1 : 3,jcol,ichain)
              BdotwPz => wbltemp(nbead_indx3*(ichain_pp-1) + 3 : nbead_indx3*ichain_pp : 3,jcol,ichain)
              SumBdotw=[sum(real(BdotwPx,kind=wp)),&
              sum(real(BdotwPy,kind=wp)),&
              sum(real(BdotwPz,kind=wp))]

              DdotFPx => DdotF(nbead_indx3*(ichain_pp-1) + 1 : nbead_indx3*ichain_pp - 2 : 3)
              DdotFPy => DdotF(nbead_indx3*(ichain_pp-1) + 2 : nbead_indx3*ichain_pp - 1 : 3) !(2 : nbeadx3-1:3)
              DdotFPz => DdotF(nbead_indx3*(ichain_pp-1) + 3 : nbead_indx3*ichain_pp : 3)
              SumDdotF=(/sum(DdotFPx),sum(DdotFPy),sum(DdotFPz)/)

              rcmP(:,ichain_pp)=rcmP(:,ichain_pp)+(Pe(iPe)*matmul(Kappareg,rcmP(:,ichain_pp))&
                +1._wp/(4*nbead_ind)*SumDdotF)*dt(iPe,idt)+coeff/nbead_ind*SumBdotw
            end do

            if (sph_flow) then
              !print *, 'before U_bead(5) = ', U_bead(5)
              call mysde%U_sph(U_seg,U_bead,q(:,ichain),rcm(:,ichain,:))
              !print *, 'after U_bead(5) = ', U_bead(5)
              ! call print_vector(U_bead(:),'in CoM calc: U_bead')
              ! call print_vector(U_seg(:),'in CoM calc: U_seg')
              ! call print_vector(q(:,ichain),'in CoM calc: q(:,ichain)')
              ! call print_vector(rcm(:,ichain,1),'in CoM calc: rcm(:,ichain,1)')
              ! call print_vector(rcm(:,ichain,2),'in CoM calc: rcm(:,ichain,2)')

              do ichain_pp=1,nchain_pp
                !call print_vector(rcmP(:,ichain_pp),'rcmP(:,ichain_pp) before')
                rcmP(1,ichain_pp) = rcmP(1,ichain_pp) + (1._wp/nbead_ind)*&
                  sum(U_bead(nbead_indx3*(ichain_pp-1) + 1 : nbead_indx3*ichain_pp - 2 : 3))*dt(iPe,idt)
                rcmP(2,ichain_pp) = rcmP(2,ichain_pp) + (1._wp/nbead_ind)*&
                  sum(U_bead(nbead_indx3*(ichain_pp-1) + 2 : nbead_indx3*ichain_pp - 1 : 3))*dt(iPe,idt)
                rcmP(3,ichain_pp) = rcmP(3,ichain_pp) + (1._wp/nbead_ind)*&
                  sum(U_bead(nbead_indx3*(ichain_pp-1) + 3 : nbead_indx3*ichain_pp : 3))*dt(iPe,idt)
                !call print_vector(rcmP(:,ichain_pp),'rcmP(:,ichain_pp) after')
              end do
            endif

            !! Blake's part
            if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
              do ichain_pp=1,nchain_pp
               rcmP(2,ichain_pp)=rcmP(2,ichain_pp)+1._wp/(4*nbead_ind)*&
                sum(divD(nbead_ind*(ichain_pp-1)+1:nbead_ind*(ichain_pp)))*dt(iPe,idt)
              end do
            end if
            !!-------------

            call print_vector(rvmrcP,'rvmrcP (1229)')
            call print_matrix(rcmP,'rcmP (1230)')
            call print_vector(rvmrcP(1:3)+rcmP(:,1),'Tether point 1')
            call print_vector(rf0(:,1,1),'rf0 (1232)')
            call print_vector(rvmrcP(10:12)+rcmP(:,2),'Tether point 2')
            call print_vector(rf0(:,1,2),'rf0 (1234)')

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
            call wall_rflc_sph(myintrn%evbw,r_sph(:,ichain),rf0(:,ichain,:))

            dr_sph_rflc = 0._wp !total amount to move sphere and tether points

            do ichain_pp=1,nchain_pp
              idx_pp_seg = (nseg_ind) * (ichain_pp -1) + 1
              idx_pp_bead = (nbead_ind) * (ichain_pp -1) + 1
              call wall_rflc(myintrn%evbw,dt(iPe,idt),itime,time,id,ichain,&
                qPx(idx_pp_seg:idx_pp_seg+(nseg_ind-1)),&
                qPy(idx_pp_seg:idx_pp_seg+(nseg_ind-1)),&
                qPz(idx_pp_seg:idx_pp_seg+(nseg_ind-1)),&
                rvmrcPx(idx_pp_bead:idx_pp_bead+(nbead_ind-1)),&
                rvmrcPy(idx_pp_bead:idx_pp_bead+(nbead_ind-1)),&
                rvmrcPz(idx_pp_bead:idx_pp_bead+(nbead_ind-1)),&
                rcmP(1,ichain_pp),rcmP(2,ichain_pp),rcmP(3,ichain_pp),&
                rf0(:,ichain,ichain_pp),r_sph(:,ichain),dr_sph_rflc(:)) !rf_in(:,ichain_pp)
                !rf0(:,jchain,ichain_pp) = rf_in(1:3,ichain_pp)
            end do

            !!!!!!!!!!!!!TYL - move sphere due to collisions
            !if (sqrt(dr_sph_rflc(1)**2+dr_sph_rflc(2)**2+dr_sph_rflc(3)**2) > .000001_wp) then
            !  call print_vector(dr_sph_rflc(:),'dr_sph_rflc = ')
            !end if



            !move the sphere
            r_sph(:,ichain) = r_sph(:,ichain) + dr_sph_rflc(:)



            !fix the tether points and center of masses
            ! do ichain_pp=1,nchain_pp
            !   !move tether points
            !   rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph_rflc(:)
            !
            !   !recalculate the center of mass
            !   rcmP(1,ichain_pp) = rcmP(1,ichain_pp) + dr_sph_rflc(1)!*(nbead_ind-1)/(nbead_ind)
            !   rcmP(2,ichain_pp) = rcmP(2,ichain_pp) + dr_sph_rflc(2)!*(nbead_ind-1)/(nbead_ind)
            !   rcmP(3,ichain_pp) = rcmP(3,ichain_pp) + dr_sph_rflc(3)!*(nbead_ind-1)/(nbead_ind)
            ! end do
            !!!!!!!!!!!!!TYL - move sphere due to collisions

            !fix the tether points, tethered springs, and center of masses
            !THIS ALGORITHM WILL CHANGE THE TETHERED SPRING
            do ichain_pp=1,nchain_pp
              idx_pp_seg = (nseg_ind) * (ichain_pp -1) + 1
              idx_pp_bead = (nbead_ind) * (ichain_pp -1) + 1

              !move tether points
              rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph_rflc(:)

              !temporarily add center of mass to rvmrc to save memory
              rvmrcPx(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPx(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) + rcmP(1,ichain_pp)
              rvmrcPy(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPy(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) + rcmP(2,ichain_pp)
              rvmrcPz(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPz(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) + rcmP(3,ichain_pp)

              !shift tethered bead
              rvmrcPx(idx_pp_bead) = rvmrcPx(idx_pp_bead) + dr_sph_rflc(1)
              rvmrcPy(idx_pp_bead) = rvmrcPy(idx_pp_bead) + dr_sph_rflc(2)
              rvmrcPz(idx_pp_bead) = rvmrcPz(idx_pp_bead) + dr_sph_rflc(3)

              !correct tethered spring (other spring vectors don't change)
              qPx(idx_pp_seg) = rvmrcPx(idx_pp_bead+1) - rvmrcPx(idx_pp_bead)
              qPy(idx_pp_seg) = rvmrcPy(idx_pp_bead+1) - rvmrcPy(idx_pp_bead)
              qPz(idx_pp_seg) = rvmrcPz(idx_pp_bead+1) - rvmrcPz(idx_pp_bead)

              !recalculate the center of mass
              rcmP(:,ichain_pp) = 0._wp
              do ibead=1,nbead_ind
                rcmP(1,ichain_pp) = rcmP(1,ichain_pp) + rvmrcPx(idx_pp_bead+ibead-1)
                rcmP(2,ichain_pp) = rcmP(2,ichain_pp) + rvmrcPy(idx_pp_bead+ibead-1)
                rcmP(3,ichain_pp) = rcmP(3,ichain_pp) + rvmrcPz(idx_pp_bead+ibead-1)
              end do
              rcmP(:,ichain_pp) = rcmP(:,ichain_pp)/nbead_ind

              !recalculate the positions relative to the center of mass
              rvmrcPx(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPx(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) - rcmP(1,ichain_pp)
              rvmrcPy(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPy(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) - rcmP(2,ichain_pp)
              rvmrcPz(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) = rvmrcPz(idx_pp_bead:idx_pp_bead+(nbead_ind-1)) - rcmP(3,ichain_pp)
            end do
            !!!!!!!!!!!!!TYL - move sphere due to collisions


          end if

          !call print_matrix(rcmP(:,:),'rcmP(:,:) after (1300)')
          !-----

        end do ! ichain loop

        !----------------------------------------------------------------
        !>>>>> data processing and outputs:
        !----------------------------------------------------------------

        !time_check3 is postprocessing
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

        !time_check4 = restart
        if (time >= time_check4) then
          time_check4=time_check4+frm_rt_rst*lambda
          ! For writing restart data to the output file
          call MPI_Gatherv(q,nsegx3*npchain,MPI_REAL_WP,qTot,q_counts,q_disps,&
           q_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          if (CoM) then
            !print *, 'before1 CoM MPI_Gatherv'
            call MPI_Gatherv(rcm,3*npchain*nchain_pp,MPI_REAL_WP,rcmT,rcm_counts,rcm_disps,&
             rcm_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            !print *, 'after1 CoM MPI_Gatherv'
          end if
          if (CoHR) then
            call MPI_Gatherv(rchr,3*npchain,MPI_REAL_WP,rchrT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (sph_move) then
            call MPI_Gatherv(r_sph,3*npchain,MPI_REAL_WP,r_sphT,rsph_counts,rsph_disps,&
             rsph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(p_sph,3*npchain,MPI_REAL_WP,p_sphT,psph_counts,psph_disps,&
             psph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(q_sph,3*npchain,MPI_REAL_WP,q_sphT,qsph_counts,qsph_disps,&
             qsph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if

          if (id == 0) then
            open(newunit=u21,file='data/q.rst.dat',status='replace')
            if (CoM) open(newunit=u22,file='data/CoM.rst.dat',status='replace')
            if (CoHR) open(newunit=u23,file='data/CoHR.rst.dat',status='replace')
            if (sph_move) then
              open(newunit=u_rsph_rst,file='data/rsph.rst.dat',status='replace')
              open(newunit=u_pqsph_rst,file='data/pqsph.rst.dat',status='replace')
            endif
            do ip = 1, p
              do ichain =1, npchain
                do iseg = 1, nseg
                  offset=3*(iseg-1)
                  write(u21,1) qTot(offset+1:offset+3,ichain,ip)
                end do
                if (CoM) then
                  do ichain_pp = 1,nchain_pp
                    write(u22,1) rcmT(1:3,ichain,ichain_pp,ip)
                  end do
                end if
                if (CoHR) write(u23,1) rchrT(1:3,ichain,ip)
                if (sph_move) then
                  write(u_rsph_rst,1) r_sphT(1:3,ichain,ip)
                  write(u_pqsph_rst,1) p_sphT(1:3,ichain,ip)
                  write(u_pqsph_rst,1) q_sphT(1:3,ichain,ip)
                endif
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
            if (sph_move) then
              write(u_rsph_rst,'(f7.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed";close(u_rsph_rst)
              write(u_pqsph_rst,'(f7.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed";close(u_pqsph_rst)
            end if
          end if ! id == 0
        end if ! mod(itime,lambda/dt)==0



        ! if (time >= time_check2) then

        !dumping data
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
            !print *, 'before2 CoM MPI_Gatherv'
            call MPI_Gatherv(rcm,3*npchain*nchain_pp,MPI_REAL_WP,rcmT,rcm_counts,rcm_disps,&
              rcm_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            !print *, 'after2 CoM MPI_Gatherv'
          end if
          if (CoHR) then
            call MPI_Gatherv(rchr,3*npchain,MPI_REAL_WP,rchrT,rc_counts,rc_disps,&
             rc_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
          end if
          if (sph_move) then
            call MPI_Gatherv(r_sph,3*npchain,MPI_REAL_WP,r_sphT,rsph_counts,rsph_disps,&
             rsph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(p_sph,3*npchain,MPI_REAL_WP,p_sphT,psph_counts,psph_disps,&
             psph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(q_sph,3*npchain,MPI_REAL_WP,q_sphT,qsph_counts,qsph_disps,&
             qsph_resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
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
                if (sph_move) then
                  open(newunit=u_rsph_flow,file='data/rsph.flow.dat',status=fstat,position='append')
                  open(newunit=u_pqsph_flow,file='data/pqsph.flow.dat',status=fstat,position='append')
                endif

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
                if (CoM) then
                  do ichain_pp = 1,nchain_pp
                    write(u26,1) rcmT(1:3,ichain,ichain_pp,ip)
                  end do
                end if
                if (CoHR) write(u27,1) rchrT(1:3,ichain,ip)
                if (sph_move) then
                  write(u_rsph_flow,1) r_sphT(1:3,ichain,ip)
                  write(u_pqsph_flow,1) p_sphT(1:3,ichain,ip)
                  write(u_pqsph_flow,1) q_sphT(1:3,ichain,ip)
                endif
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
    if (sph_move) then
      close(u_rsph_flow)
      close(u_pqsph_flow)
    endif
    if (DumpstrCalc) close(u40)
    if (cnf_srt) close(u41)
  end if

  if (id == 0) then
    deallocate(rdnt)
    if (sph_move) then
      deallocate(rdnt_sph)
      deallocate(rdnt_sph_or)
    endif
    deallocate (qTot,qxT,qyT,qzT,RTot,FphiTot)
    if (CoM) deallocate(rcmT)
    if (CoHR) deallocate(rchrT)
    if (sph_move) then
      deallocate(r_sphT)
      deallocate(p_sphT)
      deallocate(q_sphT)
    endif
  end if

  deallocate(rdn)
  deallocate(r_sph)
  deallocate(p_sph)
  deallocate(q_sph)
  if (sph_move) then
    deallocate(rdn_sph)
    deallocate(rdn_sph_or)
  endif
  deallocate(qc,qstar,Fseg,w,wbl,wbltemp,Kappa,Amat,FBr,FBrbl,Kdotq,AdotD,Amat_sph)
  if (sph_move) then
    deallocate(wbl_sph)
    deallocate(wbl_sph_or)
  endif
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
  if (CoHR) deallocate(rc_counts,rc_disps)
  if (sph_move) deallocate(rsph_counts,rsph_disps,psph_counts,psph_disps,qsph_counts,qsph_disps)
  if (CoM) deallocate(rcm_counts,rcm_disps)
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

!   ! Random numeber seeding (from H. C. Ottinger):
!   subroutine ranils(iseed)
!
!     integer,intent(in) :: iseed
!     integer,parameter :: in=2147483563,ik=40014,iq=53668,ir=12211,ntab=32
!     integer :: iv(ntab),idum,idum2,iy
!     integer :: k,j
!
!     common /ranbls/ idum,idum2,iy,iv
!
!     ! Initial seeds for two random number generators
!     idum=iseed+123456789
!     idum2=idum
!
!     ! Load the shuffle table (after 8 warm-ups)
!     do 10 j=ntab+8,1,-1
!      k=idum/iq
!      idum=ik*(idum-k*iq)-k*ir
!      if(idum < 0) idum=idum+in
!      if(j <= ntab) iv(j)=idum
!      10 continue
!      iy=iv(1)
!      return
!
!    end subroutine ranils
!
!    ! Uniform random number generator (from H. C. Ottinger):
!    real(wp) function ranuls()
!
!    integer,parameter :: in1=2147483563,ik1=40014,iq1=53668,ir1=12211,&
!    in2=2147483399,ik2=40692,iq2=52774,ir2=3791 ,&
!    ntab=32,inm1=in1-1,ndiv=1+inm1/ntab
!    real(wp),parameter :: an=1./in1
!    integer :: iv(ntab),idum,idum2,iy
!    integer :: k,j
!
!    common /ranbls/ idum,idum2,iy,iv
!
!    ! Linear congruential generator 1
!    k=idum/iq1
!    idum=ik1*(idum-k*iq1)-k*ir1
!    if(idum < 0._wp) idum=idum+in1
!
!    ! Linear congruential generator 2
!    k=idum2/iq2
!    idum2=ik2*(idum2-k*iq2)-k*ir2
!    if(idum2 < 0._wp) idum2=idum2+in2
!
!    !Shuffling and subtracting
!    j=1+iy/ndiv
!    iy=iv(j)-idum2
!    iv(j)=idum
!    if(iy < 1) iy=iy+inm1
!    ranuls=an*iy
!    return
!
!  end function ranuls
!
!  ! Gaussian random number generator (from H. C. Ottinger):
!  real(wp) function rangls()
!
!  integer :: iflag
!  real(wp) :: gauss2,x1,x2,xsq,aux
!
!  save iflag,gauss2
!  data iflag/0/
!
!  if(iflag == 0) then
!   10 continue
!
!   ! pair of uniform random numbers in [-1,1]x[-1,1]
!   x1=2*ranuls()-1
!   x2=2*ranuls()-1
!
!   ! if not in the unit circle, try again
!   xsq=x1*x1+x2*x2
!   if(xsq >= 1._wp .or. xsq == 0._wp) goto 10
!   ! pair of gaussian random numbers; return one and
!   ! save the other for next time
!   aux=sqrt(-2*log(xsq)/xsq)
!   rangls=x1*aux
!   gauss2=x2*aux
!   iflag=1
! else
!   rangls=gauss2
!   iflag=0
! endif
! return
!
! end function rangls

end subroutine dlt_bs

end module dlt_mod
