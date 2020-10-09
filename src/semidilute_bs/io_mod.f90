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
! MODULE: Input/Output
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION:
!> For reading and writing the configurational information
!
!--------------------------------------------------------------------

module io_mod

  use :: prcn_mod

  implicit none

  private :: init_conf_io  ,&
             read_conf     ,&
             read_init_conf,&
             read_dmp_conf ,&
             write_conf    ,&
             del_io


  !> A public type for reading and writing configurational info
  type conf_io

    private
    !> From input file: initial simulation status st, ext, rst
    character(len=10),public :: initmode
    !> From input file: if making animation is intended?
    logical :: MakeAnim
    !> From input file: the starting fractional segmental connectivity
    real(wp) :: qfctr(3)
    !> From input file: the starting fractional segmental connectivity
    real(wp) :: qfctr_cmbbb(3)
    !> From input file: the starting fractional segmental connectivity
    real(wp) :: qfctr_cmbar(3)
    !> mpi information for connectivity vector
    integer :: q_recvsubarray
    !> mpi information for center of mass
    integer :: rc_recvsubarray
    !> mpi information for center of mass
    integer :: cif_recvsubarray
    !> mpi information for position of the beads
    integer :: rb_recvsubarray
    !> mpi information for interparticle distance times total force
    integer :: rf_recvsubarray
    !> The handle for connectivity vector start file
    integer :: fqsthandle
    !> The handle for center of mass start file
    integer :: fcsthandle
    !> The handle for the position of the beads start file
    integer :: frbsthandle
    !> The handle for connectivity vector restart file
    integer :: fqrsthandle
    !> The handle for center of mass restart file
    integer :: fcrsthandle
    !> The handle for center of mass image flag restart file
    integer :: fcifrsthandle
    !> The handle for center of mass initial file
    integer :: fcinithandle
    !> The handle for the position of the beads restart file
    integer :: frbrsthandle
    !> The handle for connectivity vector dump file
    integer :: fqdmphandle
    !> The handle for center of mass dump file
    integer :: fcdmphandle
    !> The handle for center of mass image flag dump file
    integer :: fcifdmphandle
    !> The handle for the position of particles  dump file
    integer :: frbdmphandle
    !> The handle for the interparticle distance times total force dump file
    integer :: frfdmphandle
    !> The unit for q.rst.dat
    integer :: oldu1
    !> The unit for CoM.rst.dat
    integer :: oldu2
    !> The unit for q.equil.dat or q.final.dat
    integer :: oldu3
    !> The unit for CoM.equil.dat or CoM.final.dat
    integer :: oldu4
    !> The unit for Rb.dat
    integer :: oldu5
    !> The unit for CoM.dat
    integer :: oldu6
    !> The unit for BoxConfig.dat
    integer :: oldu7
    !> The unit for Rb.rst.dat
    integer :: oldu8
    !> The unit for Rb.equil.dat or Rb.final.dat
    integer :: oldu9

  contains

    procedure,pass(this) :: init => init_conf_io
    procedure,pass(this) :: read => read_conf
    procedure,pass(this) :: read_init => read_init_conf
    procedure,pass(this) :: read_dmp => read_dmp_conf
    procedure,pass(this) :: write => write_conf
    final :: del_io

  end type conf_io

  ! Private module variables:
  private :: DumpConf
  ! Protected module variables:
  protected :: CoMDiff

  !> From input file: if dumping configurational info is intended?
  logical,save :: DumpConf
  !> From input file: if calculation of diffusivity is intended?
  logical,save :: CoMDiff
    !> The starting bead position vector for all chains
  real(wp),allocatable,dimension(:,:),target,save :: Rbst
  ! !> The starting connectivity vector for all chains
  ! real(wp),allocatable,dimension(:,:,:),target,save :: Qst
  !> The starting connectivity vector for all chains
  real(wp),allocatable,dimension(:,:),target,save :: Qst
  !> The starting center of mass for all chains
  real(wp),allocatable,dimension(:,:,:),target,save :: rcmst
  !> The initial center of mass for all chains
  real(wp),allocatable,dimension(:,:,:),target,save :: rcminit
  !> The starting image flag of the center of mass for all chains
  integer,allocatable,dimension(:,:,:),target,save :: cmifst


contains

  !> Initializes io_mod module variables
  subroutine init_io(id,nchain,nsegx3,nbeadx3,ntotchain,ntotsegx3,ntotbeadx3,nprun)

    use,intrinsic :: iso_fortran_env
    use :: strg_mod
    use :: flow_mod, only: FlowType

    integer,intent(in) :: id,nsegx3,nbeadx3,nchain,ntotchain,ntotsegx3,ntotbeadx3,nprun
    integer :: i,j,ntokens,u1,stat,ios,il
    character(len=1024) :: line
    character(len=100) :: tokens(10)

    ! default values:
    DumpConf=.false.
    CoMDiff=.false.

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
            case ('DumpConf')
              if(tokens(j+1) == 'TRUE') then
                DumpConf=.true.
              elseif(tokens(j+1) == 'FALSE') then
                DumpConf=.false.
              else
                print *,'Incorrect Type for DumpConf.'
                stop
              end if
            case ('CoM-diff')
              if(tokens(j+1) == 'TRUE') then
                CoMDiff=.true.
              elseif(tokens(j+1) == 'FALSE') then
                CoMDiff=.false.
              else
                print '(" Inconsistent CoM-diff.")'
                stop
              end if
          end select
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

    ! allocate(Qst(nsegx3,nchain,nprun))
    ! allocate(rcmst(nchain,3,nprun))
    ! allocate(Rbst(ntotbeadx3,nprun))

    allocate(rcmst(ntotchain,3,nprun))
    allocate(Rbst(ntotbeadx3,nprun))
    allocate(Qst(ntotsegx3,nprun))

    if ((FlowType == 'Equil').and.CoMDiff) then
      ! allocate(rcminit(nchain,3,nprun))
      ! allocate(cmifst(nchain,3,nprun))

      allocate(rcminit(ntotchain,3,nprun))
      allocate(cmifst(ntotchain,3,nprun))
    end if

  end subroutine init_io

  !> Constructor for configurational io type
  ! subroutine init_conf_io(this,id,p,nchain,nsegx3,nbeadx3,MPI_REAL_WP)
  subroutine init_conf_io(this,id,p,ntotchain,ntotsegx3,ntotbeadx3,MPI_REAL_WP)

    use,intrinsic :: iso_fortran_env
    use :: flow_mod, only: FlowType
    use :: strg_mod
    use :: mpi
    !include 'mpif.h'

    class(conf_io),intent(inout) :: this
    ! integer,intent(in) :: id,p,nchain,nsegx3,nbeadx3,MPI_REAL_WP
    integer,intent(in) :: id,p,ntotchain,ntotsegx3,ntotbeadx3,MPI_REAL_WP
    integer :: ierr
    integer :: q_starts(2),q_sizes(2),q_subsizes(2)
    integer :: rc_starts(3),rc_sizes(3),rc_subsizes(3)
    integer :: rb_starts(2),rb_sizes(2),rb_subsizes(2)
    integer :: rf_starts(2),rf_sizes(2),rf_subsizes(2)
    integer :: stat1,stat2,stat3,stat4,stat5,stat6
    integer :: u1,u2,u3,u4,u5,u6,i,j,ios,ntokens,stat,il
    character(len=1024) :: line
    character(len=100) :: tokens(10)


    ! default values:
    this%qfctr=[0.7_wp,0._wp,0._wp]
    this%qfctr_cmbbb=[0.7_wp,0._wp,0._wp]
    this%qfctr_cmbar=[0.0_wp,0.7_wp,0._wp]
    this%initmode='st'
    this%MakeAnim=.false.
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
            case ('qfctr')
              call value(tokens(j+1),this%qfctr(1),ios)
              call value(tokens(j+2),this%qfctr(2),ios)
              call value(tokens(j+3),this%qfctr(3),ios)
            case ('qfctr-comb-bb')
              call value(tokens(j+1),this%qfctr_cmbbb(1),ios)
              call value(tokens(j+2),this%qfctr_cmbbb(2),ios)
              call value(tokens(j+3),this%qfctr_cmbbb(3),ios)
            case ('qfctr-comb-ar')
              call value(tokens(j+1),this%qfctr_cmbar(1),ios)
              call value(tokens(j+2),this%qfctr_cmbar(2),ios)
              call value(tokens(j+3),this%qfctr_cmbar(3),ios)
            case ('initmode')
              this%initmode=trim(adjustl(tokens(j+1)))
            case ('MakeAnim')
              if(tokens(j+1) == 'TRUE') then
                this%MakeAnim=.true.
              elseif(tokens(j+1) == 'FALSE') then
                this%MakeAnim=.false.
              else
                print '(" Inconsistent MakeAnim.")'
              end if
          end select
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)
    ! For making binary restart and dump files:
    ! q
    ! q_sizes=[nsegx3*nchain,p]
    ! q_subsizes=[nsegx3*nchain,1]
    q_sizes=[ntotsegx3,p]
    q_subsizes=[ntotsegx3,1]
    q_starts=[0,id]
    call MPI_Type_create_subarray(2,q_sizes,q_subsizes,q_starts,MPI_ORDER_FORTRAN,&
                                  MPI_REAL_WP,this%q_recvsubarray,ierr)
    call MPI_Type_commit(this%q_recvsubarray,ierr)
    call MPI_File_open(MPI_COMM_WORLD,"data/rst/q.rst.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                       MPI_INFO_NULL,this%fqrsthandle,ierr)
    if (FlowType /= 'Equil') then
      call MPI_File_open(MPI_COMM_WORLD,"data/st/q.st.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&
                         this%fqsthandle,ierr)
    end if
    ! rcm
    ! rc_sizes=[nchain,3,p]
    ! rc_subsizes=[nchain,3,1]
    rc_sizes=[ntotchain,3,p]
    rc_subsizes=[ntotchain,3,1]
    rc_starts=[0,0,id]
    call MPI_Type_create_subarray(3,rc_sizes,rc_subsizes,rc_starts,MPI_ORDER_FORTRAN,&
                                  MPI_REAL_WP,this%rc_recvsubarray,ierr)
    call MPI_Type_commit(this%rc_recvsubarray,ierr)
    call MPI_File_open(MPI_COMM_WORLD,"data/rst/CoM.rst.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                       MPI_INFO_NULL,this%fcrsthandle,ierr)
    if (FlowType /= 'Equil') then
      call MPI_File_open(MPI_COMM_WORLD,"data/st/CoM.st.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&
                         this%fcsthandle,ierr)
    end if
    ! cm image flag
    if ((FlowType == 'Equil').and.CoMDiff) then
      call MPI_Type_create_subarray(3,rc_sizes,rc_subsizes,rc_starts,MPI_ORDER_FORTRAN,&
                                    MPI_INTEGER,this%cif_recvsubarray,ierr)
      call MPI_Type_commit(this%cif_recvsubarray,ierr)
      call MPI_File_open(MPI_COMM_WORLD,"data/rst/cmif.rst.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%fcifrsthandle,ierr)
      call MPI_File_open(MPI_COMM_WORLD,"data/rst/CoM.init.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%fcinithandle,ierr)
    end if
    ! Rb
    ! rb_sizes=[nbeadx3*nchain,p]
    ! rb_subsizes=[nbeadx3*nchain,1]
    rb_sizes=[ntotbeadx3,p]
    rb_subsizes=[ntotbeadx3,1]
    rb_starts=[0,id]
    call MPI_Type_create_subarray(2,rb_sizes,rb_subsizes,rb_starts,MPI_ORDER_FORTRAN,&
                                  MPI_REAL_WP,this%rb_recvsubarray,ierr)
    call MPI_Type_commit(this%rb_recvsubarray,ierr)
    call MPI_File_open(MPI_COMM_WORLD,"data/rst/Rb.rst.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                       MPI_INFO_NULL,this%frbrsthandle,ierr)
    if (FlowType /= 'Equil') then
      call MPI_File_open(MPI_COMM_WORLD,"data/st/Rb.st.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&
                         this%frbsthandle,ierr)
    end if
    if (DumpConf) then
      ! For making binary dump files for the entire simulation period:
      ! q
      call MPI_File_open(MPI_COMM_WORLD,"data/dump/q.dmp.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%fqdmphandle,ierr)
      ! rcm
      call MPI_File_open(MPI_COMM_WORLD,"data/dump/CoM.dmp.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%fcdmphandle,ierr)
      ! cm image flag
      if ((FlowType == 'Equil').and.CoMDiff) then
        call MPI_File_open(MPI_COMM_WORLD,"data/dump/cmif.dmp.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                           MPI_INFO_NULL,this%fcifdmphandle,ierr)
      end if
      ! Rb
      call MPI_File_open(MPI_COMM_WORLD,"data/dump/Rb.dmp.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%frbdmphandle,ierr)
      ! rFphi
      rf_sizes=[4,p]
      rf_subsizes=[4,1]
      rf_starts=[0,id]
      call MPI_Type_create_subarray(2,rf_sizes,rf_subsizes,rf_starts,MPI_ORDER_FORTRAN,&
                                    MPI_REAL_WP,this%rf_recvsubarray,ierr)
      call MPI_Type_commit(this%rf_recvsubarray,ierr)
      call MPI_File_open(MPI_COMM_WORLD,"data/dump/rF.dmp.bin",MPI_MODE_RDWR+MPI_MODE_CREATE,&
                         MPI_INFO_NULL,this%frfdmphandle,ierr)
    end if ! DumpConf

    if (id == 0) then

      if (FlowType == 'Equil') then
        open(newunit=this%oldu3,file='data/rst/q.equil.dat',status='replace')
        open(newunit=this%oldu4,file='data/rst/CoM.equil.dat',status='replace')
        open(newunit=this%oldu9,file='data/rst/Rb.equil.dat',status='replace')
      else
        open(newunit=this%oldu3,iostat=stat1,file='data/rst/q.final.dat',status='replace')
        open(newunit=this%oldu4,iostat=stat2,file='data/rst/CoM.final.dat',status='replace')
        open(newunit=this%oldu9,file='data/rst/Rb.final.dat',status='replace')
      end if ! FlowType == 'Equil'

      open(newunit=this%oldu1,file='data/rst/q.rst.dat',status='replace')
      open(newunit=this%oldu2,file='data/rst/CoM.rst.dat',status='replace')
      open(newunit=this%oldu8,file='data/rst/Rb.rst.dat',status='replace')

      if (this%MakeAnim) then
        open(newunit=this%oldu5,file='data/dump/Rb.dat',status='replace',position='append')
        open(newunit=this%oldu6,file='data/dump/CoM.dat',status='replace',position='append')
        if (FlowType /= 'Equil') &
          open(newunit=this%oldu7,file='data/dump/BoxConfig.dat',status='replace',position='append')
      end if ! this%MakeAnim
    end if ! id

  end subroutine init_conf_io

!   !> Reads the configurational information for all entities
!   !! \param id The rank of the process
!   !! \param p The number of processes
!   !! \param boxsize The dimension of primary box
!   subroutine read_conf(this,id,p,boxsize,nchain,nseg,nbead,nsegx3,ntotbead,&
!                        ntotsegx3,ntotbeadx3,nprun,runrst,qmax,MPI_REAL_WP)

!     use :: flow_mod, only: FlowType
!     use :: arry_mod, only: print_vector,print_matrix
!     use :: conv_mod, only: QtoR
!     use :: mpi
!     !include 'mpif.h'

!     class(conf_io),intent(in) :: this
!     integer,intent(in) :: id,p,nseg,nbead,nsegx3,nchain,nprun,runrst
!     integer,intent(in) :: ntotbead,ntotsegx3,ntotbeadx3,MPI_REAL_WP
!     real(wp),intent(in) :: boxsize(3),qmax
!     integer ::ichain,iseg,offset,i,msec,n,time_info(8),irun,ierr,ich,os,igb
!     integer(kind=MPI_OFFSET_KIND) :: offsetMPI
!     real(wp),pointer :: rcmstP(:,:) => null()
!     integer,pointer :: cmifstP(:,:) => null()
!     real(wp),pointer :: QstPP(:) => null()
!     real(wp),pointer :: RbstP(:) => null()
!     real(wp),pointer,contiguous :: QstP(:,:) => null()
!     integer :: intvar
!     real(wp) :: realvar


!     !   %-------------------------------------------------------%
!     !   | The initial guess for Qs in case we are looking for   |
!     !   | equilibrium connector vector is arbitrary.            |
!     !   | In case we start from equilibrium or any other start  |
!     !   | condition, the Qs is sequentially read from file.     |
!     !   %-------------------------------------------------------%

!     if (FlowType == 'Equil') then
!       select case (this%initmode)

!         case ('st')

!           do irun=1, nprun
!             do ichain=1, nchain
!               do iseg=1, nseg
!                 offset=3*(iseg-1)
!                 Qst(offset+1:offset+3,ichain,irun)=[this%qfctr(1)*qmax,&
!                                                     this%qfctr(2)*qmax,&
!                                                     this%qfctr(3)*qmax]
!               end do ! iseg
!             end do ! ichain
!           end do ! irun

!           end do ! irun

!           ! For debugging:
!          rcmst(1,1:3,1)=(/-4.517_wp,2.308_wp,-3.395_wp/)
!          rcmst(2,1:3,1)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
!          rcmst(3,1:3,1)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
!          rcmst(4,1:3,1)=(/-3-0.034_wp,-3-3.173_wp,-3-1.454_wp/) ! For 3
!          ! rcmst(4,1:3,1)=(/24.766_wp,21.627_wp,23.346_wp/) ! For 3
! !          rcmst(1,1:3,2)=(/-4.517_wp,2.308_wp,-3.395_wp/)
! !          rcmst(2,1:3,2)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
! !          rcmst(3,1:3,2)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
! !          rcmst(4,1:3,2)=(/-3-0.034_wp,-3-3.173_wp,-3-1.454_wp/) ! For 3
! !
!           ! call date_and_time(values=time_info)
!           ! msec=(1000*time_info(7)+time_info(8))*((id-83)*359) ! a somewhat random integer
!           ! call random_seed(size=n) ! get the number of integers used for the seed
!           ! ! This is because we want different order of random numbers in each call
!           ! call random_seed(put=(/(i*msec,i=1,n)/)) ! give a proper seed
!           ! call random_number(rcmst) ! generate a sequence of nchain pseudo-random numbers
!           ! rcmst(:,1,:)=rcmst(:,1,:)*boxsize(1)
!           ! rcmst(:,2,:)=rcmst(:,2,:)*boxsize(2)
!           ! rcmst(:,3,:)=rcmst(:,3,:)*boxsize(3)




!           do irun=1, nprun
!             offsetMPI=nchain*3*p*sizeof(realvar)*(irun-1)
!             call MPI_File_set_view(this%fcinithandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
!                                    MPI_INFO_NULL,ierr)
!             call MPI_File_write(this%fcinithandle,rcmst(:,:,irun),nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           end do

!           if (CoMDiff) cmifst=0

!           do irun=1, nprun
!             QstP => Qst(:,:,irun)
!             QstPP(1:size(QstP)) => QstP
!             RbstP => Rbst(:,irun)
!             call QtoR(QstPP,RbstP,ntotsegx3,ntotbeadx3)
!             do igb=1, ntotbead
!               os=(igb-1)*3
!               ich=(igb-1)/nbead+1
!               RbstP(os+1)=RbstP(os+1)+rcmst(ich,1,irun)
!               RbstP(os+2)=RbstP(os+2)+rcmst(ich,2,irun)
!               RbstP(os+3)=RbstP(os+3)+rcmst(ich,3,irun)
!             end do
!           end do

!         case ('rst')

!           ! rc
!           offsetMPI=nchain*3*p*runrst*sizeof(realvar)
!           rcmstP => rcmst(:,:,runrst+1)
!           call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
!                                  "native",MPI_INFO_NULL,ierr)
!           call MPI_File_read(this%fcrsthandle,rcmstP,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           if (CoMDiff) then
!             ! cm image flag
!             offsetMPI=nchain*3*p*runrst*sizeof(intvar)
!             cmifstP => cmifst(:,:,runrst+1)
!             call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%fcifrsthandle,cmifstP,nchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!           end if
!           ! q
!           offsetMPI=ntotsegx3*p*runrst*sizeof(realvar)
!           QstP => Qst(:,:,runrst+1)
!           QstPP(1:size(QstP)) => QstP
!           call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
!                                  "native",MPI_INFO_NULL,ierr)
!           call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           ! Rb
!           offsetMPI=ntotbeadx3*p*runrst*sizeof(realvar)
!           RbstP => Rbst(:,runrst+1)
!           call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
!                                  "native",MPI_INFO_NULL,ierr)
!           call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)

!           do irun=runrst+2, nprun
!             do ichain=1, nchain
!               do iseg=1, nseg
!                 offset=3*(iseg-1)
!                 Qst(offset+1:offset+3,ichain,irun)=[this%qfctr(1)*qmax,&
!                                                     this%qfctr(2)*qmax,&
!                                                     this%qfctr(3)*qmax]
!               end do ! iseg
!             end do ! ichain
!           end do ! irun
! !         For debugging:
! !          rcmst(1,1:3,1)=(/-4.517_wp,2.308_wp,-3.395_wp/)
! !          rcmst(2,1:3,1)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
! !          rcmst(3,1:3,1)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
! !          rcmst(1,1:3,2)=(/-4.517_wp,2.308_wp,-3.395_wp/)
! !          rcmst(2,1:3,2)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
! !          rcmst(3,1:3,2)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3

!           if (runrst+1 < nprun) then

!             if (CoMDiff) cmifst(:,:,runrst+2:nprun)=0

!             call date_and_time(values=time_info)
!             msec=(1000*time_info(7)+time_info(8))*((id-83)*359) ! a somewhat random integer
!             call random_seed(size=n) ! get the number of integers used for the seed
!             ! This is because we want different order of random numbers in each call
!             call random_seed(put=(/(i*msec,i=1,n)/)) ! give a proper seed
!             call random_number(rcmst(:,:,runrst+2:nprun))
!             rcmst(:,1,runrst+2:nprun)=(rcmst(:,1,runrst+2:nprun)-0.5)*boxsize(1)
!             rcmst(:,2,runrst+2:nprun)=(rcmst(:,2,runrst+2:nprun)-0.5)*boxsize(2)
!             rcmst(:,3,runrst+2:nprun)=(rcmst(:,3,runrst+2:nprun)-0.5)*boxsize(3)
!           end if

!           do irun=runrst+2, nprun
!             offsetMPI=nchain*3*p*sizeof(realvar)*(irun-1)
!             call MPI_File_set_view(this%fcinithandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
!                                    MPI_INFO_NULL,ierr)
!             call MPI_File_write(this%fcinithandle,rcmst(:,:,irun),nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           end do

!         case ('ext')

!           do irun=1, nprun
!             ! rc
!             offsetMPI=nchain*3*p*(irun-1)*sizeof(realvar)
!             rcmstP => rcmst(:,:,irun)
!             call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%fcrsthandle,rcmstP,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             if (CoMDiff) then
!               ! cm image flag
!               offsetMPI=nchain*3*p*(irun-1)*sizeof(intvar)
!               cmifstP => cmifst(:,:,irun)
!               call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%fcifrsthandle,cmifstP,nchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!             end if
!             ! q
!             offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
!             QstP => Qst(:,:,irun)
!             QstPP(1:size(QstP)) => QstP
!             call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             ! Rb
!             offsetMPI=ntotbeadx3*p*runrst*sizeof(realvar)
!             RbstP => Rbst(:,runrst+1)
!             call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           end do ! irun

!       end select ! this%initmode

!     else ! FlowType /= Equil

!       select case (this%initmode)

!         case ('st')

!           do irun=1, nprun
!             ! rc
!             offsetMPI=nchain*3*p*(irun-1)*sizeof(realvar)
!             rcmstP => rcmst(:,:,irun)
!             call MPI_File_set_view(this%fcsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%fcsthandle,rcmstP,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             ! q
!             offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
!             QstP => Qst(:,:,irun)
!             qstPP(1:size(QstP)) => QstP
!             call MPI_File_set_view(this%fqsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%fqsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             ! Rb
!             offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
!             RbstP => Rbst(:,irun)
!             call MPI_File_set_view(this%frbsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
!                                    "native",MPI_INFO_NULL,ierr)
!             call MPI_File_read(this%frbsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!           end do ! irun

!         case ('rst','ext')
!           do irun=1, nprun
!             if (irun < runrst+2) then
!               ! rc
!               offsetMPI=nchain*3*p*(irun-1)*sizeof(realvar)
!               rcmstP => rcmst(:,:,irun)
!               call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%fcrsthandle,rcmstP,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!               ! q
!               offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
!               QstP => Qst(:,:,irun)
!               qstPP(1:size(QstP)) => QstP
!               call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!               ! Rb
!               offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
!               RbstP => Rbst(:,irun)
!               call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             else
!               ! rc
!               offsetMPI=nchain*3*p*(irun-1)*sizeof(realvar)
!               rcmstP => rcmst(:,:,irun)
!               call MPI_File_set_view(this%fcsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%fcsthandle,rcmstP,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!               ! q
!               offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
!               QstP => Qst(:,:,irun)
!               qstPP(1:size(QstP)) => QstP
!               call MPI_File_set_view(this%fqsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%fqsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!               ! Rb
!               offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
!               RbstP => Rbst(:,irun)
!               call MPI_File_set_view(this%frbsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
!                                      "native",MPI_INFO_NULL,ierr)
!               call MPI_File_read(this%frbsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!             end if
!           end do ! irun

!       end select ! this%initmode

!     end if ! FlowType

!   end subroutine read_conf


  !> Reads the configurational information for all entities
  !! \param id The rank of the process
  !! \param p The number of processes
  !! \param boxsize The dimension of primary box
  subroutine read_conf(this,id,p,boxsize,nchain,nseg,nbead,nsegx3,ntotchain,ntotbead,&
    ntotsegx3,ntotbeadx3,nprun,runrst,qmax,MPI_REAL_WP,add_cmb,nchain_cmb,nseg_cmb,&
    nseg_cmbbb,nseg_cmbar)

    use :: flow_mod, only: FlowType
    use :: arry_mod, only: print_vector,print_matrix
    use :: conv_mod, only: QtoR
    use :: mpi
    !include 'mpif.h'

    class(conf_io),intent(in) :: this
    integer,intent(in) :: id,p,nseg,nbead,nsegx3,nchain,nprun,runrst
    integer,intent(in) :: ntotchain,ntotbead,ntotsegx3,ntotbeadx3,MPI_REAL_WP
    real(wp),intent(in) :: boxsize(3),qmax
    logical,intent(in) :: add_cmb
    integer,intent(in) :: nchain_cmb,nseg_cmb,nseg_cmbbb,nseg_cmbar
    integer ::ichain,iseg,offset,i,msec,n,time_info(8),irun,ierr,ich,os,igb
    integer(kind=MPI_OFFSET_KIND) :: offsetMPI
    real(wp),pointer :: rcmstP(:,:) => null()
    integer,pointer :: cmifstP(:,:) => null()
    real(wp),pointer :: QstPP(:) => null()
    real(wp),pointer :: RbstP(:) => null()
    real(wp),pointer,contiguous :: QstP(:,:) => null()
    integer :: intvar
    real(wp) :: realvar

    integer :: offsetch


    !   %-------------------------------------------------------%
    !   | The initial guess for Qs in case we are looking for   |
    !   | equilibrium connector vector is arbitrary.            |
    !   | In case we start from equilibrium or any other start  |
    !   | condition, the Qs is sequentially read from file.     |
    !   %-------------------------------------------------------%

    if (FlowType == 'Equil') then
      select case (this%initmode)

        case ('st')

          do irun=1, nprun
            do ichain=1, nchain
              offsetch=(ichain-1)*nseg*3
              do iseg=1, nseg
                offset=offsetch+(iseg-1)*3
                Qst(offset+1:offset+3,irun)=[this%qfctr(1)*qmax,this%qfctr(2)*qmax,this%qfctr(3)*qmax]
              end do ! iseg
            end do ! ichain
            if (add_cmb) then
              do ichain=1, nchain_cmb
                offsetch=nchain*nseg*3 + (ichain-1)*nseg_cmb*3

                do iseg=1, nseg_cmb
                  offset=offsetch+(iseg-1)*3

                  if (iseg <= nseg_cmbbb) then
                    Qst(offset+1:offset+3,irun)=[this%qfctr_cmbbb(1)*qmax,this%qfctr_cmbbb(2)*qmax,this%qfctr_cmbbb(3)*qmax]
                  else
                    if (mod((iseg-nseg_cmbbb-1)/nseg_cmbar+1,2) == 0) then
                      Qst(offset+1:offset+3,irun)=-[this%qfctr_cmbar(1)*qmax,this%qfctr_cmbar(2)*qmax,this%qfctr_cmbar(3)*qmax]
                    else
                      Qst(offset+1:offset+3,irun)=[this%qfctr_cmbar(1)*qmax,this%qfctr_cmbar(2)*qmax,this%qfctr_cmbar(3)*qmax]
                    end if
                  end if

                end do ! iseg
              end do ! ichain
            endif

          end do ! irun

          ! For debugging:
        !  rcmst(1,1:3,1)=(/0.1_wp,1._wp,0.1_wp/)
        !  do ichain=2,size(rcmst,1)
        !    rcmst(ichain,1:3,1)=rcmst(1,1:3,1)+[0._wp,1._wp,0]
        !  enddo
        !  rcmst(1,1:3,1)=(/-4.517_wp,2.308_wp,-3.395_wp/)
        !  rcmst(2,1:3,1)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
        !  rcmst(3,1:3,1)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
        !  rcmst(4,1:3,1)=(/3-0.034_wp,-3-3.173_wp,-3-1.454_wp/) ! For 3
         ! rcmst(4,1:3,1)=(/24.766_wp,21.627_wp,23.346_wp/) ! For 3
!          rcmst(1,1:3,2)=(/-4.517_wp,2.308_wp,-3.395_wp/)
!          rcmst(2,1:3,2)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
!          rcmst(3,1:3,2)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
!          rcmst(4,1:3,2)=(/-3-0.034_wp,-3-3.173_wp,-3-1.454_wp/) ! For 3
!
          call date_and_time(values=time_info)
          msec=(1000*time_info(7)+time_info(8))*((id-83)*359) ! a somewhat random integer
          call random_seed(size=n) ! get the number of integers used for the seed
          ! This is because we want different order of random numbers in each call
          call random_seed(put=(/(i*msec,i=1,n)/)) ! give a proper seed
          call random_number(rcmst) ! generate a sequence of nchain pseudo-random numbers
          rcmst(:,1,:)=rcmst(:,1,:)*boxsize(1)
          rcmst(:,2,:)=rcmst(:,2,:)*boxsize(2)
          rcmst(:,3,:)=rcmst(:,3,:)*boxsize(3)

          do irun=1, nprun
            offsetMPI=ntotchain*3*p*sizeof(realvar)*(irun-1)
            call MPI_File_set_view(this%fcinithandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
                                   MPI_INFO_NULL,ierr)
            call MPI_File_write(this%fcinithandle,rcmst(:,:,irun),ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          end do

          if (CoMDiff) cmifst=0

          do irun=1, nprun
            QstPP => Qst(:,irun)
            RbstP => Rbst(:,irun)
            call QtoR(QstPP,RbstP,ntotsegx3,ntotbeadx3)

            do igb=1, nchain*nbead
              os=(igb-1)*3
              ich=(igb-1)/nbead+1
              RbstP(os+1)=RbstP(os+1)+rcmst(ich,1,irun)
              RbstP(os+2)=RbstP(os+2)+rcmst(ich,2,irun)
              RbstP(os+3)=RbstP(os+3)+rcmst(ich,3,irun)
            end do
            if (add_cmb) then
              do igb=1, nchain_cmb*(nseg_cmb+1)
                os=nchain*nbead*3+(igb-1)*3
                ich=nchain+(igb-1)/(nseg_cmb+1)+1
                RbstP(os+1)=RbstP(os+1)+rcmst(ich,1,irun)
                RbstP(os+2)=RbstP(os+2)+rcmst(ich,2,irun)
                RbstP(os+3)=RbstP(os+3)+rcmst(ich,3,irun)
              end do
            endif

            ! call print_matrix(rcmst(1:2,:,1),'rcm')
            ! ! call print_vector(QstPP,'qst')
            ! call print_vector(RbstP,'rbst')
          end do

        case ('rst')

          ! rc
          offsetMPI=ntotchain*3*p*runrst*sizeof(realvar)
          rcmstP => rcmst(:,:,runrst+1)
          call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
                                 "native",MPI_INFO_NULL,ierr)
          call MPI_File_read(this%fcrsthandle,rcmstP,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          if (CoMDiff) then
            ! cm image flag
            offsetMPI=ntotchain*3*p*runrst*sizeof(intvar)
            cmifstP => cmifst(:,:,runrst+1)
            call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%fcifrsthandle,cmifstP,ntotchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
          end if
          ! q
          offsetMPI=ntotsegx3*p*runrst*sizeof(realvar)
          QstPP => Qst(:,runrst+1)
          ! QstPP(1:size(QstP)) => QstP
          call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                                 "native",MPI_INFO_NULL,ierr)
          call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          ! Rb
          offsetMPI=ntotbeadx3*p*runrst*sizeof(realvar)
          RbstP => Rbst(:,runrst+1)
          call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
                                 "native",MPI_INFO_NULL,ierr)
          call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)

          do irun=runrst+2, nprun
            do ichain=1, nchain
              offsetch=(ichain-1)*nseg*3
              do iseg=1, nseg
                offset=offsetch+(iseg-1)*3
                Qst(offset+1:offset+3,irun)=[this%qfctr(1)*qmax,this%qfctr(2)*qmax,this%qfctr(3)*qmax]
              end do ! iseg
            end do ! ichain
            if (add_cmb) then
              do ichain=1, nchain_cmb
                offsetch=nchain*nseg*3+(ichain-1)*nseg_cmb*3

                do iseg=1, nseg_cmb
                  offset=offsetch+(iseg-1)*3

                  if (iseg <= nseg_cmbbb) then
                    Qst(offset+1:offset+3,irun)=[0.9_wp*qmax,0._wp,0._wp]
                  else
                    if (mod((iseg-nseg_cmbbb-1)/nseg_cmbar+1,2) == 0) then
                      Qst(offset+1:offset+3,irun)=[0._wp,-0.9_wp*qmax,0._wp]
                    else
                      Qst(offset+1:offset+3,irun)=[0._wp,0.9_wp*qmax,0._wp]
                    end if
                  end if

                end do ! iseg
              end do ! ichain
            endif
          end do ! irun
!         For debugging:
!          rcmst(1,1:3,1)=(/-4.517_wp,2.308_wp,-3.395_wp/)
!          rcmst(2,1:3,1)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
!          rcmst(3,1:3,1)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3
!          rcmst(1,1:3,2)=(/-4.517_wp,2.308_wp,-3.395_wp/)
!          rcmst(2,1:3,2)=(/5-0.034_wp,5-3.173_wp,5-1.454_wp/) ! For 3
!          rcmst(3,1:3,2)=(/1-0.034_wp,1-3.173_wp,1-1.454_wp/) ! For 3

          if (runrst+1 < nprun) then

            if (CoMDiff) cmifst(:,:,runrst+2:nprun)=0

            call date_and_time(values=time_info)
            msec=(1000*time_info(7)+time_info(8))*((id-83)*359) ! a somewhat random integer
            call random_seed(size=n) ! get the number of integers used for the seed
            ! This is because we want different order of random numbers in each call
            call random_seed(put=(/(i*msec,i=1,n)/)) ! give a proper seed
            call random_number(rcmst(:,:,runrst+2:nprun))
            rcmst(:,1,runrst+2:nprun)=(rcmst(:,1,runrst+2:nprun)-0.5)*boxsize(1)
            rcmst(:,2,runrst+2:nprun)=(rcmst(:,2,runrst+2:nprun)-0.5)*boxsize(2)
            rcmst(:,3,runrst+2:nprun)=(rcmst(:,3,runrst+2:nprun)-0.5)*boxsize(3)
          end if

          do irun=runrst+2, nprun
            offsetMPI=ntotchain*3*p*sizeof(realvar)*(irun-1)
            call MPI_File_set_view(this%fcinithandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
                                   MPI_INFO_NULL,ierr)
            call MPI_File_write(this%fcinithandle,rcmst(:,:,irun),ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          end do

        case ('ext')

          do irun=1, nprun
            ! rc
            offsetMPI=ntotchain*3*p*(irun-1)*sizeof(realvar)
            rcmstP => rcmst(:,:,irun)
            call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%fcrsthandle,rcmstP,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            if (CoMDiff) then
              ! cm image flag
              offsetMPI=ntotchain*3*p*(irun-1)*sizeof(intvar)
              cmifstP => cmifst(:,:,irun)
              call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%fcifrsthandle,cmifstP,ntotchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            end if
            ! q
            offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
            QstPP => Qst(:,irun)
            ! QstPP(1:size(QstP)) => QstP
            call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            ! Rb
            offsetMPI=ntotbeadx3*p*runrst*sizeof(realvar)
            RbstP => Rbst(:,runrst+1)
            call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          end do ! irun

      end select ! this%initmode

    else ! FlowType /= Equil

      select case (this%initmode)

        case ('st')

          do irun=1, nprun
            ! rc
            offsetMPI=ntotchain*3*p*(irun-1)*sizeof(realvar)
            rcmstP => rcmst(:,:,irun)
            call MPI_File_set_view(this%fcsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%fcsthandle,rcmstP,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            ! q
            offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
            QstPP => Qst(:,irun)
            ! qstPP(1:size(QstP)) => QstP
            call MPI_File_set_view(this%fqsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%fqsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            ! Rb
            offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
            RbstP => Rbst(:,irun)
            call MPI_File_set_view(this%frbsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
                                   "native",MPI_INFO_NULL,ierr)
            call MPI_File_read(this%frbsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
          end do ! irun

        case ('rst','ext')
          do irun=1, nprun
            if (irun < runrst+2) then
              ! rc
              offsetMPI=ntotchain*3*p*(irun-1)*sizeof(realvar)
              rcmstP => rcmst(:,:,irun)
              call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%fcrsthandle,rcmstP,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
              ! q
              offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
              QstPP => Qst(:,irun)
              ! qstPP(1:size(QstP)) => QstP
              call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%fqrsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
              ! Rb
              offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
              RbstP => Rbst(:,irun)
              call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%frbrsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            else
              ! rc
              offsetMPI=ntotchain*3*p*(irun-1)*sizeof(realvar)
              rcmstP => rcmst(:,:,irun)
              call MPI_File_set_view(this%fcsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%fcsthandle,rcmstP,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
              ! q
              offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
              QstPP => Qst(:,irun)
              ! qstPP(1:size(QstP)) => QstP
              call MPI_File_set_view(this%fqsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%fqsthandle,QstPP,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
              ! Rb
              offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
              RbstP => Rbst(:,irun)
              call MPI_File_set_view(this%frbsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,&
                                     "native",MPI_INFO_NULL,ierr)
              call MPI_File_read(this%frbsthandle,RbstP,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
            end if
          end do ! irun

      end select ! this%initmode

    end if ! FlowType

  end subroutine read_conf

  !> Reades the initial configurational information for all entities
  !! \param p The number of processes
  !! \param irun The run index
  !! \param rcm The center of mass for all chains
  subroutine read_init_conf(this,p,irun,ntotchain,MPI_REAL_WP)

    use :: flow_mod, only: FlowType
    use :: mpi
    !include 'mpif.h'

    class(conf_io),intent(in) :: this
    integer,intent(in) :: p,irun,ntotchain,MPI_REAL_WP
    integer ::ierr
    integer(kind=MPI_OFFSET_KIND) :: offsetMPI
    real(wp) :: realvar


    ! rc
    if ((FlowType == 'Equil').and.CoMDiff) then
      offsetMPI=ntotchain*3*p*sizeof(realvar)*(irun-1)
      call MPI_File_set_view(this%fcinithandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray&
                             &,"native",MPI_INFO_NULL,ierr)
      call MPI_File_read(this%fcinithandle,rcmst(:,:,irun),ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    end if

  end subroutine read_init_conf

  !> Reades the dumped configurational information for all entities
  !! \param p The number of processes
  !! \param irun The run index
  !! \param idmp The dumping index
  !! \param Q The connectivity vector for all chains
  !! \param rcm The center of mass for all chains
  subroutine read_dmp_conf(this,p,irun,idmp,Q,rcm,cmif,ntotchain,ntotsegx3,ndmp,MPI_REAL_WP)

    use :: force_smdlt, only: rFphi
    use :: flow_mod, only: FlowType
    use :: mpi
    !include 'mpif.h'

    class(conf_io),intent(in) :: this
    integer,intent(in) :: p,irun,idmp,ntotchain,ntotsegx3,ndmp,MPI_REAL_WP
    real(wp),intent(inout),target :: Q(:),rcm(:,:)
    integer,intent(inout) :: cmif(:,:)
    integer ::ierr
    integer(kind=MPI_OFFSET_KIND) :: offsetMPI
    integer :: intvar
    real(wp) :: realvar


    ! q
    offsetMPI=ntotsegx3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
    call MPI_File_set_view(this%fqdmphandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,&
                           "native",MPI_INFO_NULL,ierr)
    call MPI_File_read(this%fqdmphandle,Q,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    ! rc
    offsetMPI=ntotchain*3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
    call MPI_File_set_view(this%fcdmphandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray&
                           &,"native",MPI_INFO_NULL,ierr)
    call MPI_File_read(this%fcdmphandle,rcm,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    if ((FlowType == 'Equil').and.CoMDiff) then
      ! cm image flag
      offsetMPI=ntotchain*3*p*sizeof(intvar)*((irun-1)*ndmp+idmp-1)
      call MPI_File_set_view(this%fcifdmphandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray&
                             &,"native",MPI_INFO_NULL,ierr)
      call MPI_File_read(this%fcifdmphandle,cmif,ntotchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    end if
    ! rFphi
    offsetMPI=4*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
    call MPI_File_set_view(this%frfdmphandle,offsetMPI,MPI_REAL_WP,this%rf_recvsubarray&
                           &,"native",MPI_INFO_NULL,ierr)
    call MPI_File_read(this%frfdmphandle,rFphi,4,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end subroutine read_dmp_conf

!   !> Writes the configurational information for all entities
!   !! \param id The rank of the process
!   !! \param p The number of processes
!   !! \param itime The time index
!   !! \param ntime The end time index
!   !! \param irun The run index
!   !! \param idmp The dumping index
!   !! \param time Current time
!   !! \param Wi The Weissenberg number
!   !! \param dt The time step size
!   !! \param Q The connectivity vector for all chains
!   !! \param rcm The center of mass for all chains
!   !! \param Rb The position vector of beads for all chains
!   subroutine write_conf(this,id,p,itime,ntime,irun,idmp,time,Wi,dt,nchain,nbead,&
!                     nsegx3,nbeadx3,ntotsegx3,ntotbeadx3,ndmp,lambda,MPI_REAL_WP,&
!                     Q,rcm,cmif,Rb,R)

!     use :: flow_mod, only: FlowType
!     use :: arry_mod, only: print_vector,print_matrix
!     use :: force_smdlt, only: rFphi
!     use :: trsfm_mod, only: delrx_L,L1,L2
!     use :: mpi
!     !include 'mpif.h'

!     class(conf_io),intent(in) :: this
!     integer,intent(in) :: id,p,itime,ntime,irun,idmp,nsegx3,nchain
!     integer,intent(in) :: ntotsegx3,ndmp,nbead,nbeadx3,ntotbeadx3
!     integer,intent(in) :: MPI_REAL_WP
!     real(wp),intent(in) :: Wi,dt,time,lambda
!     real(wp),intent(in) :: Q(:),rcm(:,:)
!     integer,intent(in) :: cmif(:,:)
!     real(wp),intent(in) :: Rb(ntotbeadx3)
!     real(wp),intent(in),optional :: R(ntotbeadx3)
!     integer :: iseg,ibead,ichain,ierr,offsetch,offsetb,intvar
!     real(wp) :: rtpassed,realvar
!     integer(kind=MPI_OFFSET_KIND) :: offsetMPI


!     ! For dumping the configuration of the system:
!     if (DumpConf) then
!       ! q
!       offsetMPI=ntotsegx3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
!       call MPI_File_set_view(this%fqdmphandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,"native",&
!                              MPI_INFO_NULL,ierr)
!       call MPI_File_write(this%fqdmphandle,Q,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!       ! rc
!       offsetMPI=nchain*3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
!       call MPI_File_set_view(this%fcdmphandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
!                              MPI_INFO_NULL,ierr)
!       call MPI_File_write(this%fcdmphandle,rcm,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!       if ((FlowType == 'Equil').and.CoMDiff) then
!         ! cm image flag
!         offsetMPI=nchain*3*p*sizeof(intvar)*((irun-1)*ndmp+idmp-1)
!         call MPI_File_set_view(this%fcifdmphandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,"native",&
!                                MPI_INFO_NULL,ierr)
!         call MPI_File_write(this%fcifdmphandle,cmif,nchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!       end if
!       ! rFphi
!       offsetMPI=4*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
!       call MPI_File_set_view(this%frfdmphandle,offsetMPI,MPI_REAL_WP,this%rf_recvsubarray,"native",&
!                              MPI_INFO_NULL,ierr)
!       call MPI_File_write(this%frfdmphandle,rFphi,4,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)

!       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!     end if

!     ! For writing the restart files:
!     ! q
!     offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
!     call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,"native",&
!                            MPI_INFO_NULL,ierr)
!     call MPI_File_write(this%fqrsthandle,Q,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!     ! rc
!     offsetMPI=nchain*3*p*(irun-1)*sizeof(realvar)
!     call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
!                            MPI_INFO_NULL,ierr)
!     call MPI_File_write(this%fcrsthandle,rcm,nchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!     if ((FlowType == 'Equil').and.CoMDiff) then
!       ! cm image flag
!       offsetMPI=nchain*3*p*(irun-1)*sizeof(intvar)
!       call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,"native",&
!                              MPI_INFO_NULL,ierr)
!       call MPI_File_write(this%fcifrsthandle,cmif,nchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
!     end if
!     ! Rb
!     offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
!     call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,"native",&
!                            MPI_INFO_NULL,ierr)
!     call MPI_File_write(this%frbrsthandle,Rb,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
!     call MPI_Barrier(MPI_COMM_WORLD,ierr)

!     ! Providing data files to clarify the status of bin files:
!     if (itime == ntime) then

!       if (id == 0) then

!         rtpassed=time/lambda
!         write (this%oldu3,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write (this%oldu4,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write (this%oldu9,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write(this%oldu3,'(2x,a)') 'AT:'
!         write(this%oldu3,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu3,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu3,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu3,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu3,'(2x,a,1x,i8)') 'dump number:',idmp
!         write(this%oldu4,'(2x,a)') 'AT:'
!         write(this%oldu4,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu4,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu4,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu4,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu4,'(2x,a,1x,i8)') 'dump number:',idmp
!         write(this%oldu9,'(2x,a)') 'AT:'
!         write(this%oldu9,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu9,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu9,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu9,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu9,'(2x,a,1x,i8)') 'dump number:',idmp

!       end if ! id

!     else ! itime /= ntime

!       if (id == 0) then

!         rewind(this%oldu1);rewind(this%oldu2);rewind(this%oldu8)
!         rtpassed=time/lambda
!         write (this%oldu1,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write (this%oldu2,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write (this%oldu8,'(f8.3,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
!         write(this%oldu1,'(2x,a)') 'AT:'
!         write(this%oldu1,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu1,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu1,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu1,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu1,'(2x,a,1x,i8)') 'dump number:',idmp
!         write(this%oldu2,'(2x,a)') 'AT:'
!         write(this%oldu2,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu2,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu2,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu2,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu2,'(2x,a,1x,i8)') 'dump number:',idmp
!         write(this%oldu8,'(2x,a)') 'AT:'
!         write(this%oldu8,'(2x,a,1x,f10.3)') 'Wi:',Wi
!         write(this%oldu8,'(2x,a,1x,f10.3)') 'dt:',dt
!         write(this%oldu8,'(2x,a,1x,i8)') 'run number:',irun
!         write(this%oldu8,'(2x,a,1x,i8)') 'time index number:',itime
!         if (DumpConf) write(this%oldu8,'(2x,a,1x,i8)') 'dump number:',idmp

!       end if ! id

!     end if ! itime

!     if (id == 0) then
!       if (this%MakeAnim) then
!         do ichain=1, nchain
!           write(this%oldu6,'(3(f18.7,1x))') rcm(ichain,:)
!           offsetch=(ichain-1)*nbeadx3
!           do ibead=1, nbead
!             offsetb=(ibead-1)*3
! !            write(this%oldu5,'(3(f18.7,1x))') Rb(offsetch+offsetb+1:offsetch+offsetb+3)
!             write(this%oldu5,'(3(f18.7,1x))') R(offsetch+offsetb+1:offsetch+offsetb+3)+&
!                                               rcm(ichain,1:3)
!           end do ! ibead
!         end do ! ichain
!         select case (FlowType)
!           case ('PSF')
!             write(this%oldu7,'(f18.7)') delrx_L
!           case ('PEF')
!             write(this%oldu7,'(4(f18.7,1x))') L1,L2
!         end select
!       end if ! this%MakeAnim
!     end if ! id

!   end subroutine write_conf

  subroutine write_conf(this,id,p,itime,ntime,irun,idmp,time,Wi,dt,nchain,nbead,&
    nsegx3,nbeadx3,ntotchain,ntotsegx3,ntotbeadx3,ndmp,lambda,MPI_REAL_WP,Q,rcm,&
    cmif,Rb,R,add_cmb,nchain_cmb,nseg_cmb)

    use :: flow_mod, only: FlowType
    use :: arry_mod, only: print_vector,print_matrix
    use :: force_smdlt, only: rFphi
    use :: trsfm_mod, only: delrx_L,L1,L2
    use :: mpi
    !include 'mpif.h'

    class(conf_io),intent(in) :: this
    integer,intent(in) :: id,p,itime,ntime,irun,idmp,nsegx3,nchain
    integer,intent(in) :: ntotchain,ntotsegx3,ndmp,nbead,nbeadx3,ntotbeadx3
    integer,intent(in) :: MPI_REAL_WP
    real(wp),intent(in) :: Wi,dt,time,lambda
    real(wp),intent(in) :: Q(:),rcm(:,:)
    integer,intent(in) :: cmif(:,:)
    real(wp),intent(in) :: Rb(ntotbeadx3)
    real(wp),intent(in) :: R(ntotbeadx3)
    logical :: add_cmb
    integer,intent(in) :: nchain_cmb,nseg_cmb
    integer :: iseg,ibead,ichain,ierr,offsetch,offsetb,intvar
    real(wp) :: rtpassed,realvar
    integer(kind=MPI_OFFSET_KIND) :: offsetMPI


    ! For dumping the configuration of the system:
    if (DumpConf) then
      ! q
      offsetMPI=ntotsegx3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
      call MPI_File_set_view(this%fqdmphandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,"native",&
                             MPI_INFO_NULL,ierr)
      call MPI_File_write(this%fqdmphandle,Q,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
      ! rc
      offsetMPI=ntotchain*3*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
      call MPI_File_set_view(this%fcdmphandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
                             MPI_INFO_NULL,ierr)
      call MPI_File_write(this%fcdmphandle,rcm,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
      if ((FlowType == 'Equil').and.CoMDiff) then
        ! cm image flag
        offsetMPI=ntotchain*3*p*sizeof(intvar)*((irun-1)*ndmp+idmp-1)
        call MPI_File_set_view(this%fcifdmphandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,"native",&
                               MPI_INFO_NULL,ierr)
        call MPI_File_write(this%fcifdmphandle,cmif,ntotchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      end if
      ! rFphi
      offsetMPI=4*p*sizeof(realvar)*((irun-1)*ndmp+idmp-1)
      call MPI_File_set_view(this%frfdmphandle,offsetMPI,MPI_REAL_WP,this%rf_recvsubarray,"native",&
                             MPI_INFO_NULL,ierr)
      call MPI_File_write(this%frfdmphandle,rFphi,4,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end if

    ! For writing the restart files:
    ! q
    offsetMPI=ntotsegx3*p*(irun-1)*sizeof(realvar)
    call MPI_File_set_view(this%fqrsthandle,offsetMPI,MPI_REAL_WP,this%q_recvsubarray,"native",&
                           MPI_INFO_NULL,ierr)
    call MPI_File_write(this%fqrsthandle,Q,ntotsegx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    ! rc
    offsetMPI=ntotchain*3*p*(irun-1)*sizeof(realvar)
    call MPI_File_set_view(this%fcrsthandle,offsetMPI,MPI_REAL_WP,this%rc_recvsubarray,"native",&
                           MPI_INFO_NULL,ierr)
    call MPI_File_write(this%fcrsthandle,rcm,ntotchain*3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    if ((FlowType == 'Equil').and.CoMDiff) then
      ! cm image flag
      offsetMPI=ntotchain*3*p*(irun-1)*sizeof(intvar)
      call MPI_File_set_view(this%fcifrsthandle,offsetMPI,MPI_INTEGER,this%cif_recvsubarray,"native",&
                             MPI_INFO_NULL,ierr)
      call MPI_File_write(this%fcifrsthandle,cmif,ntotchain*3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    end if
    ! Rb
    offsetMPI=ntotbeadx3*p*(irun-1)*sizeof(realvar)
    call MPI_File_set_view(this%frbrsthandle,offsetMPI,MPI_REAL_WP,this%rb_recvsubarray,"native",&
                           MPI_INFO_NULL,ierr)
    call MPI_File_write(this%frbrsthandle,Rb,ntotbeadx3,MPI_REAL_WP,MPI_STATUS_IGNORE,ierr)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! Providing data files to clarify the status of bin files:
    if (itime == ntime) then

      if (id == 0) then

        rtpassed=time/lambda
        write (this%oldu3,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write (this%oldu4,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write (this%oldu9,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write(this%oldu3,'(2x,a)') 'AT:'
        write(this%oldu3,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu3,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu3,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu3,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu3,'(2x,a,1x,i8)') 'dump number:',idmp
        write(this%oldu4,'(2x,a)') 'AT:'
        write(this%oldu4,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu4,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu4,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu4,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu4,'(2x,a,1x,i8)') 'dump number:',idmp
        write(this%oldu9,'(2x,a)') 'AT:'
        write(this%oldu9,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu9,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu9,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu9,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu9,'(2x,a,1x,i8)') 'dump number:',idmp
        write(*,'(2x,a)') 'You can restart using the following:'
        write (*,'(a,f14.7)') " trst: ",rtpassed
        write(*,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(*,'(2x,a,1x,f14.7)') 'dt:',dt
        write(*,'(2x,a,1x,i8)') 'runrst:',irun-1
        write(*,'(2x,a,1x,i8)') 'trst index (second entry):',itime
        if (DumpConf) write(*,'(2x,a,1x,i8)') 'dmprst:',idmp
        write (*,*)

      end if ! id

    else ! itime /= ntime

      if (id == 0) then

        rewind(this%oldu1);rewind(this%oldu2);rewind(this%oldu8)
        rtpassed=time/lambda
        write (this%oldu1,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write (this%oldu2,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write (this%oldu8,'(f14.7,a)') rtpassed," 'Chain-Relaxation-Time(s)' Passed"
        write(this%oldu1,'(2x,a)') 'AT:'
        write(this%oldu1,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu1,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu1,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu1,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu1,'(2x,a,1x,i8)') 'dump number:',idmp
        write(this%oldu2,'(2x,a)') 'AT:'
        write(this%oldu2,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu2,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu2,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu2,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu2,'(2x,a,1x,i8)') 'dump number:',idmp
        write(this%oldu8,'(2x,a)') 'AT:'
        write(this%oldu8,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(this%oldu8,'(2x,a,1x,f14.7)') 'dt:',dt
        write(this%oldu8,'(2x,a,1x,i8)') 'run number:',irun
        write(this%oldu8,'(2x,a,1x,i8)') 'time index number:',itime
        if (DumpConf) write(this%oldu8,'(2x,a,1x,i8)') 'dump number:',idmp
        write(*,'(2x,a)') 'You can restart using the following:'
        write (*,'(a,f14.7)') " trst: ",rtpassed
        write(*,'(2x,a,1x,f14.7)') 'Wi:',Wi
        write(*,'(2x,a,1x,f14.7)') 'dt:',dt
        write(*,'(2x,a,1x,i8)') 'runrst:',irun-1
        write(*,'(2x,a,1x,i8)') 'trst index (second entry):',itime
        if (DumpConf) write(*,'(2x,a,1x,i8)') 'dmprst:',idmp
        write (*,*)

      end if ! id

    end if ! itime

! call print_vector(R,'rio')
! call print_matrix(rcm,'rcmio')

    if (id == 0) then
      if (this%MakeAnim) then
        do ichain=1, nchain
          write(this%oldu6,'(3(f18.7,1x))') rcm(ichain,:)
          offsetch=(ichain-1)*nbeadx3
          do ibead=1, nbead
            offsetb=(ibead-1)*3
            ! write(this%oldu5,'(3(f18.7,1x))') Rb(offsetch+offsetb+1:offsetch+offsetb+3)
            write(this%oldu5,'(3(f18.7,1x))') R(offsetch+offsetb+1:offsetch+offsetb+3)+&
                                              rcm(ichain,1:3)
          end do ! ibead
        end do ! ichain
        if (add_cmb) then
          do ichain=1, nchain_cmb
            write(this%oldu6,'(3(f18.7,1x))') rcm(nchain+ichain,:)
            offsetch=nchain*nbeadx3+(ichain-1)*(nseg_cmb+1)*3
            do ibead=1, nseg_cmb+1
              offsetb=(ibead-1)*3
              ! write(this%oldu5,'(3(f18.7,1x))') Rb(offsetch+offsetb+1:offsetch+offsetb+3)
              write(this%oldu5,'(3(f18.7,1x))') R(offsetch+offsetb+1:offsetch+offsetb+3)+&
              rcm(nchain+ichain,1:3)
            end do ! ibead
          end do ! ichain
        endif
        select case (FlowType)
          case ('PSF')
            write(this%oldu7,'(f18.7)') delrx_L
          case ('PEF')
            write(this%oldu7,'(4(f18.7,1x))') L1,L2
        end select
      end if ! this%MakeAnim
    end if ! id

  end subroutine write_conf

  !> Destructor for configurational io type
  subroutine del_io(this)

    use :: flow_mod, only: FlowType

    type(conf_io) :: this
    integer :: ierr

    close(this%oldu1);close(this%oldu2)
    close(this%oldu3);close(this%oldu4)
    close(this%oldu5);close(this%oldu6)
    close(this%oldu7);close(this%oldu8)
    close(this%oldu9)
    call MPI_File_close(this%fqrsthandle,ierr)
    call MPI_File_close(this%fcrsthandle,ierr)
    call MPI_File_close(this%frbrsthandle,ierr)
    if ((FlowType /= 'Equil') .and. (this%initmode == 'st')) then
      call MPI_File_close(this%fqsthandle,ierr)
      call MPI_File_close(this%fcsthandle,ierr)
      call MPI_File_close(this%frbsthandle,ierr)
    end if
    if (DumpConf) then
      call MPI_File_close(this%fqdmphandle,ierr)
      call MPI_File_close(this%fcdmphandle,ierr)
    end if
    deallocate(Qst,rcmst,Rbst)
    if ((FlowType == 'Equil').and.CoMDiff) deallocate(cmifst)

  end subroutine del_io

end module io_mod
