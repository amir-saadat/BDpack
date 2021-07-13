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
! MODULE: verlet
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: construction of verlet list
!>
!
!--------------------------------------------------------------------

module evverlet_mod

  use :: prcn_mod

  implicit none

  ! Private module procedures:
  private :: init_verlet_t,&
             init_clllst  ,&
             get_ncps     ,&
             cnstr_clllst ,&
             cnstr_nablst ,&
             del_verlet_t

  !> A public type for constructing verlet list
  type evverlet

    private
    !> The number of the cells per side
    integer :: ncps(3)
    !> Total number of cells in the cubic box
    integer :: nct
    !> The dimension of the cells in each direction
    real(wp) :: cll_sz(3)
    !> The volume of the cells
    real(wp) :: cll_vol
    !> Maximum number of beads per cell
    integer :: mbpc
    !> The array which contains the number of beads in cells
    integer,allocatable :: head(:)
    !> The array which contains the beads index in cells
    integer,allocatable :: binc(:,:)
    !> The array which contains the neighboring cell list
    integer,allocatable :: nclst(:,:)
    !> Maximum occupancy of the cells
    integer :: mocc
    !> i index in all possible interactions
    integer,allocatable :: iidx(:)
    !> j index in all possible interactions
    integer,allocatable :: jidx(:)
    !> An array for keeping track of interactions of interest
    logical,allocatable :: inside(:)
    !> The x-component of vector between beads with indices i,j
    real(wp),allocatable :: Rijx(:)
    !> The y-component of vector between beads with indices i,j
    real(wp),allocatable :: Rijy(:)
    !> The temporary y-component of vector between beads with indices i,j
    real(wp),allocatable :: Rijytmp(:)
    !> The z-component of vector between beads with indices i,j
    real(wp),allocatable :: Rijz(:)
    !> The squared distance between beads with indices i,j
    real(wp),allocatable :: Rijsq(:)
    !> Total number of interactions possible
    integer :: num_int

  contains

    procedure,pass(this) :: init => init_verlet_t
    procedure,pass(this) :: init_cll => init_clllst
    procedure,pass(this) :: cnstr_cll => cnstr_clllst
    procedure,pass(this) :: cnstr_nab => cnstr_nablst
    procedure,pass(this) :: get_ncps
    final :: del_verlet_t

  end type evverlet
    !> The array which contains the number of beads in cells
    integer,allocatable :: head(:)


  ! Private module variables:
  private :: cll_dns_ev,nnc
  ! Protected module variables:
  ! protected ::

  !> The density of particles in a cell
  real(wp),save :: cll_dns_ev
  !> Number of neighbering cells
  integer,save :: nnc
  ! !> The neighboring cells offset
  ! integer,pointer,save :: shifts_ptr(:,:)
  ! !> The coordinates for neighboring cells
  ! integer,pointer :: j_clx_ptr(:),j_cly_ptr(:),j_clz_ptr(:),j_cll_ptr(:)

contains

  !> Initializes the verlet module
  !! \param id The rank of the process
  subroutine init_evverlet(id)

    use :: arry_mod, only: print_vector
    use :: flow_mod, only: FlowType
    use :: strg_mod
    use,intrinsic :: iso_fortran_env
    use :: cmn_io_mod, only: read_input
	use :: verlet_mod, only: shifts,j_clx,j_cly,j_clz,j_cll

    integer,intent(in) :: id
#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:init_evverlet"
#endif
    call read_input('cll-dns-ev',0,cll_dns_ev,0.1_wp)
    select case (FlowType)

      case ('Equil','PSF')
        nnc=13
!         allocate(shifts(nnc,3))
!         shifts(1,:) =[ 0, 0,-1]
!         shifts(2,:) =[ 1, 0,-1]
!         shifts(3,:) =[ 1, 0, 0]
!         shifts(4,:) =[ 1, 0, 1]
!         shifts(5,:) =[-1, 1,-1]
!         shifts(6,:) =[ 0, 1,-1]
!         shifts(7,:) =[ 1, 1,-1]
!         shifts(8,:) =[-1, 1, 0]
!         shifts(9,:) =[ 0, 1, 0]
!         shifts(10,:)=[ 1, 1, 0]
!         shifts(11,:)=[-1, 1, 1]
!         shifts(12,:)=[ 0, 1, 1]
!         shifts(13,:)=[ 1, 1, 1]

      case ('PEF')
        nnc=31
!         allocate(shifts(nnc,3))
!         shifts(1,:) =[ 0, 0,-1]
!         shifts(2,:) =[ 1, 0,-1]
!         shifts(3,:) =[ 2, 0,-1]
!         shifts(4,:) =[ 3, 0,-1]
!         shifts(5,:) =[ 1, 0, 0]
!         shifts(6,:) =[ 2, 0, 0]
!         shifts(7,:) =[ 3, 0, 0]
!         shifts(8,:) =[ 1, 0, 1]
!         shifts(9,:) =[ 2, 0, 1]
!         shifts(10,:)=[ 3, 0, 1]
!         shifts(11,:)=[-3, 1,-1]
!         shifts(12,:)=[-2, 1,-1]
!         shifts(13,:)=[-1, 1,-1]
!         shifts(14,:)=[ 0, 1,-1]
!         shifts(15,:)=[ 1, 1,-1]
!         shifts(16,:)=[ 2, 1,-1]
!         shifts(17,:)=[ 3, 1,-1]
!         shifts(18,:)=[-3, 1, 0]
!         shifts(19,:)=[-2, 1, 0]
!         shifts(20,:)=[-1, 1, 0]
!         shifts(21,:)=[ 0, 1, 0]
!         shifts(22,:)=[ 1, 1, 0]
!         shifts(23,:)=[ 2, 1, 0]
!         shifts(24,:)=[ 3, 1, 0]
!         shifts(25,:)=[-3, 1, 1]
!         shifts(26,:)=[-2, 1, 1]
!         shifts(27,:)=[-1, 1, 1]
!         shifts(28,:)=[ 0, 1, 1]
!         shifts(29,:)=[ 1, 1, 1]
!         shifts(30,:)=[ 2, 1, 1]
!         shifts(31,:)=[ 3, 1, 1]
! !        this%ncps(1:2)=bs(1:2)/(sqrt(10._wp)*rc_F)
! !        this%ncps(3)=bs(3)/rc_F
    end select

!     allocate(j_clx(nnc))
!     allocate(j_cly(nnc))
!     allocate(j_clz(nnc))
!     allocate(j_cll(nnc))

  end subroutine init_evverlet

  !> Constructor for  verlet type
  !! \param rc The cutoff radius
  !! \param bs The dimension of the box
  subroutine init_verlet_t(this,rc,bs,ntotbead)

    class(evverlet),intent(inout) :: this
    real(wp),intent(in) :: rc,bs(3)
    integer,intent(in) :: ntotbead
#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:init_evverlet_t"
#endif
    this%ncps=0
    write(*,*) "init_cll"
    call this%init_cll(rc,bs,ntotbead)

  end subroutine init_verlet_t

  subroutine get_ncps(this)

    class(evverlet),intent(inout) :: this

    print *
    print '(" Initial number of cells for EV calculation: ")'
    print '(3(i10,1x))',this%ncps

  end subroutine get_ncps

  !> Initializes the cell list
  !! \param rc The cutoff radius
  !! \param bs The dimension of the box
  subroutine init_clllst(this,rc,bs,ntotbead)

    use :: verlet_mod, only: shifts,j_clx,j_cly,j_clz,j_cll
    use :: flow_mod, only: FlowType
!    use :: inp_smdlt, only: ntotbead
    use :: arry_mod, only: print_vector

    class(evverlet),intent(inout) :: this
    real(wp),intent(in) :: rc,bs(3)
    integer,intent(in) :: ntotbead
    integer :: clx,cly,clz,cll,czNxNy,cyNx,ierr
    real(wp) :: ncpsl(3)
#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:init_clllst"
#endif
    ncpsl=this%ncps

    select case (FlowType)
      case ('Equil')
        this%ncps(:)=bs(:)/rc
      case ('PSF')
        this%ncps(1)=bs(1)/(sqrt(2._wp)*rc)
        this%ncps(2:3)=bs(2:3)/rc
      case ('PEF')
        this%ncps(1)=bs(1)/(sqrt(10._wp)*rc/3)
        this%ncps(2:3)=bs(2:3)/rc
    end select
    this%cll_sz(1:3)=bs(1:3)/this%ncps(1:3)
    this%nct=this%ncps(1)*this%ncps(2)*this%ncps(3)
    this%cll_vol=bs(1)*bs(2)*bs(3)/this%nct

    this%mbpc=int(this%cll_vol*cll_dns_ev)

    if (allocated(this%binc)) deallocate(this%binc)
    allocate(this%binc(this%nct,this%mbpc))

    if (any(this%ncps /= ncpsl)) then

      if (allocated(this%head)) deallocate(this%head)
      if (allocated(this%nclst)) deallocate(this%nclst)

      allocate(this%head(this%nct))
      allocate(this%nclst(this%nct,nnc))

      do clz=0, this%ncps(3)-1
        czNxNy=clz*this%ncps(1)*this%ncps(2)
        do cly=0, this%ncps(2)-1
          cyNx=cly*this%ncps(1)
          do clx=0, this%ncps(1)-1
            cll=czNxNy+cyNx+clx+1
            ! Make sure to take only the first nnc part - later in the cnstr_nab we only need half of the neighbors since j>i always
            j_clx(1:nnc)=clx+shifts(1:nnc,1)
            j_cly(1:nnc)=cly+shifts(1:nnc,2)
            j_clz(1:nnc)=clz+shifts(1:nnc,3)
            j_clx(1:nnc)=modulo(j_clx(1:nnc),this%ncps(1))
            j_cly(1:nnc)=modulo(j_cly(1:nnc),this%ncps(2))
            j_clz(1:nnc)=modulo(j_clz(1:nnc),this%ncps(3))
            j_cll(1:nnc)=j_clz(1:nnc)*this%ncps(1)*this%ncps(2)+j_cly(1:nnc)*this%ncps(1)+j_clx(1:nnc)+1
            this%nclst(cll,:)=j_cll(1:nnc)
          end do ! clx
        end do ! cly
      end do ! clz

    end if

    this%num_int=ntotbead*nnc*this%mbpc*0.5

    if (allocated(this%iidx)) deallocate(this%iidx)
    if (allocated(this%jidx)) deallocate(this%jidx)
    if (allocated(this%inside)) deallocate(this%inside)
    if (allocated(this%Rijx)) deallocate(this%Rijx)
    if (allocated(this%Rijy)) deallocate(this%Rijy)
    if (allocated(this%Rijz)) deallocate(this%Rijz)
    if (allocated(this%Rijsq)) deallocate(this%Rijsq)

    if (FlowType == 'PEF') then
      if (allocated(this%Rijytmp)) deallocate(this%Rijytmp)
    end if
!write(*,*) "1-1"
    ! We check the memory allocation result, as num_int might be a big one
    allocate(this%iidx(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for iidx in evverlet_mod"
!write(*,*) "1-2"
    allocate(this%jidx(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for jidx in evverlet_mod"
!write(*,*) "1-3"
    allocate(this%inside(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for inside in evverlet_mod"
!write(*,*) "1-3"
    allocate(this%Rijx(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for Rijx in evverlet_mod"
!write(*,*) "1-4"
    allocate(this%Rijy(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for Rijy in evverlet_mod"
!write(*,*) "1-5"
    allocate(this%Rijz(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for Rijz in evverlet_mod"
!write(*,*) "1-6"
    allocate(this%Rijsq(this%num_int),stat=ierr)
    if ( ierr /= 0 ) stop "Memory allocation issue for Rijsq in evverlet_mod"
!write(*,*) "1-7"
    if (FlowType == 'PEF') then
      allocate(this%Rijytmp(this%num_int),stat=ierr)
     if ( ierr /= 0 ) stop "Memory allocation issue for Rijytmp in evverlet_mod"
    end if
!write(*,*) "1-8"
  end subroutine init_clllst

  !> Constructs the cell list
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  subroutine cnstr_clllst(this,Rbx,Rby,Rbz,itime,ntotbead,ntotbeadx3)

!    use :: inp_smdlt, only: ntotbead,ntotbeadx3
    use :: arry_mod, only: print_vector,print_matrix

    class(evverlet),intent(inout) :: this
    real(wp),intent(in) :: Rbx(:)
    real(wp),intent(in) :: Rby(:)
    real(wp),intent(in) :: Rbz(:)
    integer,intent(in) :: ntotbead,ntotbeadx3
    integer :: i,clx,cly,clz,cll,itime,j

#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:cnstr_clllst"
#endif
    this%head=0
    this%binc=0

    do i=1, ntotbead

      clx=int(Rbx(i)/this%cll_sz(1))
      cly=int(Rby(i)/this%cll_sz(2))
      clz=int(Rbz(i)/this%cll_sz(3))

      ! if the bead is exactly on the boundary
      if (clx == this%ncps(1)) clx=clx-1
      if (cly == this%ncps(2)) cly=cly-1
      if (clz == this%ncps(3)) clz=clz-1

      cll=clz*this%ncps(1)*this%ncps(2)+cly*this%ncps(1)+clx+1

      ! print*,'2',clx,cly,clz,cll
      ! print*, cll
	  ! print*, this%head(cll)

      this%head(cll)=this%head(cll)+1

      if (this%head(cll) >= this%mbpc) then
        print '(" Warning: cll-dns-ev is too small from evverlet_mod. ")'
        stop
      endif

      this%binc(cll,this%head(cll))=i

!      print *,'cll,head,i',cll,this%head(cll),i
!      print *,'binc',this%binc(cll,this%head(cll))

    end do
    this%mocc=maxval(this%head)

    print '(" maxocc (from the actual bead positions): ",i4," -- max bead per cell (from cell_dns_ev): ",i6)',this%mocc,this%mbpc

!if (itime==8479)then
!    call print_vector(this%head,'newhead')
!do i=1,size(this%binc,1)
!do j=1,size(this%binc,2)
!if (this%binc(i,j)/=0)then
!print*,'i,j',i,j
!print*,'binc',this%binc(i,j)
!endif
!enddo
!enddo
!    print *,'mocc',this%mocc
!endif

  end subroutine cnstr_clllst

  !> Constructs the neighbor list
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param bs the dimension of the box
  !! \param invbs the inverse of box dimensions
  !! \param nlst The neighbor list
  subroutine cnstr_nablst(this,Rbx,Rby,Rbz,rc,bs,invbs,nlst,itime,ntotbead,ntotbeadx3)

!    use :: inp_smdlt, only: ntotbead,ntotbeadx3
    use :: arry_mod, only: print_vector,print_matrix
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m,tanb,sinth,costh

    class(evverlet),intent(inout) :: this
    real(wp),intent(in),contiguous :: Rbx(:)
    real(wp),intent(in),contiguous :: Rby(:)
    real(wp),intent(in),contiguous :: Rbz(:)
    real(wp),intent(in) :: rc
    integer,intent(in) :: itime,ntotbead,ntotbeadx3
    integer,allocatable,intent(inout) :: nlst(:,:)
    integer,allocatable :: beadi_tmp(:),beadj(:),beadj_tmp(:)
    logical,allocatable :: pair(:)
    integer :: i,j,nab,idx,cll,beadi,k,intidx
    real(wp) :: bs(3),invbs(3),rcsq
#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:cnstr_nablst"
#endif
    this%iidx=0
    this%jidx=0
    allocate(beadi_tmp(this%nct))
    allocate(beadj_tmp(this%nct))
    allocate(pair(this%nct))

    ! Same-cell interactions:
    idx=1
    do i=1, this%mocc-1
      beadi_tmp=this%binc(:,i)
      do j=i+1, this%mocc
        beadj_tmp=this%binc(:,j)
        pair=beadi_tmp < beadj_tmp
        nab=count(pair)
        this%iidx(idx:(idx+nab-1))=pack(beadi_tmp,mask=pair)
        this%jidx(idx:(idx+nab-1))=pack(beadj_tmp,mask=pair)
!if(itime==6028) then
!print *,'i',i,j
!print *,'nab',nab
!call print_vector(beadi_tmp,'bi')
!call print_vector(beadj_tmp,'bj')
!call print_vector(this%iidx(idx:(idx+nab-1)),'iidx')
!call print_vector(this%jidx(idx:(idx+nab-1)),'jidx')
!end if
        idx=idx+nab
      end do
    end do
    deallocate(beadi_tmp)
    deallocate(beadj_tmp)
    deallocate(pair)

    ! Different-cell interactions:
    allocate(beadj(nnc*this%mbpc))
    allocate(beadj_tmp(nnc*this%mbpc))
    allocate(pair(nnc*this%mbpc))

    do cll=1, this%nct
      beadj=0
      beadj_tmp=0
      do j=1, nnc
        beadj_tmp((j-1)*this%mbpc+1:j*this%mbpc)=this%binc(this%nclst(cll,j),:)
      end do
      pair=beadj_tmp /= 0
      nab=count(pair)
      beadj(1:nab)=pack(beadj_tmp,mask=pair)
      do i=1, this%mbpc
        beadi=this%binc(cll,i)
        if (beadi == 0) exit
        this%iidx(idx:(idx+nab-1))=beadi
        this%jidx(idx:(idx+nab-1))=beadj(1:nab)
        idx=idx+nab
      end do ! i
    end do ! cll
    idx=idx-1

    deallocate(beadj)
    deallocate(beadj_tmp)
    deallocate(pair)

!$omp parallel default(private) shared(this,Rbx,Rby,Rbz,eps_m,sinth,costh,tanb,rc) &
!$omp shared(idx,bs,invbs,FlowType)
!$omp do simd
    do intidx=1, idx

      this%Rijx(intidx)=Rbx(this%iidx(intidx))-Rbx(this%jidx(intidx))
      this%Rijy(intidx)=Rby(this%iidx(intidx))-Rby(this%jidx(intidx))
      this%Rijz(intidx)=Rbz(this%iidx(intidx))-Rbz(this%jidx(intidx))
      ! Minimum Image Covention:
      this%Rijx(intidx)=this%Rijx(intidx)-nint(this%Rijx(intidx)*invbs(1))*bs(1)
      this%Rijy(intidx)=this%Rijy(intidx)-nint(this%Rijy(intidx)*invbs(2))*bs(2)
      this%Rijz(intidx)=this%Rijz(intidx)-nint(this%Rijz(intidx)*invbs(3))*bs(3)
      select case (FlowType)
        case ('PSF')
          this%Rijx(intidx)=this%Rijx(intidx)+eps_m*this%Rijy(intidx)
        case ('PEF')
          this%Rijytmp(intidx)=this%Rijy(intidx)
          this%Rijx(intidx)=this%Rijx(intidx)+tanb*this%Rijytmp(intidx)
          this%Rijy(intidx)=sinth*this%Rijx(intidx)+costh*this%Rijytmp(intidx)
          this%Rijx(intidx)=costh*this%Rijx(intidx)-sinth*this%Rijytmp(intidx)
      end select
      this%Rijsq(intidx)=this%Rijx(intidx)*this%Rijx(intidx) + &
                         this%Rijy(intidx)*this%Rijy(intidx) + &
                         this%Rijz(intidx)*this%Rijz(intidx)
    end do
!!$omp end do simd
    this%inside=.false.
!$omp do simd
    do intidx=1, idx
      this%inside(intidx)=this%Rijsq(intidx) <= rc**2
    end do
!!$omp end do simd
!$omp end parallel
!$ivdep
    nab=count(this%inside)
    if(allocated(nlst)) deallocate(nlst)
    allocate(nlst(nab,2))
    nlst(:,1)=pack(this%iidx,mask=this%inside)
    nlst(:,2)=pack(this%jidx,mask=this%inside)
  end subroutine cnstr_nablst


  !> Destructor for  verlet type
  subroutine del_verlet_t(this)

    type(evverlet),intent(inout) :: this
#ifdef Debuge_sequence
    write(*,*) "module:evverlet_mod:del_verlet_t"
#endif
  end subroutine del_verlet_t

  ! subroutine del_verlet()

  !   deallocate(shifts)
  !   deallocate(j_clx,j_cly,j_clz,j_cll)

  ! end subroutine del_verlet

end module evverlet_mod
