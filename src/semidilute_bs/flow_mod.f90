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
! MODULE: flow
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION:
!> Applies homogeneous flow to the entities inside the box
!--------------------------------------------------------------------

module flow_mod

  use :: prcn_mod

  implicit none

  private :: init_flow_t ,&
             apply_flow  ,&
             del_flow_t

  !> A public type for applying flow
  type flow

    private
    !> The arrays for storing sparse Kappa
    !> @{
    real(wp),allocatable :: K_vals(:)
    integer,allocatable :: K_cols(:),K_rowIdx(:)
    !> @}

  contains

    procedure,pass(this) :: init => init_flow_t
    procedure,pass(this) :: apply => apply_flow
    final :: del_flow_t

  end type flow

  protected :: FlowType

  !> The type of flow applied to the entities inside the box
  character(len=10) :: FlowType

contains

  !> Initializes flow_mod module variables
  subroutine init_flow(id)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id
    integer :: i,j,ntokens,u1,il,stat
    character(len=1024) :: line
    character(len=100) :: tokens(50)
#ifdef Debuge_sequence
    write(*,*) "module:flow_mod:init_flow"
#endif
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
          if(trim(adjustl(tokens(j))) == 'Flow-Type') then
            FlowType=trim(adjustl(tokens(j+1)))
          end if
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

  end subroutine init_flow

  !> Constructor for flow type
  !! \param
  !subroutine init_flow_t(this,nchain,nbead)
  subroutine init_flow_t(this,ntotbeadx3)
  !MB
    use :: arry_mod, only: print_vector

    class(flow),intent(inout) :: this
    !integer,intent(in) :: nchain,nbead  !MB
	integer,intent(in) :: ntotbeadx3 !MB
    integer :: maxnz,ich,os,ibx3,nbeadx3,ntotbead
#ifdef Debuge_sequence
    write(*,*) "module:flow_mod:init_flow_t"
#endif
    !nbeadx3=nbead*3  !MB
    !ntotbead=nbead*nchain
	!ntotbeadx3=ntotbead*3

    !!!!!! should be fixed for comb polymer

    ! For making K sparse (CSR):
    ! Specifying kappa based on type of flow
    select case (FlowType)
      case ('PSF')
!        !    | 0  1  0 |
!        !  k=| 0  0  0 |
!        !    | 0  0  0 |
!        maxnz=nchain*nbead
!        allocate(this%K_vals(maxnz))
!        allocate(this%K_cols(maxnz))
!        allocate(this%K_rowIdx(ntotbeadx3+1))
!        this%K_rowIdx(1)=1
!        do ich=1, nchain
!          os=(ich-1)*nbead
!          do ibx3=1, nbeadx3
!            if (mod(ibx3,3) == 1) then ! Only first row of 3x3 matrix k matters.
!              this%K_cols(os+ibx3/3+1)=os*3+ibx3+1
!             this%K_vals(os+ibx3/3+1)=1._wp
!              this%K_rowIdx(os*3+ibx3+1)=this%K_rowIdx(os*3+ibx3)+1
!            else
!              this%K_rowIdx(os*3+ibx3+1)=this%K_rowIdx(os*3+ibx3)
!            end if
!          end do
!        end do
!       !     | 0  1  0 |
        !   k=\ 0  0  0 \
        !     | 0  0  0 |
        print*, 'ntotbeadx3=', ntotbeadx3
        maxnz=ntotbeadx3/3
        print*, 'maxnz=', maxnz
        allocate(this%K_vals(maxnz))
        allocate(this%K_cols(maxnz))
        allocate(this%K_rowIdx(ntotbeadx3+1))
        this%K_rowIdx(1)=1
!        do ich=1, nchain
!          os=(ich-1)*nbead
!          do ibx3=1, nbeadx3
          do ibx3=1, ntotbeadx3
            if (mod(ibx3,3) == 1) then ! Only first row of 3x3 matrix k matters.
              !print*, "ibx3=" ,ibx3
              this%K_cols(ibx3/3+1)=ibx3+1
              this%K_vals(ibx3/3+1)=1._wp
              this%K_rowIdx(ibx3+1)=this%K_rowIdx(ibx3)+1
            else
              this%K_rowIdx(ibx3+1)=this%K_rowIdx(ibx3)
            end if
          end do
!        end do
      ! For planar extensional flow:
      case ('PEF')
!        !    | 1  0  0 |
!        !  k=| 0 -1  0 |
!        !    | 0  0  0 |
!        maxnz=nchain*2*nbead
!        allocate(this%K_vals(maxnz))
!        allocate(this%K_cols(maxnz))
!        allocate(this%K_rowIdx(ntotbeadx3+1))
!        this%K_rowIdx(1)=1
!        do ich=1, nchain
!          os=(ich-1)*nbead
!          do ibx3=1, nbeadx3
!            select case (mod(ibx3,3))
!              case (1) ! First row of kappa
!                this%K_cols((os+ibx3/3)*2+1)=os*3+ibx3
!                this%K_vals((os+ibx3/3)*2+1)=1._wp
!                this%K_rowIdx(os*3+ibx3+1)=this%K_rowIdx(os*3+ibx3)+1
!              case (2) ! Second row of kappa
!                this%K_cols((os+ibx3/3)*2+2)=os*3+ibx3
!                this%K_vals((os+ibx3/3)*2+2)=-1._wp
!                this%K_rowIdx(os*3+ibx3+1)=this%K_rowIdx(os*3+ibx3)+1
!              case (0) ! Third row of kappa
!                this%K_rowIdx(os*3+ibx3+1)=this%K_rowIdx(os*3+ibx3)
!            end select
!          end do
!        end do
		!    | 1  0  0 |
        !  k=| 0 -1  0 |
        !    | 0  0  0 |
        maxnz=2*ntotbeadx3/3
        allocate(this%K_vals(maxnz))
        allocate(this%K_cols(maxnz))
        allocate(this%K_rowIdx(ntotbeadx3+1))
        this%K_rowIdx(1)=1
!        do ich=1, nchain
!          os=(ich-1)*nbead
!          do ibx3=1, nbeadx3
           do ibx3=1, ntotbeadx3
            select case (mod(ibx3,3))
              case (1) ! First row of kappa
                this%K_cols((ibx3/3)*2+1)=ibx3
                this%K_vals((ibx3/3)*2+1)=1._wp
                this%K_rowIdx(ibx3+1)=this%K_rowIdx(ibx3)+1
              case (2) ! Second row of kappa
                this%K_cols((ibx3/3)*2+2)=ibx3
                this%K_vals((ibx3/3)*2+2)=-1._wp
                this%K_rowIdx(ibx3+1)=this%K_rowIdx(ibx3)+1
              case (0) ! Third row of kappa
                this%K_rowIdx(ibx3+1)=this%K_rowIdx(ibx3)
            end select
          end do
!        end do
       ! For uniaxial extensional flow:
       case ('UEF')

     end select

  end subroutine init_flow_t

  !> Applies flow to the entities within the box
  !! \param Pe the Peclet number
  !! \param dt the time step size
  !! \param Rb the position vector of the beads
  subroutine apply_flow(this,Pe,dt,Rb,ntotbeadx3)

    class(flow),intent(in) :: this
    real(wp),intent(in) :: Pe,dt
    real(wp),intent(inout) :: Rb(:)
    real(wp),allocatable :: KR(:)
    integer,intent(in) :: ntotbeadx3

    allocate(KR(ntotbeadx3))
#ifdef Debuge_sequence
    write(*,*) "module:flow_mod:apply_flow"
#endif

#ifdef USE_DP
   call mkl_dcsrmv('N',ntotbeadx3,ntotbeadx3,Pe*dt,'GIIF',this%K_vals,this%K_cols,&
                   this%K_rowIdx,this%K_rowIdx(2),Rb,0._wp,KR)
#elif USE_SP
   call mkl_scsrmv('N',ntotbeadx3,ntotbeadx3,Pe*dt,'GIIF',this%K_vals,this%K_cols,&
                   this%K_rowIdx,this%K_rowIdx(2),Rb,0._wp,KR)
#endif

   Rb=Rb+KR

   deallocate(KR)

  end subroutine apply_flow

  !> Destructor for trsfm type
  subroutine del_flow_t(this)

!    use :: inp_smdlt, only:

    type(flow) :: this
#ifdef Debuge_sequence
    write(*,*) "module:flow_mod:del_flow_t"
#endif
  end subroutine del_flow_t

end module flow_mod
