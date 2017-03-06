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
! MODULE: hydrodynamic interaction (HI)
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: 
!> definition and initialization of HI variables
!--------------------------------------------------------------------
module hi_mod

  use :: mpi
  use :: prcn_mod
  use :: arry_mod
  use,intrinsic :: iso_c_binding
  use :: mkl_dfti
 
  implicit none

  save

  !> @name Group1
  !! The Constants for HI long range HI calculation
  !> @{
  real(wp) :: HI_a,HI_ax2,HI_ax3,HI_ato2,HI_ato3,HI_alpha,HI_alphato2,HI_alphato3
  real(wp) :: HI_Vol,HI_c0,M1_c1,M1_c2,M1_c3,M1_c4,M1_c5,M1_c6,M1_c7,M1_c8,M1_c9
  real(wp) :: M1_c10,M1_c11,M1_c12,M1_c13,M2_c1,M2_c2,M2_c3,M2_c4,Mstar_c1,Mstar_c2
  !> The maximum value of k, Fourier vector
  real(wp) :: HI_kmax
  !> The compontents of the maximum value of k, Fourier vector
  integer :: HI_kimax(3)
  !> The maximum value of kk
  integer :: HI_kikimax
  !> pi
  real(wp),parameter :: PI=3.1415926535897958648_wp
  !> 2pi
  real(wp) :: PIx2
  !> sqrt(pi)
  real(wp) :: sqrtPI
  !> @name Group2
  !! HI input variables:
  !> @{
  real(wp),protected :: hstar,errormin,rc_D,skin_D,rlist_D,rc_Dto2
  real(wp) :: HI_M,Mstart,Minc
  integer,protected :: ncols,upfactr,p_PME,p_PMEto3
  integer :: mBlLan,mubBlLan
  integer(long) :: maxNb_list_D
  character(len=10),protected :: DecompMeth,HITens
  character(len=10),protected :: HIcalc_mode,InterpMethod
  logical,protected :: mset,Dreal_sparse_mode,kmeshset
  !> @}
  !> Number of grid points in reciprocal space used in PME method.
  integer,dimension(3) :: K_mesh
  !> @name Group3
  !! Derived from the number of grid points
  !> @{
  integer,protected :: Kto3,Kcto3,Kcto3ip
  !> @}
  !> starting value of the interation number in Lanczos algorithm
  integer :: mst
  !> @name Group3
  !! Variables used for PME:
  !> @{
  type(DFTI_DESCRIPTOR),pointer :: FFTfwDescHand => null()
  type(DFTI_DESCRIPTOR),pointer :: FFTbwDescHand => null()
  integer :: FFTStatus,r_str(0:3),c_str(0:3)
  real(wp) :: FFT_scale
  character(len=DFTI_MAX_MESSAGE_LENGTH) :: err_message
  !> @}
  !> The unit 3x3 matrix
  real(wp),dimension(3,3) :: unitDelta
  !> @name Group3
  !! Arrays used in decomposition of D, when Lanczos is utilized:
  !> @{
  real(wp),allocatable :: Ybar(:)
  real(wp),allocatable,dimension(:,:) :: aBlLan,WBlLan,VBlLan,VcntBlLan
  real(wp),allocatable,dimension(:,:) :: dV0,dV1,dVk,Dsh
  !> @}
  ! Arrays for diffusion tensor when Ewald is used:
  !> Diffusion tensor
  real(wp),allocatable,target :: Diff_tens(:,:)
  !> Coefficent tensor
  real(double),allocatable :: Coeff_tens(:,:)
  ! Arrays for different parts of DF when PME is used:
  !> Real part of diffusion tensor
  real(wp),allocatable,target :: Diff_tens_real(:,:)
  !> Reciprocal part of the diffusion tensor
  real(wp),allocatable,target :: Diff_tens_recip(:,:)
  !> Brownian noise: B.dw_bl
  real(double),allocatable,target :: dw_bltmp(:,:)
  !> Brownian noise: B.dw_bl
  real(double),pointer :: dw_bltmpP(:) => null()
  !> Random vectors
  real(wp),allocatable :: dw_bl(:,:)
  !> Force on the mesh
  real(wp),allocatable,target :: F_mesh(:,:)
  !> Force on the mesh
  real(wp),allocatable,target :: F_meshtmp(:,:)
  !> x component of the force on the mesh
  real(wp),pointer :: F_meshPx(:) => null()
  !> y component of the force on the mesh
  real(wp),pointer :: F_meshPy(:) => null()
  !> z component of the force on the mesh
  real(wp),pointer :: F_meshPz(:) => null()
  !> PME tensor P values
  real(wp),allocatable :: P_vals(:)
  !> PME tensor P columns
  integer,allocatable :: P_cols(:)
  !> PME tensor P row Index
  integer,allocatable :: P_rowInd(:)
  !> D.F self part
  real(wp),allocatable :: DF_self(:)
  !> D.F real part
  real(wp),allocatable :: DF_real(:)
  !> D.F recip part
  real(wp),allocatable,target :: DF_recip(:)
  !> D.F total
  real(wp),allocatable,target :: DF_tot(:)
  ! It is assumed that pointer attribute for Dreal_vals makes it a legitimate pointee.
  ! Note that the TARGET attribute is not allowed in derived type (In current compilers).
  ! So POINTER attribute for both pointer and pointee.
  !> values of D
  real(wp),pointer :: Dreal_vals(:)
  !> columns of D
  integer,pointer :: Dreal_cols(:)
  !> row index of D
  integer,allocatable:: Dreal_rowInd(:)
  !> @name Group4
  !! Arrays for constructing the Fourier space array, m2_alpha(k):
  !> @{
  real(wp),allocatable :: m2_vec(:),m2_vec_tmp(:)
  real(wp),allocatable,dimension(:) :: mpvecx,mpvecy,mpvecz
  real(wp),allocatable,dimension(:) :: kvecx,kvecy,kvecz
  real(wp),allocatable,dimension(:) :: mpvecxtmp,mpvecytmp,mpvecztmp
  real(wp),allocatable,dimension(:) :: kvecxtmp,kvecytmp,kvecztmp
  integer,allocatable,dimension(:) :: kiuppr,kiylowr,kizlowr
  integer,allocatable,dimension(:) :: kiupprtmp,kiylowrtmp,kizlowrtmp
  !> @}
  !> @name Group5
  !! Arrays for constructing e(ik.r) in calculation of recip. space:
  !> @{
  complex(wp),allocatable,dimension(:,:) :: eikx,eiky,eikz
  complex(wp),allocatable,dimension(:,:) :: eikxtmp,eikytmp,eikztmp
  complex(wp),allocatable,dimension(:) :: eikr
  !> @}
  !> @name Group6
  !! Arrays for constructing the Cardinal B-splines exponential interpolation:
  !> @{
  complex(wp),allocatable,dimension(:) :: b_splx2,b_sply2,b_splz2
  complex(wp),allocatable,dimension(:) :: b_splx2tmp,b_sply2tmp,b_splz2tmp
  !> @}
  !> @name Group7
  !! Arrays for constructing list_D used in PME
  !> @{
  real(wp),allocatable :: Rb0(:)
  integer,allocatable :: head_D(:),LkdLst_D(:),point_D(:)
  integer,pointer :: list_D(:),list_DP(:)
  !> @}
  !> Max value of the displacement of the beads in the box 
  real(wp) :: dispmax
  !> Cell size in each direction
  real(wp) :: CellSize_D(3)
  !> number of cells in each direction
  integer :: ncells_D(3)
  !> Total number of cells, ncells_D(1)*ncells_D(2)*ncells_D(3)
  integer :: ntotcells_D
 
contains
 
  subroutine init_hi(myrank,ntotbead,ntotbeadx3,bs)

    use :: strg_mod
    use :: iso_fortran_env
    use :: flow_mod, only: FlowType

    integer,intent(in) :: myrank,ntotbead,ntotbeadx3
    real(wp),intent(in) :: bs(3)
    integer :: u1,i,j,ios,ntokens,ierr,il,stat
    character(len=1024) :: line 
    character(len=100) :: tokens(10)

    ! default values:
    hstar=0._wp
    HITens='RPY'
    HIcalc_mode='Ewald'
    Mstart=3.5_wp;Minc=0.2_wp
    rc_D=20._wp;skin_D=0.2_wp
    InterpMethod='BSpline';p_PME=4
    K_mesh=65;kmeshset=.false.
    maxNb_list_D=500000000
    DecompMeth='Cholesky'
    ncols=1
    mBlLan=3;mubBlLan=15;mset=.false.
    errormin=real(1.e-2,kind=wp)
    upfactr=50

    open (newunit=u1,action='read',file='input.dat',status='old')

    il=1
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" hi_mod: Error reading line ",i0, " Process ID ",i0)',il,myrank
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

            case ('hstar')
              call value(tokens(j+1),hstar,ios)
            case ('HITens')
              HITens=trim(adjustl(tokens(j+1)))
            case ('HIcalc-mode')
              HIcalc_mode=trim(adjustl(tokens(j+1)))
            case ('Interp-method')
              InterpMethod=trim(adjustl(tokens(j+1)))
              call value(tokens(j+2),p_PME,ios)
            case ('Mesh-size')
              call value(tokens(j+1),K_mesh(1),ios)
              call value(tokens(j+2),K_mesh(2),ios)
              call value(tokens(j+3),K_mesh(3),ios)            
              if(tokens(j+4) == 'TRUE') then
                kmeshset=.true.
              elseif(tokens(j+4) == 'FALSE') then
                kmeshset=.false.
              else
                print '(" Incorrect kmeshset.")'
                stop
              end if
            case ('M-exp')
              call value(tokens(j+1),Mstart,ios)
              call value(tokens(j+2),Minc,ios)
            case ('rc_D')
              call value(tokens(j+1),rc_D,ios)
              call value(tokens(j+2),skin_D,ios)
            case ('Nblist-size')
              call value(tokens(j+1),maxNb_list_D,ios)
            case ('Decomp-method')
              DecompMeth=trim(adjustl(tokens(j+1)))
            case ('ncols')
              call value(tokens(j+1),ncols,ios)
            case ('m')
              call value(tokens(j+1),mBlLan,ios)
              call value(tokens(j+2),mubBlLan,ios)
              if(tokens(j+3) == 'TRUE') then
                mset=.true.
              elseif(tokens(j+3) == 'FALSE') then
                mset=.false.
              end if
            case ('errormin')
              call value(tokens(j+1),errormin,ios)
            case ('upfactr')
              call value(tokens(j+1),upfactr,ios)
          end select
        enddo
      end if
    end do ef
    close(u1)

    ! This ensures the correct setting in case h*=0:
    if (hstar == 0._wp) then
      HIcalc_mode='Ewald'
      DecompMeth='Cholesky'
    end if

    PIx2=PI*2
    sqrtPI=sqrt(PI)

    ! Initialize HI parameters:
    call setHIPar(myrank,bs)

    if (HIcalc_mode == 'Ewald') then
      allocate(Diff_tens(ntotbeadx3,ntotbeadx3))
      allocate(Coeff_tens(ntotbeadx3,ntotbeadx3))
    elseif (HIcalc_mode == 'PME') then
      call setupPME(myrank,ntotbead)
      if (Dreal_sparse_mode) then
        allocate(Dreal_vals(maxNb_list_D*9),Dreal_cols(maxNb_list_D))
        allocate(Dreal_rowInd(ntotbead+1))
        allocate(DF_tot(ntotbeadx3), &
                 DF_self(ntotbeadx3),&
                 DF_real(ntotbeadx3),&
                 DF_recip(ntotbeadx3))
      else
        allocate(Diff_tens_real(ntotbeadx3,ntotbeadx3),&
                 DF_tot(ntotbeadx3), &
                 DF_self(ntotbeadx3),&
                 DF_real(ntotbeadx3),&
                 DF_recip(ntotbeadx3))
      end if
      ! Array for saving bead positions, used in PME algorithm:
      allocate(Rb0(ntotbeadx3))
      ! Cell information used for list construction:
      ncells_D(:)=max(bs(:)/rc_D,1._wp)
      CellSize_D(1:3)=bs(1:3)/ncells_D(1:3)
      ntotcells_D=ncells_D(1)*ncells_D(2)*ncells_D(3)
      ! Arrays containing cell information:
      allocate(head_D(0:ntotcells_D-1),LkdLst_D(ntotbead))
      allocate(point_D(ntotbead),list_D(maxNb_list_D))
      
      if (myrank == 0) then
        print * 
        print '(" No. cells for diffusion tensor in real space:")'
        print '(3(i10,1x))',ncells_D(:)
      end if

    end if ! HIcalc_mode

    ! Allocation of arrays for Brownian noise:
    allocate(dw_bl(ntotbeadx3,ncols),dw_bltmp(ntotbeadx3,ncols))

    call setupKSPACE(myrank,bs,ntotbead)

    if (DecompMeth == 'Lanczos') then
      allocate(aBlLan(ntotbeadx3,ncols),WBlLan(ntotbeadx3,ncols))
      allocate(VBlLan(ntotbeadx3,mBlLan*ncols),Ybar(ntotbeadx3))
      allocate(VcntBlLan(ntotbeadx3,ncols))
      ! For calculating the average of iteration number
      mst=mBlLan
    end if
    ! For cut off checking:
    rc_Dto2=rc_D*rc_D 

    unitDelta(1,1:3)=[1._wp,0._wp,0._wp]
    unitDelta(2,1:3)=[0._wp,1._wp,0._wp]
    unitDelta(3,1:3)=[0._wp,0._wp,1._wp]

  end subroutine init_hi

  subroutine update_lst(Rb,Rbtr,itime,itrst,nchain,nbead,nbeadx3,ntotbead,bs,bo)

    real(wp),intent(in) :: Rb(:),Rbtr(:),bs(3),bo(3)
    integer,intent(in) :: itime,itrst,nchain,nbead,nbeadx3,ntotbead
    logical :: update

    if (itime == itrst+1) then
      update=.true.
    else
      ! Calculate maximum displacement since last update:
      dispmax=maxval(abs(Rb-Rb0))
      ! A conservative testing of the list skin crossing:
      dispmax=2*sqrt(3*dispmax*dispmax)
      update=dispmax > (rlist_D-rc_D)
    end if
    if (update) then
      ! Save positions for next evaluations:
      Rb0=Rb
      call cnstrlst_D(Rb,Rbtr,itime,nchain,nbead,nbeadx3,ntotbead,bs,bo)
    end if

  end subroutine update_lst
  
  subroutine cnstrlst_D(Rb,Rbtr,itime,nchain,nbead,nbeadx3,ntotbead,bs,bo)

    use :: arry_mod, only: print_vector,ResizeArray
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: delrx_L

    real(wp),intent(in) :: Rb(:),Rbtr(:),bs(3),bo(3)
    integer,intent(in) :: itime,nchain,nbead,nbeadx3,ntotbead
    integer,dimension(3) :: cell_ind,neigcell_ind_p ! should contain 0~(ncells-1) cell indices
    integer,dimension(3),target :: neigcell_ind ! should contain 0~(ncells-1) cell indices
    real(wp),dimension(3) :: rij,rbtmp,rbi,rbj
    integer :: ichain,ibead,offsetch,offset,iglobbead,jglobbead
    integer :: icell,offseti,offsetj,neigcell_ID,nlist
    integer,pointer :: neigcell_indx,neigcell_indy,neigcell_indz
    real(wp) :: rlist_Dto2,rijmagto2
    integer,parameter :: EMPTY=-1

    !-----------------------------
    ! Construction of linked-list:
    !-----------------------------
    ! Reset the array head
    head_D=EMPTY
    ! Note!!: the loops are in reverse order to make them 
    ! sorted and appropriate for sparse operations.
    do ichain=nchain, 1, -1
      offsetch=(ichain-1)*nbeadx3
      do ibead=nbead, 1, -1
        iglobbead=(ichain-1)*nbead+ibead
        offset=offsetch+(ibead-1)*3
        rbtmp=Rb(offset+1:offset+3)
        select case (FlowType)
          case ('Equil')
            cell_ind(1:3)=(rbtmp(1:3)-bo(1:3))/CellSize_D(1:3)
          case ('PSF')
            cell_ind(1)=(Rbtr(iglobbead)-bo(1))/CellSize_D(1)
            cell_ind(2:3)=(rbtmp(2:3)-bo(2:3))/CellSize_D(2:3)
        end select
        icell=cell_ind(1)*ncells_D(2)*ncells_D(3)+&
              cell_ind(2)*ncells_D(3)+cell_ind(3)
        ! Link to the previous occupant to EMPTY if you are the first
        LkdLst_D(iglobbead)=head_D(icell)
        ! The previous one goes to header
        head_D(icell)=iglobbead
      end do
    end do
    !------------------------------------------------------------
    ! Construction of Verlet neighbor-list, by using linked-list:
    !------------------------------------------------------------
    neigcell_indx => neigcell_ind(1)
    neigcell_indy => neigcell_ind(2)
    neigcell_indz => neigcell_ind(3)
    rlist_Dto2=rlist_D*rlist_D
    nlist=0
    ! For the last bead, jglobbead<iglobbead, so it is not required!
    do 20 iglobbead=1, ntotbead-1
      point_D(iglobbead)=nlist+1
      offseti=(iglobbead-1)*3
      rbi=Rb(offseti+1:offseti+3)
      select case (FlowType)
        case ('Equil')
          cell_ind(1:3)=(rbi(1:3)-bo(1:3))/CellSize_D(1:3)
        case ('PSF')
          cell_ind(1)=(Rbtr(iglobbead)-bo(1))/CellSize_D(1)
          cell_ind(2:3)=(rbi(2:3)-bo(2:3))/CellSize_D(2:3)
      end select
      ! Scan the neighbouring cells:
ncix: do neigcell_indx=cell_ind(1)-1, cell_ind(1)+1
nciy:   do neigcell_indy=cell_ind(2)-1, cell_ind(2)+1
nciz:     do neigcell_indz=cell_ind(3)-1, cell_ind(3)+1
            ! Calculate the scalar neigbor cell index:
            ! Corresponding periodic cell for neighbor cell:
            neigcell_ind_p(:)=mod(neigcell_ind(:)+ncells_D(:),ncells_D(:)) 
            neigcell_ID=neigcell_ind_p(1)*ncells_D(2)*ncells_D(3)+&
                        neigcell_ind_p(2)*ncells_D(3)+neigcell_ind_p(3)
            ! Get first bead in neighbour cell:
            jglobbead=head_D(neigcell_ID)
            do 10 while (jglobbead /= EMPTY)
              offsetj=(jglobbead-1)*3
              ! Equal is important only if we have one cell.
              if ( (iglobbead < jglobbead) .or. &
                  ((iglobbead == jglobbead).and.(any(ncells_D == 1))) ) then
                ! Calculate the distance:
                ! The third term in the RHS takes into account the correction 
                ! needed for boundary cells.
                ! We could use minimum image convension, as for using Verlet 
                ! list, we need the rc_F to be smaller than box dimension.
                rbj=Rb(offsetj+1:offsetj+3)
                select case (FlowType)
                  case ('Equil')
                    rij(1:3)=rbi(1:3)-(rbj(1:3) + &
                        floor(real(neigcell_ind(1:3))/real(ncells_D(1:3)))*bs(1:3))
                  case ('PSF')
                    rij(1)=rbi(1)-(rbj(1) + &
                           floor(real(neigcell_ind(1))/ncells_D(1))*bs(1) + &
                           floor(real(neigcell_ind(2))/ncells_D(2))*delrx_L)
                    rij(2:3)=rbi(2:3)-(rbj(2:3) + &
                             floor(real(neigcell_ind(2:3))/ncells_D(2:3))*bs(2:3))
                end select
                rijmagto2=dot_product(rij,rij)
                if (rijmagto2 <= rlist_Dto2) then
                  nlist=nlist+1
                  if (nlist == maxNb_list_D) then
                    print '(" Note!!: list_D is small and will be incremented by 10")'
                    maxNb_list_D=maxNb_list_D+10
                    print '(" Dim(list_D):",1x,i15)',maxNb_list_D
                    call ResizeArray(list_D,maxNb_list_D)
                  end if
                  list_D(nlist)=jglobbead
                end if
              end if
              jglobbead=LkdLst_D(jglobbead)
10          end do
          end do nciz
        end do nciy
      end do ncix
20  end do
    point_D(ntotbead)=nlist+1
    list_DP => list_D(1:nlist)
    ! Resizing the Dreal related arrays:
    if (maxNb_list_D > size(Dreal_cols,1)) then
      call ResizeArray(Dreal_vals,maxNb_list_D*9)
      call ResizeArray(Dreal_cols,maxNb_list_D)
    end if

  end subroutine cnstrlst_D

  subroutine del_hi(myrank,bs)

    integer,intent(in) :: myrank
    real(wp),intent(in) :: bs(3)

    deallocate(dw_bl,dw_bltmp)
    if (HIcalc_mode == 'Ewald') then
      if (hstar /= 0._wp) deallocate(Diff_tens)
    elseif (HIcalc_mode == 'PME') then
      deallocate(P_vals,P_cols,P_rowInd)
      deallocate(F_mesh) 
      if (Dreal_sparse_mode) then
        deallocate(Dreal_vals,Dreal_cols,Dreal_rowInd)
        deallocate(DF_tot,DF_self,DF_real,DF_recip)
      else
        deallocate(Diff_tens_real,DF_tot,DF_self,DF_real,DF_recip)
      end if
      ! Destroying FFT handles:
      FFTStatus=DftiFreeDescriptor(FFTfwDescHand)
      if (FFTStatus /= 0) then
        print '(" In process ID: ",i0)',myrank
        print *,'Error!!: Problem in DftiFreeDescriptor, forward.'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if
      FFTStatus=DftiFreeDescriptor(FFTbwDescHand)
      if (FFTStatus /= 0) then
        print '(" In process ID: ",i0)',myrank
        print *,'Error!!: Problem in DftiFreeDescriptor, backward.'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if       
    end if
    if (DecompMeth == 'Lanczos') then
      deallocate(aBlLan,WBlLan,VBlLan,Ybar,VcntBlLan)
    end if
    
    if ((HIcalc_mode == 'Ewald') .and. (hstar /= 0._wp)) deallocate(Coeff_tens)
    if (HIcalc_mode == 'PME') deallocate(head_D,LkdLst_D,point_D,list_D,Rb0)

    if (hstar /= 0._wp) deallocate(m2_vec)
    if (HIcalc_mode == 'PME') then
      deallocate(mpvecx,mpvecy,mpvecz)
    elseif ((HIcalc_mode == 'Ewald') .and. (hstar /= 0._wp)) then
      if ((bs(1) == bs(2)) .and. &
          (bs(2) == bs(3))) then
        deallocate(kiuppr,kiylowr,kizlowr)
      end if
    end if

  end subroutine del_hi

! This routine initializes the constants that are going to be used in HI calculator 
  subroutine setHIPar(myrank,BoxDim)
  
    use :: arry_mod, only: print_vector

    integer,intent(in) :: myrank
    real(wp),intent(in) :: BoxDim(3)
    integer :: kiki

    HI_M=Mstart ! assumtion: exp(-M^2) << 1
    ! For Rotne-Prager-Yamakawa Tensor; ewald_Beenakar_Zhou
    HI_a=sqrtPI*hstar
    HI_ax2=HI_a*2
    HI_ax3=HI_a*3
    HI_ato2=HI_a**2
    HI_ato3=HI_a**3
!    rc_D=10._wp!HI_M/HI_alpha
!    print *,'rcd',rc_d
    if (HIcalc_mode == 'PME') then
      rlist_D=rc_D+skin_D*rc_D ! User specified
      if (all(rlist_D >= BoxDim(:)/2)) then
        Dreal_sparse_mode=.false.
      else
        Dreal_sparse_mode=.true.
      end if
    end if
    HI_alpha=HI_M/rc_D ! From jain et al.
    HI_alphato2=HI_alpha**2
    HI_alphato3=HI_alpha**3
    HI_Vol=BoxDim(1)*BoxDim(2)*BoxDim(3)
    HI_c0=1 - 6*HI_a*HI_alpha/sqrtPI + 40*HI_ato3*HI_alphato3/(3*sqrtPI)
!
    M1_c1=HI_ax3/4
    M1_c2=HI_ato3/2
    M1_c3=HI_ax3*HI_alphato3
    M1_c4=9*HI_a*HI_alpha/2
    M1_c5=4*HI_ato3*HI_alpha**7
    M1_c6=20*HI_ato3*HI_alpha**5
    M1_c7=14*HI_ato3*HI_alphato3
    M1_c8=HI_ato3*HI_alpha
    M1_c9=3*M1_c2
    M1_c10=HI_ax3*HI_alpha/2
    M1_c11=16*HI_ato3*HI_alpha**5
    M1_c12=2*HI_ato3*HI_alphato3
    M1_c13=3*HI_ato3*HI_alpha
!
    M2_c1=HI_ato3/3
    M2_c2=1/(4*HI_alphato2)
    M2_c3=1/(8*HI_alpha**4)
    M2_c4=6*PI/HI_Vol
! 
    Mstar_c1=9/(32*HI_a)
    Mstar_c2=Mstar_c1/3
! 
    HI_kmax=2*HI_M**2/rc_D ! From Jain et al.
    HI_kimax(1:3)=HI_kmax*BoxDim(1:3)/PIx2
    if (myrank == 0) then
      print *
      print "(1x,a,f10.4)",'kmax:',HI_kmax
      print "(1x,a,3(i4,1x))",'kimax:',HI_kimax(1:3)
    end if
    if (HIcalc_mode == 'Ewald') then
      if ((BoxDim(1) == BoxDim(2)) .and. (BoxDim(2) == BoxDim(3))) then ! Equal box length
        HI_kikimax=floor((HI_kmax*maxval(BoxDim)/PIx2)**2)
        allocate(kiuppr(0:HI_kikimax),kiylowr(0:HI_kikimax),kizlowr(0:HI_kikimax))
        do kiki=0, HI_kikimax
          kiuppr(kiki)=int(sqrt(real(HI_kikimax-kiki,kind=wp)))
          kiylowr(kiki)=-kiuppr(kiki)
          kizlowr(kiki)=-kiuppr(kiki)
        end do
        kiylowr(0)=0
        kizlowr(0)=1
      
       end if ! BoxDim ...
!call print_vector(kiuppr,'ku')
!call print_vector(kiylowr,'kyl')
!call print_vector(kizlowr,'kzl')

!    print *,'kikim',HI_kikimax
!    call print_vector(kiuppr,'ku')
!    call print_vector(kiylowr,'k2l')
!    call print_vector(kizlowr,'k3l')
!      HI_nmax=1.0!HI_kmax*BoxDim(1)/(2*PI)

    elseif (HIcalc_mode == 'PME') then
      if (.not.kmeshset) K_mesh=2*HI_kimax+1
      if (myrank == 0) then
        print *
        print "(1x,a,3(i4,1x))",'K_mesh:',K_mesh(1:3)
        print "(1x,a,2(f10.4,1x))",'rc_D,rlist_D:',rc_D,rlist_D
        print *
        if (Dreal_sparse_mode) then
          print *, 'Note!!: The calculation of D_real is performed with Sparse Matrix Operation.'
        else
          print *, 'Note!!: The calculation of D_real is performed with Dense Matrix Operation.'
        end if
      end if
    end if ! HIcalc_mode

  end subroutine setHIPar

  subroutine updateHIpar(myrank,BoxDim,ntotbead)

    integer,intent(in) :: ntotbead,myrank
    real(wp),intent(in) :: BoxDim(3)
    integer :: kiki

    HI_M=HI_M+Minc
!   These parameter will also change:
    HI_alpha=HI_M/rc_D ! From jain et al.
    HI_alphato2=HI_alpha**2
    HI_alphato3=HI_alpha**3
    HI_c0=1 - 6*HI_a*HI_alpha/sqrtPI + 40*HI_ato3*HI_alphato3/(3*sqrtPI)
!
    M1_c3=HI_ax3*HI_alphato3
    M1_c4=9*HI_a*HI_alpha/2
    M1_c5=4*HI_ato3*HI_alpha**7
    M1_c6=20*HI_ato3*HI_alpha**5
    M1_c7=14*HI_ato3*HI_alphato3
    M1_c8=HI_ato3*HI_alpha
    M1_c10=HI_ax3*HI_alpha/2
    M1_c11=16*HI_ato3*HI_alpha**5
    M1_c12=2*HI_ato3*HI_alphato3
    M1_c13=3*HI_ato3*HI_alpha
! 
    M2_c2=1/(4*HI_alphato2)
    M2_c3=1/(8*HI_alpha**4)
!
    HI_kmax=2*HI_M**2/rc_D ! from Jain et al.
    HI_kimax(1:3)=HI_kmax*BoxDim(1:3)/PIx2
!
    if (HIcalc_mode == 'Ewald') then
      if ((BoxDim(1) == BoxDim(2)) .and. (BoxDim(2) == BoxDim(3))) then ! Equal box length
        HI_kikimax=floor((HI_kmax*maxval(BoxDim)/PIx2)**2)
        allocate(kiupprtmp(0:HI_kikimax),kiylowrtmp(0:HI_kikimax),kizlowrtmp(0:HI_kikimax))
        call move_alloc(from=kiupprtmp,to=kiuppr)
        call move_alloc(from=kiylowrtmp,to=kiylowr)
        call move_alloc(from=kizlowrtmp,to=kizlowr)
        do kiki=0, HI_kikimax
          kiuppr(kiki)=int(sqrt(real(HI_kikimax-kiki,kind=wp)))
          kiylowr(kiki)=-kiuppr(kiki)
          kizlowr(kiki)=-kiuppr(kiki)
        end do
        kiylowr(0)=0
        kizlowr(0)=1
      end if ! BoxDim ...
    end if ! HIcalc_mode
!
    print *
    print '( "M has been updated to ",f10.4," in process ID: ",i0)',HI_M,myrank
    if (HI_M > 15) then 
      print *
      print'(" M is too large, The program is terminated. Process ID: ",i0)',myrank
      stop
    else
      print "(1x,a,f10.4,1x,i0)",'kmax, process ID:',HI_kmax,myrank
      print "(1x,a,3(i4,1x),1x,i0)",'kimax, process ID:',HI_kimax(1:3),myrank
      if (HIcalc_mode == 'PME') then
        if (.not.kmeshset) K_mesh=2*HI_kimax+1
        print "(1x,a,3(i4,1x),1x,i0)",'K_mesh, process ID:',K_mesh(1:3),myrank
      end if
    end if
 
    if (HIcalc_mode == 'PME') call setupPME(myrank,ntotbead,reset=.true.)
    call setupKSPACE(myrank,BoxDim,ntotbead,reset=.true.)    
 
  end subroutine updateHIpar

  subroutine restartHIpar(myrank,BoxDim,ntotbead)

    implicit none
    real(wp) :: BoxDim(3)
    integer :: ntotbead,kiki,myrank

    HI_M=Mstart
!   These parameter will also change:
    HI_alpha=HI_M/rc_D ! From jain et al.
    HI_alphato2=HI_alpha**2
    HI_alphato3=HI_alpha**3
    HI_c0=1 - 6*HI_a*HI_alpha/sqrtPI + 40*HI_ato3*HI_alphato3/(3*sqrtPI)
!
    M1_c3=HI_ax3*HI_alphato3
    M1_c4=9*HI_a*HI_alpha/2
    M1_c5=4*HI_ato3*HI_alpha**7
    M1_c6=20*HI_ato3*HI_alpha**5
    M1_c7=14*HI_ato3*HI_alphato3
    M1_c8=HI_ato3*HI_alpha
    M1_c10=HI_ax3*HI_alpha/2
    M1_c11=16*HI_ato3*HI_alpha**5
    M1_c12=2*HI_ato3*HI_alphato3
    M1_c13=3*HI_ato3*HI_alpha
! 
    M2_c2=1/(4*HI_alphato2)
    M2_c3=1/(8*HI_alpha**4)
!
    HI_kmax=2*HI_M**2/rc_D ! from Jain et al.
    HI_kimax(1:3)=HI_kmax*BoxDim(1:3)/PIx2

    if (HIcalc_mode == 'Ewald') then
      if ((BoxDim(1) == BoxDim(2)) .and. (BoxDim(2) == BoxDim(3))) then ! Equal box length
        HI_kikimax=floor((HI_kmax*maxval(BoxDim)/PIx2)**2)
        allocate(kiupprtmp(0:HI_kikimax),kiylowrtmp(0:HI_kikimax),kizlowrtmp(0:HI_kikimax))
        call move_alloc(from=kiupprtmp,to=kiuppr)
        call move_alloc(from=kiylowrtmp,to=kiylowr)
        call move_alloc(from=kizlowrtmp,to=kizlowr)
        do kiki=0, HI_kikimax
          kiuppr(kiki)=int(sqrt(real(HI_kikimax-kiki,kind=wp)))
          kiylowr(kiki)=-kiuppr(kiki)
          kizlowr(kiki)=-kiuppr(kiki)
        end do
        kiylowr(0)=0
        kizlowr(0)=1
      end if ! BoxDim ...
    end if ! HIcalc_mode

    print *
    print '( "M has been restarted to ",f10.4," in process ID: ",i0)',HI_M,myrank
    print "(1x,a,f12.5,1x,i0)",'kmax, process ID:',HI_kmax,myrank
    print "(1x,a,3(i4,1x),1x,i0)",'kimax, process ID:',HI_kimax(1:3),myrank
    if (HIcalc_mode == 'PME') then
      if (.not.kmeshset) K_mesh=2*HI_kimax+1
      print "(1x,a,3(i4,1x),1x,i0)",'K_mesh, process ID:',K_mesh(1:3),myrank
    end if

    if (HIcalc_mode == 'PME') call setupPME(myrank,ntotbead,reset=.true.)
    call setupKSPACE(myrank,BoxDim,ntotbead,reset=.true.)
 
  end subroutine restartHIpar

  subroutine setupPME(myrank,ntotbead,reset)

    integer :: i,k,ibead,ntotbead,myrank
    logical,optional,intent(in) :: reset

    p_PMEto3=p_PME*p_PME*p_PME
    if (.not.allocated(P_vals)) allocate(P_vals(ntotbead*p_PMEto3))
    if (.not.allocated(P_cols)) allocate(P_cols(ntotbead*p_PMEto3))
    if (.not.allocated(P_rowInd)) allocate(P_rowInd(ntotbead+1))
    P_rowInd(1)=0
    do ibead=1, ntotbead
      P_rowInd(ibead+1)=P_rowInd(ibead)+p_PMEto3
    end do
    Kto3=K_mesh(1)*K_mesh(2)*K_mesh(3)
    Kcto3=2*(K_mesh(1)/2+1)*K_mesh(2)*K_mesh(3)
    if (present(reset)) then
      if (reset) then
        nullify(F_meshPx,F_meshPy,F_meshPz)
        allocate(F_meshtmp(0:Kcto3-1,3))
        call move_alloc(from=F_meshtmp,to=F_mesh)
        allocate(b_splx2tmp(0:K_mesh(1)-1),b_sply2tmp(0:K_mesh(2)-1),b_splz2tmp(0:K_mesh(3)-1))
        call move_alloc(from=b_splx2tmp,to=b_splx2)
        call move_alloc(from=b_sply2tmp,to=b_sply2)
        call move_alloc(from=b_splz2tmp,to=b_splz2)
      else
        allocate(F_mesh(0:Kcto3-1,3))
        allocate(b_splx2(0:K_mesh(1)-1),b_sply2(0:K_mesh(2)-1),b_splz2(0:K_mesh(3)-1))
      end if
    else
      allocate(F_mesh(0:Kcto3-1,3))
      allocate(b_splx2(0:K_mesh(1)-1),b_sply2(0:K_mesh(2)-1),b_splz2(0:K_mesh(3)-1))
    end if

    F_meshPx(0:Kcto3-1) => F_mesh(:,1)
    F_meshPy(0:Kcto3-1) => F_mesh(:,2)
    F_meshPz(0:Kcto3-1) => F_mesh(:,3)

    do i=0, K_mesh(1)-1
      b_splx2(i)=b_spl(i,p_PME,K_mesh(1))*conjg(b_spl(i,p_PME,K_mesh(1)))
    end do

    do i=0, K_mesh(2)-1
      b_sply2(i)=b_spl(i,p_PME,K_mesh(2))*conjg(b_spl(i,p_PME,K_mesh(2)))
    end do

    do i=0, K_mesh(3)-1
      b_splz2(i)=b_spl(i,p_PME,K_mesh(3))*conjg(b_spl(i,p_PME,K_mesh(3)))
    end do

!   For In-Place:
    r_str=[0, 1, 2*(K_mesh(1)/2+1), 2*(K_mesh(1)/2+1)*K_mesh(2)]
    c_str=[0, 1, (K_mesh(1)/2+1), (K_mesh(1)/2+1)*K_mesh(2)]
!   For Not In-Place:
!    r_str=[0, 1, K_mesh(1), K_mesh(1)*K_mesh(2)]
!    c_str=[0, 1, K_mesh(1)/2+1, (K_mesh(1)/2+1)*K_mesh(2)]

!   -----------------------
!   >>>>>> For forward FFT:
!   -----------------------
    FFTStatus=DftiCreateDescriptor(FFTfwDescHand,DFTI_DOUBLE,DFTI_REAL,3,[K_mesh(1),K_mesh(2),K_mesh(3)])
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiCreateDescriptor in init_pme, forward.")'
      print '("  Error, status = ",i0)',FFTStatus
      stop
    end if 
    FFTStatus=DftiSetValue(FFTfwDescHand,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiSetValue-CCE in init_pme, forward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if
    FFTStatus=DftiSetValue(FFTfwDescHand,DFTI_INPUT_STRIDES,r_str)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '("Error!!: Problem in DftiSetValue-output-stride in init_pme, forward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if 
    FFTStatus=DftiSetValue(FFTfwDescHand,DFTI_OUTPUT_STRIDES,c_str)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiSetValue-output-stride in init_pme, forward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if
    FFTStatus=DftiCommitDescriptor(FFTfwDescHand)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiCommitDescriptor in init_pme, forward.")'
      print '("  Error, status = ",i0)', FFTStatus
      print *, DftiErrorClass(FFTStatus,DFTI_NO_ERROR)
      print *, DftiErrorMessage(FFTStatus)
      stop
    end if
!   -----------------------
!   >>>>> For backward FFT:
!   -----------------------
    FFTStatus=DftiCreateDescriptor(FFTbwDescHand,DFTI_DOUBLE,DFTI_REAL,3,[K_mesh(1),K_mesh(2),K_mesh(3)])
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiCreateDescriptor in init_pme, backward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if    
    FFTStatus=DftiSetValue(FFTbwDescHand,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiSetValue-CCE in init_pme, backward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if
    FFTStatus=DftiSetValue(FFTbwDescHand,DFTI_INPUT_STRIDES,c_str)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiSetValue-output-stride in init_pme, backward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if 
    FFTStatus=DftiSetValue(FFTbwDescHand,DFTI_OUTPUT_STRIDES,r_str)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiSetValue-output-stride in init_pme, backward.")'
      print '("  Error, status = ",i0)', FFTStatus
      stop
    end if 
    FFTStatus=DftiCommitDescriptor(FFTbwDescHand)
    if (FFTStatus /= 0) then
      print '(" In process ID: ",i0)',myrank
      print '(" Error!!: Problem in DftiCommitDescriptor in init_pme, backward.")'
      print '("  Error, status = ",i0)', FFTStatus
      print *, DftiErrorClass(FFTStatus,DFTI_NO_ERROR)
      print *, DftiErrorMessage(FFTStatus)
      stop
    end if

  end subroutine setupPME

  subroutine setupKSPACE(myrank,BoxDim,ntotbead,reset)

    use :: flow_mod, only: FlowType
    use :: arry_mod, only: print_vector
    use :: trsfm_mod, only: eps_m

    integer,intent(in) :: myrank,ntotbead
    real(wp),intent(in) :: BoxDim(3)
    logical,optional,intent(in) :: reset
    integer :: ktot,kix,kiy,kiz,kisqmax,kiymin,kizmin,mivecx,mivecy,mivecz
    integer :: mpivecx,mpivecy,mpivecz,mtot,kiydev,kikix,kiyy,kikixy,ierr
    real(wp) :: kvec(3),k,kto2,mpvec(3)
    
    if (HIcalc_mode == 'Ewald') then       
      if (present(reset)) then
        if (reset) then
          allocate(eikxtmp(0:HI_kimax(1),ntotbead),           &
                   eikytmp(-HI_kimax(2):HI_kimax(2),ntotbead),&
                   eikztmp(-HI_kimax(3):HI_kimax(3),ntotbead))
          allocate(kvecxtmp(0:HI_kimax(1)),           &
                   kvecytmp(-HI_kimax(2):HI_kimax(2)),&
                   kvecztmp(-HI_kimax(3):HI_kimax(3)))
          allocate(m2_vec_tmp( 4*HI_kimax(1)*HI_kimax(2)*HI_kimax(3) +                                           &
                               2*(HI_kimax(1)*HI_kimax(2) + HI_kimax(2)*HI_kimax(3) + HI_kimax(1)*HI_kimax(3)) + &
                               1*(HI_kimax(1) + HI_kimax(2) + HI_kimax(3)) ))
          call move_alloc(from=kvecxtmp,to=kvecx)
          call move_alloc(from=kvecytmp,to=kvecy)
          call move_alloc(from=kvecztmp,to=kvecz)
          call move_alloc(from=m2_vec_tmp,to=m2_vec)
          call move_alloc(from=eikxtmp,to=eikx)
          call move_alloc(from=eikytmp,to=eiky)
          call move_alloc(from=eikztmp,to=eikz)
        else
          allocate(kvecx(0:HI_kimax(1)),           &
                   kvecy(-HI_kimax(2):HI_kimax(2)),&
                   kvecz(-HI_kimax(3):HI_kimax(3)))
          allocate(eikx(0:HI_kimax(1),ntotbead),           &
                   eiky(-HI_kimax(2):HI_kimax(2),ntotbead),&
                   eikz(-HI_kimax(3):HI_kimax(3),ntotbead))
          allocate(m2_vec( 4*HI_kimax(1)*HI_kimax(2)*HI_kimax(3) +                                           &
                           2*(HI_kimax(1)*HI_kimax(2) + HI_kimax(2)*HI_kimax(3) + HI_kimax(1)*HI_kimax(3)) + &
                           1*(HI_kimax(1) + HI_kimax(2) + HI_kimax(3)) ))
        end if
      else ! .not.present(reset)
        allocate(kvecx(0:HI_kimax(1)),           &
                 kvecy(-HI_kimax(2):HI_kimax(2)),&
                 kvecz(-HI_kimax(3):HI_kimax(3)))
        allocate(eikx(0:HI_kimax(1),ntotbead),           &
                 eiky(-HI_kimax(2):HI_kimax(2),ntotbead),&
                 eikz(-HI_kimax(3):HI_kimax(3),ntotbead))
        allocate(m2_vec( 4*HI_kimax(1)*HI_kimax(2)*HI_kimax(3) +                                           &
                         2*(HI_kimax(1)*HI_kimax(2) + HI_kimax(2)*HI_kimax(3) + HI_kimax(1)*HI_kimax(3)) + &
                         1*(HI_kimax(1) + HI_kimax(2) + HI_kimax(3)) ))
      end if
!     kvecx~z will be used in Diffcalc_mod. 
!     Note!!: that the formula is correct for equal or non-equal box length.
      do kix=0, HI_kimax(1)
        kvecx(kix)=PIx2/BoxDim(1)*kix
      end do
      do kiy=-HI_kimax(2), HI_kimax(2)
        kvecy(kiy)=PIx2/BoxDim(2)*kiy
      end do
      do kiz=-HI_kimax(3), HI_kimax(3)
        kvecz(kiz)=PIx2/BoxDim(3)*kiz
      end do

      if (FlowType == 'Equil') then
        if ((BoxDim(1) /= BoxDim(2)) .or. (BoxDim(2) /= BoxDim(3))) then ! Unequal box length
          ktot=0
          do kix=0, HI_kimax(1)
            kvec(1)=PIx2/BoxDim(1)*kix
!           because of inversion symmetry:
            if (kix == 0) then
              kiymin=0
            else
              kiymin=-HI_kimax(2)
            end if
            do kiy=kiymin, HI_kimax(2)
              kvec(2)=PIx2/BoxDim(2)*kiy
!             because of inversion symmetry:
              if ((kix == 0) .and. (kiy == 0)) then
                kizmin=0
              else
                kizmin=-HI_kimax(3)
              end if
              do kiz=kizmin, HI_kimax(3)
                kvec(3)=PIx2/BoxDim(3)*kiz
                kto2=dot_product(kvec,kvec)
                if ((kto2 < (HI_kmax*HI_kmax)) .and. (kto2 /= 0)) then
!                  print *,'ki',kix,kiy,kiz
                  ktot=ktot+1
                  m2_vec(ktot)=m2_alpha(kvec)
                end if
              end do ! kiz
            end do ! kiy
          end do ! kix
        else ! Equal box length
          ktot=0
kx:       do kix=0, kiuppr(0)
            kikix=kix*kix
            kvec(1)=kvecx(kix)
ky:         do kiy=kiylowr(kikix), kiuppr(kikix)
!              kiy=kiyy+kiydev
              kikixy=kikix+kiy*kiy
              kvec(2)=kvecy(kiy)
kz:           do kiz=kizlowr(kikixy), kiuppr(kikixy)
                kvec(3)=kvecz(kiz)
!                print *,'ki2',kix,kiy,kiz,kiydev
!                kvec=[kvecx(kix),kvecy(kiyy)+kvecy(kiydev),kvecz(kiz)]
                ktot=ktot+1
                m2_vec(ktot)=m2_alpha(kvec)
              end do kz
            end do ky
          end do kx
        end if ! BoxDim ...

        if (present(reset) .and. (reset)) then
          print '(" Reciprocal space set up complete in Process ID: ",i0)',myrank
          print '(" Number of reciprocal vectors: ",i8)',ktot
        else
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          if (myrank == 0) then
            print '(" Reciprocal space set up complete.")'
            print '(" Number of reciprocal vectors: ",i8)',ktot
          end if
        end if

      end if ! FlowType == 'Equil'

    elseif (HIcalc_mode == 'PME') then

!     Note!!!!!: This part should be changed in case K1,K2,K3 are different.
!      if ((K_mesh(1) /= K_mesh(2)) .or. (K_mesh(2) /= K_mesh(3))) then
!        print *,'(" Warning!!: For different K_mesh, a part in GlobalData.f90 should be changed.")'
!        stop
!      end if
      if (present(reset)) then
        if (reset) then
          allocate(mpvecxtmp(0:K_mesh(1)-1))
          allocate(mpvecytmp(0:K_mesh(2)-1))
          allocate(mpvecztmp(0:K_mesh(3)-1))
          allocate(m2_vec_tmp((K_mesh(1)/2+1)*K_mesh(2)*K_mesh(3)))
          call move_alloc(from=mpvecxtmp,to=mpvecx)
          call move_alloc(from=mpvecytmp,to=mpvecy)
          call move_alloc(from=mpvecztmp,to=mpvecz)
          call move_alloc(from=m2_vec_tmp,to=m2_vec)
        else
          allocate(mpvecx(0:K_mesh(1)-1))
          allocate(mpvecy(0:K_mesh(2)-1))
          allocate(mpvecz(0:K_mesh(3)-1))
          allocate(m2_vec((K_mesh(1)/2+1)*K_mesh(2)*K_mesh(3)))
        end if
      else          
        allocate(mpvecx(0:K_mesh(1)-1))
        allocate(mpvecy(0:K_mesh(2)-1))
        allocate(mpvecz(0:K_mesh(3)-1))
        allocate(m2_vec((K_mesh(1)/2+1)*K_mesh(2)*K_mesh(3)))
      end if
      mtot=0
mz:   do mivecz=0, K_mesh(3)-1
        if (mivecz <= K_mesh(3)/2) then
          mpivecz=mivecz
        else
          mpivecz=mivecz-K_mesh(3)
        end if
        mpvec(3)=PIx2/BoxDim(3)*mpivecz
        mpvecz(mivecz)=mpvec(3)
my:     do mivecy=0, K_mesh(2)-1
          if (mivecy <= K_mesh(2)/2) then
            mpivecy=mivecy
          else
            mpivecy=mivecy-K_mesh(2)
          end if
          mpvec(2)=PIx2/BoxDim(2)*mpivecy
          if (FlowType == 'Equil') mpvecy(mivecy)=mpvec(2)
mx:       do mivecx=0, K_mesh(1)/2
            mpivecx=mivecx
            mpvec(1)=PIx2/BoxDim(1)*mpivecx
            mpvecx(mivecx)=mpvec(1)
            if (FlowType == 'PSF') then
              mpvec(2)=mpvec(2)-mpvec(1)*eps_m
              mpvecy(mivecy)=mpvec(2)
            end if
            if (.not.all(mpvec == 0)) then
              mtot=mtot+1
              m2_vec(mtot)=b_splx2(mivecx)*b_sply2(mivecy)*b_splz2(mivecz)*m2_alpha(mpvec)
            end if
          end do mx
        end do my
      end do mz

      if (present(reset) .and. (reset)) then
        print '(" Reciprocal space set up complete in Process ID: ",i0)',myrank
        print '(" Number of reciprocal vectors: ",i8)',mtot
      else
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if (myrank == 0) then
          print '(" Reciprocal space set up complete.")'
          print '(" Number of reciprocal vectors: ",i8)',mtot
        end if
      end if

    end if ! HIcalc_mode
 
  end subroutine setupKSPACE

  function m2_alpha(kvec)
    implicit none
    real(wp) :: m2_alpha,kvec(3),kto2

    kto2=dot_product(kvec,kvec)

    m2_alpha=(HI_a-M2_c1*kto2)*(1+M2_c2*kto2+M2_c3*kto2*kto2)*(M2_c4/kto2)*exp(-M2_c2*kto2)
  end function m2_alpha

  subroutine strupdateKSPACE(BoxDim)

    use :: trsfm_mod, only: eps_m

    integer :: ktot,kix,kiy,kiz,kiydev,kikix,kiyy,kikixy,kiymin,kizmin
    real(wp) :: BoxDim(3),kvec(3),kto2
    
    if (HIcalc_mode == 'Ewald') then

      ! Unequal box length
      if ((BoxDim(1) /= BoxDim(2)) .or. (BoxDim(2) /= BoxDim(3))) then

        ktot=0
        do kix=0, HI_kimax(1)
          kiydev=nint(eps_m*kix)
          kvec(1)=kvecx(kix)
          ! because of inversion symmetry:
          if (kix == 0) then
            kiymin=0
          else
            kiymin=-HI_kimax(2)
          end if
          do kiyy=kiymin, HI_kimax(2)
!            kvec(2)=kvecy(kiyy)+kvecy(kiydev)
            kvec(2)=kvecy(kiyy)+PIx2/BoxDim(2)*kiydev
            ! because of inversion symmetry:
            if ((kix == 0) .and. (kiyy == 0)) then
              kizmin=0
            else
              kizmin=-HI_kimax(3)
            end if
            do kiz=kizmin, HI_kimax(3)
              kvec(3)=kvecz(kiz)
              kto2=dot_product(kvec,kvec)
              if ((kto2 < (HI_kmax*HI_kmax)) .and. (kto2 /= 0)) then
!                  print *,'ki',kix,kiy,kiz
                ktot=ktot+1
                m2_vec(ktot)=m2_alpha(kvec)
              end if
            end do ! kiz
          end do ! kiy
        end do ! kix

      else ! Equal box length

        ktot=0
kx:     do kix=0, kiuppr(0)
          kiydev=nint(eps_m*kix)
          kikix=kix*kix
          kvec(1)=kvecx(kix)
ky:       do kiyy=kiylowr(kikix), kiuppr(kikix)
!            kiy=kiyy+kiydev
            kikixy=kikix+kiyy*kiyy
            kvec(2)=kvecy(kiyy)+kvecy(kiydev)
kz:         do kiz=kizlowr(kikixy), kiuppr(kikixy)
              kvec(3)=kvecz(kiz)
              ktot=ktot+1
              m2_vec(ktot)=m2_alpha(kvec)
            end do kz
          end do ky
        end do kx

      end if ! BoxDim ...

    elseif (HIcalc_mode == 'PME') then
      ! Not yet considered!
    end if

  end subroutine strupdateKSPACE

  ! Calculatioan of Cardinal B-splines:
  ! Note!!: writing complex variables in the form of cmplx(..,..) results in higher accuracy 
  !  compared to exp(i..). Also having kind=wp is very important for having enough precision.
  complex(wp) function b_spl(m,p,K_mesh)
    implicit none
    integer ::m,p,K_mesh,i,k
    complex(wp) :: sum_spl
    
    sum_spl=(0.0_wp,0.0_wp)
    do k=0, p-2
      sum_spl=sum_spl+M_spl(real(k+1,kind=wp),p)*cmplx( cos(PIx2*m*k/K_mesh),sin(PIx2*m*k/K_mesh),kind=wp )
    end do ! k
    b_spl=cmplx( cos(PIx2*(p-1)*m/K_mesh),sin(PIx2*(p-1)*m/K_mesh),kind=wp )/sum_spl

  end function b_spl

  ! This function calculated the B-spline function with order p at x.
  recursive function M_spl(x,p) result(M_res)
    implicit none
    integer :: p
    real(wp) :: x,M_res
    if (p < 2) then
      stop 'Warning!!: p should be >= 2'
    end if
    if (p == 2) then
      if ((x < 0) .or. (x > 2)) then
        M_res=0.0_wp
      else
        M_res=1-abs(x-1)
      end if
    else
      M_res=x/(p-1)*M_spl(x,p-1) + (p-x)/(p-1)*M_spl(x-1,p-1)
    end if
  end function M_spl

  function W_Lag(x,p,MP) result(W_res)

    integer :: p,MP
    real(wp) :: x,W_res

    select case(p)
      case(1)

        W_res=1.0_wp

      case(2)

        select case (MP)
          case(1)
            W_res=0.5*(1-2*x)
          case(2)
            W_res=0.5*(1+2*x)
        end select ! MP
    
      case(3)

        select case (MP)
          case(1)
            W_res=(x*x-x)/2
          case(2)
            W_res=(-2*x*x+2)/2
          case(3)
            W_res=(x*x+x)/2
        end select ! MP
      
      case(4) 
   
        select case (MP)
          case(1)
            W_res=(-8*x*x*x+12*x*x+2*x-3)/48
          case(2)
            W_res=(24*x*x*x-12*x*x-54*x+27)/48
          case(3)
            W_res=(-24*x*x*x-12*x*x+54*x+27)/48
          case(4)
            W_res=(8*x*x*x+12*x*x-2*x-3)/48
        end select ! MP

      case default
 
        print *,'Warning!!: Unsupported degree of Lagrange-polynomial.'
        stop

    end select
         
  end function W_Lag

end module hi_mod
