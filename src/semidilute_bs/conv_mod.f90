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
! MODULE: Conversion
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2015
!
! DESCRIPTION: 
!> (1) Interchanging the configurational state Rb to its components
!! (2) {R,rc} <--> position vector Rb
!! (3) Q <--> Rb-rc
!! (4) Converting Q to segmental force Fseg
!--------------------------------------------------------------------

module conv_mod

  use :: prcn_mod

  implicit none

  !> A public type for configurational conversion
  type conf_conv

    private

  contains

  end type conf_conv

  ! Protected module variables:
  protected :: B_vals,B_cols,B_rowInd,Bbar_vals,&
               Bbar_cols,Bbar_rowInd

  !> @name Group1
  !! The parameters defined to construct sparse B
  !> @{
  real(wp),allocatable :: B_vals(:)
  integer,allocatable :: B_cols(:),B_rowInd(:)
  !> @}
  !> @name Group2
  !! The parameters defined to construct sparse Bbar
  !> @{
  real(wp),allocatable,save :: Bbar_vals(:)
  integer,allocatable,save :: Bbar_cols(:),Bbar_rowInd(:)
  !> @}

contains

  !> Initialization of conversion module
  subroutine init_conv(nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3)

    integer,intent(in) :: nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3
    integer :: jseg,ibead,maxnz_Bbar,ichain,offsetch1,offsetch2,offsetch3
    integer :: offsetchx2,offsetseg,offsetseg1,offsetseg2,isegx3,offsetbead1
    integer :: offsetbead2,maxnz_B,icoor,offsetcoor

    ! Constructing Bbar_tilde and B_tilde based on Bbar and B in DPL of Bird et al.:
    ! For making Bbar_tilde sparse (CSR):
    maxnz_Bbar=nchain*nsegx3*2
    allocate(Bbar_vals(maxnz_Bbar),Bbar_cols(maxnz_Bbar),Bbar_rowInd(ntotsegx3+1))
    Bbar_rowInd(1)=1
    do ichain=1, nchain
      offsetch1=(ichain-1)*nsegx3
      offsetch2=(ichain-1)*nbeadx3
      offsetchx2=offsetch1*2
      do isegx3=1, nsegx3
        offsetseg=(isegx3-1)*2
        Bbar_cols(offsetchx2+offsetseg+1)=offsetch2+isegx3
        Bbar_vals(offsetchx2+offsetseg+1)=-1._wp
        Bbar_cols(offsetchx2+offsetseg+2)=offsetch2+isegx3+3
        Bbar_vals(offsetchx2+offsetseg+2)=1._wp
        Bbar_rowInd(offsetch1+isegx3+1)=Bbar_rowInd(offsetch1+isegx3)+2
      end do
    end do
    ! For making B_tilde sparse (CSR):
    maxnz_B=nchain*nbeadx3*nseg
    allocate(B_vals(maxnz_B),B_cols(maxnz_B),B_rowInd(ntotbeadx3+1))
    B_rowInd(1)=1
    do ichain=1, nchain
      offsetch1=(ichain-1)*nbeadx3*nseg
      offsetch2=(ichain-1)*nsegx3
      offsetch3=(ichain-1)*nbeadx3
      do ibead=1, nbead
        offsetbead1=(ibead-1)*nseg*3
        offsetbead2=(ibead-1)*3
        do icoor=1, 3
          offsetcoor=(icoor-1)*nseg
          do jseg=1, nseg
            offsetseg1=offsetcoor+(jseg-1)
            offsetseg2=(jseg-1)*3
            if (ibead > jseg) then
              B_vals(offsetch1+offsetbead1+offsetseg1+1)=jseg/real(nbead,kind=wp)
              B_cols(offsetch1+offsetbead1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
            else
              B_vals(offsetch1+offsetbead1+offsetseg1+1)=-(1-jseg/real(nbead,kind=wp))
              B_cols(offsetch1+offsetbead1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
            end if
          end do ! jseg
          B_rowInd(offsetch3+offsetbead2+icoor+1)=B_rowInd(offsetch3+offsetbead2+icoor)+nseg
        end do ! icoor
      end do ! ibead
    end do ! ichain

  end subroutine init_conv

  !> Converting Q to R = (r-rcm) using sparse multiplication
  !! \param Q connectivity vectors
  !! \param R bead to center of mass vectors
  subroutine QtoR(Q,R,ntotsegx3,ntotbeadx3)

    integer,intent(in) :: ntotsegx3,ntotbeadx3
    real(wp),intent(in) :: Q(:)
    real(wp),intent(inout) :: R(:)

#ifdef USE_DP
    call mkl_dcsrmv('N',ntotbeadx3,ntotsegx3,1._wp,'GIIF',B_vals,B_cols,&
                    B_rowInd,B_rowInd(2),Q,0._wp,R)
#elif USE_SP
    call mkl_scsrmv('N',ntotbeadx3,ntotsegx3,1._wp,'GIIF',B_vals,B_cols,&
                    B_rowInd,B_rowInd(2),Q,0._wp,R)
#endif

  end subroutine QtoR

  !> Converting {R,rc} to (Rbx,Rby,Rbz)
  !! \param R bead to center of mass distance for all chains
  !! \param rcm center of mass for all chains
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  subroutine RtoRbc(R,rcm,Rbx,Rby,Rbz,nchain,nbead,ntotbead,ntotbeadx3)

!$  use :: omp_lib

    integer,intent(in) :: nchain,nbead,ntotbead,ntotbeadx3
    real(wp),intent(in) :: R(:),rcm(:,:)
    real(wp),intent(inout) :: Rbx(ntotbead),Rby(ntotbead),Rbz(ntotbead)
    integer :: igb,os,ich

!$omp parallel default(private) shared(ntotbead,nbead,Rbx,Rby,Rbz,R,rcm)
!$omp do schedule(auto)
    do igb=1, ntotbead
      os=(igb-1)*3
      ich=(igb-1)/nbead+1
      Rbx(igb)=R(os+1)+rcm(ich,1)
      Rby(igb)=R(os+2)+rcm(ich,2)
      Rbz(igb)=R(os+3)+rcm(ich,3)
    end do
!$omp end do
!$omp end parallel
    
  end subroutine RtoRbc

  !> Converting Rbc to Q
  !! \param Rb the position vector of the beads for all chains
  !! \param Q connectivity vectors
  subroutine RbctoQ(Rbx,Rby,Rbz,Q,bs,invbs,nseg,nbead,ntotseg)

    use :: trsfm_mod, only: eps_m,tanb,sinth,costh
    use :: flow_mod, only: FlowType
!$  use :: omp_lib

    integer,intent(in) :: ntotseg,nbead,nseg
    real(wp),intent(in) :: Rbx(:),Rby(:),Rbz(:)
    real(wp),intent(inout) :: Q(:)
    real(wp),intent(in) :: bs(3),invbs(3) 
    integer :: its,ich,oss,osb,is
    real(wp) :: qx,qy,qz,qytmp

!$omp parallel default(private) &
!$omp shared(ntotseg,nseg,nbead,Rbx,Rby,Rbz,bs,invbs,FlowType,eps_m,tanb,sinth,costh,Q)
!$omp do schedule(auto)
    do its=1, ntotseg
      ich=(its-1)/nseg+1
      oss=(ich-1)*nseg
      osb=(ich-1)*nbead
      is=its-oss
      qx=Rbx(osb+is+1)-Rbx(osb+is)
      qy=Rby(osb+is+1)-Rby(osb+is)
      qz=Rbz(osb+is+1)-Rbz(osb+is)
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
      Q((oss+is-1)*3+1)=qx
      Q((oss+is-1)*3+2)=qy
      Q((oss+is-1)*3+3)=qz
    end do
!$omp end do
!$omp end parallel

  end subroutine RbctoQ

  !> Converting {R,rc} to Rb
  !! \param R bead to center of mass distance for all chains
  !! \param rcm center of mass for all chains
  !! \param Rb the position vector of the beads for all chains
  subroutine RtoRb(R,rcm,Rb,nchain,nbead,ntotbead,ntotbeadx3)

!$  use :: omp_lib

    integer,intent(in) :: nchain,nbead,ntotbead,ntotbeadx3
    real(wp),intent(in) :: R(:)
    real(wp),intent(in) :: rcm(:,:)
    real(wp),intent(inout) :: Rb(:)
    integer :: igb,os,ich

!$omp parallel default(private) shared(ntotbead,nbead,Rb,R,rcm)
!$omp do schedule(auto)
    do igb=1, ntotbead
      os=(igb-1)*3
      ich=(igb-1)/nbead+1
      Rb(os+1)=R(os+1)+rcm(ich,1)
      Rb(os+2)=R(os+2)+rcm(ich,2)
      Rb(os+3)=R(os+3)+rcm(ich,3)
    end do
!$omp end do
!$omp end parallel

  end subroutine RtoRb

  !> Converting (Rx,Ry,Rz) to Rb
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  !! \param Rb the position vector of the beads for all chains
  subroutine RbctoRb(Rbx,Rby,Rbz,Rb,ntotbead)

!$  use :: omp_lib

    integer,intent(in) :: ntotbead
    real(wp),intent(in) :: Rbx(:),Rby(:),Rbz(:)
    real(wp),intent(inout) :: Rb(:)
    integer :: igb,os

!$omp parallel default(private) shared(ntotbead,Rbx,Rby,Rbz,Rb)
!$omp do schedule(auto)
    do igb=1, ntotbead
      os=(igb-1)*3
      Rb(os+1)=Rbx(igb)
      Rb(os+2)=Rby(igb)
      Rb(os+3)=Rbz(igb)
    end do
!$omp end do
!$omp end parallel

  end subroutine RbctoRb 

  !> Converting Rb to (Rx,Ry,Rz)
  !! \param Rb the position vector of the beads for all chains
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  subroutine RbtoRbc(Rb,Rbx,Rby,Rbz,ntotbead)

!$  use :: omp_lib

    integer,intent(in) :: ntotbead
    real(wp),intent(in) :: Rb(:)
    real(wp),intent(inout) :: Rbx(:),Rby(:),Rbz(:)
    integer :: igb,os

!$omp parallel default(private) shared(ntotbead,Rbx,Rby,Rbz,Rb)
!$omp do schedule(auto)
    do igb=1, ntotbead
      os=(igb-1)*3
      Rbx(igb)=Rb(os+1)
      Rby(igb)=Rb(os+2)
      Rbz(igb)=Rb(os+3)
    end do
!$omp end do
!$omp end parallel

  end subroutine RbtoRbc

  !> Converting Q to Fseg
  !! \param Q connectivity vectors
  !! \param Fseg spring forces
  subroutine QtoFseg(Q,Fseg,ForceLaw,qmax,ntotseg)

!$  use :: omp_lib

    integer,intent(in) :: ntotseg
    real(wp),intent(in) :: Q(:)
    real(wp),intent(inout),target :: Fseg(:)
    integer :: offset,iseg,itime
    real(wp) :: qtmp(3),qmag,F,qmax
    character(len=10) :: ForceLaw
    real(wp),pointer :: FsegP(:) => null()

!$omp parallel default(private) &
!$omp shared(itime,qmax,ForceLaw,ntotseg,Q,Fseg)
!$omp do schedule(auto)
   do iseg=1, ntotseg
     offset=3*(iseg-1)
     qtmp=Q(offset+1:offset+3)
     qmag=sqrt(dot_product(qtmp,qtmp))
     FsegP(1:3) => Fseg(offset+1:offset+3)
     select case (ForceLaw)
       case ('FENE')
         F = 1/(1-(qmag/qmax)**2)
       case ('WLC')
         F = 2*qmax/(3*qmag)*(0.25/((1-qmag/qmax)**2)-0.25+qmag/qmax)
       case ('ILCCP')
         F = (1-(qmag/qmax)**2/3)/(1-(qmag/qmax)**2)
       case ('Hookean')
         F = 1._wp
     end select
     FsegP=F*qtmp
   end do
!$omp end do
   nullify(FsegP)
!$omp end parallel

  end subroutine QtoFseg

  !> Deallocating conversion module arrays
  subroutine del_conv()

    deallocate(Bbar_vals,Bbar_cols,Bbar_rowInd)
    deallocate(B_vals,B_cols,B_rowInd)

  end subroutine del_conv

end module conv_mod
