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

  ! ! B for comb polymer
  ! real(wp),allocatable :: B_vals_cmb(:)
  ! integer,allocatable :: B_cols_cmb(:),B_rowInd_cmb(:)
  ! ! Bbar for comb polymer
  ! real(wp),allocatable :: Bbar_vals_cmb(:)
  ! integer,allocatable :: Bbar_cols_cmb(:),Bbar_rowInd_cmb(:)

  !> the max number of nonzero elements for B
  integer,save :: maxnz_B
  !> the max number of nonzero elements for Bbar
  integer,save :: maxnz_Bbar


contains

  !> Initialization of conversion module
  subroutine init_conv(nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3,&
    add_cmb,nchain_cmb,nseg_cmb,nseg_cmbbb,nseg_cmbar,Na,Ia)

    use :: arry_mod, only: print_vector,print_matrix

    integer,intent(in) :: nchain,nseg,nbead,nsegx3,nbeadx3,ntotsegx3,ntotbeadx3

    logical,intent(in) :: add_cmb
    integer,intent(in) :: nchain_cmb,nseg_cmb,nseg_cmbbb,nseg_cmbar,Na,Ia(:)

    integer :: jseg,ibead,ichain,offsetch1,offsetch2,offsetch3
    integer :: offsetchx2,offsetseg,offsetseg1,offsetseg2,isegx3,offsetbead1
    integer :: offsetbead2,icoor,offsetcoor


    integer :: nbead_cmb,ntotseg_cmb,iarm,kseg,kbead,ntotbead_cmb
    real(wp) :: fctr

    integer :: idx,k,mu,nu,i
    real(wp),allocatable :: Bmattest(:,:)


    if (add_cmb) then

      ! nchain_cmb=1
      ! nseg_cmb=4
      ! nbead_cmb=5
      ! nseg_cmbbb=2
      ! nseg_cmbar=2
      ! Na=1
      ! Ia=[1,2]
      ! ntotseg_cmb=4
      ! ntotbead_cmb=5

      ! nseg_cmbbb=nseg_cmb-Na*nseg_cmbar
      nbead_cmb=nseg_cmb+1
      ntotseg_cmb=nchain_cmb*nseg_cmb
      ntotbead_cmb=ntotseg_cmb+1

      !-------------------------------------
      ! Constructing Bbar_tilde_cmb and B_tilde_cmb based on Bbar and B:
      ! For making Bbar_tilde_cmb sparse (CSR):
      !-------------------------------------
      maxnz_Bbar=nchain*nseg*3*2 + nchain_cmb*nseg_cmb*3*2
      allocate(Bbar_vals(maxnz_Bbar))
      allocate(Bbar_cols(maxnz_Bbar))
      allocate(Bbar_rowInd((nchain*nseg*3+nchain_cmb*nseg_cmb*3)+1))

      Bbar_rowInd(1)=1
      do ichain=1, nchain
        offsetch1=(ichain-1)*nseg*3
        offsetch2=(ichain-1)*nbead*3
        offsetchx2=offsetch1*2
        do isegx3=1, nseg*3
          offsetseg=(isegx3-1)*2
          Bbar_cols(offsetchx2+offsetseg+1)=offsetch2+isegx3
          Bbar_vals(offsetchx2+offsetseg+1)=-1._wp
          Bbar_cols(offsetchx2+offsetseg+2)=offsetch2+isegx3+3
          Bbar_vals(offsetchx2+offsetseg+2)=1._wp

          Bbar_rowInd(offsetch1+isegx3+1)=Bbar_rowInd(offsetch1+isegx3)+2
        end do
      end do

      do ichain=1, nchain_cmb
        offsetch1=nchain*nseg*3 + (ichain-1)*nseg_cmb*3
        offsetch2=nchain*nbead*3 + (ichain-1)*nbead_cmb*3
        offsetchx2=offsetch1*2
        
        iarm=1
        icoor=0
        do isegx3=1, nseg_cmb*3

          offsetseg=(isegx3-1)*2
          
          if (isegx3 <= nseg_cmbbb*3) then            
            Bbar_cols(offsetchx2+offsetseg+1)=offsetch2+isegx3
            Bbar_vals(offsetchx2+offsetseg+1)=-1._wp
            Bbar_cols(offsetchx2+offsetseg+2)=offsetch2+isegx3+3
            Bbar_vals(offsetchx2+offsetseg+2)=1._wp
          else ! iseg > nseg_cmbbb

            if ( (isegx3+2)/3 - nseg_cmbbb - (iarm-1)*nseg_cmbar == 1) then

              Bbar_cols(offsetchx2+offsetseg+1)=offsetch2+(Ia(iarm+1)-1)*3+icoor+1
              Bbar_vals(offsetchx2+offsetseg+1)=-1._wp
              Bbar_cols(offsetchx2+offsetseg+2)=offsetch2+isegx3+3
              Bbar_vals(offsetchx2+offsetseg+2)=1._wp

              ! taking into account all 3 components associated with the first segment of the arm
              icoor=icoor+1
              if (icoor==3) then
                iarm=iarm+1
                icoor=0
              endif

            else

              Bbar_cols(offsetchx2+offsetseg+1)=offsetch2+isegx3
              Bbar_vals(offsetchx2+offsetseg+1)=-1._wp
              Bbar_cols(offsetchx2+offsetseg+2)=offsetch2+isegx3+3
              Bbar_vals(offsetchx2+offsetseg+2)=1._wp

            end if
          endif

          Bbar_rowInd(offsetch1+isegx3+1)=Bbar_rowInd(offsetch1+isegx3)+2

        end do

        
      end do

      ! call print_vector(Bbar_vals,'vals')
      ! call print_vector(Bbar_cols,'cols')
      ! call print_vector(Bbar_rowInd,'rowind')


      ! allocate(Bmattest(ntotbead_cmb*3,ntotseg_cmb*3))

      !-------------------------------------
      ! >>> For making B_tilde sparse (CSR):
      !-------------------------------------

      ! maxnz_B_cmb=nchain*nbeadx3*nseg+nchain_cmb*nbead_cmb*3*nseg_cmb
      maxnz_B=nchain*nbeadx3*nseg+nchain_cmb*nbead_cmb*3*nseg_cmb
      allocate(B_vals(maxnz_B),B_cols(maxnz_B))
      allocate(B_rowInd(nchain*nbead*3+nchain_cmb*nbead_cmb*3+1))

      ! for debugging
      ! maxnz_B_cmb=nchain_cmb*nbead_cmb*3*nseg_cmb
      ! allocate(B_vals_cmb(maxnz_B_cmb),B_cols_cmb(maxnz_B_cmb),B_rowInd_cmb(ntotbead_cmb*3+1))

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


      do ichain=1, nchain_cmb

        offsetch1=nchain*nbeadx3*nseg + (ichain-1)*nbead_cmb*3*nseg_cmb
        offsetch2=nchain*nsegx3 + (ichain-1)*nseg_cmb*3
        offsetch3=nchain*nbeadx3 + (ichain-1)*nbead_cmb*3

        ! for debugging
        ! offsetch1= (ichain-1)*nbead_cmb*3*nseg_cmb
        ! offsetch2=(ichain-1)*nseg_cmb*3
        ! offsetch3=(ichain-1)*nbead_cmb*3

        ! Constructing the elements of the first row (first bead) of B
        do icoor=1, 3
          offsetcoor=(icoor-1)*nseg_cmb
          do jseg=1, nseg_cmbbb
            offsetseg1=offsetcoor+(jseg-1)
            offsetseg2=(jseg-1)*3
            B_vals(offsetch1+offsetseg1+1)=-(nseg_cmbbb-jseg+1)/real(nbead_cmb,kind=wp)
            B_cols(offsetch1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
          end do ! jseg
          B_rowInd(offsetch3+icoor+1)=B_rowInd(offsetch3+icoor)+nseg_cmb
          
          do iarm=1, Na

            ! modifying the values of the first row, due to the segments between the arms  
            fctr=(Na-iarm+1)*nseg_cmbar/real(nbead_cmb,kind=wp)
            do jseg=Ia(iarm), Ia(iarm+1)-1
              offsetseg1=offsetcoor+(jseg-1)
              B_vals(offsetch1+offsetseg1+1)=B_vals(offsetch1+offsetseg1+1)-fctr

            enddo

            ! the elements of the first row corresponding to the arm segments
            do kseg=1, nseg_cmbar
              jseg=nseg_cmbbb+(iarm-1)*nseg_cmbar+kseg
              offsetseg1=offsetcoor+(jseg-1)
              offsetseg2=(jseg-1)*3
              B_vals(offsetch1+offsetseg1+1)=-(nseg_cmbar-kseg+1)/real(nbead_cmb,kind=wp)
              B_cols(offsetch1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
            end do

          enddo

        enddo ! icoor

        ! Constructing the rest of the rows in backbone
        do ibead=2, nseg_cmbbb+1

          offsetbead1=(ibead-1)*nseg_cmb*3
          offsetbead2=(ibead-1)*3

          do icoor=1, 3
            offsetcoor=(icoor-1)*nseg_cmb

            do jseg=1, nseg_cmb
              offsetseg1=offsetcoor+(jseg-1)
              offsetseg2=(jseg-1)*3
              B_vals(offsetch1+offsetbead1+offsetseg1+1)=B_vals(offsetch1+offsetseg1+1)
              B_cols(offsetch1+offsetbead1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
            end do ! jseg

            ! modifying the values 
            do jseg=1, ibead-1
              offsetseg1=offsetcoor+(jseg-1)
              offsetseg2=(jseg-1)*3
              B_vals(offsetch1+offsetbead1+offsetseg1+1)=B_vals(offsetch1+offsetbead1+offsetseg1+1)+1
            end do ! jseg

            B_rowInd(offsetch3+offsetbead2+icoor+1)=B_rowInd(offsetch3+offsetbead2+icoor)+nseg_cmb
          
          end do ! icoor
        
        end do ! ibead

        ! Constructing the rows for the arms
        do iarm=1, Na

          ! Constructing the rest of the rows in backbone
          do kbead=1, nseg_cmbar

            ibead=nseg_cmbbb+1+(iarm-1)*nseg_cmbar+kbead
            offsetbead1=(ibead-1)*nseg_cmb*3
            offsetbead2=(ibead-1)*3

            do icoor=1, 3
              offsetcoor=(icoor-1)*nseg_cmb

              do jseg=1, nseg_cmb
                offsetseg1=offsetcoor+(jseg-1)
                offsetseg2=(jseg-1)*3
                B_vals(offsetch1+offsetbead1+offsetseg1+1)=B_vals(offsetch1+offsetseg1+1)
                B_cols(offsetch1+offsetbead1+offsetseg1+1)=offsetch2+offsetseg2+1+(icoor-1)
              end do ! jseg

              ! modifying the values 
              do jseg=1, Ia(iarm+1)-1
                offsetseg1=offsetcoor+(jseg-1)
                offsetseg2=(jseg-1)*3
                B_vals(offsetch1+offsetbead1+offsetseg1+1)=B_vals(offsetch1+offsetbead1+offsetseg1+1)+1
              end do ! jseg

              ! modifying the values 
              do kseg=1, kbead
                jseg=nseg_cmbbb+(iarm-1)*nseg_cmbar+kseg
                offsetseg1=offsetcoor+(jseg-1)
                offsetseg2=(jseg-1)*3
                B_vals(offsetch1+offsetbead1+offsetseg1+1)=B_vals(offsetch1+offsetbead1+offsetseg1+1)+1
              end do ! jseg

              B_rowInd(offsetch3+offsetbead2+icoor+1)=B_rowInd(offsetch3+offsetbead2+icoor)+nseg_cmb

            end do ! icoor

          end do ! ibead

        enddo

      end do ! ichain

      ! call print_vector(B_vals(nchain*nbeadx3*nseg+1:maxnz_B),'vals')
      ! call print_vector(B_cols(nchain*nbeadx3*nseg+1:maxnz_B),'cols')
      ! call print_vector(B_rowInd(nchain*nbead*3:nchain*nbead*3+nchain_cmb*nbead_cmb*3+1),'rowind')

      ! Bmattest=0.0_wp
      ! ! Constructing the elements of the first row of B
      ! do k=1, nseg_cmbbb
      !  forall (i=1:3) Bmattest(i,3*(k-1)+i)=-(nseg_cmbbb-k+1)/real(nbead_cmb,kind=wp)
      ! end do
      ! do iarm=1, Na
      !  fctr=(Na-iarm+1)*nseg_cmbar/real(nbead_cmb,kind=wp)
      !  do k=Ia(iarm), Ia(iarm+1)-1
      !    forall (i=1:3)
      !      Bmattest(i,3*(k-1)+i)=Bmattest(i,3*(k-1)+i)-fctr
      !    end forall
      !  end do ! k
      !  do k=1, nseg_cmbar
      !    idx=nseg_cmbbb+(iarm-1)*nseg_cmbar+k
      !    forall (i=1:3)
      !      Bmattest(i,3*(idx-1)+i)=Bmattest(i,3*(idx-1)+i)-&
      !      (nseg_cmbar-k+1)/real(nbead_cmb,kind=wp)
      !    end forall
      !  end do ! k
      ! end do ! iarm
      ! ! Constructing the rest of the rows in backbone
      ! do nu=2, nseg_cmbbb+1
      !  forall (i=1:3) Bmattest(3*(nu-1)+i,:)=Bmattest(i,:)
      !  do k=1, nu-1
      !    forall (i=1:3)
      !      Bmattest(3*(nu-1)+i,3*(k-1)+i)=Bmattest(3*(nu-1)+i,3*(k-1)+i)+1
      !    end forall
      !  end do ! k
      ! end do ! nu
      ! ! Constructing the rows for the arms
      ! do iarm=1, Na
      !  do mu=1, nseg_cmbar
      !    nu=nseg_cmbbb+1+(iarm-1)*nseg_cmbar+mu
      !    forall (i=1:3) Bmattest(3*(nu-1)+i,:)=Bmattest(i,:)
      !    do k=1, Ia(iarm+1)-1
      !      forall (i=1:3)
      !        Bmattest(3*(nu-1)+i,3*(k-1)+i)=Bmattest(3*(nu-1)+i,3*(k-1)+i)+1
      !      end forall
      !    end do ! k
      !    do k=1, mu
      !      idx=nseg_cmbbb+(iarm-1)*nseg_cmbar+k
      !      forall (i=1:3)
      !        Bmattest(3*(nu-1)+i,3*(idx-1)+i)=Bmattest(3*(nu-1)+i,3*(idx-1)+i)+1
      !      end forall
      !    end do ! k
      !  end do ! mu
      ! end do ! iarm

      ! call print_matrix(Bmattest,'btest')

    else

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

      ! call print_vector(Bbar_vals,'vals')
      ! call print_vector(Bbar_cols,'cols')
      ! call print_vector(Bbar_rowInd,'rowind')

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

    endif ! add_cmb


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

  !> Converting Rb to Q using sparse multiplication
  !! \param Q connectivity vectors
  !! \param R bead to center of mass vectors
  subroutine RbtoQ(Rb,Q,ntotsegx3,ntotbeadx3,bs)

    integer,intent(in) :: ntotsegx3,ntotbeadx3
    real(wp),intent(in) :: Rb(:)
    real(wp),intent(inout) :: Q(:)
    real(wp),intent(in) :: bs(3)

    real(wp) :: qx,qy,qz,bsx,bsy,bsz,invbsx,invbsy,invbsz
    integer :: its

    bsx=bs(1);bsy=bs(2);bsz=bs(3)
    invbsx=1/bs(1);invbsy=1/bs(2);invbsz=1/bs(3)

#ifdef USE_DP
    call mkl_dcsrmv('N',ntotsegx3,ntotbeadx3,1._wp,'GIIF',Bbar_vals,Bbar_cols,&
      Bbar_rowInd,Bbar_rowInd(2),Rb,0._wp,Q)
#elif USE_SP
    call mkl_scsrmv('N',ntotsegx3,ntotbeadx3,1._wp,'GIIF',Bbar_vals,Bbar_cols,&
      Bbar_rowInd,Bbar_rowInd(2),Rb,0._wp,Q)
#endif

    do its=1, ntotsegx3/3

      qx=Q((its-1)*3+1)
      qy=Q((its-1)*3+2)
      qz=Q((its-1)*3+3)

      qx=qx-nint(qx*invbsx)*bsx
      qy=qy-nint(qy*invbsy)*bsy
      qz=qz-nint(qz*invbsz)*bsz
      ! select case (FlowType)
      !   case ('PSF')
      !     qx=qx+eps_m*qy
      !   case ('PEF')
      !     qytmp=qy
      !     qx=qx+tanb*qytmp
      !     qy=sinth*qx+costh*qytmp
      !     qx=costh*qx-sinth*qytmp
      ! end select

      Q((its-1)*3+1)=qx
      Q((its-1)*3+2)=qy
      Q((its-1)*3+3)=qz

    enddo

  end subroutine RbtoQ

  !> Converting {R,rc} to (Rbx,Rby,Rbz)
  !! \param R bead to center of mass distance for all chains
  !! \param rcm center of mass for all chains
  !! \param Rbx x-coordinate of the position vector
  !! \param Rby y-coordinate of the position vector
  !! \param Rbz z-coordinate of the position vector
  subroutine RtoRbc(R,rcm,Rbx,Rby,Rbz,nchain,nbead,ntotbead,ntotbeadx3,add_cmb,nchain_cmb,nseg_cmb)

!$  use :: omp_lib

    integer,intent(in) :: nchain,nbead,ntotbead,ntotbeadx3
    real(wp),intent(in) :: R(:),rcm(:,:)
    real(wp),intent(inout) :: Rbx(ntotbead),Rby(ntotbead),Rbz(ntotbead)
    logical,intent(in) :: add_cmb
    integer,intent(in) :: nchain_cmb,nseg_cmb
    integer :: igb,os,ich,os1,os2

!$omp parallel default(private) shared(nchain,ntotbead,nbead,Rbx,Rby,Rbz,R,rcm)
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

    if (add_cmb) then
!$omp parallel default(private) shared(nchain,ntotbead,nbead,Rbx,Rby,Rbz,R,rcm,add_cmb,nchain_cmb,nseg_cmb)
!$omp do schedule(auto)
      do igb=1, nchain_cmb*(nseg_cmb+1)
        os1=nchain*nbead+(igb-1)
        os2=nchain*nbead*3+(igb-1)*3
        ich=nchain+(igb-1)/(nseg_cmb+1)+1

        Rbx(os1+1)=R(os2+1)+rcm(ich,1)
        Rby(os1+1)=R(os2+2)+rcm(ich,2)
        Rbz(os1+1)=R(os2+3)+rcm(ich,3)

      end do
!$omp end do
!$omp end parallel
    endif
    
  end subroutine RtoRbc

  !> Converting Rbc to Q (only for linear chains)
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
  subroutine RtoRb(R,rcm,Rb,nchain,nbead,ntotbead,ntotbeadx3,add_cmb,nchain_cmb,nseg_cmb)

!$  use :: omp_lib

    integer,intent(in) :: nchain,nbead,ntotbead,ntotbeadx3
    real(wp),intent(in) :: R(:)
    real(wp),intent(in) :: rcm(:,:)
    real(wp),intent(inout) :: Rb(:)
    logical,intent(in) :: add_cmb
    integer,intent(in) :: nchain_cmb,nseg_cmb
    integer :: igb,os,ich

!$omp parallel default(private) shared(nchain,ntotbead,nbead,Rb,R,rcm)
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

    if (add_cmb) then
!$omp parallel default(private) shared(nchain,ntotbead,nbead,Rb,R,rcm,add_cmb,nchain_cmb,nseg_cmb)
!$omp do schedule(auto)
      do igb=1, nchain_cmb*(nseg_cmb+1)
        os=nchain*nbead*3+(igb-1)*3
        ich=nchain+(igb-1)/(nseg_cmb+1)+1
        Rb(os+1)=R(os+1)+rcm(ich,1)
        Rb(os+2)=R(os+2)+rcm(ich,2)
        Rb(os+3)=R(os+3)+rcm(ich,3)
      end do
!$omp end do
!$omp end parallel
    endif

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
