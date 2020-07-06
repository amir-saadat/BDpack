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
module chain_mod

  use :: prcn_mod
!$ use :: omp_lib

  implicit none

  type :: chain

    integer :: chain_ID
    !> The pointer to the chain section of Rb_tilde
    real(wp),pointer :: chain_Rb(:) => null()
    !> The pointer to the chain section of Rb x-component
    real(wp),pointer :: chain_Rbx(:) => null()
    !> The pointer to the chain section of Rb y-component
    real(wp),pointer :: chain_Rby(:) => null()
    !> The pointer to the chain section of Rb z-component
    real(wp),pointer :: chain_Rbz(:) => null()
    !> The pointer to the chain section of the b_img
    integer,pointer :: chain_b_img(:,:) => null()
    !> The pointer to chain section of Qdagger_tilde
    real(wp),pointer :: chain_Q(:) => null()
    !> The pointer to chain section of R_tilde
    real(wp),pointer :: chain_R(:) => null()
    !> The pointer to chain section of rcm_tilde
    real(wp),pointer :: chain_rcm(:) => null()
    !> The pointer to the chain section of the cm_img
    integer,pointer :: chain_cm_img(:) => null()
    !> The pointer to chain section of image flag
    integer,pointer :: chain_cmif(:) => null()
    !> The chain x-component of transformed position of the beads
    real(wp),pointer :: chain_Rbtrx(:) => null()
    !> The chain y-component of transformed position of the beads
    real(wp),pointer :: chain_Rbtry(:) => null()
    !> The chain x-component of transformed position of the center of mass
    real(wp),pointer :: chain_rcmtrx => null()
    !> The chain y-component of transformed position of the center of mass
    real(wp),pointer :: chain_rcmtry => null()
    !> The image flag of the center of mass of the chain
!    integer :: rcm_ImFlag(3)
  
    contains
      procedure,pass(this) :: init => init_chain
      procedure,pass(this) :: update => update_chain
      final :: del_chain

  end type chain

contains

  subroutine init_chain(this,id,nchain,nsegx3,nbead,nbeadx3,Rb,Rbx,Rby,Rbz,Q,R,rcm,cmif,Rbtr,rcmtr,&
    b_img,cm_img,nseg_cmb)

    use :: flow_mod, only: FlowType
    use :: io_mod, only: CoMDiff

    class(chain),intent(inout) :: this
    integer,intent(in) :: id,nchain,nsegx3,nbead,nbeadx3
    real(wp),intent(in),target :: Rb(:)
    real(wp),intent(in),target :: Rbx(:),Rby(:),Rbz(:)
    real(wp),intent(in),target :: Q(:)
    real(wp),intent(in),target :: R(:)
    real(wp),intent(in),target :: rcm(:,:)
    integer,intent(in),target :: cmif(:,:)
    real(wp),intent(in),target :: Rbtr(:,:)
    real(wp),intent(in),target :: rcmtr(:,:)
    integer,intent(in),target :: b_img(:,:)
    integer,intent(in),target :: cm_img(:,:)
    integer,intent(in),optional :: nseg_cmb

    integer :: offsetseg,offsetbead1,offsetbead2,nbead_cmb

    ! this%chain_ID=id
    ! this%chain_Q => Q((id-1)*nsegx3+1:(id-1)*nsegx3+nsegx3)
    ! this%chain_R => R((id-1)*nbeadx3+1:(id-1)*nbeadx3+nbeadx3)
    ! this%chain_rcm => rcm(this%chain_ID,:)
    ! if ((FlowType == 'Equil').and.CoMDiff) then
    !   this%chain_cmif => cmif(this%chain_ID,:)
    ! end if
    ! this%chain_Rb => Rb((id-1)*nbeadx3+1:(id-1)*nbeadx3+nbeadx3)
    ! this%chain_Rbx => Rbx((id-1)*nbead+1:(id-1)*nbead+nbead)
    ! this%chain_Rby => Rby((id-1)*nbead+1:(id-1)*nbead+nbead)
    ! this%chain_Rbz => Rbz((id-1)*nbead+1:(id-1)*nbead+nbead)
    ! if (FlowType == 'PEF') then
    !   this%chain_Rbtrx => Rbtr((id-1)*nbead+1:(id-1)*nbead+nbead,1)
    !   this%chain_Rbtry => Rbtr((id-1)*nbead+1:(id-1)*nbead+nbead,2)
    !   this%chain_rcmtrx => rcmtr(id,1)
    !   this%chain_rcmtry => rcmtr(id,2)
    ! end if
    ! this%chain_b_img => b_img((id-1)*nbead+1:(id-1)*nbead+nbead,:)
    ! this%chain_cm_img => cm_img(id,:)

    if (present(nseg_cmb)) then ! comb chain

      nbead_cmb=nseg_cmb+1
      offsetseg=nchain*nsegx3
      offsetbead1=nchain*nbead
      offsetbead2=nchain*nbeadx3

      this%chain_ID=nchain+id
      this%chain_Q => Q(offsetseg+(id-1)*nseg_cmb*3+1:offsetseg+(id-1)*nseg_cmb*3+nseg_cmb*3)
      this%chain_R => R(offsetbead2+(id-1)*nbead_cmb*3+1:offsetbead2+(id-1)*nbead_cmb*3+nbead_cmb*3)
      this%chain_rcm => rcm(this%chain_ID,:)
      if ((FlowType == 'Equil').and.CoMDiff) then
        this%chain_cmif => cmif(this%chain_ID,:)
      end if
      this%chain_Rb => Rb(offsetbead2+(id-1)*nbead_cmb*3+1:offsetbead2+(id-1)*nbead_cmb*3+nbead_cmb*3)
      this%chain_Rbx => Rbx(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb)
      this%chain_Rby => Rby(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb)
      this%chain_Rbz => Rbz(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb)
      if (FlowType == 'PEF') then
        this%chain_Rbtrx => Rbtr(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb,1)
        this%chain_Rbtry => Rbtr(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb,2)
        this%chain_rcmtrx => rcmtr(id,1)
        this%chain_rcmtry => rcmtr(id,2)
      end if
      this%chain_b_img => b_img(offsetbead1+(id-1)*nbead_cmb+1:offsetbead1+(id-1)*nbead_cmb+nbead_cmb,:)
      this%chain_cm_img => cm_img(id,:)

    else ! linear chain
            
      this%chain_ID=id
      this%chain_Q => Q((id-1)*nsegx3+1:(id-1)*nsegx3+nsegx3)
      this%chain_R => R((id-1)*nbeadx3+1:(id-1)*nbeadx3+nbeadx3)
      this%chain_rcm => rcm(this%chain_ID,:)
      if ((FlowType == 'Equil').and.CoMDiff) then
        this%chain_cmif => cmif(this%chain_ID,:)
      end if
      this%chain_Rb => Rb((id-1)*nbeadx3+1:(id-1)*nbeadx3+nbeadx3)
      this%chain_Rbx => Rbx((id-1)*nbead+1:(id-1)*nbead+nbead)
      this%chain_Rby => Rby((id-1)*nbead+1:(id-1)*nbead+nbead)
      this%chain_Rbz => Rbz((id-1)*nbead+1:(id-1)*nbead+nbead)
      if (FlowType == 'PEF') then
        this%chain_Rbtrx => Rbtr((id-1)*nbead+1:(id-1)*nbead+nbead,1)
        this%chain_Rbtry => Rbtr((id-1)*nbead+1:(id-1)*nbead+nbead,2)
        this%chain_rcmtrx => rcmtr(id,1)
        this%chain_rcmtry => rcmtr(id,2)
      end if
      this%chain_b_img => b_img((id-1)*nbead+1:(id-1)*nbead+nbead,:)
      this%chain_cm_img => cm_img(id,:)

    endif


  end subroutine init_chain


  !> Updating the properties of the chain
  !! \param invbs the inverse of box dimensions
  subroutine update_chain(this,bs,invbs)

    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: L1,L2
    use :: io_mod, only: CoMDiff

    class(chain),intent(inout) :: this
    real(wp),intent(in) :: bs(3),invbs(3)

    integer :: nbead
    
    nbead=size(this%chain_Rbx,1)

    this%chain_b_img(:,1)=this%chain_b_img(:,1)-this%chain_cm_img(1)
    this%chain_b_img(:,2)=this%chain_b_img(:,2)-this%chain_cm_img(2)
    this%chain_b_img(:,3)=this%chain_b_img(:,3)-this%chain_cm_img(3)
    
    select case (FlowType)
       case ('Equil')
         this%chain_rcm(1)=sum(this%chain_Rbx(:)-this%chain_b_img(:,1)*bs(1))/nbead
         this%chain_rcm(2)=sum(this%chain_Rby(:)-this%chain_b_img(:,2)*bs(2))/nbead
         this%chain_rcm(3)=sum(this%chain_Rbz(:)-this%chain_b_img(:,3)*bs(3))/nbead
         if (CoMDiff) this%chain_cmif(:)=this%chain_cmif(:)-this%chain_cm_img(:)
!         this%rcm_ImFlag=this%rcm_ImFlag-this%chain_cm_img
       case ('PSF')
         this%chain_rcm(1)=sum(this%chain_Rbx(:)-this%chain_b_img(:,1)*bs(1))/nbead
         this%chain_rcm(2)=sum(this%chain_Rby(:)-this%chain_b_img(:,2)*bs(2))/nbead
         this%chain_rcm(3)=sum(this%chain_Rbz(:)-this%chain_b_img(:,3)*bs(3))/nbead
       case ('PEF')
         this%chain_rcm(1)=sum(this%chain_Rbx(:)-this%chain_b_img(:,1)*L1(1)-this%chain_b_img(:,2)*L2(1))/nbead
         this%chain_rcm(2)=sum(this%chain_Rby(:)-this%chain_b_img(:,1)*L1(2)-this%chain_b_img(:,2)*L2(2))/nbead
         this%chain_rcm(3)=sum(this%chain_Rbz(:)-this%chain_b_img(:,3)*bs(3))/nbead
    end select

  end subroutine update_chain

  !> Destructor of the chain type
  subroutine del_chain(this)

    type(chain),intent(inout) :: this

    nullify(this%chain_Rbx,this%chain_Rby,this%chain_Rbz)
    nullify(this%chain_Rb)
    nullify(this%chain_Q,this%chain_R,this%chain_rcm)

  end subroutine del_chain

end module chain_mod
