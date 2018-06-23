!sphere sde module
!Tiras Lin
!May 2, 2017
!!!!!!!!!!!!!!!!!!!!!!!!

module sph_sde_mod

  use :: prcn_mod

  implicit none

  public :: sph_sde_t
  private

  !type containing the SDE
  type :: sph_sde_t
    real(wp) :: a_sph
  contains
    procedure,pass(this) :: init => init_sph_sde
    procedure,pass(this) :: advance => advance_sph_sde
  end type sph_sde_t

contains

  subroutine init_sph_sde(this)

    use :: cmn_io_mod, only: read_input

    !input and output arguments
    class(sph_sde_t),intent(inout) :: this

    call read_input('Sph-rad',0,this%a_sph)

  end subroutine init_sph_sde

  subroutine advance_sph_sde(this,r_sph,rf0,iPe,idt,ichain,Ftet,wbl_sph,jcol)

    !variables from other places
    use :: inp_dlt, only: hstar,nchain_pp,dt,nseg_indx3
    use :: arry_mod, only: print_vector, print_matrix

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    !input and output arguments
    class(sph_sde_t),intent(inout) :: this
    real(wp), intent(inout) :: r_sph(:,:),rf0(:,:,:)
    integer, intent(in) :: ichain,jcol,iPe,idt
    real(wp), intent(in) :: Ftet(:),wbl_sph(:,:)

    !variables used inside advance_sph_sde
    integer :: ichain_pp, offset
    real(wp),dimension(3) :: dr_sph

    dr_sph(1:3) = 0._wp

    !advance the coordinate of the sphere: tether forces
    do ichain_pp = 1,nchain_pp
      !using Fseg
      ! offset = nseg_indx3*(ichain_pp-1)
      ! dr_sph(1) = dr_sph(1) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+1)*dt(iPe,idt)
      ! dr_sph(2) = dr_sph(2) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+2)*dt(iPe,idt)
      ! dr_sph(3) = dr_sph(3) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+3)*dt(iPe,idt)

      !using Ftet
      offset = 3*(ichain_pp-1)
      dr_sph(1) = dr_sph(1) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+1)*dt(iPe,idt)
      dr_sph(2) = dr_sph(2) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+2)*dt(iPe,idt)
      dr_sph(3) = dr_sph(3) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+3)*dt(iPe,idt)
    end do
    call print_vector(dr_sph(:),'dr_sph(:)')
    !advance the coordinate of the sphere: Brownian motion
    dr_sph(1) = dr_sph(1) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(1,jcol)
    dr_sph(2) = dr_sph(2) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(2,jcol)
    dr_sph(3) = dr_sph(3) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(3,jcol)

    r_sph(:,ichain) = r_sph(:,ichain) + dr_sph(1:3)

    !update the tether points on the surface of the sphere
    do ichain_pp = 1,nchain_pp
      rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph(1:3)
    end do

  end subroutine advance_sph_sde

end module sph_sde_mod
