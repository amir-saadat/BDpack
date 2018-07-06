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

  subroutine advance_sph_sde(this,r_sph,p_sph,q_sph,rf0,iPe,idt,ichain,Fseg,wbl_sph,wbl_sph_or,jcol)

    !variables from other places
    use :: inp_dlt, only: hstar,nchain_pp,dt,nseg_indx3
    use :: arry_mod, only: print_vector, print_matrix

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    !input and output arguments
    class(sph_sde_t),intent(inout) :: this
    real(wp), intent(inout) :: r_sph(:,:),p_sph(:,:),q_sph(:,:),rf0(:,:,:)
    integer, intent(in) :: ichain,jcol,iPe,idt
    real(wp), intent(in) :: Fseg(:),wbl_sph(:,:),wbl_sph_or(:,:)

    !variables used inside advance_sph_sde
    integer :: ichain_pp, offset
    real(wp),dimension(3) :: dr_sph,dphi_sph,ax_sph
    real(wp),dimension(3) :: p_old,q_old, rf0_rel_old,rf0_rel_new
    real(wp) :: ang_sph
    real(wp),dimension(3,3) :: R

    dr_sph(1:3) = 0._wp
    dphi_sph(1:3) = 0._wp
    R(1:3,1:3) = 0._wp
    p_old(1:3)  = p_sph(:,ichain)
    q_old(1:3)  = q_sph(:,ichain)
    rf0_rel_old(1:3) = 0._wp
    rf0_rel_new(1:3) = 0._wp
    !advance the orientation of the sphere: Brownian motion
    dphi_sph(1) = dphi_sph(1) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(1,jcol)
    dphi_sph(2) = dphi_sph(2) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(2,jcol)
    dphi_sph(3) = dphi_sph(3) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(3,jcol)

    ang_sph = sqrt(dphi_sph(1)**2 + dphi_sph(2)**2 + dphi_sph(3)**2)
    ax_sph = dphi_sph(:)/ang_sph

    ! R(1,1:3) = (/1._wp,1._wp,1._wp/)
    ! R(2,1:3) = (/1._wp,1._wp,1._wp/)
    ! R(3,1:3) = (/1._wp,1._wp,1._wp/)
    ! call print_matrix(R(:,:),'R before1')
    ! R(1,1:3) = R(1,1:3) + 2._wp * (/1._wp,1._wp,1._wp/)
    ! R(2,1:3) = R(2,1:3) + (/1._wp,1._wp,1._wp/)
    ! R(3,1:3) = R(3,1:3) + (/1._wp,1._wp,1._wp/)
    ! call print_matrix(R(:,:),'R before2')

    R(1,1:3) = (/COS(ang_sph),0._wp,0._wp/)
    R(2,1:3) = (/0._wp,COS(ang_sph),0._wp/)
    R(3,1:3) = (/0._wp,0._wp,COS(ang_sph)/)

    R(1,1:3) = R(1,1:3) + SIN(ang_sph)*(/0._wp,-ax_sph(3),ax_sph(2)/)
    R(2,1:3) = R(2,1:3) + SIN(ang_sph)*(/ax_sph(3),0._wp,-ax_sph(1)/)
    R(3,1:3) = R(3,1:3) + SIN(ang_sph)*(/-ax_sph(2),ax_sph(1),0._wp/)

    R(1,1:3) = R(1,1:3) + (1-COS(ang_sph))*(/ax_sph(1)*ax_sph(1),ax_sph(1)*ax_sph(2),ax_sph(1)*ax_sph(3)/)
    R(2,1:3) = R(2,1:3) + (1-COS(ang_sph))*(/ax_sph(2)*ax_sph(1),ax_sph(2)*ax_sph(2),ax_sph(2)*ax_sph(3)/)
    R(3,1:3) = R(3,1:3) + (1-COS(ang_sph))*(/ax_sph(3)*ax_sph(1),ax_sph(3)*ax_sph(2),ax_sph(3)*ax_sph(3)/)

    !call print_matrix(R(:,:),'R = ')
    !call print_vector(p_sph(:,ichain),'p (before) = ')
    !call print_vector(q_sph(:,ichain),'q (before) = ')
    call gemv(R,p_old(:),p_sph(:,ichain),alpha=1._wp,beta=0._wp)
    call gemv(R,q_old(:),q_sph(:,ichain),alpha=1._wp,beta=0._wp)
    !call print_vector(p_sph(:,ichain),'p (after) = ')
    !call print_vector(q_sph(:,ichain),'q (after) = ')
    do ichain_pp = 1,nchain_pp
      !call print_vector(rf0(:,ichain,ichain_pp) ,'rf0(:,ichain,ichain_pp) (before)  = ')
      rf0_rel_old(1:3) = rf0(:,ichain,ichain_pp) - r_sph(:,ichain)
      !call print_vector(rf0(:,ichain,ichain_pp) ,'rf0(:,ichain,ichain_pp) (before)  = ')
      call gemv(R,rf0_rel_old(1:3),rf0_rel_new(1:3),alpha=1._wp,beta=0._wp)
      rf0(1:3,ichain,ichain_pp) = rf0_rel_new(1:3) + r_sph(1:3,ichain)
      !call print_vector(rf0(:,ichain,ichain_pp) ,'rf0(:,ichain,ichain_pp) (after)  = ')
    end do

    !advance the coordinate of the sphere: tether forces
    do ichain_pp = 1,nchain_pp
      !using Fseg
      offset = nseg_indx3*(ichain_pp-1)
      dr_sph(1) = dr_sph(1) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+1)*dt(iPe,idt)
      dr_sph(2) = dr_sph(2) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+2)*dt(iPe,idt)
      dr_sph(3) = dr_sph(3) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+3)*dt(iPe,idt)

      !using Ftet
      ! offset = 3*(ichain_pp-1)
      ! dr_sph(1) = dr_sph(1) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+1)*dt(iPe,idt)
      ! dr_sph(2) = dr_sph(2) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+2)*dt(iPe,idt)
      ! dr_sph(3) = dr_sph(3) - (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Ftet(offset+3)*dt(iPe,idt)
    end do
    !call print_vector(dr_sph(:),'dr_sph(:) (only tether forces)')
    !advance the coordinate of the sphere: Brownian motion
    dr_sph(1) = dr_sph(1) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(1,jcol)
    dr_sph(2) = dr_sph(2) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(2,jcol)
    dr_sph(3) = dr_sph(3) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(3,jcol)

    !call print_vector(dr_sph(:),'dr_sph(:) (tether forces and brownian motion)')
    r_sph(:,ichain) = r_sph(:,ichain) + dr_sph(1:3)

    !update the tether points on the surface of the sphere
    do ichain_pp = 1,nchain_pp
      rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph(1:3)
      !call print_vector(rf0(:,ichain,ichain_pp) ,'rf0(:,ichain,ichain_pp) (after translation)  = ')
    end do

  end subroutine advance_sph_sde

end module sph_sde_mod
