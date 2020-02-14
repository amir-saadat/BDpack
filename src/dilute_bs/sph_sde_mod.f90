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

  subroutine advance_sph_sde(this,r_sph,p_sph,q_sph,rf0,iPe,idt,ichain,Fseg,Fbead,wbl_sph,wbl_sph_or,jcol,Fev_sph,F_sph,Fev,Fbnd,DiffTensP,Amat,Fphi,Fphi_all,Fphi_all_temp,wbltempP1,coeff)

    !variables from other places
    use :: inp_dlt, only: hstar,nchain_pp,dt,nbead_indx3,nseg_indx3,tplgy,nbeadx3
    use :: arry_mod, only: print_vector, print_matrix

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    !input and output arguments
    class(sph_sde_t),intent(inout) :: this
    real(wp), intent(inout) :: r_sph(:,:),p_sph(:,:),q_sph(:,:),rf0(:,:,:)
    integer, intent(in) :: ichain,jcol,iPe,idt
    real(wp), intent(in) :: Fseg(:),wbl_sph(:,:),wbl_sph_or(:,:),Fev_sph(:)
    real(wp), intent(out) :: F_sph(:)
    real(wp), intent(inout) :: Fphi(:,:),Fev(:),Fbnd(:),Fbead(:),Fphi_all(:),Fphi_all_temp(:)
    real(wp), intent(in) :: DiffTensP(:,:),wbltempP1(:,:)
    real(wp), intent(in) :: Amat(:,:)

    !variables used inside advance_sph_sde
    integer :: ichain_pp, offset
    real(wp),dimension(3) :: dr_sph,dphi_sph,ax_sph
    real(wp) :: ang_sph
    real(wp),dimension(3) :: p_old,q_old, rf0_rel_old,rf0_rel_new
    real(wp),dimension(3) :: rf0_rel,rf0_rel_unit, Tseg
    real(wp),dimension(3,3) :: R
    logical :: debug_TYL
    real(wp) :: coeff


    debug_TYL = .false.

    !initializing the variables
    dr_sph(1:3) = 0._wp
    dphi_sph(1:3) = 0._wp
    R(1:3,1:3) = 0._wp
    p_old(1:3)  = p_sph(:,ichain)
    q_old(1:3)  = q_sph(:,ichain)
    rf0_rel_old(1:3) = 0._wp
    rf0_rel_new(1:3) = 0._wp
    !rf0_unit = 0._wp
    rf0_rel = 0._wp
    rf0_rel_unit = 0._wp
    Tseg = 0._wp

    !ROTATION------------------------------------------------------
    !advance the orientation of the sphere: Brownian motion
    dphi_sph(1) = dphi_sph(1) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(1,jcol)
    dphi_sph(2) = dphi_sph(2) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(2,jcol)
    dphi_sph(3) = dphi_sph(3) + sqrt(3._wp/8._wp)*sqrt(hstar*sqrtPI/(this%a_sph**3))*wbl_sph_or(3,jcol)

    !advance the orientation of the sphere: torques
    do ichain_pp = 1,nchain_pp
      offset = nseg_indx3*(ichain_pp-1)
      !rf0_unit(:) = rf0(:,ichain,ichain_pp) / (sqrt(rf0(1,ichain,ichain_pp)**2 + rf0(2,ichain,ichain_pp)**2 + rf0(3,ichain,ichain_pp)**2))
      !wait...this should be relative!

      rf0_rel(:) = rf0(:,ichain,ichain_pp) - r_sph(:,ichain)
      rf0_rel_unit(:) = rf0_rel(:) / (sqrt(rf0_rel(1)**2 + rf0_rel(2)**2 + rf0_rel(3)**2))

      Tseg(1) = rf0_rel_unit(2)*Fseg(offset+3) - rf0_rel_unit(3)*Fseg(offset+2)
      Tseg(2) = rf0_rel_unit(3)*Fseg(offset+1) - rf0_rel_unit(1)*Fseg(offset+3)
      Tseg(3) = rf0_rel_unit(1)*Fseg(offset+2) - rf0_rel_unit(2)*Fseg(offset+1)

      if (debug_TYL) then
        call print_vector(Tseg(1:3),'Sphere SDE, spring torque:')
      end if

      dphi_sph(1) = dphi_sph(1) + (3._wp/16._wp)*(hstar*sqrtPI/(this%a_sph**2))*Tseg(1)*dt(iPe,idt)
      dphi_sph(2) = dphi_sph(2) + (3._wp/16._wp)*(hstar*sqrtPI/(this%a_sph**2))*Tseg(2)*dt(iPe,idt)
      dphi_sph(3) = dphi_sph(3) + (3._wp/16._wp)*(hstar*sqrtPI/(this%a_sph**2))*Tseg(3)*dt(iPe,idt)
    enddo

    !axis angle representation of the angular displacement
    ang_sph = sqrt(dphi_sph(1)**2 + dphi_sph(2)**2 + dphi_sph(3)**2)
    ax_sph = dphi_sph(:)/ang_sph

    !ang_sph = 0.0_wp
    !ax_sph = (/1._wp,0._wp,0._wp/)

    !constructing the rotation matrix (Rodriguez formula)
    R(1,1:3) = (/COS(ang_sph),0._wp,0._wp/)
    R(2,1:3) = (/0._wp,COS(ang_sph),0._wp/)
    R(3,1:3) = (/0._wp,0._wp,COS(ang_sph)/)

    R(1,1:3) = R(1,1:3) + SIN(ang_sph)*(/0._wp,-ax_sph(3),ax_sph(2)/)
    R(2,1:3) = R(2,1:3) + SIN(ang_sph)*(/ax_sph(3),0._wp,-ax_sph(1)/)
    R(3,1:3) = R(3,1:3) + SIN(ang_sph)*(/-ax_sph(2),ax_sph(1),0._wp/)

    R(1,1:3) = R(1,1:3) + (1-COS(ang_sph))*(/ax_sph(1)*ax_sph(1),ax_sph(1)*ax_sph(2),ax_sph(1)*ax_sph(3)/)
    R(2,1:3) = R(2,1:3) + (1-COS(ang_sph))*(/ax_sph(2)*ax_sph(1),ax_sph(2)*ax_sph(2),ax_sph(2)*ax_sph(3)/)
    R(3,1:3) = R(3,1:3) + (1-COS(ang_sph))*(/ax_sph(3)*ax_sph(1),ax_sph(3)*ax_sph(2),ax_sph(3)*ax_sph(3)/)

    !Rotating the orientation of the particle and the tether points
    call gemv(R,p_old(:),p_sph(:,ichain),alpha=1._wp,beta=0._wp)
    call gemv(R,q_old(:),q_sph(:,ichain),alpha=1._wp,beta=0._wp)
    do ichain_pp = 1,nchain_pp
      rf0_rel_old(1:3) = rf0(:,ichain,ichain_pp) - r_sph(:,ichain)
      call gemv(R,rf0_rel_old(1:3),rf0_rel_new(1:3),alpha=1._wp,beta=0._wp)
      rf0(1:3,ichain,ichain_pp) = rf0_rel_new(1:3) + r_sph(1:3,ichain)
    end do


    !TRANSLATION------------------------------------------------------
    !calculation of the total force on the sphere: sum of the tethered spring forces
    F_sph = 0._wp
    F_sph = Fev_sph

    do ichain_pp = 1,nchain_pp
      offset = nseg_indx3*(ichain_pp-1)
      F_sph(1) = F_sph(1) + Fseg(offset+1)
      F_sph(2) = F_sph(2) + Fseg(offset+2)
      F_sph(3) = F_sph(3) + Fseg(offset+3)
    end do

    !fake external force for debugging
    !F_sph(2) = F_sph(2) + 400._wp


    !calculation of the spring forces on the beads
    if (tplgy == 'Linear') then
      !call gbmv(this%AmatBF,Fseg,Fbead,kl=0,m=nsegx3,alpha=-1.0_wp,trans='T')
      call gemv(Amat,Fseg,Fbead,alpha=-1._wp,trans='T')
    else
      call gemv(Amat,Fseg,Fbead,alpha=-1._wp,trans='T')
    end if

    !Sum of the forces on the beads
    Fphi(:,ichain)=Fbead+Fev+Fbnd

    !All forces on beads, and the force on the sphere at the end of the vector
    Fphi_all(1:nbeadx3) = Fphi(:,ichain)
    Fphi_all(nbeadx3+1:nbeadx3+3) = F_sph(1:3)

    !A temporary vector also holding Fphi_all.
    !In calculating the displacement of the sphere due to beads,
    !the forces on the tethered beads should not be included.
    Fphi_all_temp = Fphi_all
    do ichain_pp = 1,nchain_pp
      offset = nbead_indx3*(ichain_pp-1)
      Fphi_all_temp(offset+1:offset+3) = 0._wp
    end do

    !call print_vector(F_sph(:),'Sph SDE, F_sph:')
    !call print_vector(Fseg(:),'Sph SDE, Fseg:')
    !call print_vector(Fbead(:),'Sph SDE, Fbead:')
    !call print_vector(Fev(:),'Sph SDE, Fev:')
    !call print_vector(Fbnd(:),'Sph SDE, Fbnd:')
    !call print_vector(Fphi(:,ichain),'Sph SDE, Fphi:')
    !call print_vector(Fphi_all(:),'Sph SDE, Fphi_all:')
    !call print_vector(Fphi_all_temp(:),'Sph SDE, Fphi_all_temp:')
    !call print_vector(dr_sph,'Sph SDE, dr_sph before:')

    !!displacement of the sphere due to the force on the sphere and HI from beads
    call gemv(DiffTensP(nbeadx3+1:nbeadx3+3,:),Fphi_all_temp,dr_sph,alpha=0.25*dt(iPe,idt),beta=1._wp)

    !call print_vector(dr_sph,'Sph SDE, dr_sph after:')
    !call print_matrix(wbltempP1,'Sph SDE, wbltempP1')
    !call print_vector(wbltempP1(:,jcol),'Sph SDE, wbltempP1')

    !! displacement of the sphere due to Brownian motion 1/sqrt(2) * C*dW
    !wbltempP1 is the vector (C*dW) for all beads and the sphere
    !coeff is 1/sqrt(2)
    !print *,'coeff is ',coeff
    dr_sph(1:3) = dr_sph(1:3) + coeff*wbltempP1(nbeadx3+1:nbeadx3+3,jcol)

    !call print_vector(dr_sph,'Sph SDE, dr_sph after Brownian:')

    !moving the sphere
    r_sph(:,ichain) = r_sph(:,ichain) + dr_sph(1:3)

    !update the tether points on the surface of the sphere
    do ichain_pp = 1,nchain_pp
      rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph(1:3)
    end do

    !call print_vector(F_sph(:),'F_sph after:')
    !Translational motion of sphere moved to sde_mod
    ! !TRANSLATION------------------------------------------------------
    ! !advance the coordinate of the sphere: tether forces
    ! do ichain_pp = 1,nchain_pp
    !   !using Fseg
    !   offset = nseg_indx3*(ichain_pp-1)
    !   dr_sph(1) = dr_sph(1) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+1)*dt(iPe,idt)
    !   dr_sph(2) = dr_sph(2) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+2)*dt(iPe,idt)
    !   dr_sph(3) = dr_sph(3) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fseg(offset+3)*dt(iPe,idt)
    !
    !   if (debug_TYL) then
    !     call print_vector(Fseg(offset+1:offset+3),'Sphere SDE, spring force:')
    !   end if
    !
    ! end do
    !
    ! !advance the coordinate of the sphere: bead-sphere EV force
    ! ! dr_sph(1) = dr_sph(1) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fev_sph(1)*dt(iPe,idt)
    ! ! dr_sph(2) = dr_sph(2) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fev_sph(2)*dt(iPe,idt)
    ! ! dr_sph(3) = dr_sph(3) + (1._wp/4)*(hstar*sqrtPI/this%a_sph)*Fev_sph(3)*dt(iPe,idt)
    !
    ! !advance the coordinate of the sphere: Brownian motion
    ! dr_sph(1) = dr_sph(1) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(1,jcol)
    ! dr_sph(2) = dr_sph(2) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(2,jcol)
    ! dr_sph(3) = dr_sph(3) + (1._wp/sqrt(2._wp))*sqrt(hstar*sqrtPI/this%a_sph)*wbl_sph(3,jcol)
    ! r_sph(:,ichain) = r_sph(:,ichain) + dr_sph(1:3)
    !
    ! !update the tether points on the surface of the sphere
    ! do ichain_pp = 1,nchain_pp
    !   rf0(:,ichain,ichain_pp) = rf0(:,ichain,ichain_pp) + dr_sph(1:3)
    ! end do

  end subroutine advance_sph_sde

end module sph_sde_mod
