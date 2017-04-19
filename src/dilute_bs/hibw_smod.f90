submodule (intrn_mod:hi_smod) hibw_smod

  implicit none

contains

  module procedure init_hibw
    !add later
  end procedure init_hibw

  module procedure calc_hibw

    use :: inp_dlt, only: HITens, hstar
    use :: arry_mod, only: print_vector, print_matrix

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    integer :: osi,osj
    real(wp) :: rijmag3im, rijmag5im, rijmag7im
    real(wp), dimension(3,3) :: Sij,Pij_D,Sij_D,Hij_PD,Hji_PD
    real(wp), dimension(3,3) :: Omega_PF,Omegaij_PD,Omegaji_PD,Omega_c, Omega_W
    !real(wp), dimension(3,3) :: Omega_PD_trans

    osi=3*(i-1)
    osj=3*(j-1)

    rijmag3im = rij%mag2im*rij%magim
    rijmag5im = rij%mag2im*rijmag3im
    rijmag7im = rij%mag2im*rijmag5im
    !call print_matrix(DiffTens,'DiffTens(osi+1:osi+3,osj+1:osj+3)')
    if (HITens == 'Blake') then
      call Sij_calc(Sij,rij,rijmag3im,rijmag5im,rijmag7im)
      call Pij_D_calc(Pij_D,rij,rijmag3im,rijmag5im,rijmag7im)
      call Sij_D_calc(Sij_D,rij,rijmag3im,rijmag5im,rijmag7im)
      call H_PD_calc(Hij_PD,rij,rijmag3im,rijmag5im,rijmag7im,rij%x,rij%yim,rij%z,rij%rjy)
      call H_PD_calc(Hji_PD,rij,rijmag3im,rijmag5im,rijmag7im,-rij%x,rij%yim,-rij%z,rij%riy)

      Omega_PF = (3*sqrtPI*hstar/4)*(-Sij + 2*(rij%rjy**2)*Pij_D - 2*rij%rjy*Sij_D)
      Omegaij_PD = (3*sqrtPI*hstar/4)*Hij_PD
      Omegaji_PD = (3*sqrtPI*hstar/4)*Hji_PD
      Omega_c = (1._wp/2)*(Omegaij_PD+TRANSPOSE(Omegaji_PD))
      Omega_W = Omega_PF - (2*((sqrtPI*hstar)**2)/3)*Omega_c

      !for debugging...
      ! print *, 'x is ', rij%x
      ! print *, 'yim is ', rij%yim
      ! print *, 'y is ', rij%y
      ! print *, 'z is ', rij%z
      ! print *, 'rjy is ', rij%rjy
      ! print *, 'i is ', i
      ! print *, 'j is ', j
      ! call print_matrix(Sij,'Sij')
      ! call print_matrix(Pij_D,'Pij_D')
      ! call print_matrix(Sij_D,'Sij_D')
      ! !call print_matrix(Hij_PD,'H_PD')
      ! call print_matrix(Omega_PF,'Omega_PF')
      ! !call print_matrix(Omegaij_PD,'Omegaij_PD')
      ! call print_matrix(Omega_c,'Omega_c')
      ! call print_matrix(Omega_W,'Omega_W')

      if ((rij%rjy < hstar*sqrtPI) .OR. (rij%riy < hstar*sqrtPI) .OR. (rij%mag<2*hstar*sqrtPI))then
        DiffTens(osi+1:osi+3,osj+1:osj+3) = DiffTens(osi+1:osi+3,osj+1:osj+3)
      else
        DiffTens(osi+1:osi+3,osj+1:osj+3) = DiffTens(osi+1:osi+3,osj+1:osj+3) &
                                          + Omega_W
      end if
      !call print_matrix(DiffTens(osi+1:osi+3,osj+1:osj+3),'DiffTens(osi+1:osi+3,osj+1:osj+3)')
       !if (i == 3 .AND. j == 3)then
        !call print_matrix(DiffTens,'DiffTens(osi+1:osi+3,osj+1:osj+3)')
    !    call print_matrix(DiffTens-TRANSPOSE(DiffTens),'Difftens - transpose')
      !end if
    end if

  end procedure calc_hibw


  subroutine Sij_Calc(SijTens,rij,rijmag3im,rijmag5im,rijmag7im)

    !use :: arry_mod, only: print_vector, print_matrix


    real(wp), intent(inout) :: SijTens(:,:)
    type(dis),intent(in) :: rij
    real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im

    SijTens(1,1) = 1/rij%magim + rij%x * rij%x   / rijmag3im
    SijTens(1,2) = 0           + rij%x * rij%yim / rijmag3im
    SijTens(1,3) = 0           + rij%x * rij%z   / rijmag3im

    SijTens(2,1) = 0           + rij%x   * rij%yim / rijmag3im
    SijTens(2,2) = 1/rij%magim + rij%yim * rij%yim / rijmag3im
    SijTens(2,3) = 0           + rij%yim * rij%z   / rijmag3im

    SijTens(3,1) = 0           + rij%x   * rij%z / rijmag3im
    SijTens(3,2) = 0           + rij%yim * rij%z / rijmag3im
    SijTens(3,3) = 1/rij%magim + rij%z   * rij%z / rijmag3im
    !call print_matrix(SijTens,'Sij tensor')

  end subroutine Sij_Calc

  subroutine Pij_D_Calc(Pij_DTens,rij,rijmag3im,rijmag5im,rijmag7im)
    !use :: arry_mod, only: print_vector, print_matrix
    real(wp), intent(inout) :: Pij_DTens(:,:)
    type(dis),intent(in) :: rij
    real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im


    Pij_DTens(1,1) =   1/rijmag3im - 3*rij%x * rij%x  /rijmag5im
    Pij_DTens(1,2) = -(0           - 3*rij%x * rij%yim/rijmag5im)
    Pij_DTens(1,3) =   0           - 3*rij%x * rij%z  /rijmag5im

    Pij_DTens(2,1) =   0           - 3*rij%x   * rij%yim/rijmag5im
    Pij_DTens(2,2) = -(1/rijmag3im - 3*rij%yim * rij%yim/rijmag5im)
    Pij_DTens(2,3) =   0           - 3*rij%yim * rij%z  /rijmag5im

    Pij_DTens(3,1) =   0           - 3*rij%x   * rij%z/rijmag5im
    Pij_DTens(3,2) = -(0           - 3*rij%yim * rij%z/rijmag5im)
    Pij_DTens(3,3) =   1/rijmag3im - 3*rij%z   * rij%z/rijmag5im

    !call print_matrix(Pij_DTens,'Pij_D tensor')
  end subroutine Pij_D_Calc

  subroutine Sij_D_Calc(Sij_DTens,rij,rijmag3im,rijmag5im,rijmag7im)
    !use :: arry_mod, only: print_vector, print_matrix
    real(wp), intent(inout) :: Sij_DTens(:,:)
    type(dis),intent(in) :: rij
    real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im


    !!! please check
    call Pij_D_Calc(Sij_DTens,rij,rijmag3im,rijmag5im,rijmag7im)

    Sij_DTens = rij%yim * Sij_DTens

    Sij_DTens(1,1) = Sij_DTens(1,1)
    Sij_DTens(1,2) = Sij_DTens(1,2) - (rij%x /rijmag3im)
    Sij_DTens(1,3) = Sij_DTens(1,3)

    Sij_DTens(2,1) = Sij_DTens(2,1) + (-rij%x /rijmag3im)
    Sij_DTens(2,2) = Sij_DTens(2,2)
    Sij_DTens(2,3) = Sij_DTens(2,3) + (-rij%z /rijmag3im)

    Sij_DTens(3,1) = Sij_DTens(3,1)
    Sij_DTens(3,2) = Sij_DTens(3,2) - (rij%z /rijmag3im)
    Sij_DTens(3,3) = Sij_DTens(3,3)
    !call print_matrix(Sij_DTens,'Sij_D tensor')
  end subroutine Sij_D_Calc


  ! subroutine H_PD_Calc(HTens,rij,rijmag3im,rijmag5im,rijmag7im)
  !   !use :: arry_mod, only: print_vector, print_matrix
  !   real(wp), intent(inout) :: HTens(:,:)
  !   type(dis),intent(in) :: rij
  !   real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im
  !
  !
  !   HTens(1,1) = (1/rijmag3im) &
  !                - 3*(rij%x**2)/rijmag5im &
  !                - 6*(rij%yim)*(rij%yim - rij%ry) / rijmag5im &
  !                + 30*(rij%x**2)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im
  !
  !   HTens(1,3) = - 3* (rij%x)* (rij%z) / rijmag5im  &
  !                + 30*(rij%x)*(rij%z)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im
  !
  !   HTens(3,1) = HTens(1,3)
  !
  !   HTens(1,2) = (rij%x)*(9*(rij%yim) - 6*(rij%ry))/ rijmag5im &
  !                -30*(rij%x)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im
  !
  !   HTens(3,3) = (1/rijmag3im) &
  !                - 3*(rij%z**2)/rijmag5im &
  !                - 6*(rij%yim)*(rij%yim - rij%ry) / rijmag5im &
  !                + 30*(rij%z**2)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im
  !
  !   HTens(3,2) = (rij%z)*(9*(rij%yim) - 6*(rij%ry))/ rijmag5im &
  !                -30*(rij%z)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im
  !
  !   HTens(2,1) = -3*(rij%x) * (rij%yim - 2*rij%ry) / rijmag5im &
  !                +30*(rij%x)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im
  !
  !   HTens(2,3) = -3*(rij%z) * (rij%yim - 2*rij%ry) / rijmag5im &
  !                +30*(rij%z)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im
  !
  !   HTens(2,2) = (1/rijmag3im) &
  !                + 3*(rij%yim)*(5*rij%yim - 6*rij%ry) / rijmag5im &
  !                - 30*(rij%yim**3)*(rij%yim - rij%ry)/rijmag7im
  !
  !   !call print_matrix(HTens,'H tensor')
  ! end subroutine H_PD_Calc



  subroutine H_PD_Calc(HTens,rij,rijmag3im,rijmag5im,rijmag7im,x,yim,z,ry)
    !use :: arry_mod, only: print_vector, print_matrix
    real(wp), intent(inout) :: HTens(:,:)
    type(dis),intent(in) :: rij
    real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im
    real(wp), intent(in) :: x,yim,z,ry

    HTens(1,1) = (1/rijmag3im) &
                 - 3*(x**2)/rijmag5im &
                 - 6*(yim)*(yim - ry) / rijmag5im &
                 + 30*(x**2)*(yim)*(yim - ry)/rijmag7im

    HTens(1,3) = - 3* (x)* (z) / rijmag5im  &
                 + 30*(x)*(z)*(yim)*(yim - ry)/rijmag7im

    HTens(3,1) = HTens(1,3)

    HTens(1,2) = (x)*(9*(yim) - 6*(ry))/ rijmag5im &
                 -30*(x)*(yim**2)*(yim - ry) /rijmag7im

    HTens(3,3) = (1/rijmag3im) &
                 - 3*(z**2)/rijmag5im &
                 - 6*(yim)*(yim - ry) / rijmag5im &
                 + 30*(z**2)*(yim)*(yim - ry)/rijmag7im

    HTens(3,2) = (z)*(9*(yim) - 6*(ry))/ rijmag5im &
                 -30*(z)*(yim**2)*(yim - ry) /rijmag7im

    HTens(2,1) = -3*(x) * (yim - 2*ry) / rijmag5im &
                 +30*(x)*(yim**2)*(yim - ry) /rijmag7im

    HTens(2,3) = -3*(z) * (yim - 2*ry) / rijmag5im &
                 +30*(z)*(yim**2)*(yim - ry) /rijmag7im

    HTens(2,2) = (1/rijmag3im) &
                 + 3*(yim)*(5*yim - 6*ry) / rijmag5im &
                 - 30*(yim**3)*(yim - ry)/rijmag7im

    !call print_matrix(HTens,'H tensor')
  end subroutine H_PD_Calc


end submodule hibw_smod
