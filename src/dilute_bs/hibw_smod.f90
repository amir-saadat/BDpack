submodule (intrn_mod:hi_smod) hibw_smod

  implicit none

contains

  module procedure init_hibw
    !add later
  end procedure init_hibw

  module procedure calc_hibw

    use :: inp_dlt, only: HITens, hstar

    integer :: osi,osj
    real(wp) :: rijmag3im, rijmag5im, rijmag7im
    real(wp), dimension(3,3) :: Sij,Pij_D,Sij_D,H_PD
    real(wp), dimension(3,3) :: Omega_PF,Omega_PD,Omega_c

    osi=3*(i-1)
    osj=3*(j-1)

    rijmag3im = rij%mag2im*rij%magim
    rijmag5im = rij%mag2im*rijmag3im
    rijmag7im = rij%mag2im*rijmag5im

    if (HITens == 'Blake') then
      call Sij_calc(Sij,rij,rijmag3im,rijmag5im,rijmag7im)
      call Pij_D_calc(Pij_D,rij,rijmag3im,rijmag5im,rijmag7im)
      call Sij_D_calc(Sij_D,rij,rijmag3im,rijmag5im,rijmag7im)
      call H_PD_calc(H_PD,rij,rijmag3im,rijmag5im,rijmag7im)

      Omega_PF = (3*hstar/4)*(-Sij + 2*(rij%ry**2)*Pij_D - 2*rij%ry*Sij_D)
      Omega_PD = (3*hstar/4)*H_PD
      Omega_c = (1/2)*(Omega_PD+TRANSPOSE(Omega_PD))

      DiffTens(osi+1:osi+3,osj+1:osj+3) = DiffTens(osi+1:osi+3,osj+1:osj+3) &
                                          + Omega_PF - (2*(hstar**2)/3)*Omega_c
    end if

  end procedure calc_hibw


  subroutine Sij_Calc(SijTens,rij,rijmag3im,rijmag5im,rijmag7im)
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
  end subroutine Sij_Calc

  subroutine Pij_D_Calc(Pij_DTens,rij,rijmag3im,rijmag5im,rijmag7im)
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
  end subroutine Pij_D_Calc

  subroutine Sij_D_Calc(Sij_DTens,rij,rijmag3im,rijmag5im,rijmag7im)
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
  end subroutine Sij_D_Calc

  subroutine H_PD_Calc(HTens,rij,rijmag3im,rijmag5im,rijmag7im)
    real(wp), intent(inout) :: HTens(:,:)
    type(dis),intent(in) :: rij
    real(wp), intent(in) :: rijmag3im,rijmag5im,rijmag7im

    HTens(1,1) = (1/rijmag3im) &
                 - 3*(rij%x**2)/rijmag5im &
                 - 6*(rij%yim)*(rij%yim - rij%ry) / rijmag5im &
                 + 30*(rij%x**2)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im

    HTens(1,3) = - 3* (rij%x)* (rij%z) / rijmag5im  &
                 + 30*(rij%x)*(rij%z)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im

    HTens(3,1) = HTens(1,3)

    HTens(1,2) = (rij%x)*(9*(rij%yim) - 6*(rij%ry))/ rijmag5im &
                 -30*(rij%x)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im

    HTens(3,3) = (1/rijmag3im) &
                 - 3*(rij%z**2)/rijmag5im &
                 - 6*(rij%yim)*(rij%yim - rij%ry) / rijmag5im &
                 + 30*(rij%z**2)*(rij%yim)*(rij%yim - rij%ry)/rijmag7im

    HTens(3,2) = (rij%z)*(9*(rij%yim) - 6*(rij%ry))/ rijmag5im &
                 -30*(rij%z)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im

    HTens(2,1) = -3*(rij%x) * (rij%yim - 2*rij%ry) / rijmag5im &
                 +30*(rij%x)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im

    HTens(2,3) = -3*(rij%z) * (rij%yim - 2*rij%ry) / rijmag5im &
                 +30*(rij%z)*(rij%yim**2)*(rij%yim - rij%ry) /rijmag7im

    HTens(2,2) = (1/rijmag3im) &
                 + 3*(rij%yim)*(5*rij%yim - 6*rij%ry) / rijmag5im &
                 - 30*(rij%yim**3)*(rij%yim - rij%ry)/rijmag7im

  end subroutine H_PD_Calc

end submodule hibw_smod
