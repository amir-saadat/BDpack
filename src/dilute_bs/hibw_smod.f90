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
    integer :: digits
    integer,parameter :: longlong=selected_int_kind(18)
    integer(kind=longlong) :: int_rijmag3im,int_rijmag5im,int_rijmag7im

    integer :: osi,osj
    real(wp) :: rijmag3im, rijmag5im, rijmag7im
    real(wp) :: rijmag3im_frc, rijmag5im_frc, rijmag7im_frc
    real(wp) :: rijmag3im_frc_r, rijmag5im_frc_r, rijmag7im_frc_r
    real(wp), dimension(3,3) :: Sij,Pij_D,Sij_D,Hij_PD,Hji_PD
    real(wp), dimension(3,3) :: Omega_PF,Omegaij_PD,Omegaji_PD,Omega_c, Omega_W

    osi=3*(i-1)
    osj=3*(j-1)

    rijmag3im = rij%mag2im*rij%magim
    rijmag5im = rij%mag2im*rijmag3im
    rijmag7im = rij%mag2im*rijmag5im


    !rounding by using the method of rounding the "fraction" in base2
    ! digits = 7
    ! ! !print*,'before rij^3: ',rijmag3im
    ! ! !print*,'before rij^5: ',rijmag5im
    ! ! !print*,'before rij^7: ',rijmag7im
    ! rijmag3im_frc = FRACTION(rijmag3im)
    ! rijmag5im_frc = FRACTION(rijmag5im)
    ! rijmag7im_frc = FRACTION(rijmag7im)
    ! ! !print*,'before rij^3 frac: ',rijmag3im_frc
    ! ! !print*,'before rij^5 frac: ',rijmag5im_frc
    ! ! !print*,'before rij^7 frac: ',rijmag7im_frc
    ! int_rijmag3im = int(rijmag3im_frc*10_longlong**digits,kind=longlong)
    ! int_rijmag5im = int(rijmag5im_frc*10_longlong**digits,kind=longlong)
    ! int_rijmag7im = int(rijmag7im_frc*10_longlong**digits,kind=longlong)
    ! ! !print*,'int rij^3 frac: ',int_rijmag3im
    ! ! !print*,'int rij^5 frac: ',int_rijmag5im
    ! ! !print*,'int rij^7 frac: ',int_rijmag7im
    ! rijmag3im_frc_r = int_rijmag3im/10._wp**digits
    ! rijmag5im_frc_r = int_rijmag5im/10._wp**digits
    ! rijmag7im_frc_r = int_rijmag7im/10._wp**digits
    ! ! !print*,'rounded rij^3 frac: ',rijmag3im_frc_r
    ! ! !print*,'rounded rij^5 frac: ',rijmag5im_frc_r
    ! ! !print*,'rounded rij^7 frac: ',rijmag7im_frc_r
    ! rijmag3im = SET_EXPONENT(rijmag3im_frc_r,EXPONENT(rijmag3im))
    ! rijmag5im = SET_EXPONENT(rijmag5im_frc_r,EXPONENT(rijmag5im))
    ! rijmag7im = SET_EXPONENT(rijmag7im_frc_r,EXPONENT(rijmag7im))
    ! ! !print*,'after rij^3: ',rijmag3im
    ! ! !print*,'after rij^5: ',rijmag5im
    ! ! !print*,'after rij^7: ',rijmag7im
    ! ! !stop






    !rounding the results of rijmag3im,rijmag5im,rijmag7im
    !-------------- debug
    ! print*,'-------------------------'
    ! print*,'before rij^3: ',rijmag3im
    ! print*,'before rij^5: ',rijmag5im
    ! print*,'before rij^7: ',rijmag7im
    !-------------- debug
    ! digits = 7 !number of digits to round to
    ! int_rijmag3im=int(rijmag3im*10_longlong**digits,kind=longlong)
    ! int_rijmag5im=int(rijmag5im*10_longlong**digits,kind=longlong)
    ! int_rijmag7im=int(rijmag7im*10_longlong**digits,kind=longlong)
    ! rijmag3im=int_rijmag3im/(10._wp**digits)
    ! rijmag5im=int_rijmag5im/(10._wp**digits)
    ! rijmag7im=int_rijmag7im/(10._wp**digits)
    !-------------- debug
    ! print*,'integer rij^3: ',int_rijmag3im
    ! print*,'integer rij^5: ',int_rijmag5im
    ! print*,'integer rij^7: ',int_rijmag7im
    !
    ! print*,'after rij^3: ',rijmag3im
    ! print*,'after rij^5: ',rijmag5im
    ! print*,'after rij^7: ',rijmag7im
    ! stop
    !-------------- debug
    !done rounding

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

      DiffTens(osi+1:osi+3,osj+1:osj+3) = DiffTens(osi+1:osi+3,osj+1:osj+3) &
                                        + Omega_W
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
