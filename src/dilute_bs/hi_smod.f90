submodule (intrn_mod) hi_smod

  implicit none

  !--------------------------------------------
  !>>>> Interface to routines for hibb and hibw
  !--------------------------------------------
  interface
    !> initializes HI between beads
    !! \param this hibb object
    module subroutine init_hibb(this)
      class(hibb_t),intent(inout) :: this
    end subroutine init_hibb
    !> initializes HI between beads and the wall
    !! \param this hibw object
    module subroutine init_hibw(this)
      class(hibw_t),intent(inout) :: this
    end subroutine init_hibw
    !> calculates HI between the beads
    !! \param this hibb object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    module subroutine calc_hibb(this,i,j,rij,DiffTens)
      class(hibb_t),intent(inout) :: this
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(out) :: DiffTens(:,:)
    end subroutine calc_hibb
    !> calculates HI between the beads and the wall
    !! \param this hibw object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    module subroutine calc_hibw(this,i,j,rij,DiffTens)
      class(hibw_t),intent(inout) :: this
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: DiffTens(:,:)
    end subroutine calc_hibw


  end interface

contains

  module procedure init_hi
    call init_hibb(this%hibb)
    call init_hibw(this%hibw)
  end procedure init_hi

  ! module procedure hi_init

  !   use :: inp_dlt, only: HITens,hstar

  !   real(wp),parameter :: PI=3.1415926535897958648_wp
  !   real(wp),parameter :: sqrtPI=sqrt(PI)

  !   if (HITens == 'RPY') then
  !     ! For Rotne-Prager-Yamakawa Tensor:
  !     hi_prm%A=0.75*hstar*sqrtPI
  !     hi_prm%B=hstar**3*PI*sqrtPI/2
  !     hi_prm%C=(3.0_wp/2)*hstar**3*PI*sqrtPI
  !     hi_prm%D=2*sqrtPI*hstar
  !     hi_prm%E=9/(32*sqrtPI*hstar)
  !     hi_prm%F=3/(32*sqrtPI*hstar)
  !   elseif (HITens == 'OB') then
  !     ! For Oseen-Burgers Tensor:
  !     hi_prm%G=0.75*hstar*sqrtPI
  !   elseif (HITens == 'RegOB') then
  !     ! For Reguralized Oseen-Burgers Tensor (introduced in HCO book):
  !     hi_prm%G=0.75*hstar*sqrtPI
  !     hi_prm%O=4*PI*hstar**2/3
  !     hi_prm%P=14*PI*hstar**2/3
  !     hi_prm%R=8*PI**2*hstar**4
  !     hi_prm%S=2*PI*hstar**2
  !     hi_prm%T=hi_prm%R/3
  !   end if
  !   hi_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  ! end procedure hi_init

  module procedure calc_hi


    use :: inp_dlt, only: HITens

    integer :: osi,osj


    osj=3*(j-1)
    osi=3*(i-1)

    if (i == j) then

      DiffTens(osi+1,osj+1)=1._wp
      DiffTens(osi+1,osj+2)=0._wp
      DiffTens(osi+1,osj+3)=0._wp
      DiffTens(osi+2,osj+2)=1._wp
      DiffTens(osi+2,osj+3)=0._wp
      DiffTens(osi+3,osj+3)=1._wp

      !Tiras addition
      DiffTens(osi+2,osj+1)=0._wp
      DiffTens(osi+3,osj+1)=0._wp
      DiffTens(osi+3,osj+2)=0._wp

    else

      call calc_hibb(this%hibb,i,j,rij,DiffTens)

    endif



    ! Blake's part
    if (HITens == 'Blake') then
      call calc_hibw(this%hibw,i,j,rij,DiffTens)
    endif
    !------------
    ! Oseen solid sphere
    if (HITens == 'Osph') then
      if (i==j) then
        call calc_hibw(this%hibw,i,j,rij,DiffTens)
      else
        call calc_hibw(this%hibw,i,j,rij,DiffTens)
      endif
    endif
    !------------


  end procedure calc_hi

  module procedure calc_div

    use :: inp_dlt, only: hstar

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)
    real(wp) :: rjy_r
    !!!!!!variables that were used in rounding divD
    !integer :: digits
    !integer,parameter :: longlong=selected_int_kind(18)
    !integer(kind=longlong) :: int_rjy,int_divD

    !digits = 9

    if (j == 1) then
      divD(j)=0._wp
    else

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Note on numerical errors arising from calculation of divD:
      !T.Y. Lin, A. Saadat
      !Aug. - Dec. 2017
      !
      !While working on simulations on wall tethered polymers in shear,
      !we observed that there were numerical errors arising due to
      !calculation of divD. Specifically, while working on the default
      !(32ppn) nodes of Certainty, we noticed that results were slightly
      !different between nodes. We narrowed down the issue to the
      !calculation of divD, in particular the calculation of rjy^-2
      !and rjy^-4. To mitigate this issue, we tried combination of:
      !(1) rounding rjy first
      !(2) rounding divD at the end of the calculation
      !but neither seemed to fully surpress the issue. We ultimately
      !decided to calculate divD without rounding, because
      !(1) the errors are small
      !(2) it's unclear what inaccuracies are introduced by this
      !    artificial rounding.
      !The attempts are documented in the comments below. The necessary
      !declared variables are also commented above.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! ! Rounding rjy (making 11+ digits after decimal point zeros)
      ! ! if we want higher precision (not available in ifort 16 or less)
      ! ! int_rjy=int((real(rjy,kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! ! rjy_hp=int_rjy/(10._real128**digits)
      ! int_rjy=int(rjy*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! rjy_r=int_rjy/(10._wp**digits)
      !
      ! ! Calculating divD
      ! ! if we want higher precision (not available in ifort 16 or less)
      ! ! divD(j)=1.125*sqrtPI*hstar/rjy_hp**2 - 1.5*(sqrtPI*hstar)**3/rjy_hp**4
      ! divD(j)=1.125*sqrtPI*hstar/rjy_r**2 - 1.5*(sqrtPI*hstar)**3/rjy_r**4
      !
      !
      ! ! Rounding divD (making 11+ digits after decimal point zeros)
      ! ! if we want higher precision (not available in ifort 16 or less)
      ! ! int_divD=int((real(divD(j),kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! ! divD(j)=int_divD/(10._real128**digits)
      ! int_divD=int(divD(j)*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! divD(j)=int_divD/(10._wp**digits)

      divD(j)=1.125*sqrtPI*hstar/rjy**2 - 1.5*(sqrtPI*hstar)**3/rjy**4
      !divD(j) = 0
    end if

  end procedure calc_div

end submodule hi_smod
