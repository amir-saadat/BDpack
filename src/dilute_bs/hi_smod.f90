! submodule (intrn_mod) hi_smod
module hi_smod

  use :: prcn_mod
  use :: hibb_smod, only: hibb_t
  use :: hibw_smod, only: hibw_t

  implicit none

  type :: hi_t
    type(hibb_t) :: hibb
    type(hibw_t) :: hibw
  end type hi_t

  !--------------------------------------------
  !>>>> Interface to routines for hibb and hibw
  !--------------------------------------------
  interface
    !> initializes HI between beads
    !! \param this hibb object
    ! module subroutine init_hibb(this)
    !   class(hibb_t),intent(inout) :: this
    ! end subroutine init_hibb
    !> initializes HI between beads and the wall
    !! \param this hibw object
    ! module subroutine init_hibw(this)
    !   class(hibw_t),intent(inout) :: this
    ! end subroutine init_hibw
    !> calculates HI between the beads
    !! \param this hibb object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    ! module subroutine calc_hibb(this,i,j,rij,DiffTens)
    !   class(hibb_t),intent(inout) :: this
    !   integer,intent(in) :: i,j
    !   type(dis),intent(in) :: rij
    !   real(wp),intent(out) :: DiffTens(:,:)
    ! end subroutine calc_hibb
    !> calculates HI between the beads and the wall
    !! \param this hibw object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    ! module subroutine calc_hibw(this,i,j,rij,DiffTens)
    !   class(hibw_t),intent(inout) :: this
    !   integer,intent(in) :: i,j
    !   type(dis),intent(in) :: rij
    !   real(wp),intent(inout) :: DiffTens(:,:)
    ! end subroutine calc_hibw
  end interface

contains

  ! module procedure init_hi
  subroutine init_hi(this)

    use :: hibb_smod, only: init_hibb
    use :: hibw_smod, only: init_hibw

    class(hi_t),intent(inout) :: this

    call init_hibb(this%hibb)
    call init_hibw(this%hibw)

  ! end procedure init_hi
  end subroutine init_hi

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

  ! module procedure calc_hi
  subroutine calc_hi(this,i,j,rij,DiffTens)

    use :: cmn_tp_mod, only: dis
    use :: hibb_smod, only: calc_hibb
    use :: hibw_smod, only: calc_hibw

    use :: inp_dlt, only: HITens

    class(hi_t),intent(inout) :: this
    integer,intent(in) :: i,j
    type(dis),intent(in) :: rij
    real(wp),intent(inout) :: DiffTens(:,:)

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


  ! end procedure calc_hi
  end subroutine calc_hi

  ! module procedure calc_div
  subroutine calc_div(j,rjy,divD)

    use :: inp_dlt, only: hstar

    integer,intent(in) :: j
    real(wp),intent(in) :: rjy
    real(wp),intent(inout) :: divD(:)

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    !divD(j)=1.125*sqrtPI*hstar/rjy**2 - 1.5*(sqrtPI*hstar)**3/rjy**4


    !debugging
    if (j == 1) then
      divD(j)=0
    else
      divD(j)=1.125*sqrtPI*hstar/rjy**2 - 1.5*(sqrtPI*hstar)**3/rjy**4
    end if

    !print *, 'divD(',j, ') is: ', divD(j)
    !print *, 'rjy is', rjy

  ! end procedure calc_div
  end subroutine calc_div

! end submodule hi_smod
end module hi_smod
