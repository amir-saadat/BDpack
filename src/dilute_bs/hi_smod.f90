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


  end procedure calc_hi

  module procedure calc_div

    use :: inp_dlt, only: hstar
    use :: mpi

    ! real(real128),parameter :: PI=3.1415926535897958648_real128
    ! real(real128),parameter :: sqrtPI=sqrt(PI)
    ! real(wp) :: rjy4_r,rjy2_r!,rjy_inv
    ! real(real128) :: rjy_r,rjy_hp
    ! real(real128) :: rjy2_hp, rjy_inv_hp, a3_hp, divD_hp,a_rjy_hp
    ! integer :: digits
    ! integer,parameter :: longlong=selected_int_kind(18)
    ! integer(kind=longlong) :: int_rjy4,int_rjy2,int_rjy,int_divD
    ! integer :: ierr

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)
    real(wp) :: rjy_r
    integer :: digits
    integer,parameter :: longlong=selected_int_kind(18)
    integer(kind=longlong) :: int_rjy,int_divD

    !divD(j)=1.125*sqrtPI*hstar/rjy**2 - 1.5*(sqrtPI*hstar)**3/rjy**4

    digits = 10

    if (j == 1) then
      divD(j)=0._wp
    else

















!!!!! Uncomment this part if you want to go back to where we were
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! !rjy_inv=1/rjy
      ! ! print*,'before',rjy**4,int((rjy**4)*10**15,kind=longlong)
      ! ! int_rjy=int((rjy**4)*10**15,kind=longlong)
      ! ! rjy_r=int_rjy/10._wp**15
      ! ! print*,'after',rjy_r
      ! ! stop


      ! ! print*,'before',rjy**4,(rjy**4)*int(10m,kind=longlong)**15
      ! ! int_rjy=(rjy**4)*int(10,kind=longlong)**15
      ! ! rjy_r=int_rjy/10._wp**15
      ! ! print*,'after',rjy_r




      ! int_rjy4=int((rjy**4)*10_longlong**digits,kind=longlong) ! integer of rjy^4: cast to longlong integer
      ! int_rjy2=int((rjy**2)*10_longlong**digits,kind=longlong) ! integer of rjy^2: cast to longlong integer
      ! int_rjy=int((real(rjy,kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer





      ! rjy4_r=int_rjy4/(10._wp**digits)
      ! rjy2_r=int_rjy2/(10._wp**digits)
      ! rjy_r=int_rjy/(10._real128**digits)

      ! if ( (id==1 .and. itime==1) .or. &
      !      (id==1 .and. itime==19) ) then
      !   if (j<=3) then
      !     print*,'rjy',rjy_r
      !   endif
      ! endif

      ! !-------------- debug
      ! ! print*,'-------------------------'
      ! ! print*,'before r^4: ',rjy**4
      ! ! print*,'before r^2: ',rjy**2
      ! !
      ! ! print*,'integer r^4: ',int_rjy4
      ! ! print*,'integer r^2: ',int_rjy2
      ! !
      ! ! print*,'after r^4: ',rjy4_r
      ! ! print*,'after r^2: ',rjy2_r
      ! !stop
      ! !-------------- debug

      ! !original
      ! !divD(j)=1.125*sqrtPI*hstar/rjy**2 - 1.5*(sqrtPI*hstar)**3/rjy**4

      ! !sketchy rounding
      ! !divD(j)=1.125_wp*sqrtPI*hstar/rjy2_r - 1.5_wp*(sqrtPI*hstar)**3/rjy4_r

      ! !factoring out rjy**2
      ! rjy_inv_hp = 1._real128/rjy_r
      ! if ( (id==1 .and. itime==1) .or. &
      !      (id==1 .and. itime==190) ) then
      !   if (j<=3) then
      !     print*,'rjy_inv',rjy_inv_hp
      !     ! print*,'sqrt',PI,sqrtPI
      !   endif
      ! endif
      ! a3_hp = sqrtPI*hstar*sqrtPI*hstar*sqrtPI*hstar

      ! a_rjy_hp = sqrtPI*.0795_real128*rjy_inv_hp
      ! !print*,sqrtPI
      ! !print*,real(hstar,kind=real128)
      ! !print*,'useless?   ',.795_real128
      ! !divD(j)= rjy_inv_hp*rjy_inv_hp*(1.125_real128*sqrtPI*hstar - 1.5_real128*a3_hp*rjy_inv_hp*rjy_inv_hp)
      ! !divD_hp = 1.125_real128*sqrtPI*hstar - 1.5_real128*a3_hp!*rjy_inv_hp*rjy_inv_hp
      ! !divD_hp = - 1.5_real128*a3_hp*rjy_inv_hp*rjy_inv_hp
      ! !divD_hp=divD_hp*rjy_inv_hp
      ! !divD_hp=divD_hp*rjy_inv_hp


      ! divD_hp = - 1.5_real128*a_rjy_hp*a_rjy_hp*a_rjy_hp*rjy_inv_hp
      ! ! if (j==2) then
      ! !   print*,divD_hp
      ! !
      ! !
      ! !   !write( *, * ) 'Press Enter to continue'
      ! !   !read( *, * )
      ! !   !call MPI_Barrier(MPI_COMM_WORLD,ierr)
      ! ! endif

      ! ! print*,'divD_hp1',j,divD_hp
      ! divD(j)=divD(j)+divD_hp
      ! int_divD=int((real(divD(j),kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! divD(j)=int_divD/(10._real128**digits)
      ! ! divD(j)=divD(j)*1000
      ! ! print*,'divD_hp2',j,divD_hp
      ! !divD(j)=divD(j)*rjy_inv_hp

      ! !high precision
      ! !rjy2_hp = rjy**2
      ! !divD(j)= (1/rjy2_hp)*( 1.125*sqrtPI*hstar - 1.5*(sqrtPI*hstar)**3/rjy2_hp)

      ! if ((id == 1) .and. (j==10) .and. (itime<=20)) then
      ! print*,'divD',divD(j)
      ! endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








      ! Rounding rjy (making 11+ digits after decimal point zeros)
      ! if we want higher precision (not available in ifort 16 or less)
      ! int_rjy=int((real(rjy,kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! rjy_hp=int_rjy/(10._real128**digits)
      int_rjy=int(rjy*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      rjy_r=int_rjy/(10._wp**digits)

      ! Calculating divD
      ! if we want higher precision (not available in ifort 16 or less)
      ! divD(j)=1.125*sqrtPI*hstar/rjy_hp**2 - 1.5*(sqrtPI*hstar)**3/rjy_hp**4
      divD(j)=1.125*sqrtPI*hstar/rjy_r**2 - 1.5*(sqrtPI*hstar)**3/rjy_r**4


      ! Rounding divD (making 11+ digits after decimal point zeros)
      ! if we want higher precision (not available in ifort 16 or less)
      ! int_divD=int((real(divD(j),kind=real128))*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      ! divD(j)=int_divD/(10._real128**digits)
      int_divD=int(divD(j)*10_longlong**digits,kind=longlong) ! integer of rjy: cast to longlong integer
      divD(j)=int_divD/(10._wp**digits)

    end if

  end procedure calc_div

end submodule hi_smod
