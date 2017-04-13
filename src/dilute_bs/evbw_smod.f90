submodule (intrn_mod) evbw_smod

  implicit none

  ! type :: evbw
  !   ! For Cubic
  !   real(wp) :: delw
  !   real(wp) :: prf
  !   real(wp) :: rmagmin
  !   ! For Reflc-bc
  !   real(wp) :: a
  !   integer,allocatable :: w_coll(:)
  !   integer,allocatable :: ia_time(:,:)
  ! contains
  !   procedure,pass(this) :: init => evbw_init
  !   procedure,pass(this) :: updt => evbw_updt
  !   procedure,pass(this) :: 
  ! end type

contains

  module procedure init_evbw

  end procedure init_evbw

  module procedure evbw_init
    
    use :: inp_dlt, only: nseg,nbead,EV_bw,Aw,N_Ks,qmax
    use :: cmn_io_mod, only: read_input

    ! Bead-wall excluded volume interaction
    select case (EV_bw)
      case ('Cubic')
!        evbw_prm%delw=0.5_wp*sqrt( (nseg**2-1._wp)/(2._wp*nseg) )
        evbw_prm%delw=0.5_wp*sqrt( 3.0 )
        evbw_prm%prf=Aw*N_Ks/(3*qmax*evbw_prm%delw**2)
        evbw_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|
      case ('Rflc_bc')
        call read_input('Bead-rad',0,evbw_prm%a)

        allocate(evbw_prm%w_coll(2:nbead))
        allocate(evbw_prm%ia_time(2:nbead,10000)) !! should be fixed
        ! Initializing the variables
        evbw_prm%w_coll=0
        evbw_prm%ia_time=1
    end  select

  end subroutine evbw_init

  module procedure evbwcalc

    use :: inp_dlt, only: EV_bw

    integer :: osi

    osi=3*(i-1)

    if (EV_bw == 'Cubic') then

      if (ry <= evbw_prm%delw ) then
        Fev(osi+2)=Fev(osi+2)+3*evbw_prm%prf*(ry-evbw_prm%delw)**2
      end if

    elseif (EV_bw == 'Gaussian') then
    end if ! EV_bw

  end procedure evbwcalc

  module procedure wall_rflc

    use :: inp_dlt, only: nbead,qmax 
    use :: arry_mod, only: print_vector

    integer :: ib

    ! To save memory, rcmy is added to Rvy to get rvy
    Ry=Ry+rcmy
    evbw_prm%w_coll=evbw_prm%w_coll-floor( Ry(2:nbead)/qmax )
    ! do ib=2, nbead
    !   evbw_prm%ia_time(ib,evbw_prm%w_coll(ib)) = &
    !   evbw_prm%ia_time(ib,evbw_prm%w_coll(ib))+1+floor( Ry(ib)/qmax )
    ! end do

    !Ry=abs(Ry)-2*evbw_prm%a*floor( Ry/qmax )
    rcmy=Ry(1)
    do ib=2, nbead
      if (Ry(ib) < evbw_prm%a) then
        Ry(ib)=abs(Ry(ib))!+2*evbw_prm%a
      end if
      rcmy=rcmy+Ry(ib)
    end do
    rcmy=rcmy/nbead
    Ry=Ry-rcmy

  end procedure wall_rflc

  module procedure del_evbw

    use :: inp_dlt, only: EV_bw

    select case (EV_bw)
      case ('Cubic')
      case ('Rflc_bc')
        deallocate(evbw_prm%w_coll)
        deallocate(evbw_prm%ia_time)
    end  select

  end procedure del_evbw


end submodule
