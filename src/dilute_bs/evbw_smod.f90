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
  character(len=99),parameter :: fmt="(1x,i3,1x,i7)"

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
        ! allocate(evbw_prm%ia_time(2:nbead,10000)) !! should be fixed
        ! Initializing the variables
        evbw_prm%w_coll=0
        evbw_prm%ia_time=1
        if (id == 0) then
          allocate(evbw_prm%w_coll_t(2:nbead))
          open(newunit=evbw_prm%u0,file='data/w_coll.dat',status='replace',position='append')
          write(evbw_prm%u0,*) "# bead index, Total number of collisions #"
          write(evbw_prm%u0,*) "# -------------------------------------- #"
        end if
        ! write(fnme,"(A,i0.2,'.dat')") 'data/ia_time',id

        ! allocate(evbw_prm%ia_time_t(2:nbead))

        ! open(newunit=uarm(iarm),file=trim(adjustl(fnme)),&
        !      status='replace',position='append')
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
    !evbw_prm%w_coll=evbw_prm%w_coll-floor( Ry(2:nbead)/qmax )
    ! do ib=2, nbead
    !   evbw_prm%ia_time(ib,evbw_prm%w_coll(ib)) = &
    !   evbw_prm%ia_time(ib,evbw_prm%w_coll(ib))+1+floor( Ry(ib)/qmax )
    ! end do

    !Ry=abs(Ry)-2*evbw_prm%a*floor( Ry/qmax )
    rcmy=Ry(1)
    do ib=2, nbead
      if (Ry(ib) < evbw_prm%a) then
        evbw_prm%w_coll(ib)=evbw_prm%w_coll(ib)+1
        Ry(ib)=abs(Ry(ib))+2*evbw_prm%a
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
        if (id == 0) then
          deallocate(evbw_prm%w_coll_t)
        endif
    end  select

  end procedure del_evbw

  module procedure print_wcll

    use :: mpi
    use :: inp_dlt, only: nbead
    use :: arry_mod, only: print_vector

    integer :: ib,ierr

    call MPI_Reduce(evbw_prm%w_coll,evbw_prm%w_coll_t,nbead-1,&
                    MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    ! call MPI_Reduce(evbw_prm%w_coll,w_cll_tot,nbead-1,MPI_REAL_WP,&
    !   MPI_SUM,0,MPI_COMM_WORLD,ierr)
    ! if (id==0) then
    ! call print_vector(evbw_prm%w_coll,'collisionsid0')
    ! else
    ! call print_vector(evbw_prm%w_coll,'collisionsid1')
    ! endif
    if (id == 0) then
      write(evbw_prm%u0,*) 'TIME',time
      do ib=2,nbead
        write(evbw_prm%u0,fmt) ib,evbw_prm%w_coll_t(ib)
      enddo
      ! call print_vector(evbw_prm%w_coll_t,'total collisions')
    endif

  end procedure print_wcll

end submodule
