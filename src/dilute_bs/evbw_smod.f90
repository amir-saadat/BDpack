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
  character(len=99),parameter :: fmt3xi="(1x,i3,1x,i3,1x,i7)"

contains

  module procedure init_evbw

  end procedure init_evbw

  module procedure evbw_init

    use :: inp_dlt, only: nseg,nbead,EV_bw,Aw,N_Ks,qmax,ntime,npchain,nchain
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

        allocate(evbw_prm%w_coll(2:nbead,npchain))
        allocate(evbw_prm%ia_time(2:nbead,500,npchain))

        ! Initializing the variables

        evbw_prm%w_coll=0
        evbw_prm%ia_time=1


        if (id == 0) then
          allocate(evbw_prm%w_coll_t(2:nbead,npchain))
          allocate(evbw_prm%ia_time_t(2:nbead,500,npchain))
          open(newunit=evbw_prm%u_wc,file='data/w_coll.dat',status='replace',position='append')
          write(evbw_prm%u_wc,*) "# chain index, bead index, Total number of collisions #"
          write(evbw_prm%u_wc,*) "# --------------------------------------------------- #"
          open(newunit=evbw_prm%u_ia,file='data/ia_time.dat',status='replace',position='append')
          write(evbw_prm%u_ia,*) "# chain index, bead index, Inter-arrival time unit #"
          write(evbw_prm%u_ia,*) "# ------------------------------------------------ #"
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

    use :: mpi
    use :: inp_dlt, only: nbead,qmax,tplgy,npchain,lambda
    use :: arry_mod, only: print_vector

    integer :: ib,ierr,sz,sz_t
    integer,allocatable :: ia_tmp(:,:,:)

    if ((it == 1)) then
      evbw_prm%w_coll(:,ich)=0
      evbw_prm%ia_time(:,:,ich)=1
    endif

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


        if (evbw_prm%ia_time(ib,evbw_prm%w_coll(ib,ich)+1,ich) > int(lambda/dt)) then

          evbw_prm%w_coll(ib,ich)=evbw_prm%w_coll(ib,ich)+1

        else

          evbw_prm%ia_time(ib,evbw_prm%w_coll(ib,ich)+1,ich) = 1

        endif


        Ry(ib)=abs(Ry(ib))+2*evbw_prm%a
        select case (tplgy)
          case ('Linear')
            qy(ib-1)=Ry(ib)-Ry(ib-1)
            if (ib < nbead) &
              qy(ib)=Ry(ib+1)-Ry(ib)
          case ('Comb')
        end select

      else


        evbw_prm%ia_time(ib,evbw_prm%w_coll(ib,ich)+1,ich) = &
        evbw_prm%ia_time(ib,evbw_prm%w_coll(ib,ich)+1,ich) + 1


      end if


      sz=size(evbw_prm%ia_time,dim=2)

      if ( evbw_prm%w_coll(ib,ich)+1 > sz ) then
        print '(" Geometric resizing of ia_time array... in rank: ",i5)',id
        sz=2*size(evbw_prm%ia_time,dim=2)
        allocate(ia_tmp(2:nbead,sz,npchain))
        ia_tmp=1
        ia_tmp(:,1:sz/2,:)=evbw_prm%ia_time(:,:,:)
        call move_alloc(from=ia_tmp,to=evbw_prm%ia_time)
      endif

      call MPI_Reduce(sz,sz_t,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      if (id == 0) then
        if ( sz_t > size(evbw_prm%ia_time_t,dim=2) ) then
          print '(" Geometric resizing of ia_time_t array: ",i5)',id
          allocate(ia_tmp(2:nbead,sz_t,npchain))
          ! ia_tmp=1
          ! ia_tmp(1:sz_t/2)=evbw_prm%ia_time_t
          call move_alloc(from=ia_tmp,to=evbw_prm%ia_time_t)
        endif
      endif


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
          deallocate(evbw_prm%ia_time_t)
        endif
    end  select

  end procedure del_evbw

  module procedure print_wcll

    use :: mpi
    use :: inp_dlt, only: nbead,npchain,ntime
    use :: arry_mod, only: print_vector

    integer :: ich,ib,iwc,osch,ierr,ncount_wc,ncount_ia,iproc,tag



    ncount_wc=(nbead-1)*npchain
    ncount_ia=(nbead-1)*size(evbw_prm%ia_time,dim=2)*npchain

    if (id == 0) then

      write(evbw_prm%u_wc,'(" TIME: ",f9.2)') time
      write(evbw_prm%u_ia,'(" TIME: ",f9.2)') time

      do ich=1, npchain
        do ib=2,nbead
          write(evbw_prm%u_wc,fmt3xi) ich,ib,evbw_prm%w_coll(ib,ich)
          do iwc=1, evbw_prm%w_coll(ib,ich)
            write(evbw_prm%u_ia,fmt3xi) ich,ib,evbw_prm%ia_time(ib,iwc,ich)
          enddo
        enddo
      enddo

      ncount_ia=(nbead-1)*size(evbw_prm%ia_time_t,dim=2)*npchain

      do iproc=1, nproc-1

        tag=1100+iproc

        call MPI_Recv(evbw_prm%w_coll_t,ncount_wc,MPI_INTEGER,iproc,tag,&
          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        call MPI_Recv(evbw_prm%ia_time_t,ncount_ia,MPI_INTEGER,iproc,tag,&
          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        osch=iproc*npchain

        do ich=1, npchain
          do ib=2,nbead
            write(evbw_prm%u_wc,fmt3xi) osch+ich,ib,evbw_prm%w_coll_t(ib,ich)
            do iwc=1, evbw_prm%w_coll_t(ib,ich)
              write(evbw_prm%u_ia,fmt3xi) osch+ich,ib,evbw_prm%ia_time_t(ib,iwc,ich)
            enddo
          enddo
        enddo

      enddo

    else

      tag=1100+id

      call MPI_Send(evbw_prm%w_coll,ncount_wc,MPI_INTEGER,0,tag,&
        MPI_COMM_WORLD,ierr)

      call MPI_Send(evbw_prm%ia_time,ncount_ia,MPI_INTEGER,0,tag,&
        MPI_COMM_WORLD,ierr)

    endif

    ! wait untill receiving all values
    call MPI_Barrier(MPI_COMM_WORLD,ierr)



    ! call MPI_Reduce(evbw_prm%w_coll,evbw_prm%w_coll_t,nbead-1,&
    !                   MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)



    ! call MPI_Reduce(evbw_prm%w_coll,w_cll_tot,nbead-1,MPI_REAL_WP,&
    !   MPI_SUM,0,MPI_COMM_WORLD,ierr)
    ! if (id==0) then
    ! call print_vector(evbw_prm%w_coll,'collisionsid0')
    ! else
    ! call print_vector(evbw_prm%w_coll,'collisionsid1')
    ! endif
    ! if (id == 0) then
    !   write(evbw_prm%u_wc,'(" TIME: ",f9.2)') time
    !   do ib=2,nbead
    !     write(evbw_prm%u_wc,fmtii) ib,evbw_prm%w_coll_t(ib)
    !   enddo
      ! call print_vector(evbw_prm%w_coll_t,'total collisions')
    ! endif

  end procedure print_wcll


end submodule
