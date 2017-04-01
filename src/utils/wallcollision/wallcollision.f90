!-------------------------------------------------------------------------
!Tiras Y. Lin
!
!March 24, 2017
!Program to count the number of times each bead contacts the wall and
!determines PDFs of the times between contacts and times spent at the wall
!-------------------------------------------------------------------------

program wallcollision

  implicit none

  !----------------------------------------
  !         variable declarations
  !----------------------------------------
  integer,parameter :: dp=selected_real_kind(p=15,r=307)

  integer :: narg,iarg,nchain,nseg,nbead,ntotseg,ichain,iseg
  integer :: iskip,nskip
  integer :: nseg_bb
  integer :: itime, ntime
  integer :: totcollisions
  integer :: iwall_idx,iaway_idx
  integer :: nbin,ibin
  real(dp),allocatable :: midbin(:,:)
  real(dp), allocatable :: Prob_t_wall(:,:),Prob_t_away(:,:)
  character(len=20) :: arg,buffer,dmpFile
  logical :: dmpExist
  real(dp) :: b,delw
  real(dp),allocatable :: q(:,:,:,:)
  real(dp),allocatable,dimension(:,:) :: Qeex,Qeey,Qeez
  real(dp),allocatable,dimension(:,:,:) :: rx,ry,rz
  logical, allocatable,dimension(:,:,:) :: nearwall
  integer,allocatable,target :: collisions(:)
  integer,allocatable,target :: t_wall(:,:), t_away(:,:)
  integer,allocatable :: t_wall_curr(:,:), t_away_curr(:,:), t_curr(:,:)
  integer,allocatable :: wall_idx(:),away_idx(:)
  real(dp),allocatable,dimension(:) :: rx_mean,ry_mean,rz_mean

  !----------------------------------------
  !         formats
  !----------------------------------------
1 format(3(e14.7,1x))
!2 format(f9.4,1x,f16.9)
2 format(f9.4,1x,e14.7)
3 format(i,1x,e14.7)
4 format(i,1x,i)
5 format(i)
6 format(i,1x,f11.6,1x,f11.6,1x,f11.6)
7 format(f11.6,1x,f11.6,1x,f11.6)

  !----------------------------------------
  !         reading in inputs
  !----------------------------------------
  narg=command_argument_count()

  do iarg=1, narg
    call get_command_argument(iarg,arg)
    select case (adjustl(arg))
      case ('--file')
        call get_command_argument(iarg+1,dmpFile)
        inquire(file=dmpFile,exist=dmpExist)!check if it exist
        if(.not.dmpExist)then
          print '(" File ",a,"not found")',dmpFile
          stop
        else
          print '(" File Name: ",a)',dmpFile
        end if
      case ('--nCh')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nchain
        print '(" No. chains: ",i10)',nchain
      case ('--nSeg')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nseg
        nbead=nseg+1
        print '(" No. segments: ",i10)',nseg
        ntotseg=nchain*nseg
      case ('--nSkip')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nskip
        print '(" Total number of lines to be skipped: ",i10)',nskip
      case ('--b')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') b
        print '(" b: ",f10.3)',b
      case ('--delw')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') delw
        print '(" Min distance from wall (collision) ",f10.3)',delw
      case ('--nTime')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') ntime
        print '(" No. chains: ",i10)',ntime
      case ('--nbin')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nbin
        print '(" number of bins ",i10)',nbin
      case ('--help')
        call print_help()
        stop
    end select
  end do ! iarg

  !----------------------------------------
  !         variable allocations
  !----------------------------------------
  allocate(q(3,nseg,nchain,ntime))
  allocate(collisions(nseg))
  allocate(Qeex(nchain,ntime),Qeey(nchain,ntime),Qeez(nchain,ntime))
  allocate(rx(nseg,nchain,ntime),ry(nseg,nchain,ntime),rz(nseg,nchain,ntime))
  allocate(nearwall(nseg,nchain,ntime))
  allocate(t_wall(nseg,ntime*nchain),t_away(nseg,ntime*nchain))
  allocate(t_curr(nseg,nchain))
  allocate(wall_idx(nseg),away_idx(nseg))
  allocate(midbin(nseg,nbin))
  allocate(Prob_t_wall(nseg,nbin),Prob_t_away(nseg,nbin))
  allocate(rx_mean(nseg),ry_mean(nseg),rz_mean(nseg))

  !----------------------------------------
  !         opening files
  !----------------------------------------
  open (unit=1,file=adjustl(dmpFile),status='old')
  open (unit=2,file='wallcoll_bw_pdf.dat',status='unknown')
  open (unit=3,file='wallcoll_t_wall.dat',status='unknown')
  open (unit=4,file='wallcoll_t_away.dat',status='unknown')
  open (unit=5,file='wallcoll_r_mean.dat',status='unknown')
  open (unit=6,file='wallcoll_endbead.dat',status='unknown')

  !----------------------------------------
  !  data initialization & read in data
  !----------------------------------------
  !recall that a=0, sets all elements of a to 0
  nseg_bb=nseg
  Qeex=0._dp;Qeey=0._dp;Qeez=0._dp
  rx=0._dp;ry=0._dp;rz=0._dp
  collisions = 0
  nearwall = .FALSE.
  totcollisions = 0
  t_wall = 0; t_away = 0
  wall_idx = 0; away_idx = 0
  t_curr = 1 !initialized to 1, because the first time step will happen regardless.
  Prob_t_wall = 0; Prob_t_away = 0
  rx_mean = 0; ry_mean= 0; rz_mean = 0

  !read in files
  ! do iskip=1, nskip
  !   read(1,*)
  ! end do
  do itime=1, ntime
    do ichain=1, nchain
      do iseg=1, nseg
        read(1,*) q(1:3,iseg,ichain,itime)
        if (iseg <= nseg_bb) then
          Qeex(ichain,itime)=Qeex(ichain,itime)+q(1,iseg,ichain,itime)
          Qeey(ichain,itime)=Qeey(ichain,itime)+q(2,iseg,ichain,itime)
          Qeez(ichain,itime)=Qeez(ichain,itime)+q(3,iseg,ichain,itime)
          rx(iseg,ichain,itime) = Qeex(ichain,itime)
          ry(iseg,ichain,itime) = Qeey(ichain,itime)
          rz(iseg,ichain,itime) = Qeez(ichain,itime)
        end if
      end do !segs
    end do !chains
  end do !time

  do iseg=1,nseg
    do ichain=1, nchain
      do itime=1, ntime
        rx_mean(iseg) = rx_mean(iseg) + rx(iseg,ichain,itime)
        ry_mean(iseg) = ry_mean(iseg) + ry(iseg,ichain,itime)
        rz_mean(iseg) = rz_mean(iseg) + rz(iseg,ichain,itime)

      end do
    end do
    rx_mean(iseg) = rx_mean(iseg)/(nchain*ntime)
    ry_mean(iseg) = ry_mean(iseg)/(nchain*ntime)
    rz_mean(iseg) = rz_mean(iseg)/(nchain*ntime)
  end do

  !----------------------------------------
  !    consider the first time step
  !----------------------------------------
  itime = 1
  do ichain=1, nchain
    do iseg=1, nseg

      !if it's now near the wall
      if (ry(iseg,ichain,itime)<=delw) then
        collisions(iseg) = collisions(iseg)+1
        nearwall(iseg,ichain,itime) = .TRUE.
      elseif (ry(iseg,ichain,itime)>delw) then
        nearwall(iseg,ichain,itime) = .FALSE.
      end if

    end do
  end do

  !----------------------------------------
  !  count collisions/times for times>1
  !----------------------------------------
  do itime=2, ntime
    do ichain=1, nchain
      do iseg=1, nseg

        !if it's now near the wall, and it wasn't previously
        if (ry(iseg,ichain,itime)<=delw .AND. nearwall(iseg,ichain,itime-1)==.FALSE.) then
          collisions(iseg) = collisions(iseg) + 1
          nearwall(iseg,ichain,itime) = .TRUE.

          away_idx(iseg) = away_idx(iseg) + 1
          t_away(iseg,away_idx(iseg)) = t_curr(iseg,ichain)
          t_curr(iseg,ichain) = 1

        !it's not near the wall, and it wasn't previously
        elseif (ry(iseg,ichain,itime)>delw .AND. nearwall(iseg,ichain,itime-1)==.FALSE.) then
          nearwall(iseg,ichain,itime) = .FALSE.
          t_curr(iseg,ichain) = t_curr(iseg,ichain) + 1

        !it's near the wall, but was already here the last time step
        elseif (ry(iseg,ichain,itime)<=delw .AND. nearwall(iseg,ichain,itime-1)==.TRUE.) then
          nearwall(iseg,ichain,itime) = .TRUE.
          t_curr(iseg,ichain) = t_curr(iseg,ichain) + 1

        !the bead just left the wall
        elseif (ry(iseg,ichain,itime)>delw .AND. nearwall(iseg,ichain,itime-1)==.TRUE.) then
          nearwall(iseg,ichain,itime) = .FALSE.

          wall_idx(iseg) = wall_idx(iseg) + 1
          t_wall(iseg,wall_idx(iseg)) = t_curr(iseg,ichain)
          t_curr(iseg,ichain) = 1

        end if
      end do
    end do
  end do

  !count up the total collisions
  do iseg = 1,nseg
    totcollisions = totcollisions + collisions(iseg)
  end do

  !note that wall_idx(iseg), is exactly the number of elements in t_wall(iseg,:). i.e. wall_idx(iseg)+1 and so on will be 0's.

  !----------------------------------------
  !    calculate pdfs and write output
  !----------------------------------------
  !collision pdf
  do iseg = 1,nseg
    write(2,3) iseg, (float(collisions(iseg))/totcollisions)
  end do

  !mean polymer shape
  do iseg = 1,nseg
    write(5,6) iseg, rx_mean(iseg), ry_mean(iseg), rz_mean(iseg)
  end do

  !mean polymer shape
  do itime = 1,ntime
    do ichain = 1,nchain
      write(6,7) rx(nseg,ichain,itime), ry(nseg,ichain,itime), rz(nseg,ichain,itime)
    end do
  end do

  !t_wall for end bead
  call pdf(t_wall,wall_idx,nseg,nbin,Prob_t_wall,midbin)
  do iseg =1, nseg
    write(3,5) iseg
    do ibin = 1,nbin
      write(3,2) midbin(iseg,ibin), Prob_t_wall(iseg,ibin)
    end do
  end do

  !t_away for end bead
  call pdf(t_away,away_idx,nseg,nbin,Prob_t_away,midbin)
  do iseg =1, nseg
    write(4,5) iseg
    do ibin = 1,nbin
      write(4,2) midbin(iseg,ibin), Prob_t_away(iseg,ibin)
    end do
  end do

  !----------------------------------------
  !     deallocating
  !----------------------------------------
  deallocate(q)
  deallocate(collisions)
  deallocate(Qeex,Qeey,Qeez)
  deallocate(rx,ry,rz)
  deallocate(nearwall)
  deallocate(t_wall,t_away)
  deallocate(t_curr)
  deallocate(wall_idx,away_idx)
  deallocate(midbin)
  deallocate(Prob_t_wall,Prob_t_away)
  deallocate(rx_mean,ry_mean,rz_mean)

contains

!subroutine to calculate the pdf
  subroutine pdf(times,times_idx,nseg,nbin,Prob,midbin)
    implicit none
    !arguments to the subroutine
    integer, intent(in):: times(:,:)
    integer, intent(in) :: times_idx(:)
    integer, intent(in) :: nseg, nbin
    real(dp), intent(out) :: Prob(:,:)
    real(dp),intent(inout) :: midbin(:,:)

    !declarations for use only inside this subroutine
    integer :: iseg, ibin, idx
    real(dp),allocatable :: dbin(:)
    integer, allocatable :: alpha(:),beta(:)
    real(dp), allocatable :: sum(:)

    allocate(alpha(nseg),beta(nseg))
    allocate(dbin(nseg))
    allocate(sum(nseg))

    sum = 0._dp

    do iseg = 1,nseg
      alpha(iseg) = minval(times(iseg,1:times_idx(iseg)))
      beta(iseg) = maxval(times(iseg,1:times_idx(iseg)))
      dbin(iseg) = float(beta(iseg)-alpha(iseg))/nbin
    end do

    do iseg = 1,nseg
      do ibin=1, nbin
        midbin(iseg,ibin)=alpha(iseg)+(dbin(iseg)/2)+(ibin-1)*dbin(iseg)
      end do
    end do

    do iseg = 1,nseg
      do idx = 1,times_idx(iseg)
        do ibin = 1,nbin
          if (ibin /= nbin) then
            if (times(iseg,idx)>=(alpha(iseg) + dbin(iseg)*(ibin-1)) .and. times(iseg,idx)<(alpha(iseg) + dbin(iseg)*(ibin))) then
              Prob(iseg,ibin) = Prob(iseg,ibin) + 1._dp/(dbin(iseg)*times_idx(iseg))
            end if
          else
            !if (times(iseg,idx)>=(alpha(iseg) + dbin(iseg)*(ibin-1)) .and. times(iseg,idx)<=(alpha(iseg) + dbin(iseg)*(ibin))) then
            if (times(iseg,idx)>=(alpha(iseg) + dbin(iseg)*(ibin-1)) .and. times(iseg,idx)<=  beta(iseg)) then
            !if (times(iseg,idx)>=(alpha(iseg) + dbin(iseg)*(ibin-1))) then
              Prob(iseg,ibin) = Prob(iseg,ibin) + 1._dp/(dbin(iseg)*times_idx(iseg))
            end if
          end if
        end do
      end do
    end do


    do iseg = 1,nseg
      do ibin = 1,nbin
        sum(iseg) = sum(iseg) + Prob(iseg,ibin)* dbin(iseg)
      end do
      if (abs(sum(iseg)-1)>1.e-10) then
        print *, 'Sum of probs for bead', iseg, ' is >1E-10'
      end if
    end do

    !deallocate the variables used in this subroutine
    deallocate(alpha, beta)
    deallocate(dbin)
    deallocate(sum)

  end subroutine pdf

  subroutine print_help()
    print '(a)', 'usage: [OPTIONS]'
    print '(a)', ''
    print '(a)', 'wallcollision options:'
    print '(a)', ''
    print '(a)', ' --help        print usage information and exit'
    print '(a)', ' --file        name of the file containing segmental q'
    print '(a)', ' --nCh         total No. chains'
    print '(a)', ' --nSeg        No. segments in a chain'
    print '(a)', ' --nSkip       Total No. lines to be skipped'
    print '(a)', ' --b           squared maximum of segment length'
    print '(a)', ' --delw        the distance from the wall where F_ev acts'
    print '(a)', ' --nTime       number of time steps'
    print '(a)', ' --nbin        number of bins for the pdf'
  end subroutine print_help

end program wallcollision
