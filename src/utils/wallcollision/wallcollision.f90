!-------------------------------------------------------------------------
!Tiras Y. Lin
!
!March 24, 2017
!Program to count the number of times each bead contacts the wall and
!determines PDFs of the times between contacts and times spent at the wall
!-------------------------------------------------------------------------

program wallcollision

  use :: cmn_io_mod, only: read_input

  implicit none

  !----------------------------------------
  !         variable declarations
  !----------------------------------------
  integer,parameter :: dp=selected_real_kind(p=15,r=307)

  integer :: nchain,nseg,nbead,ichain,iseg
  integer :: iskip,nskip
  integer :: itime, ntime
  integer :: icoll
  integer :: totcollisions,totcollisions_all
  integer :: nbin,ibin
  real(dp),allocatable :: midbin(:,:)
  !character(len=20) :: arg,buffer,dmpFile_R,dmpFile_C,dmpFile_ia,dmpFile_coll,dmpFile_collall
  character(len=20) :: RdmpFile,rcdmpFile,iadmpFile,wcdmpFile,wcalldmpFile
  character(len=20) :: temp
  real(dp) :: lambda
  real(dp),allocatable :: rrel(:,:,:,:)
  real(dp),allocatable :: CoM(:,:,:)
  real(dp),allocatable,dimension(:,:,:) :: rx,ry,rz
  integer,allocatable,target :: collisions(:),collisions_all(:)
  integer,allocatable,target :: collisions_chain(:,:),collisions_all_chain(:,:)
  integer,allocatable,target :: ia_times(:,:)
  integer,allocatable,target :: ia_times_idx(:)
  real(dp), allocatable :: Prob_ia_times(:,:)
  real(dp),allocatable,dimension(:) :: rx_mean,ry_mean,rz_mean
  real(dp) :: temp1,temp2,time

  !----------------------------------------
  !         formats
  !----------------------------------------
1 format(3(e14.7,1x))
2 format(f13.3,1x,e14.7)
3 format(i,1x,e14.7)
4 format(i,1x,i)
5 format(i)
6 format(i,1x,f11.6,1x,f11.6,1x,f11.6)
7 format(f11.6,1x,f11.6,1x,f11.6)

  ! !----------------------------------------
  ! !         reading in inputs
  ! !----------------------------------------
  ! narg=command_argument_count()
  !
  ! do iarg=1, narg
  !   call get_command_argument(iarg,arg)
  !   select case (adjustl(arg))
  !     case ('--fileR')
  !       call get_command_argument(iarg+1,dmpFile_R)
  !       inquire(file=dmpFile_R,exist=dmpExist)!check if it exist
  !       if(.not.dmpExist)then
  !         print '(" File ",a,"not found")',dmpFile_R
  !         stop
  !       else
  !         print '(" Relative position file Name: ",a)',dmpFile_R
  !       end if
  !     case ('--fileC')
  !       call get_command_argument(iarg+1,dmpFile_C)
  !       inquire(file=dmpFile_C,exist=dmpExist)!check if it exist
  !       if(.not.dmpExist)then
  !         print '(" File ",a,"not found")',dmpFile_C
  !         stop
  !       else
  !         print '(" CoM file Name: ",a)',dmpFile_C
  !       end if
  !     case ('--fileia')
  !       call get_command_argument(iarg+1,dmpFile_ia)
  !       inquire(file=dmpFile_ia,exist=dmpExist)!check if it exist
  !       if(.not.dmpExist)then
  !         print '(" File ",a,"not found")',dmpFile_ia
  !         stop
  !       else
  !         print '(" ia times file Name: ",a)',dmpFile_ia
  !       end if
  !     case ('--filecoll')
  !       call get_command_argument(iarg+1,dmpFile_coll)
  !       inquire(file=dmpFile_coll,exist=dmpExist)!check if it exist
  !       if(.not.dmpExist)then
  !         print '(" File ",a,"not found")',dmpFile_coll
  !         stop
  !       else
  !         print '(" wallcollisions file Name: ",a)',dmpFile_coll
  !       end if
  !     case ('--filecollall')
  !       call get_command_argument(iarg+1,dmpFile_collall)
  !       inquire(file=dmpFile_collall,exist=dmpExist)!check if it exist
  !       if(.not.dmpExist)then
  !         print '(" File ",a,"not found")',dmpFile_collall
  !         stop
  !       else
  !         print '(" wallcollisions (all) file Name: ",a)',dmpFile_collall
  !       end if
  !     case ('--nCh')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(i10)') nchain
  !       print '(" No. chains: ",i10)',nchain
  !     case ('--nSeg')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(i10)') nseg
  !       nbead=nseg+1
  !       print '(" No. segments: ",i10)',nseg
  !       ntotseg=nchain*nseg
  !     case ('--nSkip')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(i10)') nskip
  !       print '(" Total number of lines to be skipped: ",i10)',nskip
  !     case ('--b')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(f10.3)') b
  !       print '(" b: ",f10.3)',b
  !     !case ('--delw')
  !     !  call get_command_argument(iarg+1,buffer)
  !     !  read(buffer,'(f10.3)') delw
  !     !  print '(" Min distance from wall (collision) ",f10.3)',delw
  !     case ('--nTime')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(i10)') ntime
  !       print '(" No. time steps: ",i10)',ntime
  !     case ('--nbin')
  !       call get_command_argument(iarg+1,buffer)
  !       read(buffer,'(i10)') nbin
  !       print '(" number of bins ",i10)',nbin
  !     case ('--help')
  !       call print_help()
  !       stop
  !   end select
  ! end do ! iarg

  !----------------------------------------
  !         read inputs from input file
  !----------------------------------------
  print '("")'
  print '(" Reading input file...")'
  print '(" ---------------------")'

  call read_input('RdmpFile', 0,RdmpFile, def='R.flow.dat')
  call read_input('rcdmpFile',0,rcdmpFile,def='CoM.flow.dat')
  call read_input('iadmpFile', 0,iadmpFile, def='ia_time.dat')
  call read_input('wcdmpFile',0,wcdmpFile,def='w_coll.dat')
  call read_input('wcalldmpFile',0,wcalldmpFile,def='w_coll_all.dat')
  call read_input('nchain',0,nchain)
  call read_input('nseg', 0,nseg )
  call read_input('lambda', 0,lambda)
  call read_input('ntime',0,ntime)
  call read_input('nbin',0,nbin)

  print '(" ---------------------")'
  print '("")'

  !----------------------------------------
  !         variable allocations
  !----------------------------------------
  allocate(rrel(3,nseg,nchain,ntime))
  allocate(CoM(3,nchain,ntime))
  allocate(collisions(nseg))
  allocate(collisions_chain(nseg,nchain))
  allocate(collisions_all(nseg))
  allocate(collisions_all_chain(nseg,nchain))
  allocate(rx(nseg,nchain,ntime),ry(nseg,nchain,ntime),rz(nseg,nchain,ntime))
  allocate(midbin(nseg,nbin))
  allocate(Prob_ia_times(nseg,nbin))
  allocate(rx_mean(nseg),ry_mean(nseg),rz_mean(nseg))

  !----------------------------------------
  !         opening files
  !----------------------------------------
  open (unit=1,file=adjustl(RdmpFile),status='old')
  open (unit=2,file='wc_numcoll.dat',status='unknown')
  !open (unit=3,file='wallcoll_t_wall.dat',status='unknown')
  open (unit=4,file='wc_ia_time.dat',status='unknown')
  open (unit=5,file='wc_r_mean.dat',status='unknown')
  open (unit=6,file='wc_endbead.dat',status='unknown')
  open (unit=7,file=adjustl(rcdmpFile),status='old')
  open (unit=8,file='wc_allbead.dat',status='unknown')
  open (unit=9,file=adjustl(iadmpFile),status='old')
  open (unit=10,file=adjustl(wcdmpFile),status='old')
  open (unit=11,file=adjustl(wcalldmpFile),status='old')
  open (unit=12,file='wc_collfreq.dat',status='unknown')


  !----------------------------------------
  !  data initialization & read in data
  !----------------------------------------
  rx=0._dp;ry=0._dp;rz=0._dp
  collisions = 0
  collisions_chain = 0
  collisions_all = 0
  collisions_all_chain = 0
  totcollisions = 0
  totcollisions_all = 0
  Prob_ia_times = 0
  rx_mean = 0; ry_mean= 0; rz_mean = 0

  !read in files
  ! do iskip=1, nskip
  !   read(1,*)
  ! end do
  do itime=1, ntime
    do ichain=1, nchain
      read(7,*) CoM(1:3,ichain,itime)
      read(1,*) !don't need coordinate of the tethered bead
      do iseg=1, nseg
        read(1,*) rrel(1:3,iseg,ichain,itime)
        rx(iseg,ichain,itime) = CoM(1,ichain,itime)+ rrel(1,iseg,ichain,itime)
        ry(iseg,ichain,itime) = CoM(2,ichain,itime)+ rrel(2,iseg,ichain,itime)
        rz(iseg,ichain,itime) = CoM(3,ichain,itime)+ rrel(3,iseg,ichain,itime)
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

  read(10,*) !number of lines to skip
  read(10,*) !number of lines to skip
  read(10,*) !number of lines to skip
  read(11,*)
  read(11,*)
  read(11,*) temp, temp, temp, time
  print *, 'total time is', time
  do ichain=1,nchain
    do iseg=1,nseg
      read(10,*) temp1,temp2, collisions_chain(iseg,ichain)
      read(11,*) temp1,temp2, collisions_all_chain(iseg,ichain)
    end do
  end do

  do iseg=1,nseg
    do ichain=1,nchain
      collisions(iseg) = collisions(iseg) + collisions_chain(iseg,ichain)
      collisions_all(iseg) = collisions_all(iseg) + collisions_all_chain(iseg,ichain)
    end do
  end do

  !count up the total collisions
  do iseg = 1,nseg
    print *, 'Number of collisions for bead ',iseg,' is: ', collisions(iseg)
    totcollisions = totcollisions + collisions(iseg)
    print *, 'Number of collisions_all for bead ',iseg,' is: ', collisions_all(iseg)
    totcollisions_all = totcollisions_all + collisions_all(iseg)
  end do

  !allocate(ia_times_chain(nseg,nchain))
  allocate(ia_times(nseg,maxval(collisions(:))))
  allocate(ia_times_idx(nseg))
  ia_times_idx=1
  read(9,*) !number of lines to skip
  read(9,*) !number of lines to skip
  read(9,*) !number of lines to skip
  do ichain=1,nchain
    do iseg=1,nseg
      do icoll=1,collisions_chain(iseg,ichain)
        read(9,*) temp1, temp2, ia_times(iseg,ia_times_idx(iseg))
        !print *, 'iseg is: ',iseg, 'ia_time is: ',ia_times(iseg,ia_times_idx(iseg))
        ia_times_idx(iseg) = ia_times_idx(iseg) +1
      end do
    end do
  end do

  !----------------------------------------
  !    calculate pdfs and write output
  !----------------------------------------
  !collision pdf
  do iseg = 1,nseg
    write(2,3) iseg, collisions(iseg)
  end do

  !collision frequency
  do iseg = 1,nseg
    write(12,3) iseg, (float(collisions_all(iseg))/(time/lambda)/nchain)
  end do

  !mean polymer shape
  do iseg = 1,nseg
    write(5,6) iseg, rx_mean(iseg), ry_mean(iseg), rz_mean(iseg)
  end do

  !coordinates of the end bead
  do itime = 1,ntime
    do ichain = 1,nchain
      write(6,7) rx(nseg,ichain,itime), ry(nseg,ichain,itime), rz(nseg,ichain,itime)
    end do
  end do

  !coordinates of all beads
  do itime = 1,ntime
    do ichain = 1,nchain
      do iseg = 1,nseg
        write(8,7) rx(iseg,ichain,itime), ry(iseg,ichain,itime), rz(iseg,ichain,itime)
      end do
    end do
  end do

  !ia_times for all beads
  call pdf(ia_times,collisions,nseg,nbin,Prob_ia_times,midbin)
  do iseg =1, nseg
    write(4,5) iseg
    do ibin = 1,nbin
      write(4,2) midbin(iseg,ibin), Prob_ia_times(iseg,ibin)
    end do
  end do

  !----------------------------------------
  !     deallocating
  !----------------------------------------
  deallocate(rrel)
  deallocate(CoM)
  deallocate(collisions)
  deallocate(collisions_all)
  deallocate(collisions_chain)
  deallocate(collisions_all_chain)
  deallocate(rx,ry,rz)
  deallocate(midbin)
  deallocate(Prob_ia_times)
  deallocate(rx_mean,ry_mean,rz_mean)
  deallocate(ia_times)
  deallocate(ia_times_idx)

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
      if (times_idx(iseg) > 0) then
        alpha(iseg) = minval(times(iseg,1:(times_idx(iseg))))
        print *, 'alpha is', alpha(iseg)
        beta(iseg) = maxval(times(iseg,1:times_idx(iseg)))
        print *, 'beta is', beta(iseg)
        dbin(iseg) = float(beta(iseg)-alpha(iseg))/nbin
      else
        alpha(iseg) = 0
        beta(iseg) = 0
        dbin(iseg) = 0
      end if
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
            if (times(iseg,idx)>=(alpha(iseg) + dbin(iseg)*(ibin-1)) .and. times(iseg,idx)<=  beta(iseg)) then
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
        print *, 'Sum of probs for bead', iseg, ' is not 1'
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
    print '(a)', ' --fileR        name of the file containing R coordinate'
    print '(a)', ' --fileC        name of the file containing CoM coordinate'
    print '(a)', ' --fileia        name of the file containing IA times'
    print '(a)', ' --filecoll        name of the file containing collisions'
    print '(a)', ' --nCh         total No. chains'
    print '(a)', ' --nSeg        No. segments in a chain'
    print '(a)', ' --nSkip       Total No. lines to be skipped'
    print '(a)', ' --b           squared maximum of segment length'
    !print '(a)', ' --delw        the distance from the wall where F_ev acts'
    print '(a)', ' --nTime       number of time steps'
    print '(a)', ' --nbin        number of bins for the pdf'
  end subroutine print_help

end program wallcollision
