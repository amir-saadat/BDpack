!Tiras Lin
!Program to count the number of times each bead contacts the wall
!March 23, 2017

program wallcollision
  
  implicit none
  
  !declare variables
  integer,parameter :: dp=selected_real_kind(p=15,r=307)
  
  integer :: narg,iarg,nchain,nseg,nbead,ntotseg,ichain,iseg
  integer :: iskip,nskip
  integer :: nseg_bb
  integer :: wall_dir
  integer :: itime, ntime
  integer :: totcollisions
  integer :: iwall_idx,iaway_idx
  character(len=20) :: arg,buffer,dmpFile
  logical :: dmpExist
  real(dp) :: b,delw
  real(dp),allocatable :: q(:,:,:,:)
  real(dp),allocatable,dimension(:,:) :: Qeex,Qeey,Qeez
  real(dp),allocatable,dimension(:,:,:) :: rx,ry,rz
  logical, allocatable,dimension(:,:,:) :: nearwall
  integer,allocatable,target :: collisions(:)
  integer,allocatable,target :: t_wall(:,:), t_away(:,:)
  integer,allocatable ::  t_wall_curr(:,:), t_away_curr(:,:)
  integer,allocatable :: wall_idx(:),away_idx(:)
  
  !formats for saving files
1 format(3(e14.7,1x))
2 format(f9.4,1x,e14.7)
3 format(i,1x,e14.7)
4 format(i,1x,i)

  !set default parameters
  wall_dir = 2 !the wall normal is pointing in the y direction
  
  !read in inputs
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
        print '(" Distance from wall to start consider collision: ",f10.3)',delw
      case ('--nTime')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') ntime
        print '(" No. chains: ",i10)',ntime
      case ('--help')
        call print_help()
        stop
    end select
  end do ! iarg
  
  !allocate variables
  allocate(q(3,nseg,nchain,ntime)) 
  allocate(collisions(nseg))
  allocate(Qeex(nchain,ntime),Qeey(nchain,ntime),Qeez(nchain,ntime))
  allocate(rx(nseg,nchain,ntime),ry(nseg,nchain,ntime),rz(nseg,nchain,ntime))
  allocate(nearwall(nseg,nchain,ntime))
  allocate(t_wall(nseg,ntime*nchain),t_away(nseg,ntime*nchain)) 
  allocate(t_wall_curr(nseg,nchain), t_away_curr(nseg,nchain))
  allocate(wall_idx(nseg),away_idx(nseg))
  
  !open the files for reading and writing
  open (unit=1,file=adjustl(dmpFile),status='old')
  open (unit=2,file='beadwall_pdf.dat',status='unknown')
  open (unit=3,file='t_wall.dat',status='unknown')
  
  !start parsing the data file and store into end to end and bead vectors
  nseg_bb=nseg
  Qeex=0._dp;Qeey=0._dp;Qeez=0._dp
  rx=0._dp;ry=0._dp;rz=0._dp
  collisions = 0
  nearwall = .FALSE.
  totcollisions = 0
  t_wall = 0; t_away = 0
  
  !do iskip=1, nskip
  !  read(1,*)
  !end do
  
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
        else
        end if
      end do !segs
    end do !chains
  end do !time
  
  !set the initial condition for nearwall
  itime = 1
  do ichain=1, nchain
    do iseg=1, nseg
      
      !if it's now near the wall
      if (ry(iseg,ichain,itime)<=delw) then
        collisions(iseg) = collisions(iseg)+1
        nearwall(iseg,ichain,itime) = .TRUE.
      end if
      
    end do 
  end do 
  
  
  !process the data and count the collisions
  do itime=2, ntime
    do ichain=1, nchain
      do iseg=1, nseg
        
        !if it's now near the wall, and it wasn't previously
        if (ry(iseg,ichain,itime)<=delw .AND. nearwall(iseg,ichain,itime-1)==.FALSE.) then
          collisions(iseg) = collisions(iseg)+1
          nearwall(iseg,ichain,itime) = .TRUE.
          
        !do nothing: it's not near the wall, and it wasn't previously
        elseif (ry(iseg,ichain,itime)>delw .AND. nearwall(iseg,ichain,itime-1)==.FALSE.) then
          nearwall(iseg,ichain,itime) = .FALSE.
        !do nothing: it's near the wall, but was already here the last time step
        elseif (ry(iseg,ichain,itime)<=delw .AND. nearwall(iseg,ichain,itime-1)==.TRUE.) then
          nearwall(iseg,ichain,itime) = .TRUE.
        !the bead just left the wall
        elseif (ry(iseg,ichain,itime)>delw .AND. nearwall(iseg,ichain,itime-1)==.TRUE.) then
          nearwall(iseg,ichain,itime) = .FALSE.

        end if
      end do 
    end do 
  end do 

  !count up the total collisions
  do iseg = 1,nseg
    totcollisions = totcollisions + collisions(iseg)
  end do
  
  !now, count up the time each bead spends near the wall, and the time each bead
  !spend away from the wall. Store these in t_wall and t_away
  t_wall_curr = 0; t_away_curr = 0
  wall_idx = 0; away_idx = 0
  
  do itime=1, ntime-1
    do ichain=1, nchain
      do iseg=1, nseg
      
        if (nearwall(iseg,ichain,itime) == .TRUE.) then
          t_wall_curr(iseg,ichain) = t_wall_curr(iseg,ichain) + 1
          
          if (nearwall(iseg,ichain,itime+1) == .FALSE.) then
            wall_idx(iseg) = wall_idx(iseg) + 1
            t_wall(iseg,wall_idx(iseg)) = t_wall_curr(iseg,ichain)
            t_wall_curr(iseg,ichain) = 0
          end if
          
        elseif (nearwall(iseg,ichain,itime) == .FALSE.) then
          t_away_curr(iseg,ichain) = t_away_curr(iseg,ichain) + 1
          
          if (nearwall(iseg,ichain,itime+1) == .TRUE.) then
            away_idx(iseg) = away_idx(iseg) + 1
            t_away(iseg,away_idx(iseg)) = t_away_curr(iseg,ichain)
            t_away_curr(iseg,ichain) = 0
          end if
          
        end if
        
      end do 
    end do 
  end do 
  
  !note that wall_idx(iseg), is exactly the number of elements in t_wall(iseg,:)
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! writing output
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !collision pdf
  do iseg = 1,nseg
    write(2,3) iseg, (float(collisions(iseg))/totcollisions)
  end do

  !t_away
  iseg = 40
  do iaway_idx = 1, away_idx(iseg)
    write(3,4) iaway_idx, t_away(iseg,iaway_idx)
  end do


  !deallocate
  deallocate(q) 
  deallocate(collisions)
  deallocate(Qeex,Qeey,Qeez)
  deallocate(rx,ry,rz)
  deallocate(nearwall)
  
contains

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
  end subroutine print_help
  
end program wallcollision