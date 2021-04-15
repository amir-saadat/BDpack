!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is a free software: you can redistribute it and/or modify      |
!|  it under the terms of the GNU General Public License as published by  |
!|  the Free Software Foundation, either version 3 of the License, or     |
!|  (at your option) any later version.                                   |
!|                                                                        |
!|  BDpack is distributed in the hope that it will be useful,             |
!|  but WITHOUT ANY WARRANTY; without even the implied warranty of        |
!|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
!|  GNU General Public License for more details.                          |
!|                                                                        |
!|  You should have received a copy of the GNU General Public License     |
!|  along with BDpack.  If not, see <http://www.gnu.org/licenses/>.       |
!%------------------------------------------------------------------------%
module smdlt_mod

  use :: prcn_mod

  implicit none

contains

  subroutine smdlt_bs(p,id)

    use :: box_mod
    use :: io_mod
    use :: arry_mod, only: logspace,linspace,print_vector,print_matrix
    use :: tmng_mod, only: tick,tock,et_whole,init_tmng,reportTiming,doTiming
    use :: pp_smdlt, only: init_pp,data_time_init,data_run_init,CorrFcn,del_pp,CorFun,kchk
    use :: conv_mod, only: init_conv,QtoR,QtoFseg,del_conv
    use :: flow_mod, only: FlowType
    use :: sprforce_mod, only: ForceLaw,qmx
    use :: hi_mod, only: ncols,dw_bl,hstar,strupdateKSPACE
    use :: mpi
    !include 'mpif.h'

    integer,intent(in) :: p,id
    integer :: offset,offsettot,iPe,idt,itime,ibead,ierr
    integer :: ichain,icol,jcol,kbead,kchain,kcol,iseed
    integer :: irun,tgap,idmp,dmpmx
    integer(long) :: count0
    real(wp) :: rtpassed,sqrtdt,wx,wy,wz,eps,time
    ! real(wp),allocatable :: rdn(:,:,:) ! random array
    real(wp),allocatable :: rdn(:,:) ! random array
    real(wp),parameter :: c1=14.14858378_wp,c2=1.21569221_wp
    type(box) :: MainBox
    integer :: nWi,itrst,nrun,runrst,ndmp,dmprst,nprun,ndt
    real(wp) :: Wii,Wif,tend,tss,trst,dti,dtf
    real(wp) :: wx1,wx2,wy1,wy2,wz1,wz2
    character(len=10) :: WiSpacing,dtSpacing
    real(wp),allocatable :: dt(:),Wi(:),Pe(:)
    integer,allocatable :: ntime(:)
    logical :: DumpConf
#ifdef Debuge_sequence
	write(*,*) "module:smdlt_mod:smdlt_bs"
#endif
    !----------------------------------------------
    !>>> Initialization of SDE in Semi-dilute BD:
    !----------------------------------------------

    call init()

    if (id == 0) then
      if (mod(nrun,p) /= 0) then
        print *
        print '(" Error: No. processes should be a multiple of No. Boxes.")'
        stop
      end if
      write (*,*)
        write (*,*) "%------------------------------------------------------------%"
        write (*,*) "| ***Start of BDpack program to perform Brownian dynamics*** |"
        write (*,*) "|           simulation for semidilute solution               |"
        write (*,*) "%------------------------------------------------------------%"
    endif

    ! Initialization of random number generator by all processes:
!    if (id==0) then
   iseed=657483726
!    elseif (id==1) then
!    iseed=63726
!    end if
    ! iseed=seedgen(id)
    print '(" Process Rank: ",i0," Random Number Seed: ",i0)',id,iseed
    call ranils(iseed)

    !----------------------------------------------
    !>>> Initialization of the modules:
    !----------------------------------------------
    call init_box(id,nprun)
    call init_tmng(id)
    call init_pp(id,nrun)

    ! Instantiation of the box:
    call MainBox%init(id,p,nprun,runrst)
    !----------------------------------------------
    !>>> Time integration of SDE:
    !----------------------------------------------
    if (tend==0) then
	  write(*,*) "Tend is not specified"
	  tend=1
	end if
    ! Calculating Peclet number:
    Pe(:)=Wi(:)/lambda
    ! Calculating total number of iterations:
    ntime(:)=ceiling(tend*lambda/dt(:))

    ! Allocating total and local random arrays
    ! allocate(rdn(nbeadx3,ncols,nchain))
    allocate(rdn(ntotbeadx3,ncols))

  ! Loop over Pe numbers:
  do iPe=1, nWi
    ! Loop over dts:
    do idt=1, ndt
      tgap=ceiling(tend*lambda/(ndmp*dt(idt)))

      if (id == 0) then
        print *
        print '(" Time index and value between restarts: ",i,1x,f14.7)', tgap, tgap*dt(idt)
      end if

      ! Initializing for run averaging:
      call data_run_init(id,ntime(idt),tgap,ndmp,nprun,ntotchain,ntotbeadx3)
      ! Loop over run:
      do irun=runrst+1, nprun

        if (id == 0) then
          write (*,*)
          write (*,*) "%--------------------------------------------%"
          write (*,*) "| Start of time integration in all processes |"
          write (*,*) "%--------------------------------------------%"
          write (*,'(7x,a)') 'Wi            dt           process run'
          write (*,'(7x,a)') '---------------------------------------'
          write (*,'(f14.7,1x,e10.2,1x,i6)') Wi(iPe),dt(idt),irun
          write (*,*)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if (doTiming) call tick(count0)
        if (irun == runrst+1) then
          time=itrst*dt(idt)
          idmp=dmprst
        else
          itrst=0
          time=0._wp
          idmp=0
        end if
        sqrtdt=sqrt(dt(idt))
        call data_time_init(id)
        do itime=itrst+1, ntime(idt)
          ! constructing a block of random numbers
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then


            ! do ichain=1, nchain
            !   do icol=1, ncols
            !     do ibead=1, 3*nbead
            !       rdn(ibead,icol,ichain)=ranuls()-0.5
            !     end do
            !   end do
            ! end do

            do icol=1, ncols
              do ibead=1, 3*ntotbead
                rdn(ibead,icol)=ranuls()-0.5
              end do
            end do

            ! do ichain=1, nchain
            !   do icol=1, ncols
            !     do ibead=1, nbead
            !       wx=ranuls()-0.5;wx=sqrtdt*wx*(c1*wx**2+c2)
            !       wy=ranuls()-0.5;wy=sqrtdt*wy*(c1*wy**2+c2)
            !       wz=ranuls()-0.5;wz=sqrtdt*wz*(c1*wz**2+c2)
            !       offsettot=(ichain-1)*nbeadx3+3*(ibead-1)
            !       dw_bl(offsettot+1:offsettot+3,icol)=[wx,wy,wz]
            !     end do
            !   end do
            ! end do

          end if
          ! Time passed based on time step and strain based on flow strength:
          time=time+dt(idt)
          if (FlowType /= 'Equil') then
            eps=Pe(iPe)*time
            ! The modified strain because of the periodic feature of the box:
            if (FlowType == 'PSF') then
              if (hstar /= 0._wp) call strupdateKSPACE(MainBox%size)
            elseif (FlowType == 'PEF') then
			  print*, "PEF Not implemented"
            end if ! FlowType
          end if ! FlowType
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then
            jcol=1
          end if
          if ((mod(itime,ceiling(tend*lambda/(100*dt(idt)))) == 0) .and. (id == 0)) then
            rtpassed=time/lambda

            print '(" >>> Time steps and Simulation Time(steps*dt) Passed: ", i, f10.5)',itime,time
            print '(" >>> Chain-Relaxation(s) [Simulation Time/Tou_R] Passed: ", f10.5)',rtpassed
          end if
!print*,'itime',itime,id
          ! Constructing the random vector,dW, for the whole Box:
          if ((mod(itime,ncols) == 1) .or. (ncols == 1)) then

            ! do kchain=1, nchain
            !   do kcol=1, ncols
            !     do kbead=1, nbead
            !       offset=3*(kbead-1)
            !       wx=rdn(offset+1,kcol,kchain);wx=sqrtdt*wx*(c1*wx**2+c2)
            !       wy=rdn(offset+2,kcol,kchain);wy=sqrtdt*wy*(c1*wy**2+c2)
            !       wz=rdn(offset+3,kcol,kchain);wz=sqrtdt*wz*(c1*wz**2+c2)
            !       offsettot=(kchain-1)*nbeadx3+offset
            !       dw_bl(offsettot+1:offsettot+3,kcol)=[wx,wy,wz]
            !     end do
            !   end do
            ! end do

            do kcol=1, ncols
              do kbead=1, ntotbead
                offset=3*(kbead-1)
                wx=rdn(offset+1,kcol);wx=sqrtdt*wx*(c1*wx**2+c2)
                wy=rdn(offset+2,kcol);wy=sqrtdt*wy*(c1*wy**2+c2)
                wz=rdn(offset+3,kcol);wz=sqrtdt*wz*(c1*wz**2+c2)
                dw_bl(offset+1:offset+3,kcol)=[wx,wy,wz]
              end do
            end do

          end if

! print*,'bbb',size(MainBox%Boxhi_d%P_vals)

          ! Box advancement:
          call MainBox%move(itime,ntime(idt),irun,Pe(iPe),dt(idt),jcol,id,eps,itrst)

          ! Data Processing:
          ! To be done at each lambda/?*dt iteration
          ! 1/(dt) and lambda/dt is the itime for segment and chain relaxation time.
          if ( (mod(itime,tgap) == 0) .or. (itime == ntime(idt)) ) then

            print *
            print '(" >>> Dumping restart files at time (Simulation Time(steps*dt)): ",f14.7)',time
			print '(" >>> Dumping restart files at Chain-Relaxation(s) [Simulation Time/Tou_R]: ", f10.5)',rtpassed

!          if ( ((time-time_check3) >= -1.d-10) .or. (itime == ntime(iPe,idt)) ) then
            idmp=idmp+1
            if ((irun == runrst+1) .and. (itime == ntime(idt))) then
              dmpmx=idmp
              dmprst=0
            end if

            call MainBox%write(id,p,itime,ntime(idt),irun,idmp,time,Wi(iPe),dt(idt),ndmp)

          end if ! mod(itime,..)

          jcol=jcol+1 ! col in the block of random number columns

        end do ! time loop

        if (id == 0) then
          if (doTiming) then
            et_whole=tock(count0)
            call reportTiming(nchain,nbead)
          end if
        end if
        ! resetting kchk
        if (irun /= nprun) kchk=0
      end do ! run loop

      if (DumpConf) then
        do irun=1, nprun
          call data_time_init(id)
          call MainBox%read_init(p,irun)
          time=0._wp
          do idmp=1, dmpmx
            ! reading the data:
            call MainBox%read_dmp(p,irun,idmp,ndmp)
            ! Q to R:
!            call QtoR(MainBox%Qdagger_tilde,MainBox%R_tilde)
!****!           Q to Fseg:
!****            call QtoFseg(Qdagger_tilde,Fseg_tilde,ForceLaw,qmx)

            if (idmp /= dmpmx) then
              itime=idmp*tgap
            else
              itime=ntime(idt)
            end if
            time=itime*dt(idt)

            call MainBox%calc_mf(id,irun,itime,time,Wi(iPe),Pe(iPe),dt(idt),tgap,&
                                 ntime(idt),tss,trst,nrun,nprun,tend)
          end do ! idmp
          ! resetting kchk
          if (irun /= nprun) kchk=0
        end do ! irun
        if (CorFun) call CorrFcn(dt(idt),Wi(iPe),tgap,id,p,nchain,nchain_cmb,nbead,tend,&
                                 lambda,nprun,MPI_REAL_WP)
!                        CorrFcn(dt,    Wi,    tgap,myrank,p,nchain,nchain_cmb,nbead,tend,lambda,nprun,MPI_REAL_WP)
      end if ! DumpConf
    end do ! dt loop
  end do ! Pe loop

    !----------------------------------------------
    !>>> Termination of the modules:
    !----------------------------------------------
    call del_box(id)
    call del_pp(id)
    deallocate(rdn)
    call del()

contains

  subroutine init()

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer :: j,ntokens,u1,il,stat,ios
    character(len=1024) :: line
    character(len=100) :: tokens(50)
#ifdef Debuge_sequence
	write(*,*) "module:smdlt_mod:init"
#endif
    ! default setting:
    nWi=1;Wii=0._wp;Wif=0._wp;WiSpacing='Linear'
    tend=10._wp;tss=5._wp;trst=0._wp;itrst=0
    ndt=1;dti=0.01_wp;dtf=0.01_wp;dtSpacing='Linear'
    nrun=1;runrst=0
    DumpConf=.false.;ndmp=50;dmprst=0

    open (newunit=u1,action='read',file='input.dat',status='old')
    il=1
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" smdlt_mod: Error reading line ",i0, " Process ID ",i0)',il,id
        stop
      elseif (line(1:1) == '#') then
        il=il+1
        cycle ef ! commented line
      else
        il=il+1
      end if
      call parse(line,': ',tokens,ntokens)
      if (ntokens > 0) then
        do j=1,ntokens
          select case (trim(adjustl(tokens(j))))
            case ('nWi')
              call value(tokens(j+1),nWi,ios)
            case ('Wi')
              call value(tokens(j+1),Wii,ios)
              call value(tokens(j+2),Wif,ios)
              WiSpacing=trim(adjustl(tokens(j+3)))
            case ('tend')
              call value(tokens(j+1),tend,ios)
            case ('tss')
              call value(tokens(j+1),tss,ios)
            case ('trst')
              call value(tokens(j+1),trst,ios)
              call value(tokens(j+2),itrst,ios)
            case ('ndt')
              call value(tokens(j+1),ndt,ios)
            case ('dt')
              call value(tokens(j+1),dti,ios)
              call value(tokens(j+2),dtf,ios)
              dtSpacing=trim(adjustl(tokens(j+3)))
            case ('nrun')
              call value(tokens(j+1),nrun,ios)
            case ('runrst')
              call value(tokens(j+1),runrst,ios)
            case ('DumpConf')
              if(tokens(j+1) == 'TRUE') then
                DumpConf=.true.
              elseif(tokens(j+1) == 'FALSE') then
                DumpConf=.false.
              else
                print *,'Incorrect Type for DumpConf.'
              end if
            case ('ndmp')
              call value(tokens(j+1),ndmp,ios)
            case ('dmprst')
              call value(tokens(j+1),dmprst,ios)
          end select
        enddo
      end if
    end do ef
    close(u1)

    ! Wi number setting
    allocate(Wi(nWi),Pe(nWi))
    if (WiSpacing == 'Log') then
      call logspace(Wii,Wif,Wi)
    elseif (WiSpacing == 'Linear') then
      call linspace(Wii,Wif,Wi)
    end if
    ! Pe number setting
!    Pe(:)=Wi(:)/lambda
!    if (id == 0) then
!      write(*,'(7x,a)') 'Wi       Pe       dt'
!      write(*,'(7x,a)') '--------------------'
!    end if
    ! dt setting as a function of Wi
    allocate(dt(ndt),ntime(ndt))
    if (dtSpacing == 'Log') then
      call logspace(dti,dtf,dt)
    elseif (dtSpacing == 'Linear') then
      call linspace(dti,dtf,dt)
    end if
!    ntime(:)=ceiling(tend*lambda/dt(:))
    ! This gives the number of independent runs for each process:
    nprun=int(nrun/p)

  end subroutine init

  subroutine del()
#ifdef Debuge_sequence
	write(*,*) "module:smdlt_mod:del"
#endif
    deallocate(dt,Wi,Pe,ntime)

  end subroutine

  integer  function seedG()

    real(wp) :: seedtmp
    integer :: s,i,msec,n,time_info(8)
    integer(long), allocatable :: seed(:)
    
    call date_and_time(values=time_info)
    msec=(60000*time_info(6)+1000*time_info(7)+time_info(8))  !*((myrank-83)*359) ! a random integer
    call random_seed(size=n) ! get the number of integers used for the seed
    ! This is because we want different order of random numbers in each call
    allocate(seed(n))
    do i=1,n
      seed(i)=i*msec
    end do
    call random_seed(put=seed) ! give a proper seed
    call random_number(seedtmp) ! generate a sequence of nchain pseudo
    seedG=floor(300000000*seedtmp)
    deallocate(seed)
  end function seedG
  
  function seedgen(myrank)

    integer(long) :: seedgen
    real(wp) :: seedtmp
    integer,intent(in) :: myrank
    integer :: s,i,msec,n,time_info(8)

     call date_and_time(values=time_info)
     msec=(1000*time_info(7)+time_info(8))*((myrank-83)*359) ! a random integer
     call random_seed(size=n) ! get the number of integers used for the seed
     ! This is because we want different order of random numbers in each call
     call random_seed(put=(/(i*msec,i=1,n)/)) ! give a proper seed
     call random_number(seedtmp) ! generate a sequence of nchain pseudo
     seedgen=floor(2000000000*seedtmp)

  end function seedgen

  ! Random numeber seeding (from H. C. Ottinger):
  subroutine ranils(iseed)

    integer,intent(in) :: iseed
    integer,parameter :: in=2147483563,ik=40014,iq=53668,ir=12211,ntab=32
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Initial seeds for two random number generators
	!idum=iseed+123456789
	idum = 500000000+seedG()+123456789  !MB
    idum2=idum

    ! Load the shuffle table (after 8 warm-ups)
    do 10 j=ntab+8,1,-1
       k=idum/iq
       idum=ik*(idum-k*iq)-k*ir
       if(idum < 0) idum=idum+in
          if(j <= ntab) iv(j)=idum
    10 continue
    iy=iv(1)
    return

  end subroutine ranils

  ! Uniform random number generator (from H. C. Ottinger):
  real(wp) function ranuls()

    integer,parameter :: in1=2147483563,ik1=40014,iq1=53668,ir1=12211,&
                         in2=2147483399,ik2=40692,iq2=52774,ir2=3791 ,&
                         ntab=32,inm1=in1-1,ndiv=1+inm1/ntab
    real(wp),parameter :: an=1./in1
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Linear congruential generator 1
    k=idum/iq1
    idum=ik1*(idum-k*iq1)-k*ir1
    if(idum < 0._wp) idum=idum+in1

    ! Linear congruential generator 2
    k=idum2/iq2
    idum2=ik2*(idum2-k*iq2)-k*ir2
    if(idum2 < 0._wp) idum2=idum2+in2

    !Shuffling and subtracting
    j=1+iy/ndiv
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy < 1) iy=iy+inm1
    ranuls=an*iy
    return

  end function ranuls

  ! Gaussian random number generator (from H. C. Ottinger):
  real(wp) function rangls()

    integer :: iflag
    real(wp) :: gauss2,x1,x2,xsq,aux

    save iflag,gauss2
    data iflag/0/

    if(iflag == 0) then
    10 continue

    ! pair of uniform random numbers in [-1,1]x[-1,1]
    x1=2*ranuls()-1
    x2=2*ranuls()-1

    ! if not in the unit circle, try again
    xsq=x1*x1+x2*x2
    if(xsq >= 1._wp .or. xsq == 0._wp) goto 10
      ! pair of gaussian random numbers; return one and
      ! save the other for next time
      aux=sqrt(-2*log(xsq)/xsq)
      rangls=x1*aux
      gauss2=x2*aux
      iflag=1
    else
      rangls=gauss2
      iflag=0
    endif
    return

  end function rangls

  end subroutine smdlt_bs

end module smdlt_mod
