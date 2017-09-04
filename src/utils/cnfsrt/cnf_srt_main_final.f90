
!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2017:                                            |
!|  Fluid Mechanics Lab (Shaqfeh Group)                                   |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@stanford.edu>                       |
!|  Author:    Tiras Y. Lin  <tiraslin@stanford.edu>                      |
!|  Advisor:   Eric S. G. Shaqfeh <esgs@stanford.edu>                     |
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

! The Program is written for sorting different configurations
! originally written by Vydia Venkataramani, modified by Amir Saadat
program cnfsrt

  use :: iso_fortran_env
  use :: prcn_mod, only: wp,dot
  use :: cmn_io_mod, only: read_input

  implicit none


  real(wp),allocatable,target :: Rx(:,:,:),Ry(:,:,:),Rz(:,:,:)
  real(wp),pointer :: RPx(:),RPy(:),RPz(:)
  real(wp),allocatable :: x(:),y(:),z(:),r(:)
  integer,allocatable :: ib(:),ib1(:),cnf_tp(:,:),icount3(:)
  real(wp) :: iconfig(6), pconfig(6)
  real(wp) :: rmax,xi,yi,zi,xj,yj,zj,dist,xmin,ymin,zmin
  real(wp) :: xmax,ymax,zmax,xr,yr,zr,xcap,ycap,zcap,rmin
  real(wp) :: xproj,yproj,zproj,xpres,xnext,dia,radius
  real(wp) :: bead_left,bead_right,avg,avgsum
  integer :: i,ntype,maxlength,ichain,ibead,jbead,ilen
  integer :: min_idx,max_idx,maxbright,maxbright1,low,last,info,info1
  integer :: k,ij,j,ibcount,sumbright,lcount,left,markerc,ik
  integer :: lend,marker1,marker2,markerl,markerr,icount
  integer :: u_R,u_tp,u_tp_int,u_tp_ave,u_traj_x,u_traj_r,u_tr_prob
  integer :: iter,ndmp,tplgy,stat,dmp_f,dmp_l,idmp

  !------------------------ input
  integer :: nchain,nseg,nseg_bb,nbead,nbead_bb,itraj
  real(wp) :: bin_sz ! number of pixels in a bin, used for averaging brightness
  real(wp) :: npixel ! number of pixels in unit length
  real(wp) :: qmax
  real(wp) :: dmp_fr_f,dmp_fr_l ! min and max fraction of total number of dumps
                                ! to use for evaluating the average configuration
  character(len=20) :: RdmpFile,tplgy_inp
  logical :: bead_overlap ! whether brightness due to bead overlap is required
  !------------------------
  integer,parameter :: FOLD=1,HALF_DUMB=2,KINK=3,DUMB=4,COIL=5,EXT=6
  integer,parameter :: LINEAR=1,COMB=2
  real(wp) :: a,b,c
  !------------------------


  ! Reading input parameters:

  call read_input('RdmpFile', 0, RdmpFile, def='R.flow.dat'  )
  call read_input('Topology', 0, tplgy_inp, def='Linear'  )
  call read_input('nchain',0, nchain)
  call read_input('nseg', 0, nseg )
  call read_input('nseg_bb', 0, nseg_bb, def=nseg )
  call read_input('qmax', 0, qmax )
  call read_input('traj_idx', 0, itraj )
  call read_input('npixel', 0, npixel, def=1._wp )
  call read_input('bin_sz', 0, bin_sz, def=5._wp )
  call read_input('min_dmp_fr', 0, dmp_fr_f, def=0._wp )
  call read_input('max_dmp_fr', 0, dmp_fr_l, def=0.75_wp )
  call read_input('bead_overlap', 0, bead_overlap, def=.false.)

  print '("")'
  print '(" Reading input file...")'
  print '(" ---------------------")'

  print '(" Dump File = ",a)',RdmpFile
  print '(" Topology = ",a)',tplgy_inp
  print '(" nchain = ",i)',nchain
  print '(" nseg = ",i)',nseg
  print '(" nseg_bb = ",i)',nseg_bb    
  print '(" qmax = ",f15.5)',qmax
  print '(" itraj = ",i)',itraj
  print '(" npixel = ",f15.5)',npixel
  print '(" bin_sz = ",f15.5)',bin_sz
  print '(" min dmp_fr = ",f15.5)',dmp_fr_f
  print '(" max dmp_fr = ",f15.5)',dmp_fr_l
  print '(" bead overlap = ",L5)',bead_overlap

  print '(" ---------------------")'
  print '("")'

  nbead=nseg+1
  nbead_bb=nseg_bb+1
  ! Allocations
  allocate(x(nbead_bb),y(nbead_bb),z(nbead_bb),r(nbead_bb),cnf_tp(nchain,6))
  maxlength=int(qmax*nseg_bb*npixel)+1
  print '(" maximum length (should be >1) = ",i5)',maxlength
  allocate(ib(maxlength),ib1(maxlength),icount3(maxlength))

  ! initialization
  select case (tplgy_inp)
  case('Linear')
    tplgy=LINEAR
  case('Comb')
    tplgy=COMB
  endselect

  do i = 1, 6
    iconfig(i) = 0
    pconfig(i) = 0._wp
  enddo

  ! calculate ndmp

  do iter= 1, 2
    open (newunit=u_R,action='read',file=adjustl(RdmpFile),status='old')
    ndmp=1
ef: do
      do ichain=1, nchain
        do ibead=1, nbead
          if (iter == 1) then            
            read(u_R,*,iostat=stat) a,b,c
            ! print*,a,b,c
          else
            if (ibead <= nbead_bb) then              
              read(u_R,*,iostat=stat) Rx(ibead,ichain,ndmp),Ry(ibead,ichain,ndmp),Rz(ibead,ichain,ndmp)
            else
              read(u_R,*,iostat=stat)
            endif
          endif
          if (stat == iostat_end) then
            ! if (ichain < nchain .or. ibead < nbead) then
            !   print '(" Error reading file, ichain, ibead: ",2(i,1x))',ichain,ibead
            !   print*,'iter',iter,ndmp
            !   stop
            ! endif
            if (iter == 1) then
              allocate(Rx(nbead_bb,nchain,ndmp))
              allocate(Ry(nbead_bb,nchain,ndmp))
              allocate(Rz(nbead_bb,nchain,ndmp))
            endif
            exit ef ! end of file
          endif
        enddo
      enddo
      ndmp=ndmp+1
    enddo ef
    close(u_R)
  enddo
  ndmp=ndmp-1

  ! opening files
  open (newunit=u_tp,file='Type_char.dat',status='replace')
  open (newunit=u_tp_int,file='Type_int.dat',status='replace')
  open (newunit=u_tp_ave,file='Type_ave.dat',status='replace')
  open (newunit=u_traj_x,file='Traj_X.dat',status='replace')
  open (newunit=u_traj_r,file='Traj_r.dat',status='replace')
  open (newunit=u_tr_prob,file='trProb.dat',status='replace')

  ! initialization
  cnf_tp=0
  dmp_f=min(nint(dmp_fr_f*ndmp),ndmp)
  dmp_l=min(nint(dmp_fr_l*ndmp),ndmp)
  write (u_tp_ave,'(" #----- Most probable configuarion -----#")')


  do idmp=1, ndmp

    ! initialization
    iconfig=0
    ! write (u_tp,'(" #----- DUMP STEP: ",i," ------#")') idmp
    ! write (u_tp_int,'(" #----- DUMP STEP: ",i," ------#")') idmp
    ! write (u_traj,'(" #----- DUMP STEP: ",i," ------#")') idmp
    ! write (u_tr_prob,'(" #----- DUMP STEP: ",i," ------#")') idmp

    do ichain=1, nchain

      RPx => Rx(:,ichain,idmp)
      RPy => Ry(:,ichain,idmp)
      RPz => Rz(:,ichain,idmp)

      ! initialization
      ntype=0
      ib=0
      ib1=0

      ! Calculating the largest distance between the beads
      rmax=0._wp
      do ibead=1, nbead_bb-1
        xi=RPx(ibead)
        yi=RPy(ibead)
        zi=RPz(ibead)
        do jbead=1, nbead_bb
          xj=RPx(jbead)
          yj=RPy(jbead)
          zj=RPz(jbead)
          dist=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

          if (dist > rmax) then
            min_idx=ibead
            max_idx=jbead
            rmax=dist
          end if
        end do ! jbead
      end do ! ibead

      xmin=RPx(min_idx)
      ymin=RPy(min_idx)
      zmin=RPz(min_idx)

      xmax=RPx(max_idx)
      ymax=RPy(max_idx)
      zmax=RPz(max_idx)

      do ibead=1, nbead_bb
        x(ibead)=RPx(ibead)-xmin
        y(ibead)=RPy(ibead)-ymin
        z(ibead)=RPz(ibead)-zmin
        if (ichain.eq.itraj) then
          write (u_traj_x,*) x(ibead),y(ibead),z(ibead)
        endif
      end do
      !call print_vector(x,'x2')
      !call print_vector(y,'y2')
      !call print_vector(z,'z2')

      !**************************************************************
      ! Calculating the unit vector along the longest length
      !**************************************************************
      xr=xmax-xmin
      yr=ymax-ymin
      zr=zmax-zmin

      xcap=xr/rmax
      ycap=yr/rmax
      zcap=zr/rmax

      !**************************************************************
      ! Finding the projection of each bead onto the longest length
      ! essentially the problem is one dimensional after this point
      !**************************************************************
      rmin=0._wp
      do ibead=1, nbead_bb
        r(ibead)=0._wp
        xproj=x(ibead)*xcap
        yproj=y(ibead)*ycap
        zproj=z(ibead)*zcap
        r(ibead)=(xproj+yproj+zproj)*npixel
        if (r(ibead) < rmin) then
          rmin=r(ibead)
        end if
      end do

      !**************************************************************
      ! Repositioning the beads such that the minimum r = 0
      !**************************************************************
      do ibead=1, nbead_bb
        r(ibead)=r(ibead)-rmin
      end do

      !**************************************************************
      ! Accounting for bond overlaps
      !**************************************************************
      maxbright=int(rmax*npixel)+1
      ! print '(" maximum brightness (should be >1) = ",i)',maxbright
      !print *,'rmax',rmax
      !print *,'maxbright',maxbright
      if (r(1) > r(2)) then     
        low=int(r(2))+1
        last=int(r(1))+1
        info=0
      else
        low=int(r(1))+1 
        last=int(r(2))+1
        info=1
      end if

      ! from low x to high x (to enter the loop)
      do ilen=low, last
        ib(ilen)=ib(ilen)+1
      end do

      do ibead=2, nbead_bb-1
        xpres=r(ibead)
        xnext=r(ibead+1)
        if (xpres < xnext) then
          info1 = 0
        else
          info1 = 1
        end if

        ! 1 is added to prevent double counting if there is no overlap
        if ((info == 0).and.(info1 == 0)) then
          low=int(xpres)+2
          last=int(xnext)+1
        end if

        if ((info == 0).and.(info1 == 1)) then
          low=int(xnext)+1
          last=int(xpres)+1
        end if

        ! 1 is subtracted from last to prevent double counting
        if ((info == 1).and.(info1 == 1)) then
          low=int(xnext)+1
          last=int(xpres)
        end if

        if ((info == 1).and.(info1 == 0)) then
          low=int(xpres)+1
          last=int(xnext)+1
        end if

        do ilen=low, last
          ib(ilen)=ib(ilen)+1
        end do
        info=info1
      end do

      !**************************************************************
      ! Accounting for bead overlap
      !**************************************************************

      if (bead_overlap) then
        ! the diameter is the equilibrium length of the spring
        ! maybe maximum extension would make more sense
        ! not sure, should be checked later
        dia=1._wp*npixel ! parameter set based on your assumption
        radius=dia/2

        do i=1, maxbright
          icount3(i)=0._wp
        end do

        ! if the spring are collapsed compared their equilibrium state
        ! then add brightness index
        do ibead=1, nbead_bb
          bead_left=r(ibead)-radius
          bead_right=r(ibead)+radius
          if (bead_left < 0._wp) then
            bead_left = 0._wp
          end if

          if (bead_right > real(maxbright,kind=wp)) then
            bead_right=real(maxbright,kind=wp)
          end if
          low=int(bead_left)+1
          last=int(bead_right)+1
          do ilen=low, last
            ib(ilen)=ib(ilen)+1
            icount3(ilen)=1
          end do
        end do

        ! subtract 1 which is added for all beads
        do ilen=1, maxbright
          if (icount3(ilen) == 1) then
            ib(ilen)=ib(ilen)-1
          end if
        end do
      endif

      !**************************************************************
      ! Averaging over configurations
      !**************************************************************
      !call print_vector(ib(1:maxbright),'ib1')

      ! every bin_sz units are binned together
      if (mod(maxbright,int(bin_sz)) == 0) then
        maxbright1=maxbright/bin_sz
      else
        maxbright1=int(maxbright/bin_sz)+1
      end if

      k=0
      ij=1
      do i=1, maxbright1        
        j=1 ! for each bin start from 1
        ibcount=0
        do while ((j <= bin_sz).and.(ij <= maxbright))
          ibcount=ibcount+ib(ij)
          j=j+1
          ij=ij+1
        end do
        ! calculate the average of ib at each bin
        if (j == bin_sz) then
          avg=real(ibcount,kind=wp)/bin_sz
        else
          avg=real(ibcount,kind=wp)/(j-1)
        end if

        ! round it up if avg greater than half
        if ((real(int(avg),kind=wp)+0.5_wp) < avg) then
          ib1(i)=int(avg)+1
        else
          ib1(i)=int(avg)
        end if
        k=k+1
      end do

      ! replace ib data with ib1 (averaged over bins)
      do i=1, maxbright
        ib(i)=0
      end do
      maxbright=maxbright1
      do i=1, maxbright
        ib(i)=ib1(i)
      end do

      !call print_vector(ib(1:maxbright),'ib2')
      !**************************************************************
      ! Checking for the minimum resolution
      !**************************************************************

      ! it is an obvious coil if ...
      if (rmax <= 0.05_wp*qmax*nseg) then
        iconfig(5)=iconfig(5)+1
        ntype=5
      else

        !**************************************************************
        ! **Assigning for 2/3, 2/5 folds, coil and half dumbell**
        !**************************************************************

        ! we have brightness of 1 in one side and greater than 1 on the other
        if (((ib(1) == 1).and.(ib(maxbright) > 1)).or. &
          ((ib(maxbright) == 1).and.(ib(1) > 1))) then

          if ((ib(1) == 1).and.(ib(maxbright) > 1)) then
            j=1
            i=maxbright
          end if
          if ((ib(maxbright) == 1).and.(ib(1) > 1)) then
            j=0
            i=1
          end if

          icount = 0 ! total number of pixels that have ib>1
          sumbright = 0 ! sum of the brightness of all pixels that have ib>1
          do while (ib(i) > 1)
            icount=icount+1
            sumbright=sumbright+ib(i)
            if (j == 1) then ! based on the condition defined above
              i=i-1
            else
              i=i+1
            end if  
          end do
          avgsum=sumbright/(1._wp*icount)
          lcount=0 ! total number of pixels with ib=1
          if (j == 1) then
            k=1
          else
            k=maxbright
          end if
          do while (ib(k) == 1)
            lcount = lcount + 1
            if (j == 1) then
              k=k+1
            else
              k=k-1
            end if
          end do
          if (j == 1) then
            left=k
            last=i
          else
            left=i
            last=k 
          end if

          markerc=0
          do ik=left, last
            if (ib(ik) > (int(avgsum))) then
              markerc=1
              lend=ik
            end if
          end do

          if (markerc == 1) then ! kink
            iconfig(3)=iconfig(3)+1
            ntype=3
          else
            icount=maxbright-lcount
            if ((0.25_wp <= ((1._wp*icount)/maxbright))) then ! fold
              iconfig(1)=iconfig(1)+1
              ntype=1
            end if
            if (((1._wp*icount)/maxbright) < 0.25_wp) then ! half dumbell
              iconfig(2)=iconfig(2)+1
              ntype=2
            end if
          end if
        end if

        !**************************************************************
        ! **Assigning for kink or fully extended**
        !**************************************************************

        if ((ib(1) == 1).and.(ib(maxbright) == 1)) then
          marker1=0
          do i=2, (maxbright-1)
            if (ib(i) == 1) then
              marker1=marker1+1
            end if
          end do
          if (marker1 == (maxbright-2)) then ! extended
            iconfig(6)=iconfig(6)+1
            ntype=6
          else ! kink
            iconfig(3)=iconfig(3)+1
            ntype=3
          end if
        end if

        !**************************************************************
        ! **Assigning for coil or dumbell**
        !**************************************************************

        if ((ib(1) > 1).and.(ib(maxbright) > 1)) then
          marker2=0
          markerl=2
          do while ((ib(markerl) > 1).and.(markerl < maxbright))
            markerl=markerl+1
          end do
          markerr=maxbright-1
          do while ((ib(markerr) > 1).and.(markerr > 1))
            markerr=markerr-1
          end do
          if ((markerl == (maxbright)).and.(markerr == 1)) then ! coil
            iconfig(5)=iconfig(5)+1
            ntype=5
          else ! 
            do i=markerl, markerr
              if (ib(i) > 1) then
                marker2=1
              end if
            end do
            if (marker2 /= 1) then ! dumbell
              iconfig(4)=iconfig(4)+1
              ntype=4
            else ! coil
              iconfig(5)=iconfig(5)+1
              ntype=5
            end if
          end if           
        end if

      end if

      if (idmp >= dmp_fr_f .and. idmp <= dmp_fr_l) then
        cnf_tp(ichain,ntype)=cnf_tp(ichain,ntype)+1
      endif

      select case (ntype)
      case(FOLD)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' fold'
        write (u_tp_int,*) int(ntype)
      case(HALF_DUMB)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' half dumbbell'
        write (u_tp_int,*) int(ntype)
      case(KINK)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' kink'
        write (u_tp_int,*) int(ntype)
      case(DUMB)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' dumbbell'
        write (u_tp_int,*) int(ntype)
      case(COIL)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' coil'
        write (u_tp_int,*) int(ntype)
      case(EXT)
        write (u_tp,'(" chain: ",i," is ",a)') ichain, ' fully extended'
        write (u_tp_int,*) int(ntype)
      end select

      if (ichain.eq.itraj) then
        do ibead = 1, nbead_bb
          write (u_traj_r,*) r(ibead)
        end do
      endif


    end do ! ichain

    do i = 1, 6
      pconfig(i) = iconfig(i)/real(nchain,kind=wp)
      write(u_tr_prob,'(1x,6(f10.5,1x))'),pconfig(i)
    enddo

  enddo ! idmp

  print*,'cnf',cnf_tp
  write(u_tp_ave,'(i)'),maxloc(cnf_tp,2)

  close(u_tp)
  close(u_tp_int)
  close(u_tp_ave)
  close(u_traj_x)
  close(u_traj_r)
  close(u_tr_prob)
  deallocate(x,y,z,r)
  deallocate(Rx,Ry,Rz)
  deallocate(cnf_tp)
  deallocate(ib,ib1)


end program cnfsrt
