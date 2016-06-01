!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is free software: you can redistribute it and/or modify        |
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

! The Program is written for generating initial configuration
program cnfgen

  implicit none

  integer,parameter :: dp=selected_real_kind(p=15,r=307)

  integer :: nchain,nseg,nbead,narg,iarg
  character(len=20) :: arg,buffer,output
  real(dp) :: hstar,a
  real(dp) :: qSegSqAve
  real(dp) :: ratio ! the ratio of segment where overlap is not allowed
  real(dp),parameter :: PI=3.1415926535897958648d0

  real(dp),allocatable,target :: q(:,:,:)
  real(dp),pointer :: qtemp2(:)
  real(dp),allocatable :: qtemp1(:,:),qtemp3(:,:)
  real(dp),allocatable :: qx(:),qy(:),qz(:)
  real(dp) :: x,qetoeAve,qetoex,qetoey,qetoez,sqqetoe,qetoe,rijmag
  integer :: iseed,ichain,jbead,icount,ind1,ind2,ind3,ibead

  ! Formats
1 format(3(e14.7,1x))

  ! Defaults:
  hstar=0.25
  qSegSqAve=3._dp
  ratio=1._dp

  ! Check if arguments are found
  narg=command_argument_count()

  do iarg=1, narg
    call get_command_argument(iarg,arg)
    select case (adjustl(arg))
      case ('--file')
        call get_command_argument(iarg+1,output)
      case ('--nCh')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nchain
      case ('--nSeg')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nseg
        nbead=nseg+1
      case ('--hs')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') hstar
      case ('--q')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') qSegSqAve
      case ('--r')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') ratio
      case ('--help')
        call print_help()
        stop
    end select
  end do ! iarg
  print '(" File Name: ",a)',output
  print '(" No. chains: ",i10)',nchain
  print '(" No. segments: ",i10)',nseg
  print '(" h*: ",f10.3)',hstar
  print '(" <q_seg^2>: ",f10.3)',qSegSqAve
  print '(" <No overlap ratio>: ",f10.3)',ratio

  a=sqrt(PI)*hstar

  ! Allocation:
  allocate(q(3,nseg,nchain))
  allocate(qtemp1(3,nseg),qtemp3(3,nseg))
  allocate(qx(nseg),qy(nseg),qz(nseg))

  ! output files for the connecting vecotrs of chain segments.
  open (unit=1,file=output,status='unknown')

  ! initialization of random number generator
  iseed=657483726
  call ranils(iseed)
  call random_seed(iseed)

  do ichain=1, nchain
    qtemp1=0.0d0;qtemp3=0.0d0;jbead=2;icount=0
    do while (jbead <= nbead)
      qtemp2 => q(1:3,jbead-1,ichain)
2     call random_number(x)
      ind1=1+mod(int(x*1.d4),3)
      ind2=1+mod((int(x*1.d4)+1),3)
      ind3=1+mod((int(x*1.d4)+2),3)
      qtemp2(ind1)=sqrt(qSegSqAve/3)*(2*ranuls()-1)
      qtemp2(ind2)=sqrt(qSegSqAve/3)*(2*ranuls()-1)
      qtemp2(ind3)=2*ranuls()-1
      
      ! As we want |qSeg|=<qSeg>equil. For qtemp2(3), the random sign of rangls() is used:
      qtemp2(ind3)=qtemp2(ind3)/abs(qtemp2(ind3))*sqrt(qSegSqAve-qtemp2(ind1)**2-&
                   qtemp2(ind2)**2)
      ! qtemp3(:,ibead) is basically r(ibead->jbead), but as the move may be rejected, we 
      ! differentiate it with qtemp1(:,ibead), whose components are the same as qtemp3 but
      ! fixed.                                                                            
      do ibead=1, jbead-1
        qtemp3(1:3,ibead)=qtemp3(1:3,ibead)+qtemp2(1:3)
        rijmag=sqrt(qtemp3(1,ibead)**2+qtemp3(2,ibead)**2+qtemp3(3,ibead)**2)
        if (rijmag/(2*a) <= ratio) then
          qtemp3=qtemp1
          icount=icount+1
          if (icount >= 1000) then
            print '(" icount exceeded its maximum allowed value:",i5)',icount
            stop
          end if
          goto 2
        end if
      end do ! ibead
      qtemp1=qtemp3
      write(1,1) qtemp2(1:3)
      jbead=jbead+1
    end do ! while
  end do ! ichain

  qetoeAve=0._dp
  do ichain=1, nchain      
    qx=q(1,:,ichain);qy=q(2,:,ichain);qz=q(3,:,ichain) 
    qetoex=sum(qx);qetoey=sum(qy);qetoez=sum(qz)
    sqqetoe=qetoex**2+qetoey**2+qetoez**2
    qetoe=sqrt(sqqetoe)
    qetoeAve=qetoeAve+qetoe
  end do
  qetoeAve=qetoeAve/nchain
  print '(" <Qee>: ",f10.3)',qetoeAve

contains

  ! Random numeber seeding (from H. C. Ottinger):
  subroutine ranils(iseed)
  
    integer,intent(in) :: iseed
    integer,parameter :: in=2147483563,ik=40014,iq=53668,ir=12211,ntab=32
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Initial seeds for two random number generators
    idum=iseed+123456789
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
  real(dp) function ranuls()
  
    integer,parameter :: in1=2147483563,ik1=40014,iq1=53668,ir1=12211,&
                         in2=2147483399,ik2=40692,iq2=52774,ir2=3791 ,&
                         ntab=32,inm1=in1-1,ndiv=1+inm1/ntab
    real(dp),parameter :: an=1./in1
    integer :: iv(ntab),idum,idum2,iy
    integer :: k,j

    common /ranbls/ idum,idum2,iy,iv

    ! Linear congruential generator 1
    k=idum/iq1
    idum=ik1*(idum-k*iq1)-k*ir1
    if(idum < 0._dp) idum=idum+in1

    ! Linear congruential generator 2
    k=idum2/iq2
    idum2=ik2*(idum2-k*iq2)-k*ir2
    if(idum2 < 0._dp) idum2=idum2+in2

    !Shuffling and subtracting
    j=1+iy/ndiv
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy < 1) iy=iy+inm1
    ranuls=an*iy
    return
   
  end function ranuls  

  subroutine print_help()
    print '(a)', 'usage: bin/cnfgen [OPTIONS]'
    print '(a)', ''
    print '(a)', 'cnfgen options:'
    print '(a)', ''
    print '(a)', ' --help        print usage information and exit'
    print '(a)', ' --file        name of the file for output'
    print '(a)', ' --nchain      total No. chains'
    print '(a)', ' --nseg        No. segments in a chain'
    print '(a)', ' --hs          the hydrodynamic diameter of a bead; def: 0.25'
    print '(a)', ' --r           the segment fraction length w/o overlap. def: 1'
  end subroutine print_help

end program cnfgen

