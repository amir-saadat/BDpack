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

! The Program is written for calculating the PDF of an ensemble of chains.
program pdf

  implicit none

  integer,parameter :: dp=selected_real_kind(p=15,r=307)

  integer :: narg,iarg,nchain,nseg,nbead,ntotseg,k,ichain,iseg
  ! number of bins for specifying the range of |Qseg| nad |Qee|
  integer :: nbinSeg,nbinCh,ibinCh,ibinSeg,itotseg
  character(len=20) :: arg,buffer,dmpFile,yaxis,xaxis
  logical :: dmpExist
  real(dp) :: b,Qeemax,qSegmax
  real(dp),allocatable :: q(:,:,:)
  real(dp),allocatable,dimension(:) :: QeeSq,Qeex,Qeey,Qeez
  real(dp),allocatable,dimension(:) :: qmagSegSq
  real(dp),allocatable,target :: Qee(:),QeeRel(:)
  real(dp),allocatable,target :: qmagSeg(:),qmagSegRel(:)
  real(dp),pointer :: QeeP(:),qmagSegP(:)
  ! The probability arrays
  real(dp),allocatable :: ProbSeg(:),ProbCh(:)
  ! The value of middle distances array
  real(dp),allocatable :: midSeg(:),midCh(:)
  real(dp) :: qmagSegSqAve,qmagSegAve,QeeAve,QeeSqAve,alfaCh,betaCh
  real(dp) :: alfaSeg,betaSeg,dCh,dSeg,x,Tol,sum

  ! Formats
1 format(3(e14.7,1x))
2 format(f9.4,1x,e14.7)

  ! Defaults:
  nbinSeg=50;nbinCh=10
  yaxis='pdf'
  xaxis='abs'

  ! Check if arguments are found
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
      case ('--b')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(f10.3)') b
        print '(" b: ",f10.3)',b
        Qeemax=sqrt(b)*nseg;qSegmax=sqrt(b)
      case ('--nbinCh')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nbinCh
      case ('--nbinSeg')
        call get_command_argument(iarg+1,buffer)
        read(buffer,'(i10)') nbinSeg
      case ('--yaxis')
        call get_command_argument(iarg+1,yaxis)
      case ('--xaxis')
        call get_command_argument(iarg+1,xaxis)
      case ('--help')
        call print_help()
        stop
    end select
  end do ! iarg
  ! Allocation:
  allocate(q(3,nseg,nchain))
  allocate(Qee(nchain),QeeRel(nchain),QeeSq(nchain))
  allocate(Qeex(nchain),Qeey(nchain),Qeez(nchain))
  allocate(qmagSeg(ntotseg),qmagSegRel(ntotseg),qmagSegSq(ntotseg))
  allocate(ProbSeg(nbinSeg),ProbCh(nbinCh))
  allocate(midSeg(nbinSeg),midCh(nbinCh))
        
  ! input files for the connecting vecotrs of chain segments.
  open (unit=1,file=adjustl(dmpFile),status='old')
  open (unit=2,file='PQee.dat',status='unknown')
  open (unit=3,file='PQspr.dat',status='unknown')
  ! Initial positions 
  Qee=0._dp;Qeex=0._dp;Qeey=0._dp;Qeez=0._dp
  do ichain=1, nchain
    do iseg=1, nseg
      read(1,*) q(1:3,iseg,ichain)
      Qeex(ichain)=Qeex(ichain)+q(1,iseg,ichain) 
      Qeey(ichain)=Qeey(ichain)+q(2,iseg,ichain) 
      Qeez(ichain)=Qeez(ichain)+q(3,iseg,ichain) 
    end do
  end do

  ! Magnitude of Q
  k=1
  qmagSegSqAve=0._dp;QeeSqAve=0._dp
  do ichain=1, nchain
    QeeSq(ichain)=Qeex(ichain)**2+Qeey(ichain)**2+Qeez(ichain)**2
    QeeSqAve=QeeSqAve+QeeSq(ichain)
    Qee(ichain)=sqrt(QeeSq(ichain))
    do iseg=1, nseg
      qmagSegSq(k)=q(1,iseg,ichain)**2+q(2,iseg,ichain)**2+q(3,iseg,ichain)**2
      qmagSeg(k)=sqrt(qmagSegSq(k))
      qmagSegSqAve=qmagSegSqAve+qmagSegSq(k)
      k=k+1
    end do
  end do
  QeeRel=Qee/Qeemax
  qmagSegRel=qmagSeg/qSegmax

  QeeSqAve=QeeSqAve/nchain
  QeeAve=sqrt(QeeSqAve)
  qmagSegSqAve=qmagSegSqAve/(k-1)
  qmagSegAve=sqrt(qmagSegSqAve)
  if (xaxis == 'abs') then
!    alfaCh=0._dp*Qeemax
!    betaCh=1._dp*Qeemax
!    alfaSeg=0._dp*qSegmax
!    betaSeg=1._dp*qSegmax
    alfaCh=minval(Qee)
    betaCh=maxval(Qee)
    alfaSeg=minval(qmagSeg)
    betaSeg=maxval(qmagSeg)
    QeeP => Qee
    qmagSegP => qmagSeg
  elseif (xaxis == 'norm') then
!    alfaCh=minval(Qee)/Qeemax
!    betaCh=maxval(Qee)/Qeemax
!    alfaSeg=minval(qmagSeg)/qSegmax
!    betaSeg=maxval(qmagSeg)/qSegmax
    alfaCh=0._dp*Qeemax/Qeemax
    betaCh=1._dp*Qeemax/Qeemax
    alfaSeg=0._dp*qSegmax/qSegmax
    betaSeg=1._dp*qSegmax/qSegmax
    QeeP => QeeRel
    qmagSegP => qmagSegRel
  end if

  print *
  print *,'<Qee^2>:',QeeSqAve
  print *,'<Qseg^2>^1/2:',qmagSegAve
  print *,'Alfa and Beta for chain:',alfaCh,betaCh
  print *,'Alfa and Beta for segments',alfaSeg,betaSeg
  print *,'Total No. segments:',k-1

  ! finding the probability of each of intervals between alfaCh and betaCh
  dCh=(betaCh-alfaCh)/nbinCh
  dSeg=(betaSeg-alfaSeg)/nbinSeg
  print *,'bin length for chain and segment discretization:',dCh,dSeg
  ! the middle points of intervals        
  do ibinCh=1, nbinCh
    midCh(ibinCh)=alfaCh+(dCh/2)+(ibinCh-1)*dCh
  end do 
  do ibinSeg=1, nbinSeg
    midSeg(ibinSeg)=alfaSeg+(dSeg/2)+(ibinSeg-1)*dSeg
  end do 

  ProbCh=0._dp;ProbSeg=0._dp
  x=alfaCh;k=1;Tol=real(1.e-14,kind=dp);sum=0._dp
3 if (x < betaCh-Tol) then
    do ichain=1, nchain
      if (k /= nbinCh) then
        if ((QeeP(ichain) >= x) .and. (QeeP(ichain) < (x+dCh))) then
          if (yaxis == 'pdf') then
            ProbCh(k)=ProbCh(k)+1/(dCh*nchain)
          elseif (yaxis == 'prob') then
            ProbCh(k)=ProbCh(k)+1._dp/nchain
          end if
        end if
      else
        if ((QeeP(ichain) >= x) .and. (QeeP(ichain) <= betaCh)) then
          if (yaxis == 'pdf') then
            ProbCh(k)=ProbCh(k)+1/(dCh*nchain)
          elseif (yaxis == 'prob') then
            ProbCh(k)=ProbCh(k)+1._dp/nchain
          end if
        end if      
      endif
    end do
    write(2,2) midCh(k),ProbCh(k)
    if (yaxis == 'pdf') then
      sum=sum+ProbCh(k)*dCh
    elseif (yaxis == 'prob') then
      sum=sum+ProbCh(k)
    end if
    x=x+dCh;k=k+1
    goto 3
  endif
  print *
  print *, 'No. bins in the range [Qeemin Qeemax]',k-1
  print *, 'the summation of the probabilities, should be ~1',sum

  x=alfaSeg;k=1;Tol=real(1.e-14,kind=dp);sum=0._dp
4 if (x < betaSeg-Tol) then
    do itotseg=1, ntotseg
      if (k /= nbinSeg) then
        if ((qmagSegP(itotseg) >= x) .and. (qmagSegP(itotseg) < (x+dSeg))) then
          if (yaxis == 'pdf') then
            ProbSeg(k)=ProbSeg(k)+1/(dSeg*ntotseg)
          elseif (yaxis == 'prob') then
            ProbSeg(k)=ProbSeg(k)+1._dp/ntotseg
          end if
        end if
      else
        if ((qmagSegP(itotseg) >= x) .and. (qmagSegP(itotseg) <= betaSeg)) then
          if (yaxis == 'pdf') then
            ProbSeg(k)=ProbSeg(k)+1/(dSeg*ntotseg)
          elseif (yaxis == 'prob') then
            ProbSeg(k)=ProbSeg(k)+1._dp/ntotseg
          end if
        endif      
      end if
    end do
    write(3,2) midSeg(k),ProbSeg(k)
    if (yaxis == 'pdf') then
      sum=sum+ProbSeg(k)*dSeg
    elseif (yaxis == 'prob') then
      sum=sum+ProbSeg(k)
    end if
    x=x+dSeg;k=k+1
    goto 4
  endif
  print *
  print *, 'No. bins in the range [qSegmin qSegmax]',k-1
  print *, 'the summation over the probabilities, should be ~1',sum

  deallocate(q)
  deallocate(Qee,QeeRel,QeeSq)
  deallocate(Qeex,Qeey,Qeez)
  deallocate(qmagSeg,qmagSegRel,qmagSegSq)
  deallocate(ProbSeg,ProbCh)
  deallocate(midSeg,midCh)

contains

  subroutine print_help()
    print '(a)', 'usage: bin/pdf [OPTIONS]'
    print '(a)', ''
    print '(a)', 'pdf options:'
    print '(a)', ''
    print '(a)', ' --help        print usage information and exit'
    print '(a)', ' --file        name of the file containing segmental q'
    print '(a)', ' --nchain      total No. chains'
    print '(a)', ' --nseg        No. segments in a chain'
    print '(a)', ' --b           squared maximum of segment length'
    print '(a)', ' --nbinCh      No. bins for chain pdf. def: 10'
    print '(a)', ' --nbinSeg     No. bins for segment pdf. def: 50'
    print '(a)', ' --yaxis       yaxis printing mode: prob or pdf. def: pdf'
    print '(a)', ' --xaxis       yaxis printing mode: abs or norm. def: abs'
  end subroutine print_help

end program pdf
