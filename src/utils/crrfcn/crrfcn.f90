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

! The Program is written for calculating the correlation function.
program crrfcn

  use :: iso_fortran_env
  use :: prcn_mod, only: wp,dot
  use :: cmn_io_mod, only: read_input

  implicit none

  integer :: nchain,nbead
  ! integer :: nseg_bb,nseg_ar,Na
  ! character(len=20) :: tplgy
  integer :: un_R,un_rc,un_Ct,stat
  integer :: ndmp,idmp,ich,ib
  integer :: t,t0,tt0,t_idx,tt0_idmp
  integer :: tend,tgap,trun,tcrr,t_mx,tt0_mx,tt0_mn
  real(wp) :: cf_tmp,Wi,dt,lambda
  character(len=20) :: cf_tp,RdmpFile,rcdmpFile

  real(wp),allocatable :: R(:,:,:,:),rc(:,:,:),rb(:,:,:,:)
  real(wp),allocatable,target :: x(:,:)
  real(wp),allocatable :: CF0(:),CFt(:),CF0t(:)
  ! real(wp),allocatable :: CF0t(:)
  integer,allocatable :: norm(:)
  real(wp) :: xbar,fctr
  character(len=99),parameter :: fmtfe2f="(f8.2,1x,e11.3,1x,2(f20.7,2x))"
  character(len=99),parameter :: fmtfe3f="(f8.2,1x,e11.3,1x,3(f14.7,2x))"

  ! Reading input parameters:

  print '("")'
  print '(" Reading input file...")'
  print '(" ---------------------")'

  call read_input('RdmpFile', 0,RdmpFile, def='R.flow.dat'  )
  call read_input('rcdmpFile',0,rcdmpFile,def='CoM.flow.dat')
  call read_input('cf-type',  0,cf_tp,    def='xx'          )
  call read_input('nchain',0,nchain)
  call read_input('nbead', 0,nbead )
  call read_input('tcrr',  0,tcrr  )
  call read_input('tgap',  0,tgap  )
  call read_input('dt',    0,dt    )
  call read_input('Wi',    0,Wi    )

  print '(" ---------------------")'
  print '("")'

  ! select case (tplgy)
  !   case ('Linear')
  !     nseg_bb=nseg
  !   case ('Comb')
  !     nseg_bb=nseg-Na*nseg_arx
  ! end select


  ! Opening the files:


  select case (cf_tp)
  case ('xx')

    open (newunit=un_Ct,file='Ct_xx.dat',status='replace')
    write(un_Ct,*) " # time, <Ree_xx(t).Ree_xx(0)>.<Ree_xx(0).Ree_xx(0)> #"
    write(un_Ct,*) " #---------------------------------------------------#" 

  end select

  ! calculate ndmp
  open (newunit=un_R,action='read',file=adjustl(RdmpFile),status='old')
  open (newunit=un_rc,action='read',file=adjustl(rcdmpFile),status='old')
  ndmp=0

ef: do
      ndmp=ndmp+1
      do ich=1, nchain
        read(un_rc,*,iostat=stat)
        if (stat == iostat_end) then
          exit ef ! end of file
        endif
        do ib=1, nbead
          read(un_R,*,iostat=stat)
          if (stat == iostat_end) then
            print '(" Error: RdmpFile is not consistent with rcdmpFile.")'
            stop
          endif
        enddo
      enddo
    enddo ef
  close(un_R)
  close(un_rc)

ndmp=ndmp-1
print *,ndmp

  ! Allocation

  allocate(rc(3,nchain,ndmp))
  allocate(R(3,nbead,nchain,ndmp))
  allocate(rb(3,nbead,nchain,ndmp))

  open (newunit=un_R,action='read',file=adjustl(RdmpFile),status='old')
  open (newunit=un_rc,action='read',file=adjustl(rcdmpFile),status='old')

  do idmp=1, ndmp

    do ich=1, nchain

      read(un_rc,*) rc(1:3,ich,idmp)

      do ib=1, nbead

        read(un_R,*) R(1:3,ib,ich,idmp)

        rb(1:3,ib,ich,idmp)=R(1:3,ib,ich,idmp)+rc(1:3,ich,idmp)
        ! print*,'rb',rb(1:3,ib,ich,idmp)

      enddo
    enddo

  end do
  close(un_R)
  close(un_rc)

  select case (cf_tp)

  case ('xx')



    allocate(x(nchain,ndmp))
    xbar=0._wp

    do idmp=1,ndmp

      do ich=1,nchain

        x(ich,idmp)=rb(1,nbead,ich,idmp)-rb(1,1,ich,idmp)
        xbar=xbar+x(ich,idmp)
        ! print*,'cf',cf_vrvlx(ich,idmp)

      enddo

    enddo

    xbar=xbar/(ndmp*nchain)

  case ('xy')


  end select

  print *,'xbar',xbar


  !!!??? ndmp
  ! we have tgap
  !tgap=ceiling(dmpFrq*lambda/dt)

  ! Total number of time steps in simulation
  trun=ndmp*tgap

  ! The number of points for representation of CorFun (Here it is practically ntime).
  ! tcor=ceiling(tend*lambda/dt)

  ! The maximum distance between two correlated values.
  t_mx=min(trun,tgap+tcrr)-tgap

  ! allocate(CF(0:ttmpmx/tgap))
  ! allocate(sdCF(0:ttmpmx/tgap))

  allocate(CF0t( -t_mx/tgap : t_mx/tgap ))
  allocate(norm( -t_mx/tgap : t_mx/tgap ))
  allocate(CF0(nchain))
  allocate(CFt(nchain))

  ! CF=0._wp;sdCF=0._wp;
  ! CF00=0._wp
  norm=0
  CF0t=0._wp

  print *,'xbar',xbar

  select case (cf_tp)

    case ('xx')

      do idmp=1, ndmp
        t0=idmp*tgap ! TAU0
        ! CF0_tmp => x(:,idmp) ! cf(TAU0)
        CF0=x(:,idmp) - xbar ! cf(TAU0)
        tt0_mx=min(trun,t0+tcrr) ! ( TAU+TAU0 )max
        tt0_mn=max(tgap,t0-tcrr) ! ( TAU+TAU0 )min
        ! do tt0=t0, tt0_mx, tgap ! TAU+TAU0
        ! print *,'tt0_mn',tt0_mn
        ! print *,'tt0_mx',tt0_mx
        ! print *,'ndmp',ndmp
        do tt0=tt0_mn, tt0_mx, tgap ! TAU+TAU0
          tt0_idmp=tt0/tgap
          t=tt0-t0
          t_idx=t/tgap
          ! CFt => x(:,tt0_idmp)
          CFt=x(:,tt0_idmp) - xbar ! cf(TAU+TAU0)
          cf_tmp=dot(CF0,CFt)
          CF0t(t_idx)=CF0t(t_idx)+cf_tmp ! Correlation function
          norm(t_idx)=norm(t_idx)+1
          ! print *,'idmp',idmp
          ! print *,'tt0_idmp',tt0_idmp
          ! print *,'t_idx',t_idx
          ! print *,'cf_tmp',cf_tmp,CF0t(t_idx),norm(t_idx)
        end do ! tt0
      end do ! idmp

      do t_idx=-t_mx/tgap, t_mx/tgap
        CF0t(t_idx)=CF0t(t_idx)/norm(t_idx)
        if (t_idx < 0) then
          fctr=CF0t(0)/norm(0)
        else
          fctr=CF0t(0)
        endif
        write(un_Ct,fmtfe2f) Wi,dt,t_idx*tgap*dt,CF0t(t_idx)/fctr
      end do

    end select

 ! do t_idx=-t_mx/tgap, t_mx/tgap
 !    select case (cf_tp)
 !    case ('xx')
 !      ! print *,'t_idx',t_idx
 !      ! print *,'norm',norm(t_idx)
 !      CF0t(t_idx)=CF0t(t_idx)/norm(t_idx)
 !      print *,'t_idx',t_idx
 !      print *,'cf_tmp2',CF0t(t_idx),norm(t_idx),CF0t(0)
 !      write(un_Ct,fmtfe2f) Wi,dt,t_idx*tgap*dt,CF0t(t_idx)/CF0t(0)
 !    end select
 !  end do

  ! Deallocations:

  deallocate(norm)
  deallocate(CF0t)
  deallocate(CF0)
  deallocate(CFt)

end program crrfcn
