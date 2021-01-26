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
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!
! MODULE: Timing
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, Dec 2016
!
! DESCRIPTION: 
!> Data for calculating the time spent in each procedure
!--------------------------------------------------------------------
module tmng_mod

  use :: prcn_mod

  implicit none
  save

  !> Total time
  real(wp) :: et_whole
  real(wp) :: et_CF
  real(wp) :: et_vlt
  real(wp) :: et_HI
  real(wp) :: et_DT
  real(wp) :: et_DEC
  real(wp) :: et_FFT
  real(wp) :: et_IFFT
  real(wp) :: et_SPR
  real(wp) :: et_INT
  real(wp) :: et_INF
  real(wp) :: et_PME
  real(wp) :: et_R
  real(wp) :: et_K
  real(wp) :: et_PR
  real(wp) :: et_QR
  real(wp) :: et_CM
  real(wp) :: et_EW
  real(wp) :: et_EIKR
  real(wp) :: et_DCR
  real(wp) :: et_DCK

  !> cuFFT time
  real(wp) :: et_cufft
  real(wp) :: et_mklfft
  real(wp) :: et_cusparse
  real(wp) :: et_mklsparse
  real(wp) :: et_cuinfl

  !> Number of times HI routine is called
  integer :: HIcount
  !> Number of times PME routine is called
  integer :: PMEcount
  
  !> If the timing report is intended
  logical,protected :: doTiming

contains

  subroutine tick(t)

    integer(long),intent(out) :: t
     
    call system_clock(t)

  end subroutine tick

  real(wp) function tock(t)
   
    integer(long),intent(in) :: t
    integer(long) :: now,count_rate

    call system_clock(now,count_rate)
    tock=real(now-t,kind=wp)/real(count_rate,kind=wp)
  
  end function tock

  !> Initializes the data in tmng_mod
  subroutine init_tmng(id)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env

    integer,intent(in) :: id
    character(len=1024) :: line
    character(len=100) :: tokens(50)
    integer :: i,j,ntokens,u1,il,stat

    ! Default values:
    doTiming=.false.

    open (newunit=u1,action='read',file='input.dat',status='old')
    il=1
ef: do
      read(u1,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" io_mod: Error reading line ", i0, " Process ID ", i0)', il,id
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
            case ('Timer-rep')
              if(tokens(j+1) == 'TRUE') then
                doTiming=.true.
              elseif(tokens(j+1) == 'FALSE') then
                doTiming=.false.
              else
                print '(" Inconsistent Timing-rep.")'
                stop
              end if
          end select
        end do ! j
      end if ! ntokens
    end do ef
    close(u1)

  end subroutine init_tmng

  subroutine reportTiming(nchain,nbead)

    use :: hi_mod, only: HIcalc_mode

    integer,intent(in) :: nchain,nbead
    integer :: u1
    
    open (newunit=u1,file='data/Timing.dat',status='replace',position='append')
    write(u1,*) "########## The report for timing in second #############"
    write(u1,*) "# nchain, nbead #"
    write(u1,*) "# ------------- #"
    write(u1,'(1x,2(i5,1x))') nchain, nbead
    write(u1,*) "# ------------- #"
    write(u1,'(1x,a,1x,f18.6)') 'Time spent in whole simulation:',et_whole
    write(u1,*) "# ------------- #"
    write(u1,*) 'Detail time spent in Box_move():'
    write(u1,'(1x,a,1x,f18.6)') 'Calculating verlet list:',et_vlt
    write(u1,'(1x,a,1x,f18.6)') 'Calculating conservative forces:',et_CF
    write(u1,'(1x,a,1x,f18.6)') 'Calculating hydrodynamic interaction:',et_HI
    write(u1,'(1x,a,1x,f18.6)') 'Predictor step:',et_PR
    write(u1,'(1x,a,1x,f18.6)') 'chain_move():',et_CM
    write(u1,'(1x,a,1x,f18.6)') 'Q->R:',et_QR
    write(u1,*) "# ------------- #"
    write(u1,'(1x,a,1x,f18.6,1x,a,i10)') 'Total time spent in calcHI():',et_HI,'Number of calls:',HIcount
    write(u1,'(1x,2(a,1x,f18.6,1x))') 'Calculation of diffucion tensor:',et_DT,'Per call:',et_DT/HIcount
    write(u1,'(1x,2(a,1x,f18.6,1x))') 'Calculation of Brownian noise:',et_DEC,'Per call:',et_DEC/HIcount
    write(u1,*) "# ------------- #"
    select case (HIcalc_mode)
      case ('Direct')
        write(u1,*) 'Detail time spent in EWALD:'
        write(u1,'(1x,a,1x,f18.6)') 'Total EWALD  time:',et_EW
        write(u1,'(1x,a,1x,f18.6)') 'Calculating eIKR:',et_EIKR
        write(u1,'(1x,a,1x,f18.6)') 'Real Part:',et_R
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part:',et_K
      case ('PME')
        write(u1,*) 'Detail time spent in construction of PME related vector tensors:'
        write(u1,'(1x,a,1x,f18.6)') 'Real Part:',et_DCR
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part:',et_DCK
        write(u1,*) 'Detail time spent in PME:'
        write(u1,'(1x,a,1x,f18.6)') 'Total PME time:',et_PME
        write(u1,'(1x,a,1x,f18.6)') 'Real Part:',et_R
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part:',et_K
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part; Spreading:',et_SPR
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part; FFT:',et_FFT
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part; Influence:',et_INF
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part; IFFT:',et_IFFT
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part; Interpolation:',et_INT
        write(u1,*) "# ------------- #"
        write(u1,'(1x,a,1x,i10)') 'Number of times PME called:',PMEcount
        write(u1,*) 'Detail time spent on PME per call:'
        write(u1,'(1x,a,1x,f12.6)') 'Total PME time:',et_PME/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Real Part:',et_R/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part:',et_K/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part; Spreading:',et_SPR/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part; FFT:',et_FFT/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part; Influence:',et_INF/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part; IFFT:',et_IFFT/PMEcount
        write(u1,'(1x,a,1x,f12.6)') 'Recip Part; Interpolation:',et_INT/PMEcount


        write(u1,'(1x,a,1x,f18.6)') 'Recip Part on CPU; FFT:',et_mklfft/PMEcount
        write(u1,'(1x,a,1x,f18.6)') 'Recip Part on GPU; FFT:',et_cuFFT/PMEcount


    end select

    close(u1)
    
  end subroutine reportTiming


end module tmng_mod
