!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016: University of Tennessee-Knoxville          |
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
module pp_mod

  use :: prcn_mod
  
  implicit none
  
  save
  
  real(wp),private :: tauxx,tauxy,tauyy,tauzz,tauxxyy,tauyyzz
  real(wp),private :: sdxx,sdxy,sdyy,sdzz,sdxxyy,sdyyzz
  real(wp),private :: qetoeAve,sdqetoeAve,sqqetoeAve,sdsqqetoeAve
  real(wp),private :: sqqsprAve,sdsqqsprAve

  real(wp),private :: tAvesqqsprAveTot,tAvesqqetoeAveTot,tAveqetoeAveTot
  real(wp),private :: tAvetauxyTot,tAvetauxxyyTot,tAvetauyyzzTot
  real(wp),private :: tAvesdsqqsprAveTot,tAvesdsqqetoeAveTot,tAvesdqetoeAveTot
  real(wp),private :: tAvesdxyTot,tAvesdxxyyTot,tAvesdyyzzTot
  real(wp),private :: sdtAvesqqsprAveTot
  real(wp),private :: sdtAvesqqetoeAveTot,sdtAveqetoeAveTot
  real(wp),private :: sdtAvetauxyTot,sdtAvetauxxyyTot,sdtAvetauyyzzTot
  real(wp),allocatable,private :: qee_ar(:),sdqee_ar(:)
  real(wp),allocatable,private :: sqqee_ar(:),sdsqqee_ar(:)
  real(wp),allocatable,private :: qee_art(:),sdqee_art(:)
  real(wp),allocatable,private :: sqqee_art(:),sdsqqee_art(:)
  real(wp),allocatable,private :: tAvqee_art(:),tAvsqqee_art(:)
  real(wp),allocatable,private :: tAvsdqee_art(:),tAvsdsqqee_art(:)
  real(wp),allocatable,private :: sdtAvqee_art(:),sdtAvsqqee_art(:)

  real(wp),private :: tAveDcmAveTot,tAvesdDcmAveTot,sdtAveDcmAveTot
  real(wp),private :: tAveDchrAveTot,tAvesdDchrAveTot,sdtAveDchrAveTot
  real(wp),private :: tAveRgSqAveTot,tAvesdRgSqAveTot,sdtAveRgSqAveTot
  real(wp),private :: tAveAspherAveTot,tAvesdAspherAveTot,sdtAveAspherAveTot
  real(wp),private :: tAvecosTh,tAvesdcosTh,sdtAvecosTh
  real(wp),private :: tAvemAveTot,tAvesdmAveTot,sdtAvemAveTot
  real(wp),private :: tAveLAveTot,tAvesdLAveTot,sdtAveLAveTot

  integer,private :: jcount,kcount
  integer :: ue

contains

  subroutine pp_init()

    use :: inp_mod, only: StrCalc,initmode,tplgy,Na,CoM,CoHR,RgCalc,cosThCalc,&
                          AveIterCalc,DecompMeth,iflow

    integer :: iarm
    character(len=1024) :: fnme1,fnme2

    ue=33 ! The largest unit excluding arms of the comb polymer
    if (StrCalc) then
      if (initmode == 'rst') then
        open (unit=31,file='data/trQsprSq.dat',status='unknown',position='append')
        open (unit=33,file='data/trQeeSq.dat',status='unknown',position='append')
        open (unit=7,file='data/trQeeRel.dat',status='unknown',position='append')
      else
        open (unit=31,file='data/trQsprSq.dat',status='replace',position='append')
        open (unit=33,file='data/trQeeSq.dat',status='replace',position='append')
        open (unit=7,file='data/trQeeRel.dat',status='replace',position='append')
      end if
      write(31,*) "# Wi, dt, Time, <Qspr.Qspr>, sd<Qspr.Qspr> #"
      write(31,*) "# ---------------------------------------- #"
      write(33,*) "# Wi, dt, Time, <Qee.Qee>, sd<Qee.Qee> #"
      write(33,*) "# ------------------------------------ #"
      write(7,*) "# Wi, dt, time, <|Qee|>/|Qee|_max, sd<Qee>/|Qee|_max #" 
      write(7,*) "# -------------------------------------------------- #"
      if (tplgy == 'Comb') then
        do iarm=1, Na
          write(fnme1,"(A,i0.2,'.dat')") 'data/trQeeSqArm',iarm
          write(fnme2,"(A,i0.2,'.dat')") 'data/trQeeRelArm',iarm
          if (initmode == 'rst') then
            open (ue+iarm,file=trim(adjustl(fnme1)),&
                  status='unknown',position='append')
            open (ue+Na+iarm,file=trim(adjustl(fnme2)),&
                  status='unknown',position='append')
          else
            open (ue+iarm,file=trim(adjustl(fnme1)),&
                  status='replace',position='append')
            open (ue+Na+iarm,file=trim(adjustl(fnme2)),&
                  status='replace',position='append')
          end if
          write(ue+iarm,*) "# Wi, dt, Time, <Qee.Qee>, sd<Qee.Qee> #"
          write(ue+iarm,*) "# ------------------------------------ #"
          write(ue+Na+iarm,*) "# Wi, dt, time, <|Qee|>/|Qee|_max, sd<Qee>/|Qee|_max #" 
          write(ue+Na+iarm,*) "# -------------------------------------------------- #"
        end do
      end if
      if (initmode == 'rst') then
        open (unit=8,file='data/EtavsTime.dat',status='unknown',position='append')
      else
        open (unit=8,file='data/EtavsTime.dat',status='replace',position='append')
      end if 
      write(8,*) "# Wi, dt, Time, Etap=<Taupxy>/Pe, sd<Taupxy>/Pe #"
      write(8,*) "# --------------------------------------------- #"
      open (unit=30,file='data/QeeSq.dat',status='unknown',position='append')
      write(30,*) "# Wi, dt, <<Qee^2>t, <sd<Qee^2>>t, sd<<Qee^2>>t #"
      write(30,*) "# --------------------------------------------- #"
      open (unit=32,file='data/QsprSq.dat',status='unknown',position='append')
      write(32,*) "# Wi, dt, <<Qspr^2>t, <sd<Qspr^2>>t, sd<<Qspr^2>>t #"
      write(32,*) "# ------------------------------------------------ #"
      open (unit=9,file='data/Qee.dat',status='unknown',position='append')
      write(9,*) "# Wi, dt, <<|Qee|>/|Qee|_mx>t, <sd<|Qee|>/|Qee|_mx>t, sd<<|Qee|>/|Qee|_mx>t #"
      write(9,*) "# ------------------------------------------------------------------------- #"
      if (tplgy == 'Comb') then
        do iarm=1, Na
          write(fnme1,"(A,i0.2,'.dat')") 'data/QeeSqArm',iarm
          write(fnme2,"(A,i0.2,'.dat')") 'data/QeeRelArm',iarm
          open (ue+2*Na+iarm,file=trim(adjustl(fnme1)),&
                status='unknown',position='append')
          open (ue+3*Na+iarm,file=trim(adjustl(fnme2)),&
                status='unknown',position='append')
          write(ue+2*Na+iarm,*) "# Wi, dt, <<Qee^2>>t, sd<Qee.Qee> #"
          write(ue+2*Na+iarm,*) "# ------------------------------------ #"
          write(ue+3*Na+iarm,*) "# Wi, dt, time, <|Qee|>/|Qee|_max, sd<Qee>/|Qee|_max #" 
          write(ue+3*Na+iarm,*) "# -------------------------------------------------- #"
        end do
      end if
      open (unit=10,file='data/EtavsWi.dat',status='replace',position='append')
      write(10,*) "# Wi, dt, <Etap=-<Taupxy>/Pe>t, <sd<Taupxy>/Pe>t, sd<Etap>t #"
      write(10,*) "# --------------------------------------------------------- #"
      open (unit=11,file='data/Psi1vsWi.dat',status='replace',position='append') 
      write(11,*) "# Wi, dt, <Psi1=-<Taupxx-Taupyy>/Pe^2>t, <sd<Taipxx-Taupyy>/Pe^2>t, sd<Psi1>t #"
      write(11,*) "# --------------------------------------------------------------------------- #"
      open (unit=12,file='data/Psi2vsWi.dat',status='replace',position='append') 
      write(12,*) "# Wi, dt, <Psi2=-<Taupyy-Taupzz>/Pe^2>t, <sd<Taipyy-Taupzz>/Pe^2>t, sd<Psi2>t #"
      write(12,*) "# --------------------------------------------------------------------------- #"
      if (iflow >= 3) then ! For Elongational Flow
        open (unit=13,file='data/EtaElongvsWi.dat',status='unknown',position='append') 
        write(13,*) "# Wi, dt, <Eta_el=-<Taupxx-Taupyy>/Pe>t, <sd<Taipxx-Taupyy>/Pe>t, sd<Eta_el>t #"
        write(13,*) "# --------------------------------------------------------------------------- #"
        if (initmode == 'rst') then
          open (unit=14,file='data/EtaElongvsEps.dat',status='unknown',position='append') 
        else
          open (unit=14,file='data/EtaElongvsEps.dat',status='replace',position='append') 
        end if 
        write(14,*) "# Eps(Hencky Strain), dt, Eta_elong=-<Taupxx-Taupyy>/Pe, sd<Taipxx-Taupyy>/Pe #"
        write(14,*) "# --------------------------------------------------------------------------- #"
      end if
    end if ! StrCalc
    if (CoM) then
      open (unit=15,file='data/Dcm.dat',status='unknown',position='append')
      write(15,*) "# nbead, dt, <Dcm>t, <sd(Dcm)>t, sd<Dcm>t #"
      write(15,*) "# --------------------------------------- #"            
    end if
    if (CoHR) then
      open (unit=16,file='data/Dchr.dat',status='unknown',position='append')
      write(16,*) "# nbead, dt, <Dchr>t, <sd(Dchr)>t, sd<Dchr>t #"
      write(16,*) "# ------------------------------------------ #"            
    end if
    if (RgCalc) then
      open (unit=17,file='data/RgSq.dat',status='unknown',position='append')
      write(17,*) "# nbead, dt, <Rg2>t, <sd(Rg2)>t, sd<Rg2>t #"
      write(17,*) "# --------------------------------------- #"            
      open (unit=18,file='data/Asphericity.dat',status='unknown',position='append')
      write(18,*) "# nbead, dt, <Asphericity>t, <sd(Asphericity)>t, sd<Asphericity>t #"
      write(18,*) "# --------------------------------------------------------------- #"            
    end if
    if (cosThCalc) then
      if (initmode == 'rst') then
        open (unit=28,file='data/cosThvsTime.dat',status='unknown',position='append')
      else
        open (unit=28,file='data/cosThvsTime.dat',status='replace',position='append')
      end if 
      write(28,*) "# Wi, dt, Time, <cosTh>, sd<cosTh> #"
      write(28,*) "# -------------------------------- #"
      open (unit=29,file='data/cosTh.dat',status='unknown',position='append')
      write(29,*) "# nbead, dt, <cosTh>t, <sd(cosTh)>t, sd<cosTh>t #"
      write(29,*) "# --------------------------------------------- #"            
    end if
 
    if (AveIterCalc) then
      if (DecompMeth == 'Lanczos') then         
        open (unit=19,file='data/tAvemAveTot.dat',status='replace',position='append')
        write(19,*) "# nbead, dt, <m>t, <sd(m)>t, sd<m>t #"
        write(19,*) "# ------------------------------------------- #"            
        if (initmode == 'rst') then
          open (unit=20,file='data/AvemvsTime.dat',status='unknown',position='append') 
        else
          open (unit=20,file='data/AvemvsTime.dat',status='replace',position='append') 
        end if 
        write(20,*) "# iflow, Wi, dt, Time, <m>, sd<m> #"
        write(20,*) "# ------------------------------- #"
      elseif (DecompMeth == 'Chebyshev') then
        open (unit=19,file='data/tAveLAveTot.dat',status='replace',position='append')
        write(19,*) "# nbead, dt, <m>t, <sd(m)>t, sd<m>t #"
        write(19,*) "# ------------------------------------------- #"            
        if (initmode == 'rst') then
          open (unit=20,file='data/AveLvsTime.dat',status='unknown',position='append') 
        else
          open (unit=20,file='data/AveLvsTime.dat',status='replace',position='append') 
        end if 
        write(20,*) "# iflow, Wi, dt, Time, <L>, sd<L> #"
        write(20,*) "# ------------------------------- #"
      end if ! DecompMeth
    end if



  end subroutine pp_init

  subroutine pp_init_tm(id)

    use :: inp_mod, only: StrCalc,tplgy,Na,CoM,CoHR,RgCalc,cosThCalc,&
                          AveIterCalc,DecompMeth

    integer,intent(in) :: id
    integer :: iarm

    jcount=0
    if (StrCalc) then
      if (tplgy == 'Comb') then
        if (.not.allocated(qee_ar)) allocate(qee_ar(Na))
        if (.not.allocated(sdqee_ar)) allocate(sdqee_ar(Na))
        if (.not.allocated(sqqee_ar)) allocate(sqqee_ar(Na))
        if (.not.allocated(sdsqqee_ar)) allocate(sdsqqee_ar(Na))
      end if
      if (id == 0) then
        tAvesqqsprAveTot=0._wp;tAvesqqetoeAveTot=0._wp;tAveqetoeAveTot=0._wp
        tAvetauxyTot=0._wp;tAvetauxxyyTot=0._wp;tAvetauyyzzTot=0._wp
        tAvesdsqqsprAveTot=0._wp;tAvesdsqqetoeAveTot=0._wp;tAvesdqetoeAveTot=0._wp
        tAvesdxyTot=0._wp;tAvesdxxyyTot=0._wp;tAvesdyyzzTot=0._wp
        sdtAvesqqetoeAveTot=0._wp;sdtAveqetoeAveTot=0._wp;sdtAvetauxyTot=0._wp
        sdtAvetauxxyyTot=0._wp;sdtAvetauyyzzTot=0._wp
        if (tplgy == 'Comb') then
          if (.not.allocated(qee_art)) allocate(qee_art(Na))
          if (.not.allocated(sdqee_art)) allocate(sdqee_art(Na))
          if (.not.allocated(sqqee_art)) allocate(sqqee_art(Na))
          if (.not.allocated(sdsqqee_art)) allocate(sdsqqee_art(Na))
  
          if (.not.allocated(tAvqee_art)) allocate(tAvqee_art(Na))
          if (.not.allocated(tAvsdqee_art)) allocate(tAvsdqee_art(Na))
          if (.not.allocated(sdtAvqee_art)) allocate(sdtAvqee_art(Na))
  
          if (.not.allocated(tAvsqqee_art)) allocate(tAvsqqee_art(Na))
          if (.not.allocated(tAvsdsqqee_art)) allocate(tAvsdsqqee_art(Na))
          if (.not.allocated(sdtAvsqqee_art)) allocate(sdtAvsqqee_art(Na))
          do iarm=1, Na
            tAvsqqee_art=0._wp;tAvqee_art=0._wp
            tAvsdsqqee_art=0._wp;tAvsdqee_art=0._wp
            sdtAvsqqee_art=0._wp;sdtAvqee_art=0._wp
          end do
        end if ! tplgy == Comb
        if (CoM) then
          tAveDcmAveTot=0._wp;tAvesdDcmAveTot=0._wp;sdtAveDcmAveTot=0._wp
        end if
        if (CoHR) then
          tAveDchrAveTot=0._wp;tAvesdDchrAveTot=0._wp;sdtAveDchrAveTot=0._wp
        end if
        if (RgCalc) then
          tAveRgSqAveTot=0._wp;tAvesdRgSqAveTot=0._wp;sdtAveRgSqAveTot=0._wp
          tAveAspherAveTot=0._wp;tAvesdAspherAveTot=0._wp;sdtAveAspherAveTot=0._wp
        end if
        if (cosThCalc) then
          tAvecosTh=0._wp;tAvesdcosTh=0._wp;sdtAvecosTh=0._wp
        end if
        if (AveIterCalc) then
          kcount=0
          if (DecompMeth == 'Lanczos') then
            tAvemAveTot=0._wp;tAvesdmAveTot=0._wp;sdtAvemAveTot=0._wp
          elseif (DecompMeth == 'Chebyshev') then
            tAveLAveTot=0._wp;tAvesdLAveTot=0._wp;sdtAveLAveTot=0._wp
          end if
        end if
      end if ! id == 0
    end if ! StrCalc

  end subroutine pp_init_tm

 
  subroutine data_prcs(id,itime,time,idt,iPe,q,rvmrc,Fphi,rcm,rcmstart,rchr,&
                       nseg_bb,mch,Lch,MPI_REAL_WP)

    use :: mpi
    use :: inp_mod, only: npchain,nchain,nbead,nseg,tss,trst,lambda,Wi,dt,&
                          ntime,qmax,nseg_ar,Pe,ntotang,iflow,StrCalc,tplg&
                          &y,Na,cosThCalc,CoM,CoHR,RgCalc,nbeadx3,AveIterC&
                          &alc,DecompMeth,initmode

    integer,intent(in) :: id,itime,nseg_bb,MPI_REAL_WP,idt,iPe
    real(wp),intent(in) :: time
    real(wp),intent(in) :: Fphi(:,:),rcmstart(:,:)
    real(wp),intent(in),target :: q(:,:),rvmrc(:,:)
    real(wp),intent(in),target :: rcm(:,:),rchr(:,:)
    integer,intent(in) :: mch(:),Lch(:)
    
    real(wp) :: tauxxTot,tauxyTot,tauyyTot,tauzzTot,tauxxyyTot,tauyyzzTot
    real(wp) :: sdxxTot,sdxyTot,sdyyTot,sdzzTot,sdxxyyTot,sdyyzzTot
    real(wp) :: qetoeAveTot,sdqetoeAveTot,sqqetoeAveTot,sdsqqetoeAveTot
    real(wp) :: sqqsprAveTot,sdsqqsprAveTot
    real(wp) :: cosThAve,sdcosThAve,qi(3),qj(3),qi_mag,qj_mag
    real(wp) :: DcmAve,sdDcmAve,rcmdiff(3),rchrdiff(3),DchrAve,sdDchrAve
    real(wp) :: cosThAveTot,sdcosThAveTot
    real(wp) :: DcmAveTot,sdDcmAveTot,DchrAveTot,sdDchrAveTot
    real(wp) :: RgSqAve,sdRgSqAve,AspherAve,sdAspherAve
    real(wp) :: RgSqAveTot,sdRgSqAveTot,AspherAveTot,sdAspherAveTot
    real(wp) :: RgSqTens(3,3),traceRgSqTens,RgSqTensEVbar,RgSqTensEV(3)
    real(wp) :: RgSqEVdiff(3),Aspher,cosTh
    real(wp) :: mAve,sdmAve,mAveTot,sdmAveTot
    real(wp) :: LAve,sdLAve,LAveTot,sdLAveTot
    real(wp),pointer :: qP(:),rcmP(:),rchrP(:),rvmrcP(:),RPi(:),RPj(:)
    integer :: offseti,offsetj,ichain,ibead,jchain,i,j,ierr,info,iarm
 
    ! Foramats used:
1   format(6(f8.5,1x))
2   format(f8.2,1x,e11.3,1x,f14.7,2x,f20.8,2x,f14.7)
3   format(f8.2,1x,e11.3,1x,2(f20.7,2x))
4   format(i4,1x,e11.3,1x,2(f14.7,1x),2(f7.2,2x))
5   format(f8.2,1x,e11.3,1x,3(f14.7,2x))
7   format(i4,1x,e11.3,1x,3(f14.7,2x))
8   format(i4,1x,e11.3,1x,2(f14.7,2x))

    if (StrCalc) then
      call StressCalc(q,rvmrc,Fphi,nseg_bb)
      ! Reduction of terms for stress and relative extension
      call MPI_Reduce(sqqsprAve,sqqsprAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sqqetoeAve,sqqetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(qetoeAve,qetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (tplgy == 'Comb') then
        call MPI_Reduce(sqqee_ar,sqqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(qee_ar,qee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if              
      call MPI_Reduce(tauxx,tauxxTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tauxy,tauxyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tauyy,tauyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tauzz,tauzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tauxxyy,tauxxyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tauyyzz,tauyyzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      ! Reduction of terms for standard deviation of stress and relative extension
      call MPI_Reduce(sdsqqsprAve,sdsqqsprAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdsqqetoeAve,sdsqqetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdqetoeAve,sdqetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (tplgy == 'Comb') then
        call MPI_Reduce(sdsqqee_ar,sdsqqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdqee_ar,sdqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if              
      call MPI_Reduce(sdxx,sdxxTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdxy,sdxyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdyy,sdyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdzz,sdzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdxxyy,sdxxyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdyyzz,sdyyzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    end if ! StrCalc
    !----------------------------------------------------------------------------!
    ! <cos(theta)>: The ensemble average of consec. angle between any 2 segments.!
    !----------------------------------------------------------------------------!
    if (cosThCalc) then
      cosThAve=0._wp;sdcosThAve=0._wp
      do ichain=1, npchain
        do ibead=2, nbead-1
          offseti=(ibead-2)*3
          offsetj=(ibead-1)*3
          qi(1:3)=q(offseti+1:offseti+3,ichain)
          qj(1:3)=q(offsetj+1:offsetj+3,ichain)
          qi_mag=sqrt(dot_product(qi,qi))
          qj_mag=sqrt(dot_product(qj,qj))
!                cosTh=dot_product(qi,qj)/(qi_mag*qj_mag)
          cosTh=dot_product(qi,qj)
          cosThAve=cosThAve+cosTh
          sdcosThAve=sdcosThAve+cosTh*cosTh
        end do ! ibead
      end do ! ichain
      call MPI_Reduce(cosThAve,cosThAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdcosThAve,sdcosThAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    end if
    if (time >= tss*lambda) then

      ! Center of mass diffusion calculation
      if (CoM) then
        DcmAve=0.0_wp;sdDcmAve=0.0_wp
        do jchain=1, npchain
          rcmP => rcm(:,jchain)
          rcmdiff=rcmP-rcmstart(:,jchain)                
!                DcmAve=DcmAve+dot(rcmP,rcmP)
!                sdDcmAve=sdDcmAve+(dot(rcmP,rcmP))**2                       
          DcmAve=DcmAve+dot(rcmdiff,rcmdiff)
          sdDcmAve=sdDcmAve+(dot(rcmdiff,rcmdiff))**2                       
        end do
        call MPI_Reduce(DcmAve,DcmAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdDcmAve,sdDcmAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      ! Center of hydrodynamic resistance diffusion calculation
      if (CoHR) then
        DchrAve=0._wp;sdDchrAve=0._wp
        do jchain=1, npchain
          rchrP => rchr(:,jchain)
          rchrdiff=rchrP-rcmstart(:,jchain)
!                DchrAve=DchrAve+dot(rchrP,rchrP)                      
!                sdDchrAve=sdDchrAve+(dot(rchrP,rchrP))**2                       
          DchrAve=DchrAve+dot(rchrdiff,rchrdiff)                      
          sdDchrAve=sdDchrAve+(dot(rchrdiff,rchrdiff))**2
        end do
        call MPI_Reduce(DchrAve,DchrAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdDchrAve,sdDchrAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
      if (RgCalc) then
        ! R=rv-rc=Bmat.q
        RgSqAve=0._wp;sdRgSqAve=0._wp;AspherAve=0._wp;sdAspherAve=0._wp
        do jchain=1, npchain
          qP => q(:,jchain)
          rvmrcP => rvmrc(:,jchain)
!                call gemv(Bmat,qP,rvmrcP)
          do j=1, 3
            RPj => rvmrcP(j:nbeadx3-(3-j):3)
            do i=1, 3
              RPi => rvmrcP(i:nbeadx3-(3-i):3)
              RgSqTens(i,j)=1._wp/nbead*dot(RPi,RPj)
            end do
          end do
          traceRgSqTens=0._wp
          do i=1, 3
            traceRgSqTens=traceRgSqTens+RgSqTens(i,i)
          end do
          RgSqTensEVbar=traceRgSqTens/3
          call syev(RgSqTens,RgSqTensEV,jobz='N',info=info)
          if (info /= 0) then
            write(*,*) 'Unsuccessful eigenvalue computation For Rg^2 Tensor in main'
            write(*,'(a,1x,i3)') 'info:',info
            stop
          end if
          RgSqEVdiff(:)=RgSqTensEV(:)-RgSqTensEVbar
          Aspher=1._wp/6/RgSqTensEVbar**2*dot(RgSqEVdiff,RgSqEVdiff)
          RgSqAve=RgSqAve+traceRgSqTens
          sdRgSqAve=sdRgSqAve+traceRgSqTens*traceRgSqTens
          AspherAve=AspherAve+Aspher
          sdAspherAve=sdAspherAve+Aspher*Aspher
        end do
        call MPI_Reduce(RgSqAve,RgSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdRgSqAve,sdRgSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(AspherAve,AspherAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdAspherAve,sdAspherAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,&
                        ierr)
      end if ! RgCalc

    end if ! time >= tss*lambda
  
    if (AveIterCalc) then
      if (DecompMeth == 'Lanczos') then
        mAve=0._wp;sdmAve=0._wp
        do jchain=1, npchain
          mAve=mAve+real(mch(jchain),kind=wp)
          sdmAve=sdmAve+(real(mch(jchain),kind=wp))**2                       
        end do
        call MPI_Reduce(mAve,mAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdmAve,sdmAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      elseif (DecompMeth == 'Chebyshev') then
        LAve=0._wp;sdLAve=0._wp
        do jchain=1, npchain
          LAve=LAve+Lch(jchain)                     
          sdLAve=sdLAve+Lch(jchain)**2                       
        end do
        call MPI_Reduce(LAve,LAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdLAve,sdLAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
    end if ! AveIterCalc

    if (id == 0) then

      if (StrCalc) then          
        sdsqqsprAveTot=sqrt(abs(sdsqqsprAveTot-sqqsprAveTot**2)/(nchain-1))
        sdsqqetoeAveTot=sqrt(abs(sdsqqetoeAveTot-sqqetoeAveTot**2)/(nchain-1))
        sdqetoeAveTot=sqrt(abs(sdqetoeAveTot-qetoeAveTot**2)/(nchain-1))
        sdxxTot=sqrt(abs(sdxxTot-tauxxTot**2)/(nchain-1))
        sdxyTot=sqrt(abs(sdxyTot-tauxyTot**2)/(nchain-1))
        sdxxyyTot=sqrt(abs(sdxxyyTot-tauxxyyTot**2)/(nchain-1))
        sdyyzzTot=sqrt(abs(sdyyzzTot-tauyyzzTot**2)/(nchain-1))
        if (initmode == 'rst') then
          write(31,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqsprAveTot,sdsqqsprAveTot
          write(33,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqetoeAveTot,sdsqqetoeAveTot
          select case (tplgy)
            case ('Linear')
              write(7,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qetoeAveTot/(qmax*nseg),&
                         sdqetoeAveTot/(qmax*nseg)
            case ('Comb')
              write(7,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qetoeAveTot/(qmax*nseg_bb),&
                         sdqetoeAveTot/(qmax*nseg_bb)
              do iarm=1, Na
                sdsqqee_art(iarm)=sqrt(abs(sdsqqee_art(iarm)-sqqee_art(iarm)**2)/(nchain-1))
                sdqee_art(iarm)=sqrt(abs(sdqee_art(iarm)-qee_art(iarm)**2)/(nchain-1))
                write(ue+iarm,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqee_art(iarm),&
                                 sdsqqee_art(iarm)
                write(ue+Na+iarm,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qee_art(iarm)/&
                                 (qmax*nseg_ar),sdqee_art(iarm)/(qmax*nseg_ar)
              end do
          end select
        else ! initmode == 'st'
          write(31,2) Wi(iPe),dt(iPe,idt),time,sqqsprAveTot,sdsqqsprAveTot
          write(33,2) Wi(iPe),dt(iPe,idt),time,sqqetoeAveTot,sdsqqetoeAveTot
          select case (tplgy)
            case ('Linear')
              write(7,2) Wi(iPe),dt(iPe,idt),time,qetoeAveTot/(qmax*nseg),&
                         sdqetoeAveTot/(qmax*nseg)
            case ('Comb')
              write(7,2) Wi(iPe),dt(iPe,idt),time,qetoeAveTot/(qmax*nseg_bb),&
                         sdqetoeAveTot/(qmax*nseg_bb)
              do iarm=1, Na
                sdsqqee_art(iarm)=sqrt(abs(sdsqqee_art(iarm)-sqqee_art(iarm)**2)/(nchain-1))
                sdqee_art(iarm)=sqrt(abs(sdqee_art(iarm)-qee_art(iarm)**2)/(nchain-1))
                write(ue+iarm,2) Wi(iPe),dt(iPe,idt),time,sqqee_art(iarm),sdsqqee_art(iarm)
                write(ue+Na+iarm,2) Wi(iPe),dt(iPe,idt),time,qee_art(iarm)/(qmax*nseg_ar),&
                                    sdqee_art(iarm)/(qmax*nseg_ar)
              end do
          end select
        end if
        if (iflow /= 1) then
          if (iflow < 3) then
            if (initmode == 'rst') then
              write(8,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)
            else
              write(8,2) Wi(iPe),dt(iPe,idt),time,-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)
            end if
          else
            if (initmode == 'rst') then 
              write(14,3) Pe(iPe)*(time+trst*lambda),dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe)
            else
              write(14,3) Pe(iPe)*time,dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe)
            end if
          end if
        end if
      end if ! StrCalc
      if (cosThCalc) then
        cosThAveTot=cosThAveTot/ntotang
        sdcosThAveTot=sdcosThAveTot/ntotang
        sdcosThAveTot=sqrt(abs(sdcosThAveTot-cosThAveTot*cosThAveTot)/(ntotang-1))
        if (initmode == 'rst') then
          write(28,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),cosThAveTot,sdcosThAveTot
        else
          write(28,2) Wi(iPe),dt(iPe,idt),time,cosThAveTot,sdcosThAveTot
        end if
      end if ! cosThCalc
      if (AveIterCalc) then
        if (DecompMeth == 'Lanczos') then
          mAveTot=mAveTot/nchain
          sdmAveTot=sdmAveTot/nchain
          sdmAveTot=sqrt(abs(sdmAveTot-mAveTot**2)/(nchain-1))
          if (initmode == 'rst') then
            write(20,4) iflow,Wi(iPe),dt(iPe,idt),(time+trst*lambda),mAveTot,sdmAveTot
          else
            write(20,4) iflow,Wi(iPe),dt(iPe,idt),time,mAveTot,sdmAveTot
          end if
          kcount=kcount+1
          tAvemAveTot=tAvemAveTot+mAveTot
          tAvesdmAveTot=tAvesdmAveTot+sdmAveTot
          sdtAvemAveTot=sdtAvemAveTot+mAveTot*mAveTot
        elseif (DecompMeth == 'Chebyshev') then          
          kcount=kcount+1
          LAveTot=LAveTot/nchain
          sdLAveTot=sdLAveTot/nchain
          sdLAveTot=sqrt(abs(sdLAveTot-LAveTot**2)/(nchain-1))
          if (initmode == 'rst') then
            write(20,4) iflow,Wi(iPe),dt(iPe,idt),(time+trst*lambda),LAveTot,sdLAveTot
          else
            write(20,4) iflow,Wi(iPe),dt(iPe,idt),time,LAveTot,sdLAveTot
          end if
          tAveLAveTot=tAveLAveTot+LAveTot
          tAvesdLAveTot=tAvesdLAveTot+sdLAveTot
          sdtAveLAveTot=sdtAveLAveTot+LAveTot*LAveTot
        end if
      end if ! AveIterCalc
      ! After ?*'chain-largest (or) characteristic-relaxation-time' average over time:
      if (time >= tss*lambda) then
        jcount=jcount+1
        if (StrCalc) then                      
          ! Terms for stress and relative extension
          tAvesqqsprAveTot=tAvesqqsprAveTot+sqqsprAveTot
          tAvesqqetoeAveTot=tAvesqqetoeAveTot+sqqetoeAveTot
          tAveqetoeAveTot=tAveqetoeAveTot+qetoeAveTot
          tAvetauxyTot=tAvetauxyTot+tauxyTot
          tAvetauxxyyTot=tAvetauxxyyTot+tauxxyyTot
          tAvetauyyzzTot=tAvetauyyzzTot+tauyyzzTot
          ! Terms for standard deviation of stress and relative extension
          tAvesdsqqsprAveTot=tAvesdsqqsprAveTot+sdsqqsprAveTot
          tAvesdsqqetoeAveTot=tAvesdsqqetoeAveTot+sdsqqetoeAveTot
          tAvesdqetoeAveTot=tAvesdqetoeAveTot+sdqetoeAveTot
          tAvesdxyTot=tAvesdxyTot+sdxyTot
          tAvesdxxyyTot=tAvesdxxyyTot+sdxxyyTot
          tAvesdyyzzTot=tAvesdyyzzTot+sdyyzzTot
          ! Standard deviation of terms for stress and relative extension
          sdtAvesqqsprAveTot=sdtAvesqqsprAveTot+sqqsprAveTot*sqqsprAveTot
          sdtAvesqqetoeAveTot=sdtAvesqqetoeAveTot+sqqetoeAveTot*sqqetoeAveTot
          sdtAveqetoeAveTot=sdtAveqetoeAveTot+qetoeAveTot*qetoeAveTot
          sdtAvetauxyTot=sdtAvetauxyTot+tauxyTot*tauxyTot
          sdtAvetauxxyyTot=sdtAvetauxxyyTot+tauxxyyTot*tauxxyyTot
          sdtAvetauyyzzTot=sdtAvetauyyzzTot+tauyyzzTot*tauyyzzTot
          if (tplgy == 'Comb') then
            do iarm=1, Na
              tAvsdsqqee_art(iarm)=tAvsdsqqee_art(iarm)+sdsqqee_art(iarm)
              tAvsdqee_art(iarm)=tAvsdqee_art(iarm)+sdqee_art(iarm)
              sdtAvsqqee_art(iarm)=sdtAvsqqee_art(iarm)+sqqee_art(iarm)**2
              sdtAvqee_art(iarm)=sdtAvqee_art(iarm)+qee_art(iarm)**2
              tAvsqqee_art(iarm)=tAvsqqee_art(iarm)+sqqee_art(iarm)
              tAvqee_art(iarm)=tAvqee_art(iarm)+qee_art(iarm)
            end do
          end if
        end if
        if (CoM) then
          if (initmode == 'rst') then
            DcmAveTot=DcmAveTot/(6*(time+trst*lambda)*nchain)
          else
            DcmAveTot=DcmAveTot/(6*time*nchain)
          end if
          tAveDcmAveTot=tAveDcmAveTot+DcmAveTot
          if (initmode == 'rst') then                
            sdDcmAveTot=sdDcmAveTot/(36*(time+trst*lambda)**2*nchain)
          else
            sdDcmAveTot=sdDcmAveTot/(36*time**2*nchain)
          end if
          sdDcmAveTot=sqrt(abs(sdDcmAveTot-DcmAveTot**2)/(nchain-1))
          tAvesdDcmAveTot=tAvesdDcmAveTot+sdDcmAveTot
          sdtAveDcmAveTot=sdtAveDcmAveTot+DcmAveTot*DcmAveTot
        end if
        if (CoHR) then
          if (initmode == 'rst') then
            DchrAveTot=DchrAveTot/(6*(time+trst*lambda)*nchain)
          else
            DchrAveTot=DchrAveTot/(6*time*nchain)
          end if
          tAveDchrAveTot=tAveDchrAveTot+DchrAveTot
          if (initmode == 'rst') then
            sdDchrAveTot=sdDchrAveTot/(36*(time+trst*lambda)**2*nchain)
          else
            sdDchrAveTot=sdDchrAveTot/(36*time**2*nchain)
          end if
          sdDchrAveTot=sqrt(abs(sdDchrAveTot-DchrAveTot**2)/(nchain-1))
          tAvesdDchrAveTot=tAvesdDchrAveTot+sdDchrAveTot
          sdtAveDchrAveTot=sdtAveDchrAveTot+DchrAveTot*DchrAveTot
        end if
        if (RgCalc) then
          RgSqAveTot=RgSqAveTot/nchain
          AspherAveTot=AspherAveTot/nchain
          tAveRgSqAveTot=tAveRgSqAveTot+RgSqAveTot
          tAveAspherAveTot=tAveAspherAveTot+AspherAveTot
          sdRgSqAveTot=sdRgSqAveTot/nchain
          sdAspherAveTot=sdAspherAveTot/nchain
          sdRgSqAveTot=sqrt(abs(sdRgSqAveTot-RgSqAveTot**2)/(nchain-1))
          sdAspherAveTot=sqrt(abs(sdAspherAveTot-AspherAveTot**2)/(nchain-1))
          tAvesdRgSqAveTot=tAvesdRgSqAveTot+sdRgSqAveTot
          tAvesdAspherAveTot=tAvesdAspherAveTot+sdAspherAveTot
          sdtAveRgSqAveTot=sdtAveRgSqAveTot+RgSqAveTot*RgSqAveTot
          sdtAveAspherAveTot=sdtAveAspherAveTot+AspherAveTot*AspherAveTot
        end if ! RgCalc
        if (cosThCalc) then
          ! For cosTh:
          tAvecosTh=tAvecosTh+cosThAveTot
          tAvesdcosTh=tAvesdcosTh+sdcosThAveTot
          sdtAvecosTh=sdtAvecosTh+cosThAveTot*cosThAveTot
        end if ! cosThCalc
        
      end if
      ! Finalize the time averages in the last iteration 
      if (itime == ntime(iPe,idt)) then
        if (jcount /= 0) then
          if (StrCalc) then          
            ! Terms for stress and relative extension
            tAvesqqsprAveTot=tAvesqqsprAveTot/jcount
            tAvesqqetoeAveTot=tAvesqqetoeAveTot/jcount
            tAveqetoeAveTot=tAveqetoeAveTot/jcount
            tAvetauxyTot=tAvetauxyTot/jcount
            tAvetauxxyyTot=tAvetauxxyyTot/jcount
            tAvetauyyzzTot=tAvetauyyzzTot/jcount
            ! Terms for standard deviation of stress and relative extension
            tAvesdsqqsprAveTot=tAvesdsqqsprAveTot/jcount
            tAvesdsqqetoeAveTot=tAvesdsqqetoeAveTot/jcount
            tAvesdqetoeAveTot=tAvesdqetoeAveTot/jcount
            tAvesdxyTot=tAvesdxyTot/jcount
            tAvesdxxyyTot=tAvesdxxyyTot/jcount
            tAvesdyyzzTot=tAvesdyyzzTot/jcount
            ! Standard deviation of terms for stress and relative extension
            sdtAvesqqsprAveTot=sdtAvesqqsprAveTot/jcount
            sdtAvesqqetoeAveTot=sdtAvesqqetoeAveTot/jcount
            sdtAveqetoeAveTot=sdtAveqetoeAveTot/jcount
            sdtAvetauxyTot=sdtAvetauxyTot/jcount
            sdtAvetauxxyyTot=sdtAvetauxxyyTot/jcount
            sdtAvetauyyzzTot=sdtAvetauyyzzTot/jcount
            sdtAvesqqsprAveTot=sqrt(abs(sdtAvesqqsprAveTot-tAvesqqsprAveTot**2)/(jcount-1))
            sdtAvesqqetoeAveTot=sqrt(abs(sdtAvesqqetoeAveTot-tAvesqqetoeAveTot**2)/(jcount-1))
            sdtAveqetoeAveTot=sqrt(abs(sdtAveqetoeAveTot-tAveqetoeAveTot**2)/(jcount-1))
            if (tplgy == 'Comb') then
              do iarm=1, Na
                tAvsdsqqee_art(iarm)=tAvsdsqqee_art(iarm)/jcount
                tAvsdqee_art(iarm)=tAvsdqee_art(iarm)/jcount
                tAvsqqee_art(iarm)=tAvsqqee_art(iarm)/jcount
                tAvqee_art(iarm)=tAvqee_art(iarm)/jcount
                sdtAvsqqee_art(iarm)=sdtAvsqqee_art(iarm)/jcount
                sdtAvqee_art(iarm)=sdtAvqee_art(iarm)/jcount
                sdtAvsqqee_art(iarm)=sqrt(abs(sdtAvsqqee_art(iarm)-tAvsqqee_art(iarm)**2)/(jcount-1))
                sdtAvqee_art(iarm)=sqrt(abs(sdtAvqee_art(iarm)-tAvqee_art(iarm)**2)/(jcount-1))
              end do
            end if
            sdtAvetauxyTot=sqrt(abs(sdtAvetauxyTot-tAvetauxyTot**2)/(jcount-1))
            sdtAvetauxxyyTot=sqrt(abs(sdtAvetauxxyyTot-tAvetauxxyyTot**2)/(jcount-1))
            sdtAvetauyyzzTot=sqrt(abs(sdtAvetauyyzzTot-tAvetauyyzzTot**2)/(jcount-1))
            write(32,5) Wi(iPe),dt(iPe,idt),tAvesqqsprAveTot,tAvesdsqqsprAveTot,sdtAvesqqsprAveTot
            write(30,5) Wi(iPe),dt(iPe,idt),tAvesqqetoeAveTot,tAvesdsqqetoeAveTot,sdtAvesqqetoeAveTot
            select case (tplgy)
              case ('Linear')
                write(9,5) Wi(iPe),dt(iPe,idt),tAveqetoeAveTot/(qmax*nseg),tAvesdqetoeAveTot/&
                          (qmax*nseg),sdtAveqetoeAveTot/(qmax*nseg)
              case ('Comb')
                write(9,5) Wi(iPe),dt(iPe,idt),tAveqetoeAveTot/(qmax*nseg_bb),tAvesdqetoeAveTot/&
                          (qmax*nseg_bb),sdtAveqetoeAveTot/(qmax*nseg_bb)
                do iarm=1, Na
                  write(ue+2*Na+iarm,5) Wi(iPe),dt(iPe,idt),tAvsqqee_art(iarm),tAvsdsqqee_art(iarm),&
                                        sdtAvsqqee_art(iarm)
                  write(ue+3*Na+iarm,5) Wi(iPe),dt(iPe,idt),tAvqee_art(iarm)/(qmax*nseg_ar),&
                                        tAvsdqee_art(iarm)/(qmax*nseg_ar),sdtAvqee_art(iarm)/&
                                        (qmax*nseg_ar)
                end do
            end select
            if (iflow /= 1) then
              write(10,5) Wi(iPe),dt(iPe,idt),-tAvetauxyTot/Pe(iPe),tAvesdxyTot/Pe(iPe),&
                          sdtAvetauxyTot/Pe(iPe)
              write(11,5) Wi(iPe),dt(iPe,idt),-tAvetauxxyyTot/Pe(iPe)**2,tAvesdxxyyTot/Pe(iPe)**2,&
                          sdtAvetauxxyyTot/Pe(iPe)**2
              write(12,5) Wi(iPe),dt(iPe,idt),-tAvetauyyzzTot/Pe(iPe)**2,tAvesdyyzzTot/Pe(iPe)**2,&
                          sdtAvetauyyzzTot/Pe(iPe)**2
              if (iflow >= 3) then ! For Elongational Flow
                write(13,5) Wi(iPe),dt(iPe,idt),-tAvetauxxyyTot/Pe(iPe),tAvesdxxyyTot/Pe(iPe),&
                            sdtAvetauxxyyTot/Pe(iPe) 
              end if
            end if ! iflow /= 1
          end if ! StrCalc
          if (CoM) then
            tAveDcmAveTot=tAveDcmAveTot/jcount
            tAvesdDcmAveTot=tAvesdDcmAveTot/jcount
            sdtAveDcmAveTot=sdtAveDcmAveTot/jcount
            sdtAveDcmAveTot=sqrt(abs(sdtAveDcmAveTot-tAveDcmAveTot**2)/(jcount-1))
            write(15,7) nbead,dt(iPe,idt),tAveDcmAveTot,tAvesdDcmAveTot,sdtAveDcmAveTot 
          end if
          if (CoHR) then
            tAveDchrAveTot=tAveDchrAveTot/jcount
            tAvesdDchrAveTot=tAvesdDchrAveTot/jcount
            sdtAveDchrAveTot=sdtAveDchrAveTot/jcount
            sdtAveDchrAveTot=sqrt(abs(sdtAveDchrAveTot-tAveDchrAveTot**2)/(jcount-1))
            write(16,7) nbead,dt(iPe,idt),tAveDchrAveTot,tAvesdDchrAveTot,sdtAveDchrAveTot
          end if
          if (RgCalc) then
            tAveRgSqAveTot=tAveRgSqAveTot/jcount
            tAveAspherAveTot=tAveAspherAveTot/jcount
            tAvesdRgSqAveTot=tAvesdRgSqAveTot/jcount
            tAvesdAspherAveTot=tAvesdAspherAveTot/jcount
            sdtAveRgSqAveTot=sdtAveRgSqAveTot/jcount
            sdtAveAspherAveTot=sdtAveAspherAveTot/jcount
            sdtAveRgSqAveTot=sqrt(abs(sdtAveRgSqAveTot-tAveRgSqAveTot**2)/(jcount-1))
            sdtAveAspherAveTot=sqrt(abs(sdtAveAspherAveTot-tAveAspherAveTot**2)/(jcount-1))
            write(17,7) nbead,dt(iPe,idt),tAveRgSqAveTot,tAvesdRgSqAveTot,sdtAveRgSqAveTot 
            write(18,7) nbead,dt(iPe,idt),tAveAspherAveTot,tAvesdAspherAveTot,sdtAveAspherAveTot
          end if ! RgCalc
          if (cosThCalc) then
            ! For cosTh:
            tAvecosTh=tAvecosTh/jcount
            tAvesdcosTh=tAvesdcosTh/jcount
            sdtAvecosTh=sdtAvecosTh/jcount
            sdtAvecosTh=sqrt(abs(sdtAvecosTh-tAvecosTh*tAvecosTh)/(jcount-1))
            write(29,7) nbead,dt(iPe,idt),tAvecosTh,tAvesdcosTh,sdtAvecosTh
          end if ! cosThCalc             
        else ! so, jcount=0
          if (StrCalc) then          
            write(9,3) Wi(iPe),dt(iPe,idt),qetoeAveTot/(qmax*nseg),sdqetoeAveTot/(qmax*nseg)
            if (iflow /= 1) then
              write(10,3) Wi(iPe),dt(iPe,idt),-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)           
              write(11,3) Wi(iPe),dt(iPe,idt),-tauxxyyTot/(Pe(iPe)**2),sdxxyyTot/(Pe(iPe)**2)
              write(12,3) Wi(iPe),dt(iPe,idt),-tauyyzzTot/(Pe(iPe)**2),sdyyzzTot/(Pe(iPe)**2)
              if (iflow >= 3) then ! For Elongational Flow
                write(13,3) Wi(iPe),dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe) 
              end if
            end if ! iflow /= 1
          end if ! StrCalc
          if (CoM) write(15,8) nbead,dt(iPe,idt),DcmAveTot,sdDcmAveTot 
          if (CoHR) write(16,8) nbead,dt(iPe,idt),DchrAveTot,sdDchrAveTot
          if (RgCalc) then
            write(17,8) nbead,dt(iPe,idt),RgSqAveTot,sdRgSqAveTot
            write(18,8) nbead,dt(iPe,idt),AspherAveTot,sdAspherAveTot
          end if              
        end if ! jcount /= 0
        if (AveIterCalc) then
          if (kcount /= 0) then
            if (DecompMeth == 'Lanczos') then
              tAvemAveTot=tAvemAveTot/kcount
              tAvesdmAveTot=tAvesdmAveTot/kcount
              sdtAvemAveTot=sdtAvemAveTot/kcount
              sdtAvemAveTot=sqrt(abs(sdtAvemAveTot-tAvemAveTot**2)/(kcount-1))
              write(19,7) nbead,dt(iPe,idt),tAvemAveTot,tAvesdmAveTot,sdtAvemAveTot
            elseif (DecompMeth == 'Chebyshev') then
              tAveLAveTot=tAveLAveTot/kcount
              tAvesdLAveTot=tAvesdLAveTot/kcount
              sdtAveLAveTot=sdtAveLAveTot/kcount
              sdtAveLAveTot=sqrt(abs(sdtAveLAveTot-tAveLAveTot**2)/(kcount-1))
              write(19,7) nbead,dt(iPe,idt),tAveLAveTot,tAvesdLAveTot,sdtAveLAveTot
            end if
          end if ! kcount /= 0                      
        end if ! AveIterCalc
      end if ! itime == ntime

    end if ! id == 0

  end subroutine

!-------------------------------------------------------------------------!
! Stress is calculated for each chain by summing up the contribution      !
! of all segments. Then all chain stresses in each process will be        !
! added up.                                                               !
! Note that standard deviation is calculated based on the formula for     !
! sample standard deviation of the sampled mean:                          !
! sigma_mean = 1/sqrt(N-1) * sigma ;                                      !
! where sigma is the standard deviation of the distribution:              !
! sigma = sqrt(ave(X^2)-(ave(x))^2)                                       !
! Therefore, we will calculate tauxx and sdxx in each process and finally !
! they will be summed up to give total stress and also appropriate term   !
! for getting the deviation.                                              !
!-------------------------------------------------------------------------!
  subroutine StressCalc(q,rvmrc,Fphi,nseg_bb)
  
    use :: inp_mod, only: nseg,npchain,nchain,ForceLaw,tplgy,nseg_ar
  
    integer,intent(in) :: nseg_bb
    real(wp),intent(in),target :: q(:,:),rvmrc(:,:),Fphi(:,:)
    real(wp),dimension(:),pointer :: qP,qPx,qPy,qPz,RPx,RPy,RPz
    real(wp),dimension(:),pointer :: FphiPx,FphiPy,FphiPz
    real(wp) :: sqqeetmp,qeetmp,qetoex,qetoey,qetoez
    real(wp) :: sqqspr,sqqetoe,qetoe,txx,txy,tyy,tzz,txxyy,tyyzz
    real(wp),allocatable :: qee_artmp(:,:)
    integer :: ichain,iarm,os
  
    sqqetoeAve=0._wp;qetoeAve=0._wp
    sdsqqetoeAve=0._wp;sdqetoeAve=0._wp
    tauxx=0._wp;tauxy=0._wp;tauyy=0._wp
    tauzz=0._wp;tauxxyy=0._wp;tauyyzz=0._wp
    sdxx=0._wp;sdxy=0._wp;sdyy=0._wp
    sdzz=0._wp;sdxxyy=0._wp;sdyyzz=0._wp
    if (tplgy == 'Comb') then
      allocate(qee_artmp(3,size(qee_ar)))
      sqqee_ar=0._wp;qee_ar=0._wp
      sdsqqee_ar=0._wp;sdqee_ar=0._wp
    end if
    do ichain=1, npchain
      qP => q(:,ichain)
!      call gemv(Bmat,qP,R) ! R=Bmat.q !
      qPx => q(1:3*nseg-2:3,ichain)
      qPy => q(2:3*nseg-1:3,ichain) 
      qPz => q(3:3*nseg:3,ichain)
      RPx => rvmrc(1:3*(nseg+1)-2:3,ichain)
      RPy => rvmrc(2:3*(nseg+1)-1:3,ichain) 
      RPz => rvmrc(3:3*(nseg+1):3,ichain)
      select case (tplgy)
        case ('Linear')
          qetoex=sum(qPx)
          qetoey=sum(qPy)
          qetoez=sum(qPz)
        case ('Comb')
          qetoex=sum(qPx(1:nseg_bb))
          qetoey=sum(qPy(1:nseg_bb))
          qetoez=sum(qPz(1:nseg_bb))
          do iarm=1, size(qee_ar)
            os=nseg_bb+(iarm-1)*nseg_ar
            qee_artmp(1,iarm)=sum(qPx(os+1:os+nseg_ar))
            qee_artmp(2,iarm)=sum(qPy(os+1:os+nseg_ar))
            qee_artmp(3,iarm)=sum(qPz(os+1:os+nseg_ar))
            sqqeetmp=dot_product(qee_artmp(:,iarm),qee_artmp(:,iarm))
            qeetmp=sqrt(sqqeetmp)
            sqqee_ar(iarm)=sqqee_ar(iarm)+sqqeetmp
            qee_ar(iarm)=qee_ar(iarm)+qeetmp
            sdsqqee_ar(iarm)=sdsqqee_ar(iarm)+sqqeetmp*sqqeetmp
            sdqee_ar(iarm)=sdqee_ar(iarm)+qeetmp*qeetmp
          end do
      end select
      sqqspr=dot_product(q(:,ichain),q(:,ichain))
      sqqetoe=qetoex**2+qetoey**2+qetoez**2
      qetoe=sqrt(sqqetoe)
!     qetoe=qetoex
      sqqsprAve=sqqsprAve+sqqspr
      sqqetoeAve=sqqetoeAve+sqqetoe
      qetoeAve=qetoeAve+qetoe
      sdsqqsprAve=sdsqqsprAve+sqqspr*sqqspr
      sdsqqetoeAve=sdsqqetoeAve+sqqetoe*sqqetoe
      sdqetoeAve=sdqetoeAve+qetoe*qetoe
      FphiPx => Fphi(1:3*(nseg+1)-2:3,ichain)
      FphiPy => Fphi(2:3*(nseg+1)-1:3,ichain) 
      FphiPz => Fphi(3:3*(nseg+1):3,ichain)
  
      txx=dot(RPx,FphiPx)+nseg
      txy=dot(RPx,FphiPy)
      tyy=dot(RPy,FphiPy)+nseg
      tzz=dot(RPz,FphiPz)+nseg
      txxyy=txx-tyy
      tyyzz=tyy-tzz
  
      tauxx=tauxx+txx
      tauxy=tauxy+txy 
      tauyy=tauyy+tyy
      tauzz=tauzz+tzz
      tauxxyy=tauxxyy+txxyy
      tauyyzz=tauyyzz+tyyzz
      sdxx=sdxx+txx*txx
      sdxy=sdxy+txy*txy
      sdyy=sdyy+tyy*tyy
      sdzz=sdzz+tzz*tzz
      sdxxyy=sdxxyy+txxyy*txxyy
      sdyyzz=sdyyzz+tyyzz*tyyzz
    end do! chain loop
    
    sqqsprAve=sqqsprAve/(nchain*nseg)
    sqqetoeAve=sqqetoeAve/nchain
    qetoeAve=qetoeAve/nchain
    sdsqqsprAve=sdsqqsprAve/nchain
    sdsqqetoeAve=sdsqqetoeAve/nchain
    sdqetoeAve=sdqetoeAve/nchain
    tauxx = tauxx/nchain
    tauxy = tauxy/nchain
    tauyy = tauyy/nchain
    tauzz = tauzz/nchain
    tauxxyy = tauxxyy/nchain
    tauyyzz = tauyyzz/nchain
    sdxx = sdxx/nchain
    sdxy = sdxy/nchain
    sdyy = sdyy/nchain
    sdzz = sdzz/nchain
    sdxxyy = sdxxyy/nchain
    sdyyzz = sdyyzz/nchain
  
    if (tplgy == 'Comb') then
      do iarm=1, size(qee_ar)
        sqqee_ar(iarm)=sqqee_ar(iarm)/nchain
        qee_ar(iarm)=qee_ar(iarm)/nchain
        sdsqqee_ar(iarm)=sdsqqee_ar(iarm)/nchain
        sdqee_ar(iarm)=sdqee_ar(iarm)/nchain
      end do
      deallocate(qee_artmp)
    end if
   
  end subroutine StressCalc

end module pp_mod
