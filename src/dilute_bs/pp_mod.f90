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
module pp_mod

  use :: prcn_mod
  
  implicit none
  
  save
  
  real(wp),private :: tauxx,tauxy,tauyy,tauzz,tauxxyy,tauyyzz
  real(wp),private :: sdxx,sdxy,sdyy,sdzz,sdxxyy,sdyyzz
  real(wp),private :: qetoeAve,sdqetoeAve,sqqetoeAve,sdsqqetoeAve
  real(wp),private :: sqqsprAve,sdsqqsprAve
  real(wp),private :: XAve,sdXAve,XSqAve,sdXSqAve
  real(wp),private :: XbbAve,sdXbbAve,XbbSqAve,sdXbbSqAve

  real(wp),private :: tAvesqqsprAveTot,tAvesqqetoeAveTot,tAveqetoeAveTot
  real(wp),private :: tAvesdsqqsprAveTot,tAvesdsqqetoeAveTot,tAvesdqetoeAveTot
  real(wp),private :: sdtAvesqqsprAveTot
  real(wp),private :: sdtAvesqqetoeAveTot,sdtAveqetoeAveTot
  real(wp),private :: tAveXAveTot,tAvesdXAveTot,sdtAveXAveTot
  real(wp),private :: tAveXSqAveTot,tAvesdXSqAveTot,sdtAveXSqAveTot
  real(wp),private :: tAveXbbAveTot,tAvesdXbbAveTot,sdtAveXbbAveTot
  real(wp),private :: tAveXbbSqAveTot,tAvesdXbbSqAveTot,sdtAveXbbSqAveTot
  real(wp),private :: tAvetauxyTot,tAvetauxxyyTot,tAvetauyyzzTot
  real(wp),private :: tAvesdxyTot,tAvesdxxyyTot,tAvesdyyzzTot
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

  ! File units
  integer,private :: u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20
  integer,private :: u28,u29,u30,u31,u32,u33,u35,u36,u37,u38,u39,u40,u41,u42
  integer,allocatable,private :: uarm(:)
  integer,allocatable,private :: uch(:),uch1(:)

  ! Variables in conf_sort
  real(wp),allocatable,private :: x(:),y(:),z(:),r(:)
  integer,allocatable,private :: ib(:),ib1(:),icount3(:)

contains

  subroutine pp_init(id)

    use :: mpi
    use :: inp_dlt, only: StrCalc,initmode,tplgy,Na,CoM,CoHR,RgCalc,cosThCalc,&
                          AveIterCalc,DecompMeth,iflow,npchain,indvlext

    integer,intent(in) :: id
    integer :: iarm,ichain,ierr
    character(len=1024) :: fnme1,fnme2

!    ue=39 ! The largest unit excluding arms of the comb polymer
    if (id == 0) then
      if (StrCalc) then
        if (initmode == 'rst') then
          open(newunit=u31,file='data/trQsprSq.dat',status='unknown',position='append')
          open(newunit=u33,file='data/trQeeSq.dat',status='unknown',position='append')
          open(newunit=u7,file='data/trQeeRel.dat',status='unknown',position='append')
          open(newunit=u35,file='data/trXRel.dat',status='unknown',position='append')
          open(newunit=u36,file='data/trXSq.dat',status='unknown',position='append')
        else
          open(newunit=u31,file='data/trQsprSq.dat',status='replace',position='append')
          open(newunit=u33,file='data/trQeeSq.dat',status='replace',position='append')
          open(newunit=u7,file='data/trQeeRel.dat',status='replace',position='append')
          open(newunit=u35,file='data/trXRel.dat',status='replace',position='append')
          open(newunit=u36,file='data/trXSq.dat',status='replace',position='append')
        end if
        write(u31,*) "# Wi, dt, Time, <Qspr.Qspr>, sd<Qspr.Qspr> #"
        write(u31,*) "# ---------------------------------------- #"
        write(u33,*) "# Wi, dt, Time, <Qee.Qee>, sd<Qee.Qee> #"
        write(u33,*) "# ------------------------------------ #"
        write(u7,*) "# Wi, dt, time, <|Qee|>/Lc, sd<Qee>/Lc #" 
        write(u7,*) "# ------------------------------------ #"
        write(u35,*) "# Wi, dt, Time, <X>/Lc, sd<X>/Lc #" 
        write(u35,*) "# ------------------------------ #"
        write(u36,*) "# Wi, dt, Time, <X.X>, sd<X.X> #"
        write(u36,*) "# ---------------------------- #"
        if (tplgy == 'Comb') then
          allocate(uarm(4*Na))
          do iarm=1, Na
            write(fnme1,"(A,i0.2,'.dat')") 'data/trQeeSqArm',iarm
            write(fnme2,"(A,i0.2,'.dat')") 'data/trQeeRelArm',iarm
            if (initmode == 'rst') then
              open(newunit=uarm(iarm),file=trim(adjustl(fnme1)),&
                   status='unknown',position='append')
              open(newunit=uarm(Na+iarm),file=trim(adjustl(fnme2)),&
                   status='unknown',position='append')
            else
              open(newunit=uarm(iarm),file=trim(adjustl(fnme1)),&
                   status='replace',position='append')
              open(newunit=uarm(Na+iarm),file=trim(adjustl(fnme2)),&
                   status='replace',position='append')
            end if
            write(uarm(iarm),*) "# Wi, dt, Time, <Qee.Qee>, sd<Qee.Qee> #"
            write(uarm(iarm),*) "# ------------------------------------ #"
            write(uarm(Na+iarm),*) "# Wi, dt, time, <|Qee|>/|Qee|_max, sd<Qee>/|Qee|_max #" 
            write(uarm(Na+iarm),*) "# -------------------------------------------------- #"
          end do
          if (initmode == 'rst') then
            open(newunit=u39,file='data/trXbbRel.dat',status='unknown',position='append')
            open(newunit=u40,file='data/trXbbSq.dat',status='unknown',position='append')
          else
            open(newunit=u39,file='data/trXbbRel.dat',status='replace',position='append')
            open(newunit=u40,file='data/trXbbSq.dat',status='replace',position='append')
          end if
          write(u39,*) "# Wi, dt, Time, <Xbb>/Lc, sd<Xbb>/Lc #" 
          write(u39,*) "# ---------------------------------- #"
          write(u40,*) "# Wi, dt, Time, <Xbb.Xbb>, sd<Xbb.Xbb> #"
          write(u40,*) "# ------------------------------------ #"
        end if
        if (iflow /= 1) then
          if (initmode == 'rst') then
            open(newunit=u8,file='data/EtavsTime.dat',status='unknown',position='append')
          else
            open(newunit=u8,file='data/EtavsTime.dat',status='replace',position='append')
          end if 
          write(u8,*) "# Wi, dt, Time, Etap=<Taupxy>/Pe, sd<Taupxy>/Pe #"
          write(u8,*) "# --------------------------------------------- #"
        end if
        open (newunit=u30,file='data/QeeSq.dat',status='unknown',position='append')
        write(u30,*) "# Wi, dt, <<Qee^2>>t, <sd<Qee^2>>t, sd<<Qee^2>>t #"
        write(u30,*) "# --------------------------------------------- #"
        open (newunit=u32,file='data/QsprSq.dat',status='unknown',position='append')
        write(u32,*) "# Wi, dt, <<Qspr^2>>t, <sd<Qspr^2>>t, sd<<Qspr^2>>t #"
        write(u32,*) "# ------------------------------------------------ #"
        open (newunit=u9,file='data/QeeRel.dat',status='unknown',position='append')
        write(u9,*) "# Wi, dt, <<|Qee|>/Lc>t, <sd<|Qee|>/Lc>t, sd<<|Qee|>/Lc>t #"
        write(u9,*) "# ------------------------------------------------------- #"
        open (newunit=u37,file='data/XRel.dat',status='unknown',position='append')
        write(u37,*) "# Wi, dt, <<X>/Lc>t, <sd<X>/Lc>t, sd<<X>/Lc>t #"
        write(u37,*) "# ------------------------------------------- #"
        open (newunit=u38,file='data/XSq.dat',status='unknown',position='append')
        write(u38,*) "# Wi, dt, <<X.X>>t, <sd<X.X>>t, sd<<X.X>>t #"
        write(u38,*) "# ---------------------------------------- #"
        if (tplgy == 'Comb') then
          do iarm=1, Na
            write(fnme1,"(A,i0.2,'.dat')") 'data/QeeSqArm',iarm
            write(fnme2,"(A,i0.2,'.dat')") 'data/QeeRelArm',iarm
            open(newunit=uarm(2*Na+iarm),file=trim(adjustl(fnme1)),&
                 status='unknown',position='append')
            open(newunit=uarm(3*Na+iarm),file=trim(adjustl(fnme2)),&
                  status='unknown',position='append')
            write(uarm(2*Na+iarm),*) "# Wi, dt, <<Qee^2>>t, sd<Qee.Qee> #"
            write(uarm(2*Na+iarm),*) "# ------------------------------------ #"
            write(uarm(3*Na+iarm),*) "# Wi, dt, time, <|Qee|>/|Qee|_max, sd<Qee>/|Qee|_max #" 
            write(uarm(3*Na+iarm),*) "# -------------------------------------------------- #"
          end do
          open (newunit=u41,file='data/XbbRel.dat',status='unknown',position='append')
          write(u41,*) "# Wi, dt, <<Xbb>/Lc>t, <sd<Xbb>/Lc>t, sd<<Xbb>/Lc>t #"
          write(u41,*) "# ------------------------------------------------- #"
          open (newunit=u42,file='data/XbbSq.dat',status='unknown',position='append')
          write(u42,*) "# Wi, dt, <<Xbb.Xbb>>t, <sd<Xbb.Xbb>>t, sd<<Xbb.Xbb>>t #"
          write(u42,*) "# ---------------------------------------------------- #"
        end if
        if (iflow /= 1) then
          open(newunit=u10,file='data/EtavsWi.dat',status='replace',position='append')
          write(u10,*) "# Wi, dt, <Etap=-<Taupxy>/Pe>t, <sd<Taupxy>/Pe>t, sd<Etap>t #"
          write(u10,*) "# --------------------------------------------------------- #"
          open(newunit=u11,file='data/Psi1vsWi.dat',status='replace',position='append') 
          write(u11,*) "# Wi, dt, <Psi1=-<Taupxx-Taupyy>/Pe^2>t, <sd<Taipxx-Taupyy>/Pe^2>t, sd<Psi1>t #"
          write(u11,*) "# --------------------------------------------------------------------------- #"
          open(newunit=u12,file='data/Psi2vsWi.dat',status='replace',position='append') 
          write(u12,*) "# Wi, dt, <Psi2=-<Taupyy-Taupzz>/Pe^2>t, <sd<Taipyy-Taupzz>/Pe^2>t, sd<Psi2>t #"
          write(u12,*) "# --------------------------------------------------------------------------- #"
        end if
        if (iflow >= 3) then ! For Elongational Flow
          open(newunit=u13,file='data/EtaElongvsWi.dat',status='unknown',position='append') 
          write(u13,*) "# Wi, dt, <Eta_el=-<Taupxx-Taupyy>/Pe>t, <sd<Taipxx-Taupyy>/Pe>t, sd<Eta_el>t #"
          write(u13,*) "# --------------------------------------------------------------------------- #"
          if (initmode == 'rst') then
            open(newunit=u14,file='data/EtaElongvsEps.dat',status='unknown',position='append') 
          else
            open(newunit=u14,file='data/EtaElongvsEps.dat',status='replace',position='append') 
          end if 
          write(u14,*) "# Eps(Hencky Strain), dt, Eta_elong=-<Taupxx-Taupyy>/Pe, sd<Taipxx-Taupyy>/Pe #"
          write(u14,*) "# --------------------------------------------------------------------------- #"
        end if
      end if ! StrCalc
      if (CoM) then
        open(newunit=u15,file='data/Dcm.dat',status='unknown',position='append')
        write(u15,*) "# nbead, dt, <Dcm>t, <sd(Dcm)>t, sd<Dcm>t #"
        write(u15,*) "# --------------------------------------- #"            
      end if
      if (CoHR) then
        open(newunit=u16,file='data/Dchr.dat',status='unknown',position='append')
        write(u16,*) "# nbead, dt, <Dchr>t, <sd(Dchr)>t, sd<Dchr>t #"
        write(u16,*) "# ------------------------------------------ #"            
      end if
      if (RgCalc) then
        open(newunit=u17,file='data/RgSq.dat',status='unknown',position='append')
        write(u17,*) "# nbead, dt, <Rg2>t, <sd(Rg2)>t, sd<Rg2>t #"
        write(u17,*) "# --------------------------------------- #"            
        open (newunit=u18,file='data/Asphericity.dat',status='unknown',position='append')
        write(u18,*) "# nbead, dt, <Asphericity>t, <sd(Asphericity)>t, sd<Asphericity>t #"
        write(u18,*) "# --------------------------------------------------------------- #"            
      end if
      if (cosThCalc) then
        if (initmode == 'rst') then
          open(newunit=u28,file='data/cosThvsTime.dat',status='unknown',position='append')
        else
          open(newunit=u28,file='data/cosThvsTime.dat',status='replace',position='append')
        end if 
        write(u28,*) "# Wi, dt, Time, <cosTh>, sd<cosTh> #"
        write(u28,*) "# -------------------------------- #"
        open(newunit=u29,file='data/cosTh.dat',status='unknown',position='append')
        write(u29,*) "# nbead, dt, <cosTh>t, <sd(cosTh)>t, sd<cosTh>t #"
        write(u29,*) "# --------------------------------------------- #"            
      end if
   
      if (AveIterCalc) then
        if (DecompMeth == 'Lanczos') then         
          open(newunit=u19,file='data/tAvemAveTot.dat',status='replace',position='append')
          write(u19,*) "# nbead, dt, <m>t, <sd(m)>t, sd<m>t #"
          write(u19,*) "# ------------------------------------------- #"            
          if (initmode == 'rst') then
            open(newunit=u20,file='data/AvemvsTime.dat',status='unknown',position='append') 
          else
            open(newunit=u20,file='data/AvemvsTime.dat',status='replace',position='append') 
          end if 
          write(u20,*) "# iflow, Wi, dt, Time, <m>, sd<m> #"
          write(u20,*) "# ------------------------------- #"
        elseif (DecompMeth == 'Chebyshev') then
          open(newunit=u19,file='data/tAveLAveTot.dat',status='replace',position='append')
          write(u19,*) "# nbead, dt, <m>t, <sd(m)>t, sd<m>t #"
          write(u19,*) "# ------------------------------------------- #"            
          if (initmode == 'rst') then
            open(newunit=u20,file='data/AveLvsTime.dat',status='unknown',position='append') 
          else
            open(newunit=u20,file='data/AveLvsTime.dat',status='replace',position='append') 
          end if 
          write(u20,*) "# iflow, Wi, dt, Time, <L>, sd<L> #"
          write(u20,*) "# ------------------------------- #"
        end if ! DecompMeth
      end if
    end if
    if (StrCalc .and. indvlext) then
      allocate(uch(npchain))
      do ichain=1, npchain
        write(fnme1,"(A,i0.4,'.dat')") 'data/trXRelCh',id*npchain+ichain
        if (initmode == 'rst') then
          open(newunit=uch(ichain),file=trim(adjustl(fnme1)),status='unknown',&
               position='append')
        else
          open(newunit=uch(ichain),file=trim(adjustl(fnme1)),status='replace',&
               position='append')
        end if
        write(uch(ichain),*) "# Wi, dt, Time, <X>/Lc, sd<X>/Lc #"
        write(uch(ichain),*) "# ------------------------------ #"
      end do
      if (tplgy == 'Comb') then
        allocate(uch1(npchain))
        do ichain=1, npchain
          write(fnme1,"(A,i0.4,'.dat')") 'data/trXbbRelCh',id*npchain+ichain
          if (initmode == 'rst') then
            open(newunit=uch1(ichain),file=trim(adjustl(fnme1)),status='unknown',&
                 position='append')
          else
            open(newunit=uch1(ichain),file=trim(adjustl(fnme1)),status='replace',&
                 position='append')
          end if
          write(uch1(ichain),*) "# Wi, dt, Time, <Xbb>/Lc, sd<Xbb>/Lc #"
          write(uch1(ichain),*) "# ---------------------------------- #"
        end do
      end if
    end if

  end subroutine pp_init

  subroutine pp_init_tm(id)

    use :: inp_dlt, only: StrCalc,tplgy,Na,CoM,CoHR,RgCalc,cosThCalc,&
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
        tAveXAveTot=0._wp;tAveXSqAveTot=0._wp
        tAvetauxyTot=0._wp;tAvetauxxyyTot=0._wp;tAvetauyyzzTot=0._wp
        tAvesdsqqsprAveTot=0._wp;tAvesdsqqetoeAveTot=0._wp;tAvesdqetoeAveTot=0._wp
        tAvesdXAveTot=0._wp;tAvesdXSqAveTot=0._wp
        tAvesdxyTot=0._wp;tAvesdxxyyTot=0._wp;tAvesdyyzzTot=0._wp
        sdtAvesqqetoeAveTot=0._wp;sdtAveqetoeAveTot=0._wp;sdtAvetauxyTot=0._wp
        sdtAvetauxxyyTot=0._wp;sdtAvetauyyzzTot=0._wp
        sdtAveXAveTot=0._wp;sdtAveXSqAveTot=0._wp
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
          tAveXbbAveTot=0._wp;tAveXbbSqAveTot=0._wp
          tAvesdXbbAveTot=0._wp;tAvesdXbbSqAveTot=0._wp
          sdtAveXbbAveTot=0._wp;sdtAveXbbSqAveTot=0._wp
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
    use :: inp_dlt, only: npchain,nchain,nbead,nseg,tss,trst,lambda,Wi,dt,&
                          ntime,qmax,nseg_ar,Pe,ntotang,iflow,StrCalc,tplg&
                          &y,Na,cosThCalc,CoM,CoHR,RgCalc,nbeadx3,AveIterC&
                          &alc,DecompMeth,initmode,cosmode

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
    real(wp) :: XAveTot,sdXAveTot,XSqAveTot,sdXSqAveTot
    real(wp) :: XbbAveTot,sdXbbAveTot,XbbSqAveTot,sdXbbSqAveTot
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
5   format(f8.2,1x,e11.3,1x,3(f20.7,2x))
7   format(i4,1x,e11.3,1x,3(f14.7,2x))
8   format(i4,1x,e11.3,1x,2(f14.7,2x))

    if (StrCalc) then
      call conf_anlzr(q,rvmrc,Fphi,nseg_bb,id,time,idt,iPe)
      ! Reduction of terms for stress and relative extension
      call MPI_Reduce(sqqsprAve,sqqsprAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sqqetoeAve,sqqetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(qetoeAve,qetoeAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(XAve,XAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(XSqAve,XSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (tplgy == 'Comb') then
        call MPI_Reduce(sqqee_ar,sqqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(qee_ar,qee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(XbbAve,XbbAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(XbbSqAve,XbbSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
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
      call MPI_Reduce(sdXAve,sdXAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(sdXSqAve,sdXSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if (tplgy == 'Comb') then
        call MPI_Reduce(sdsqqee_ar,sdsqqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdqee_ar,sdqee_art,Na,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdXbbAve,sdXbbAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdXbbSqAve,sdXbbSqAveTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
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
          if (cosmode == 'dot') then
            cosTh=dot_product(qi,qj)
          else 
            cosTh=dot_product(qi,qj)/(qi_mag*qj_mag)
          end if
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
          ! DcmAve=DcmAve+dot(rcmP,rcmP)
          ! sdDcmAve=sdDcmAve+(dot(rcmP,rcmP))**2                       
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
          ! DchrAve=DchrAve+dot(rchrP,rchrP)                      
          ! sdDchrAve=sdDchrAve+(dot(rchrP,rchrP))**2                       
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
          ! call gemv(Bmat,qP,rvmrcP)
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
        sdXAveTot=sqrt(abs(sdXAveTot-XAveTot**2)/(nchain-1))
        sdXSqAveTot=sqrt(abs(sdXSqAveTot-XSqAveTot**2)/(nchain-1))
        sdxxTot=sqrt(abs(sdxxTot-tauxxTot**2)/(nchain-1))
        sdxyTot=sqrt(abs(sdxyTot-tauxyTot**2)/(nchain-1))
        sdxxyyTot=sqrt(abs(sdxxyyTot-tauxxyyTot**2)/(nchain-1))
        sdyyzzTot=sqrt(abs(sdyyzzTot-tauyyzzTot**2)/(nchain-1))
        if (initmode == 'rst') then
          write(u31,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqsprAveTot,sdsqqsprAveTot
          write(u33,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqetoeAveTot,sdsqqetoeAveTot
          write(u36,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),XSqAveTot,sdXSqAveTot
          select case (tplgy)
            case ('Linear')
              write(u7,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qetoeAveTot/(qmax*nseg),&
                         sdqetoeAveTot/(qmax*nseg)
              write(u35,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),XAveTot/(qmax*nseg),&
                          sdXAveTot/(qmax*nseg)
            case ('Comb')
              write(u7,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qetoeAveTot/(qmax*nseg_bb),&
                         sdqetoeAveTot/(qmax*nseg_bb)
              write(u35,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),XAveTot/(qmax*nseg_bb),&
                          sdXAveTot/(qmax*nseg_bb)
              do iarm=1, Na
                sdsqqee_art(iarm)=sqrt(abs(sdsqqee_art(iarm)-sqqee_art(iarm)**2)/(nchain-1))
                sdqee_art(iarm)=sqrt(abs(sdqee_art(iarm)-qee_art(iarm)**2)/(nchain-1))
                write(uarm(iarm),2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),sqqee_art(iarm),&
                                 sdsqqee_art(iarm)
                write(uarm(Na+iarm),2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),qee_art(iarm)/&
                                 (qmax*nseg_ar),sdqee_art(iarm)/(qmax*nseg_ar)
              end do
              sdXbbAveTot=sqrt(abs(sdXbbAveTot-XbbAveTot**2)/(nchain-1))
              sdXbbSqAveTot=sqrt(abs(sdXbbSqAveTot-XbbSqAveTot**2)/(nchain-1))
              write(u39,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),XbbAveTot/(qmax*nseg_bb),&
                           sdXbbAveTot/(qmax*nseg_bb)
              write(u40,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),XbbSqAveTot,sdXbbSqAveTot
          end select
        else ! initmode == 'st'
          write(u31,2) Wi(iPe),dt(iPe,idt),time,sqqsprAveTot,sdsqqsprAveTot
          write(u33,2) Wi(iPe),dt(iPe,idt),time,sqqetoeAveTot,sdsqqetoeAveTot
          write(u36,2) Wi(iPe),dt(iPe,idt),time,XSqAveTot,sdXSqAveTot
          select case (tplgy)
            case ('Linear')
              write(u7,2) Wi(iPe),dt(iPe,idt),time,qetoeAveTot/(qmax*nseg),&
                         sdqetoeAveTot/(qmax*nseg)
              write(u35,2) Wi(iPe),dt(iPe,idt),time,XAveTot/(qmax*nseg),sdXAveTot/(qmax*nseg)
            case ('Comb')
              write(u7,2) Wi(iPe),dt(iPe,idt),time,qetoeAveTot/(qmax*nseg_bb),&
                         sdqetoeAveTot/(qmax*nseg_bb)
              write(u35,2) Wi(iPe),dt(iPe,idt),time,XAveTot/(qmax*nseg_bb),sdXAveTot/(qmax*nseg_bb)
              do iarm=1, Na
                sdsqqee_art(iarm)=sqrt(abs(sdsqqee_art(iarm)-sqqee_art(iarm)**2)/(nchain-1))
                sdqee_art(iarm)=sqrt(abs(sdqee_art(iarm)-qee_art(iarm)**2)/(nchain-1))
                write(uarm(iarm),2) Wi(iPe),dt(iPe,idt),time,sqqee_art(iarm),sdsqqee_art(iarm)
                write(uarm(Na+iarm),2) Wi(iPe),dt(iPe,idt),time,qee_art(iarm)/(qmax*nseg_ar),&
                                    sdqee_art(iarm)/(qmax*nseg_ar)
              end do
              sdXbbAveTot=sqrt(abs(sdXbbAveTot-XbbAveTot**2)/(nchain-1))
              sdXbbSqAveTot=sqrt(abs(sdXbbSqAveTot-XbbSqAveTot**2)/(nchain-1))
              write(u39,2) Wi(iPe),dt(iPe,idt),time,XbbAveTot/(qmax*nseg_bb),sdXbbAveTot/(qmax*nseg_bb)
              write(u40,2) Wi(iPe),dt(iPe,idt),time,XbbSqAveTot,sdXbbSqAveTot
          end select
        end if
        if (iflow /= 1) then
          if (iflow < 3) then
            if (initmode == 'rst') then
              write(u8,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)
            else
              write(u8,2) Wi(iPe),dt(iPe,idt),time,-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)
            end if
          else
            if (initmode == 'rst') then 
              write(u14,3) Pe(iPe)*(time+trst*lambda),dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe)
            else
              write(u14,3) Pe(iPe)*time,dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe)
            end if
          end if
        end if
      end if ! StrCalc
      if (cosThCalc) then
        cosThAveTot=cosThAveTot/ntotang
        sdcosThAveTot=sdcosThAveTot/ntotang
        sdcosThAveTot=sqrt(abs(sdcosThAveTot-cosThAveTot*cosThAveTot)/(ntotang-1))
        if (initmode == 'rst') then
          write(u28,2) Wi(iPe),dt(iPe,idt),(time+trst*lambda),cosThAveTot,sdcosThAveTot
        else
          write(u28,2) Wi(iPe),dt(iPe,idt),time,cosThAveTot,sdcosThAveTot
        end if
      end if ! cosThCalc
      if (AveIterCalc) then
        if (DecompMeth == 'Lanczos') then
          mAveTot=mAveTot/nchain
          sdmAveTot=sdmAveTot/nchain
          sdmAveTot=sqrt(abs(sdmAveTot-mAveTot**2)/(nchain-1))
          if (initmode == 'rst') then
            write(u20,4) iflow,Wi(iPe),dt(iPe,idt),(time+trst*lambda),mAveTot,sdmAveTot
          else
            write(u20,4) iflow,Wi(iPe),dt(iPe,idt),time,mAveTot,sdmAveTot
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
            write(u20,4) iflow,Wi(iPe),dt(iPe,idt),(time+trst*lambda),LAveTot,sdLAveTot
          else
            write(u20,4) iflow,Wi(iPe),dt(iPe,idt),time,LAveTot,sdLAveTot
          end if
          tAveLAveTot=tAveLAveTot+LAveTot
          tAvesdLAveTot=tAvesdLAveTot+sdLAveTot
          sdtAveLAveTot=sdtAveLAveTot+LAveTot*LAveTot
        end if
      end if ! AveIterCalc
      if (time >= tss*lambda) then
        jcount=jcount+1
        if (StrCalc) then                      
          ! Terms for stress and relative extension
          tAvesqqsprAveTot=tAvesqqsprAveTot+sqqsprAveTot
          tAvesqqetoeAveTot=tAvesqqetoeAveTot+sqqetoeAveTot
          tAveqetoeAveTot=tAveqetoeAveTot+qetoeAveTot
          tAveXAveTot=tAveXAveTot+XAveTot
          tAveXSqAveTot=tAveXSqAveTot+XSqAveTot
          tAvetauxyTot=tAvetauxyTot+tauxyTot
          tAvetauxxyyTot=tAvetauxxyyTot+tauxxyyTot
          tAvetauyyzzTot=tAvetauyyzzTot+tauyyzzTot
          ! Terms for standard deviation of stress and relative extension
          tAvesdsqqsprAveTot=tAvesdsqqsprAveTot+sdsqqsprAveTot
          tAvesdsqqetoeAveTot=tAvesdsqqetoeAveTot+sdsqqetoeAveTot
          tAvesdqetoeAveTot=tAvesdqetoeAveTot+sdqetoeAveTot
          tAvesdXAveTot=tAvesdXAveTot+sdXAveTot
          tAvesdXSqAveTot=tAvesdXSqAveTot+sdXSqAveTot
          tAvesdxyTot=tAvesdxyTot+sdxyTot
          tAvesdxxyyTot=tAvesdxxyyTot+sdxxyyTot
          tAvesdyyzzTot=tAvesdyyzzTot+sdyyzzTot
          ! Standard deviation of terms for stress and relative extension
          sdtAvesqqsprAveTot=sdtAvesqqsprAveTot+sqqsprAveTot*sqqsprAveTot
          sdtAvesqqetoeAveTot=sdtAvesqqetoeAveTot+sqqetoeAveTot*sqqetoeAveTot
          sdtAveqetoeAveTot=sdtAveqetoeAveTot+qetoeAveTot*qetoeAveTot
          sdtAveXAveTot=sdtAveXAveTot+XAveTot*XAveTot
          sdtAveXSqAveTot=sdtAveXSqAveTot+XSqAveTot*XSqAveTot
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
            tAveXbbAveTot=tAveXbbAveTot+XbbAveTot
            tAveXbbSqAveTot=tAveXbbSqAveTot+XbbSqAveTot
            tAvesdXbbAveTot=tAvesdXbbAveTot+sdXbbAveTot
            tAvesdXbbSqAveTot=tAvesdXbbSqAveTot+sdXbbSqAveTot
            sdtAveXbbAveTot=sdtAveXbbAveTot+XbbAveTot*XbbAveTot
            sdtAveXbbSqAveTot=sdtAveXbbSqAveTot+XbbSqAveTot*XbbSqAveTot
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
            tAveXAveTot=tAveXAveTot/jcount
            tAveXSqAveTot=tAveXSqAveTot/jcount
            tAvetauxyTot=tAvetauxyTot/jcount
            tAvetauxxyyTot=tAvetauxxyyTot/jcount
            tAvetauyyzzTot=tAvetauyyzzTot/jcount
            ! Terms for standard deviation of stress and relative extension
            tAvesdsqqsprAveTot=tAvesdsqqsprAveTot/jcount
            tAvesdsqqetoeAveTot=tAvesdsqqetoeAveTot/jcount
            tAvesdqetoeAveTot=tAvesdqetoeAveTot/jcount
            tAvesdXAveTot=tAvesdXAveTot/jcount
            tAvesdXSqAveTot=tAvesdXSqAveTot/jcount
            tAvesdxyTot=tAvesdxyTot/jcount
            tAvesdxxyyTot=tAvesdxxyyTot/jcount
            tAvesdyyzzTot=tAvesdyyzzTot/jcount
            ! Standard deviation of terms for stress and relative extension
            sdtAvesqqsprAveTot=sdtAvesqqsprAveTot/jcount
            sdtAvesqqetoeAveTot=sdtAvesqqetoeAveTot/jcount
            sdtAveqetoeAveTot=sdtAveqetoeAveTot/jcount
            sdtAveXAveTot=sdtAveXAveTot/jcount
            sdtAveXSqAveTot=sdtAveXSqAveTot/jcount
            sdtAvetauxyTot=sdtAvetauxyTot/jcount
            sdtAvetauxxyyTot=sdtAvetauxxyyTot/jcount
            sdtAvetauyyzzTot=sdtAvetauyyzzTot/jcount
            sdtAvesqqsprAveTot=sqrt(abs(sdtAvesqqsprAveTot-tAvesqqsprAveTot**2)/(jcount-1))
            sdtAvesqqetoeAveTot=sqrt(abs(sdtAvesqqetoeAveTot-tAvesqqetoeAveTot**2)/(jcount-1))
            sdtAveqetoeAveTot=sqrt(abs(sdtAveqetoeAveTot-tAveqetoeAveTot**2)/(jcount-1))
            sdtAveXAveTot=sqrt(abs(sdtAveXAveTot-tAveXAveTot**2)/(jcount-1))
            sdtAveXSqAveTot=sqrt(abs(sdtAveXSqAveTot-tAveXSqAveTot**2)/(jcount-1))
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
              tAveXbbAveTot=tAveXbbAveTot/jcount
              tAveXbbSqAveTot=tAveXbbSqAveTot/jcount
              tAvesdXbbAveTot=tAvesdXbbAveTot/jcount
              tAvesdXbbSqAveTot=tAvesdXbbSqAveTot/jcount
              sdtAveXbbAveTot=sdtAveXbbAveTot/jcount
              sdtAveXbbSqAveTot=sdtAveXbbSqAveTot/jcount
              sdtAveXbbAveTot=sqrt(abs(sdtAveXbbAveTot-tAveXbbAveTot**2)/(jcount-1))
              sdtAveXbbSqAveTot=sqrt(abs(sdtAveXbbSqAveTot-tAveXbbSqAveTot**2)/(jcount-1))              
            end if
            sdtAvetauxyTot=sqrt(abs(sdtAvetauxyTot-tAvetauxyTot**2)/(jcount-1))
            sdtAvetauxxyyTot=sqrt(abs(sdtAvetauxxyyTot-tAvetauxxyyTot**2)/(jcount-1))
            sdtAvetauyyzzTot=sqrt(abs(sdtAvetauyyzzTot-tAvetauyyzzTot**2)/(jcount-1))
            write(u32,5) Wi(iPe),dt(iPe,idt),tAvesqqsprAveTot,tAvesdsqqsprAveTot,sdtAvesqqsprAveTot
            write(u30,5) Wi(iPe),dt(iPe,idt),tAvesqqetoeAveTot,tAvesdsqqetoeAveTot,sdtAvesqqetoeAveTot
            write(u38,5) Wi(iPe),dt(iPe,idt),tAveXSqAveTot,tAvesdXSqAveTot,sdtAveXSqAveTot
            select case (tplgy)
              case ('Linear')
                write(u9,5) Wi(iPe),dt(iPe,idt),tAveqetoeAveTot/(qmax*nseg),tAvesdqetoeAveTot/&
                          (qmax*nseg),sdtAveqetoeAveTot/(qmax*nseg)
                write(u37,5) Wi(iPe),dt(iPe,idt),tAveXAveTot/(qmax*nseg),tAvesdXAveTot/(qmax*nseg),&
                            sdtAveXAveTot/(qmax*nseg)
              case ('Comb')
                write(u9,5) Wi(iPe),dt(iPe,idt),tAveqetoeAveTot/(qmax*nseg_bb),tAvesdqetoeAveTot/&
                          (qmax*nseg_bb),sdtAveqetoeAveTot/(qmax*nseg_bb)
                write(u37,5) Wi(iPe),dt(iPe,idt),tAveXAveTot/(qmax*nseg_bb),tAvesdXAveTot/(qmax*nseg_bb),&
                            sdtAveXAveTot/(qmax*nseg_bb)
                do iarm=1, Na
                  write(uarm(2*Na+iarm),5) Wi(iPe),dt(iPe,idt),tAvsqqee_art(iarm),tAvsdsqqee_art(iarm),&
                                        sdtAvsqqee_art(iarm)
                  write(uarm(3*Na+iarm),5) Wi(iPe),dt(iPe,idt),tAvqee_art(iarm)/(qmax*nseg_ar),&
                                        tAvsdqee_art(iarm)/(qmax*nseg_ar),sdtAvqee_art(iarm)/&
                                        (qmax*nseg_ar)
                end do
                write(u41,5) Wi(iPe),dt(iPe,idt),tAveXbbAveTot/(qmax*nseg_bb),tAvesdXbbAveTot/(qmax*nseg_bb),&
                            sdtAveXbbAveTot/(qmax*nseg_bb)
                write(u42,5) Wi(iPe),dt(iPe,idt),tAveXbbSqAveTot,tAvesdXbbSqAveTot,sdtAveXbbSqAveTot
            end select
            if (iflow /= 1) then
              write(u10,5) Wi(iPe),dt(iPe,idt),-tAvetauxyTot/Pe(iPe),tAvesdxyTot/Pe(iPe),&
                          sdtAvetauxyTot/Pe(iPe)
              write(u11,5) Wi(iPe),dt(iPe,idt),-tAvetauxxyyTot/Pe(iPe)**2,tAvesdxxyyTot/Pe(iPe)**2,&
                          sdtAvetauxxyyTot/Pe(iPe)**2
              write(u12,5) Wi(iPe),dt(iPe,idt),-tAvetauyyzzTot/Pe(iPe)**2,tAvesdyyzzTot/Pe(iPe)**2,&
                          sdtAvetauyyzzTot/Pe(iPe)**2
              if (iflow >= 3) then ! For Elongational Flow
                write(u13,5) Wi(iPe),dt(iPe,idt),-tAvetauxxyyTot/Pe(iPe),tAvesdxxyyTot/Pe(iPe),&
                            sdtAvetauxxyyTot/Pe(iPe) 
              end if
            end if ! iflow /= 1
          end if ! StrCalc
          if (CoM) then
            tAveDcmAveTot=tAveDcmAveTot/jcount
            tAvesdDcmAveTot=tAvesdDcmAveTot/jcount
            sdtAveDcmAveTot=sdtAveDcmAveTot/jcount
            sdtAveDcmAveTot=sqrt(abs(sdtAveDcmAveTot-tAveDcmAveTot**2)/(jcount-1))
            write(u15,7) nbead,dt(iPe,idt),tAveDcmAveTot,tAvesdDcmAveTot,sdtAveDcmAveTot 
          end if
          if (CoHR) then
            tAveDchrAveTot=tAveDchrAveTot/jcount
            tAvesdDchrAveTot=tAvesdDchrAveTot/jcount
            sdtAveDchrAveTot=sdtAveDchrAveTot/jcount
            sdtAveDchrAveTot=sqrt(abs(sdtAveDchrAveTot-tAveDchrAveTot**2)/(jcount-1))
            write(u16,7) nbead,dt(iPe,idt),tAveDchrAveTot,tAvesdDchrAveTot,sdtAveDchrAveTot
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
            write(u17,7) nbead,dt(iPe,idt),tAveRgSqAveTot,tAvesdRgSqAveTot,sdtAveRgSqAveTot 
            write(u18,7) nbead,dt(iPe,idt),tAveAspherAveTot,tAvesdAspherAveTot,sdtAveAspherAveTot
          end if ! RgCalc
          if (cosThCalc) then
            ! For cosTh:
            tAvecosTh=tAvecosTh/jcount
            tAvesdcosTh=tAvesdcosTh/jcount
            sdtAvecosTh=sdtAvecosTh/jcount
            sdtAvecosTh=sqrt(abs(sdtAvecosTh-tAvecosTh*tAvecosTh)/(jcount-1))
            write(u29,7) nbead,dt(iPe,idt),tAvecosTh,tAvesdcosTh,sdtAvecosTh
          end if ! cosThCalc             
        else ! so, jcount=0
          if (StrCalc) then          
            write(u9,3) Wi(iPe),dt(iPe,idt),qetoeAveTot/(qmax*nseg),sdqetoeAveTot/(qmax*nseg)
            if (iflow /= 1) then
              write(u10,3) Wi(iPe),dt(iPe,idt),-tauxyTot/Pe(iPe),sdxyTot/Pe(iPe)           
              write(u11,3) Wi(iPe),dt(iPe,idt),-tauxxyyTot/(Pe(iPe)**2),sdxxyyTot/(Pe(iPe)**2)
              write(u12,3) Wi(iPe),dt(iPe,idt),-tauyyzzTot/(Pe(iPe)**2),sdyyzzTot/(Pe(iPe)**2)
              if (iflow >= 3) then ! For Elongational Flow
                write(u13,3) Wi(iPe),dt(iPe,idt),-tauxxyyTot/Pe(iPe),sdxxyyTot/Pe(iPe) 
              end if
            end if ! iflow /= 1
          end if ! StrCalc
          if (CoM) write(u15,8) nbead,dt(iPe,idt),DcmAveTot,sdDcmAveTot 
          if (CoHR) write(u16,8) nbead,dt(iPe,idt),DchrAveTot,sdDchrAveTot
          if (RgCalc) then
            write(u17,8) nbead,dt(iPe,idt),RgSqAveTot,sdRgSqAveTot
            write(u18,8) nbead,dt(iPe,idt),AspherAveTot,sdAspherAveTot
          end if              
        end if ! jcount /= 0
        if (AveIterCalc) then
          if (kcount /= 0) then
            if (DecompMeth == 'Lanczos') then
              tAvemAveTot=tAvemAveTot/kcount
              tAvesdmAveTot=tAvesdmAveTot/kcount
              sdtAvemAveTot=sdtAvemAveTot/kcount
              sdtAvemAveTot=sqrt(abs(sdtAvemAveTot-tAvemAveTot**2)/(kcount-1))
              write(u19,7) nbead,dt(iPe,idt),tAvemAveTot,tAvesdmAveTot,sdtAvemAveTot
            elseif (DecompMeth == 'Chebyshev') then
              tAveLAveTot=tAveLAveTot/kcount
              tAvesdLAveTot=tAvesdLAveTot/kcount
              sdtAveLAveTot=sdtAveLAveTot/kcount
              sdtAveLAveTot=sqrt(abs(sdtAveLAveTot-tAveLAveTot**2)/(kcount-1))
              write(u19,7) nbead,dt(iPe,idt),tAveLAveTot,tAvesdLAveTot,sdtAveLAveTot
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
  subroutine conf_anlzr(q,rvmrc,Fphi,nseg_bb,id,time,idt,iPe)
  
    use :: inp_dlt, only: nseg,npchain,nchain,ForceLaw,tplgy,nseg_ar,applFext,&
                          Wi,dt,trst,lambda,qmax,initmode,indvlext
  
    integer,intent(in) :: nseg_bb,id,idt,iPe
    real(wp),intent(in) :: time
    real(wp),intent(in),target :: q(:,:),rvmrc(:,:),Fphi(:,:)
    real(wp),dimension(:),pointer :: qP,qPx,qPy,qPz,RPx,RPy,RPz
    real(wp),dimension(:),pointer :: FphiPx,FphiPy,FphiPz
    real(wp) :: sqqeetmp,qeetmp,qetoex,qetoey,qetoez,X,Xbb
    real(wp) :: sqqspr,sqqetoe,qetoe,txx,txy,tyy,tzz,txxyy,tyyzz
    real(wp),allocatable :: qee_artmp(:,:)
    integer :: ichain,iarm,os
  
1   format(f8.2,1x,e11.3,1x,f14.7,2x,f20.8)

    sqqetoeAve=0._wp;qetoeAve=0._wp
    sqqsprAve=0._wp;XAve=0._wp;XSqAve=0._wp
    sdsqqetoeAve=0._wp;sdqetoeAve=0._wp
    sdsqqsprAve=0._wp;sdXAve=0._wp;sdXSqAve=0._wp
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
          Xbb=maxval(RPx(1:nseg_bb+1))-minval(RPx(1:nseg_bb+1))
          if (indvlext) then
            if (initmode == 'rst') then
              write(uch1(ichain),1) Wi(iPe),dt(iPe,idt),(time+trst*lambda),Xbb/(qmax*nseg_bb)
            else
              write(uch1(ichain),1) Wi(iPe),dt(iPe,idt),time,Xbb/(qmax*nseg_bb)
            end if
          end if
          XbbAve=XbbAve+Xbb
          XbbSqAve=XbbSqAve+Xbb**2
          sdXbbAve=sdXbbAve+Xbb**2
          sdXbbSqAve=sdXbbSqAve+Xbb**4
      end select
      sqqspr=dot_product(q(:,ichain),q(:,ichain))
      sqqetoe=qetoex**2+qetoey**2+qetoez**2
      if (applFext) then
        qetoe=qetoex
      else
        qetoe=sqrt(sqqetoe)
      end if
      X=maxval(RPx)-minval(RPx)
      if (indvlext) then
        if (initmode == 'rst') then
          write(uch(ichain),1) Wi(iPe),dt(iPe,idt),(time+trst*lambda),X/(qmax*nseg_bb)
        else
          write(uch(ichain),1) Wi(iPe),dt(iPe,idt),time,X/(qmax*nseg_bb)
        end if
      end if
      sqqsprAve=sqqsprAve+sqqspr
      sqqetoeAve=sqqetoeAve+sqqetoe
      qetoeAve=qetoeAve+qetoe
      XAve=XAve+X
      XSqAve=XSqAve+X**2
      sdsqqsprAve=sdsqqsprAve+sqqspr*sqqspr
      sdsqqetoeAve=sdsqqetoeAve+sqqetoe*sqqetoe
      sdqetoeAve=sdqetoeAve+qetoe*qetoe
      sdXAve=sdXAve+X**2
      sdXSqAve=sdXSqAve+X**4
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
    XAve=XAve/nchain
    XSqAve=XSqAve/nchain
    sdsqqsprAve=sdsqqsprAve/(nchain*nseg)
    sdsqqetoeAve=sdsqqetoeAve/nchain
    sdqetoeAve=sdqetoeAve/nchain
    sdXAve=sdXAve/nchain
    sdXSqAve=sdXSqAve/nchain
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
      XbbAve=XbbAve/nchain
      XbbSqAve=XbbSqAve/nchain
      sdXbbAve=sdXbbAve/nchain
      sdXbbSqAve=sdXbbSqAve/nchain      
    end if
   
  end subroutine conf_anlzr

  ! This routine is written by Vydia Venkataramani and modified by Amir Saadat.
  ! It is for sorting the ensemble to different configurational classes:
  ! Description of diffrent configuration types:
  ! Type 1   :   Folds(2/3 or 2/5) 
  ! Type 2   :   Half-Dumbells 
  ! Type 3   :   Kinks
  ! Type 4   :   Dumbells 
  ! Type 5   :   Coils
  ! Type 6   :   Extended
  subroutine conf_sort(q,rvmrc,nseg,nbead,cnf_tp)

    use :: inp_dlt, only: npchain,qmax,residx
    use :: arry_mod, only: print_vector,print_matrix

    integer,intent(inout) :: nseg,nbead,cnf_tp(:)
    real(wp),intent(in),target :: q(:,:),rvmrc(:,:)
    real(wp),pointer :: qPx(:),qPy(:),qPz(:)
    real(wp),pointer :: RPx(:),RPy(:),RPz(:)
!    real(wp) :: pconfig(6)
!    integer :: iconfig(6)
    real(wp) :: rmax,xi,yi,zi,xj,yj,zj,dist,xmin,ymin,zmin
    real(wp) :: xmax,ymax,zmax,xr,yr,zr,xcap,ycap,zcap,rmin
    real(wp) :: xproj,yproj,zproj,xpres,xnext,dia,radius
    real(wp) :: bead_left,bead_right,avg,avgsum
    integer :: i,ntype,maxlength,ichain,ibead,jbead,ilen
    integer :: min,max,maxbright,maxbright1,low,last,info,info1
    integer :: k,ij,j,ibcount,sumbright,lcount,left,markerc,ik
    integer :: lend,marker1,marker2,markerl,markerr,icount

    if (.not.allocated(x)) allocate(x(nbead))
    if (.not.allocated(y)) allocate(y(nbead))
    if (.not.allocated(z)) allocate(z(nbead))
    if (.not.allocated(r)) allocate(r(nbead))
    maxlength=int(qmax*nseg*residx)+1
    if (.not.allocated(ib))  allocate(ib(maxlength))
    if (.not.allocated(ib1)) allocate(ib1(maxlength))
    if (.not.allocated(icount3)) allocate(icount3(maxlength))

    do i=1, 6
!      iconfig(i)=0
!      pconfig(i)=0._wp
    end do

    do ichain=1, npchain
      qPx => q(1:3*nseg-2:3,ichain)
      qPy => q(2:3*nseg-1:3,ichain) 
      qPz => q(3:3*nseg:3,ichain)
      RPx => rvmrc(1:3*nbead-2:3,ichain)
      RPy => rvmrc(2:3*nbead-1:3,ichain) 
      RPz => rvmrc(3:3*nbead:3,ichain)
      ntype=0
      ! Calculating the position of each bead
!      x(1)=0._wp
!      y(1)=0._wp
!      z(1)=0._wp
!      do ibead=2, nbead
!        x(ibead)=x(ibead-1)+qPx(ibead-1)
!        y(ibead)=y(ibead-1)+qPy(ibead-1)
!        z(ibead)=z(ibead-1)+qPz(ibead-1)
!      end do          

      ib=0
      ib1=0
!      ! Calculating the largest distance between the beads
!      rmax=0._wp
!      do ibead=1, nbead-1
!        xi=x(ibead)
!        yi=y(ibead)
!        zi=z(ibead)
!        do jbead=1, nbead
!          xj=x(jbead)
!          yj=y(jbead)
!          zj=z(jbead)
!          dist=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
!
!          if (dist > rmax) then
!            min=ibead
!            max=jbead
!            rmax=dist
!          end if
!        end do ! jbead
!      end do ! ibead
!
!      xmin=x(min)
!      ymin=y(min)
!      zmin=z(min)
!
!      xmax=x(max)
!      ymax=y(max)
!      zmax=z(max)
!
!      ! Moving beads w.r.t min
!      do ibead=1, nbead
!         x(ibead)=x(ibead)-xmin
!         y(ibead)=y(ibead)-ymin
!         z(ibead)=z(ibead)-zmin
!      end do

      ! Calculating the largest distance between the beads
      rmax=0._wp
      do ibead=1, nbead-1
        xi=RPx(ibead)
        yi=RPy(ibead)
        zi=RPz(ibead)
        do jbead=1, nbead
          xj=RPx(jbead)
          yj=RPy(jbead)
          zj=RPz(jbead)
          dist=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

          if (dist > rmax) then
            min=ibead
            max=jbead
            rmax=dist
          end if
        end do ! jbead
      end do ! ibead

      xmin=RPx(min)
      ymin=RPy(min)
      zmin=RPz(min)

      xmax=RPx(max)
      ymax=RPy(max)
      zmax=RPz(max)

      do ibead=1, nbead
        x(ibead)=RPx(ibead)-xmin
        y(ibead)=RPy(ibead)-ymin
        z(ibead)=RPz(ibead)-zmin
      end do
!call print_vector(x,'x2')
!call print_vector(y,'y2')
!call print_vector(z,'z2')

      ! Calculating the unit vector along the longest length
      xr=xmax-xmin
      yr=ymax-ymin
      zr=zmax-zmin

      xcap=xr/rmax
      ycap=yr/rmax
      zcap=zr/rmax

      ! Finding the projection of each bead onto the longest length
      rmin=0._wp
      do ibead=1, nbead
        r(ibead)=0._wp
        xproj=x(ibead)*xcap
        yproj=y(ibead)*ycap
        zproj=z(ibead)*zcap
        r(ibead)=(xproj+yproj+zproj)*residx
        if (r(ibead) < rmin) then
          rmin=r(ibead)
        end if
      end do

      ! Repositioning the beads such that the minimum r = 0
      do ibead=1, nbead
        r(ibead)=r(ibead)-rmin
      end do

      ! Accounting for bond overlaps
      maxbright=int(rmax*residx)+1
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

      do ilen=low, last
        ib(ilen)=ib(ilen)+1
      end do

      do ibead=2, nbead-1
        xpres=r(ibead)
        xnext=r(ibead+1)
        if (xpres < xnext) then
          info1 = 0
        else
          info1 = 1
        end if

        if ((info == 0).and.(info1 == 0)) then
          low=int(xpres)+2
          last=int(xnext)+1
        end if

        if ((info == 0).and.(info1 == 1)) then
          low=int(xnext)+1
          last=int(xpres)+1
        end if

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

      ! Accounting for bead overlap

      dia=1._wp*residx ! parameter set based on your assumption
!      dia=qmax*residx ! parameter set based on your assumption
      radius=dia/2
           
      do i=1, maxbright
        icount3(i)=0._wp
      end do
      
      do ibead=1, nbead
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

      do ilen=1, maxbright
        if (icount3(ilen) == 1) then
          ib(ilen)=ib(ilen)-1
        end if
      end do

      ! Averaging over configurations
!call print_vector(ib(1:maxbright),'ib1')

      if (mod(maxbright,residx) == 0) then
        maxbright1=maxbright/residx
      else
        maxbright1=int(maxbright/residx)+1
      end if

      k=0
      ij=1
      do i=1, maxbright1
        j=1
        ibcount=0
        do while ((j <= residx).and.(ij <= maxbright))
          ibcount=ibcount+ib(ij)
          j=j+1
          ij=ij+1
        end do
        if (j == residx) then
          avg=real(ibcount,kind=wp)/residx
        else
          avg=real(ibcount,kind=wp)/(j-1)
        end if

        if ((real(int(avg),kind=wp)+0.5_wp) < avg) THEN
          ib1(i)=int(avg)+1
        else
          ib1(i)=int(avg)
        end if
        k=k+1
      end do

      do i=1, maxbright
        ib(i)=0
      end do
      maxbright=maxbright1
      do i=1, maxbright
        ib(i)=ib1(i)
      end do

!call print_vector(ib(1:maxbright),'ib2')
      ! Checking for the minimum resolution

!      if (rmax <= 15._wp) then
      if (rmax <= 0.05_wp*qmax*nseg) then
!        iconfig(5)=iconfig(5)+1
        ntype=5
      else

        ! **Assigning for 2/3, 2/5 folds, coil and half dumbell**

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
           
          icount = 0   
          sumbright = 0
          do while (ib(i) > 1)
            icount=icount+1
            sumbright=sumbright+ib(i)
            if (j == 1) then
              i=i-1
            else
              i=i+1
            end if  
          end do
          avgsum=sumbright/(1._wp*icount)
          lcount=0
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

          if (markerc == 1) then
!            iconfig(3)=iconfig(3)+1
            ntype=3
          else
            icount=maxbright-lcount
            if ((0.25_wp <= ((1._wp*icount)/maxbright))) then
!              iconfig(1)=iconfig(1)+1
              ntype=1
            end if
            if (((1._wp*icount)/maxbright) < 0.25_wp) then
!              iconfig(2)=iconfig(2)+1
              ntype=2
            end if
          end if
        end if

        ! **Assigning for kink or fully extended**
              
        if ((ib(1) == 1).and.(ib(maxbright) == 1)) THEN
          marker1=0
          do i=2, (maxbright-1)
            if (ib(i) == 1) then
              marker1=marker1+1
            end if
          end do
          if (marker1 == (maxbright-2)) then
!            iconfig(6)=iconfig(6)+1
            ntype=6
          else
!            iconfig(3)=iconfig(3)+1
            ntype=3
          end if
        end if

        ! **Assigning for coil or dumbell**
              
        if ((ib(1) > 1).and.(ib(maxbright) > 1)) THEN
          marker2=0
          markerl=2
          do while ((ib(markerl) > 1).and.(markerl < maxbright))
            markerl=markerl+1
          end do
          markerr=maxbright-1
          do while ((ib(markerr) > 1).and.(markerr > 1))
            markerr=markerr-1
          end do
          if ((markerl == (maxbright)).and.(markerr == 1)) then
!            iconfig(5)=iconfig(5)+1
            ntype=5
          else
            do i=markerl, markerr
              if (ib(i) > 1) then
                marker2=1
              end if
            end do
            if (marker2 /= 1) then
!              iconfig(4)=iconfig(4)+1
              ntype=4
            else
!              iconfig(5)=iconfig(5)+1
              ntype=5
            end if
          end if           
        end if
              
      end if

      cnf_tp(ichain)=ntype
                
    end do ! ichain

!    do i=1, 6
!      pconfig(i)=real(iconfig(i),kind=wp)/nchain
!    end do

  end subroutine conf_sort

  subroutine del_pp(id)

    use :: inp_dlt, only: StrCalc,iflow,tplgy,Na,CoM,CoHR,RgCalc,AveIterCalc,&
                          cosThCalc,nchain,indvlext

    integer,intent(in) :: id
    integer :: iarm,ichain

    if (id == 0) then
      if (StrCalc) then
        close(u31);close(u33);close(u7);close(u35);close(u36)
        close(u30);close(u9);close(u32);close(u37);close(u38)
        if (iflow /= 1) then
          close(u8);close(u10);close(u11);close(u12)
          if (iflow >= 3) close(u13);close(u14) 
        end if
        select case (tplgy)
          case ('Linear')
          case ('Comb')
            do iarm=1, 4*Na
              close(uarm(iarm))
            end do
            deallocate(uarm)
            close(u39);close(u40);close(u41);close(u42)
        end select
      end if
      if (CoM) close(u15)
      if (CoHR) close(u16)
      if (RgCalc) close(u17);close(u18)
      if (AveIterCalc) close(u19);close(u20)
      if (cosThCalc) close(u28);close(u29)
    end if
    if (StrCalc.and.indvlext) then
      deallocate(uch)
      if (tplgy == 'Comb') deallocate(uch1)
    end if

  end subroutine del_pp

end module pp_mod
