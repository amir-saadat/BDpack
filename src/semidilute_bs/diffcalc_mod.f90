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
! MODULE: diffusion calculator
!
!> @author
!> Amir Saadat, The University of Tennessee-Knoxville, June 2014
!
! DESCRIPTION:
!> Calculates the diffusion tensor on CPU
!--------------------------------------------------------------------
module diffcalc_mod

  use :: prcn_mod
  use :: types
  use :: hi_mod, only: HI_a,HI_ax3,HI_ato3,HI_alpha,HI_alphato2,HI_alphato3,rc_D,HI_Vol,HI_c0,M1_c1,M1_c2,&
   M1_c3,M1_c4,M1_c5,M1_c6,M1_c7,M1_c8,M1_c9,M1_c10,M1_c11,M1_c12,M1_c13,M2_c1,M2_c2,M2_c3,M2_c4,Mstar_c1,&
   Mstar_c2,HI_kmax,hstar,unitDelta,Diff_tens,Dreal_vals,Dreal_cols,Dreal_rowInd,DF_self,DF_real,DF_recip,&
   Diff_tens_real,Diff_tens_recip,HIcalc_mode,Dreal_sparse_mode,point_D,list_DP,p_PME,maxNb_list_D,K_mesh,&
        F_mesh,Kto3,Kcto3,FFTStatus,FFTfwDescHand,FFTbwDescHand,r_str,c_str,b_splx2,M_spl,b_sply2,b_splz2,&
 InterpMethod,W_Lag,eikx,eiky,eikz,eikr,HI_kimax,m2_vec,F_meshPx,F_meshPy,F_meshPz,P_vals,P_cols,P_rowInd,&
 p_PMEto3,mpvecx,mpvecy,mpvecz,HI_ax2,rc_Dto2,kvecx,kvecy,kvecz,kiuppr,kiylowr,kizlowr
  use :: tmng_mod, only: PMEcount,et_PME,et_FFT,et_IFFT,et_SPR,et_INT,et_INF,et_R,et_K,et_EIKR,et_EW,tick,&
                         tock,et_DCR,et_DCK,doTiming
!$ use :: omp_lib

  implicit none

contains

  subroutine calcDiffTens_cpu(Rb,ntotbead,boxsize,boxorigin,itime)

    use :: arry_mod, only: print_vector,print_matrix
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m,delrx_L
    use :: hi_mod, only: PIx2,sqrtPI


    integer,intent(in) :: ntotbead
    real(wp),intent(in) :: Rb(:)
    real(wp),dimension(:,:),pointer :: Diff_tensP => null(),Diff_tens_recipP => null()
    integer :: ichain,jchain,ibead,ibeadub,jbead,iglob,jglob,offsetjch,offsetich
    integer :: iglobbead,jglobbead,iboxbead,kix,kiy,kiz
    integer :: kiymin,kizmin,ktot,tid,kiydev,kiyy,kikix,kikixy
    integer,parameter :: FORWARD=1,BACKWARD=-1
    real(wp),dimension(3) :: rij,rijn0,rijhat,kvec,ri,rin0,rjn0,rijhat0,boxsize,boxorigin,khat
    integer,dimension(3) :: orient,kivec,nmax
    integer,target :: nvec(3)
    integer,pointer :: n1,n2,n3
    real(wp) :: r,r2,kmag,kto2,rrhat(6),invk2,t0,t1,Coeff,kkhat(6),invr2
    integer(long) :: count0,count1
    integer :: itime

    nvec=0;

    if (HIcalc_mode == 'Ewald') then

      if (doTiming) call tick(count0)
      !------------------------------------------------------
      !>>> Constructing EXP(ik.r) for all beads and kvectors:
      !------------------------------------------------------
      if (doTiming) call tick(count1)
      do iboxbead=1, ntotbead
        iglob=(iboxbead-1)*3
        ri(1:3)=Rb(iglob+1:iglob+3)
        ! Calculating kix, kiy, kiz = 0, -1, and 1 explicitly:
        eikx(0,iboxbead)=(1.0_wp,0.0_wp)
        eiky(0,iboxbead)=(1.0_wp,0.0_wp)
        eikz(0,iboxbead)=(1.0_wp,0.0_wp)
        ! Note that having kind=wp is really important for having enough precision.
        eikx(1,iboxbead)=cmplx( cos(PIx2/boxsize(1)*ri(1)) , &
                                sin(PIx2/boxsize(1)*ri(1)),kind=wp)
        eiky(1,iboxbead)=cmplx( cos(PIx2/boxsize(2)*ri(2)) , &
                                sin(PIx2/boxsize(2)*ri(2)),kind=wp)
        eikz(1,iboxbead)=cmplx( cos(PIx2/boxsize(3)*ri(3)) , &
                                sin(PIx2/boxsize(3)*ri(3)),kind=wp)

        eiky(-1,iboxbead)=conjg(eiky(1,iboxbead))
        eikz(-1,iboxbead)=conjg(eikz(1,iboxbead))

        ! Calculating remaining kix, kiy, and kiz by recurrence:
        do kix=2, HI_kimax(1)
          eikx(kix,iboxbead)=eikx(kix-1,iboxbead)*eikx(1,iboxbead)
        end do
        do kiy=2, HI_kimax(2)
          eiky(kiy,iboxbead)=eiky(kiy-1,iboxbead)*eiky(1,iboxbead)
          eiky(-kiy,iboxbead)=conjg(eiky(kiy,iboxbead))
        end do
        do kiz=2, HI_kimax(3)
          eikz(kiz,iboxbead)=eikz(kiz-1,iboxbead)*eikz(1,iboxbead)
          eikz(-kiz,iboxbead)=conjg(eikz(kiz,iboxbead))
        end do
      end do ! iboxbead
      if (doTiming) et_EIKR=et_EIKR+tock(count1)

      ! loop over all beads:
!$omp parallel default(private) &
!$omp shared(ntotbead,Rb,Diff_tens,HI_c0,boxsize,m2_vec,PIx2,HI_ax2,doTiming,FlowType) &
!$omp shared(kvecx,kvecy,kvecz,eikx,eiky,eikz,unitDelta,HI_kmax,HI_kimax,rc_Dto2,et_R,et_K) &
!$omp shared(kiuppr,kiylowr,kizlowr,eps_m,delrx_L,rc_D,itime)
      ! OpenMP preparation:
     tid=omp_get_thread_num()
!$omp do schedule(auto)
      do jglobbead=1, ntotbead
        jglob=(jglobbead-1)*3
        do iglobbead=1, jglobbead
          iglob=(iglobbead-1)*3
          ! Calculation of the pair particle distance:
          rin0(1:3)=Rb(iglob+1:iglob+3)
          rjn0(1:3)=Rb(jglob+1:jglob+3)
          rijn0=rin0-rjn0
          ! Calculating Dij,PBC:
          Diff_tensP => Diff_tens(iglob+1:iglob+3,jglob+1:jglob+3)
          !------------------------------------------------
          !>>> Calculate the self part of diffusion matrix:
          !------------------------------------------------
          if (iglobbead == jglobbead) then
            Diff_tensP=HI_c0*unitDelta ! Correction for the diagonal terms of D
          else
            Diff_tensP=0._wp
          end if
          !-----------------------------------
          !>>> Calculate the real space sum:
          !-----------------------------------
         if (tid == 0) then ! Only on master
           if (doTiming) t0=omp_get_wtime()
         end if
          ! The summation over real space (loop over primary and all periodic image boxes):
          nmax=ceiling(rc_D/boxsize)
          n1 => nvec(1)
          n2 => nvec(2)
          n3 => nvec(3)
n1lp:     do n1=-nmax(1),nmax(1)
n2lp:       do n2=-nmax(2),nmax(2)
n3lp:         do n3=-nmax(3),nmax(3)
                if ((iglobbead == jglobbead).and.(all(nvec == 0))) cycle n3lp
                select case (FlowType)
                  case ('Equil')
                    rij=rijn0+nvec*boxsize
                  case ('PSF')
                    rij(1)=rijn0(1)+nvec(1)*boxsize(1)+nvec(2)*delrx_L
                    rij(2:3)=rijn0(2:3)+nvec(2:3)*boxsize(2:3)
                end select
                r2=dot_product(rij,rij)
                if (r2-rc_Dto2 <= 1.e-6) then
                  invr2=1/r2
                  rrhat(1)=rij(1)*rij(1)*invr2
                  rrhat(2)=rij(1)*rij(2)*invr2
                  rrhat(3)=rij(1)*rij(3)*invr2
                  rrhat(4)=rij(2)*rij(2)*invr2
                  rrhat(5)=rij(2)*rij(3)*invr2
                  rrhat(6)=rij(3)*rij(3)*invr2
                  r=sqrt(r2)
                  ! Adding the M1 contribution:
                  Diff_tensP=Diff_tensP+M1(r,r2,rrhat)
                  ! Overlaps in Primary box:
                  if (r < (HI_ax2)) Diff_tensP=Diff_tensP+Mstar(r,r2,rrhat)
                end if
              end do n3lp
            end do n2lp
          end do n1lp

         if (tid == 0) then ! Only master
           if (doTiming) then
             t1=omp_get_wtime()
             et_R=et_R+(t1-t0)
           end if
         end if
          !---------------------------------------
          !>>> Calculate the reciprocal space sum:
          !---------------------------------------
         if (tid == 0) then ! Only master
           if (doTiming) t0=omp_get_wtime()
         end if

          if ((boxsize(1) /= boxsize(2)) .or. (boxsize(2) /= boxsize(3))) then ! Unequal box length

            ktot=0
kxneb:      do kix=0, HI_kimax(1)
              if (FlowType == 'PSF') kiydev=nint(eps_m*kix)
              kvec(1)=kvecx(kix)
              ! because of inversion symmetry:
              if (kix == 0) then
                kiymin=0
              else
                kiymin=-HI_kimax(2)
              end if
kyneb:        do kiyy=kiymin, HI_kimax(2)
                select case (FlowType)
                  case ('Equil')
                    kvec(2)=kvecy(kiyy)
                  case ('PSF')
                    kvec(2)=kvecy(kiyy)+PIx2/boxsize(2)*kiydev
                end select
                ! because of inversion symmetry:
                if ((kix == 0) .and. (kiyy == 0)) then
                  kizmin=0
                else
                  kizmin=-HI_kimax(3)
                end if
kzneb:          do kiz=kizmin, HI_kimax(3)
                  kvec(3)=kvecz(kiz)
                  kto2=dot_product(kvec,kvec);invk2=1/kto2
                  if ((kto2 < (HI_kmax*HI_kmax)) .and. (kto2 /= 0)) then
                    kkhat(1)=kvec(1)*kvec(1)*invk2
                    kkhat(2)=kvec(1)*kvec(2)*invk2
                    kkhat(3)=kvec(1)*kvec(3)*invk2
                    kkhat(4)=kvec(2)*kvec(2)*invk2
                    kkhat(5)=kvec(2)*kvec(3)*invk2
                    kkhat(6)=kvec(3)*kvec(3)*invk2
                    ktot=ktot+1
                    select case (FlowType)
                      case ('Equil')
                        Coeff=2*m2_vec(ktot)* &
                              eikx(kix ,jglobbead)*conjg(eikx(kix ,iglobbead)) * &
                              eiky(kiyy,jglobbead)*conjg(eiky(kiyy,iglobbead)) * &
                              eikz(kiz ,jglobbead)*conjg(eikz(kiz ,iglobbead))
                      case ('PSF')
                        Coeff=2*m2_vec(ktot)* &
                              eikx(kix   ,jglobbead)*conjg(eikx(kix   ,iglobbead)) * &
                              eiky(kiyy  ,jglobbead)*conjg(eiky(kiyy  ,iglobbead)) * &
                              eiky(kiydev,jglobbead)*conjg(eiky(kiydev,iglobbead)) * &
                              eikz(kiz   ,jglobbead)*conjg(eikz(kiz   ,iglobbead))
                    end select
                    Diff_tensP(1,1)=Diff_tensP(1,1)+Coeff*(1._wp-kkhat(1))
                    Diff_tensP(1,2)=Diff_tensP(1,2)-Coeff*kkhat(2);Diff_tensP(2,1)=Diff_tensP(1,2)
                    Diff_tensP(1,3)=Diff_tensP(1,3)-Coeff*kkhat(3);Diff_tensP(3,1)=Diff_tensP(1,3)
                    Diff_tensP(2,2)=Diff_tensP(2,2)+Coeff*(1._wp-kkhat(4))
                    Diff_tensP(2,3)=Diff_tensP(2,3)-Coeff*kkhat(5);Diff_tensP(3,2)=Diff_tensP(2,3)
                    Diff_tensP(3,3)=Diff_tensP(3,3)+Coeff*(1._wp-kkhat(6))
                  end if
                end do kzneb
              end do kyneb
            end do kxneb

          else ! Equal box length

            ktot=0
kxeb:       do kix=0, kiuppr(0)
              if (FlowType == 'PSF') kiydev=nint(eps_m*kix)
              kikix=kix*kix
              kvec(1)=kvecx(kix)
kyeb:         do kiyy=kiylowr(kikix), kiuppr(kikix)
                kikixy=kikix+kiyy*kiyy
                select case (FlowType)
                  case ('Equil')
                    kvec(2)=kvecy(kiyy)
                  case ('PSF')
                    kvec(2)=kvecy(kiyy)+kvecy(kiydev)
                end select
kzeb:           do kiz=kizlowr(kikixy), kiuppr(kikixy)
                  kvec(3)=kvecz(kiz)
                  kto2=dot_product(kvec,kvec);invk2=1/kto2
                  kkhat(1)=kvec(1)*kvec(1)*invk2
                  kkhat(2)=kvec(1)*kvec(2)*invk2
                  kkhat(3)=kvec(1)*kvec(3)*invk2
                  kkhat(4)=kvec(2)*kvec(2)*invk2
                  kkhat(5)=kvec(2)*kvec(3)*invk2
                  kkhat(6)=kvec(3)*kvec(3)*invk2
                  ktot=ktot+1
                  select case (FlowType)
                    case ('Equil')
                      Coeff=2*m2_vec(ktot)* &
                            eikx(kix,jglobbead) *conjg(eikx(kix ,iglobbead)) * &
                            eiky(kiyy,jglobbead)*conjg(eiky(kiyy,iglobbead)) * &
                            eikz(kiz,jglobbead) *conjg(eikz(kiz ,iglobbead))
                    case ('PSF')
                      Coeff=2*m2_vec(ktot)* &
                            eikx(kix   ,jglobbead)*conjg(eikx(kix   ,iglobbead)) * &
                            eiky(kiyy  ,jglobbead)*conjg(eiky(kiyy  ,iglobbead)) * &
                            eiky(kiydev,jglobbead)*conjg(eiky(kiydev,iglobbead)) * &
                            eikz(kiz   ,jglobbead)*conjg(eikz(kiz   ,iglobbead))
                  end select
                  Diff_tensP(1,1)=Diff_tensP(1,1)+Coeff*(1._wp-kkhat(1))
                  Diff_tensP(1,2)=Diff_tensP(1,2)-Coeff*kkhat(2);Diff_tensP(2,1)=Diff_tensP(1,2)
                  Diff_tensP(1,3)=Diff_tensP(1,3)-Coeff*kkhat(3);Diff_tensP(3,1)=Diff_tensP(1,3)
                  Diff_tensP(2,2)=Diff_tensP(2,2)+Coeff*(1._wp-kkhat(4))
                  Diff_tensP(2,3)=Diff_tensP(2,3)-Coeff*kkhat(5);Diff_tensP(3,2)=Diff_tensP(2,3)
                  Diff_tensP(3,3)=Diff_tensP(3,3)+Coeff*(1._wp-kkhat(6))
                end do kzeb
              end do kyeb
            end do kxeb
!
          end if ! boxsize

         if (tid == 0) then ! Only master
           if (doTiming) then
             t1=omp_get_wtime()
             et_K=et_K+(t1-t0)
           end if
         end if

        end do !iglobbead
      end do ! jglobbead
!$omp end do
!$omp end parallel

      if (doTiming) et_EW=et_EW+tock(count0)

    elseif (HIcalc_mode == 'PME') then
      if (doTiming) call tick(count0)
      call calcDiff_real(Rb,ntotbead,boxsize,boxorigin) ! Calculation of real part of D.
      if (doTiming) then
        et_DCR=et_DCR+tock(count0)
        call tick(count0)
      end if
      call calcDiff_recip(Rb,ntotbead,boxsize,boxorigin) ! Calculation of reciprocal part of D.
      if (doTiming) et_DCK=et_DCK+tock(count0)

    end if ! HIcalc_mode

  end subroutine calcDiffTens_cpu


  function M1(r,rto2,rrhat)

    use :: hi_mod, only: sqrtPI

    real(wp) :: M1(3,3),m1_1,m1_2,r,rto2,rto3,rto4,rrhat(6)

    rto3=rto2*r; rto4=rto2*rto2
    m1_1=(erfc(HI_alpha*r)*(M1_c1/r+M1_c2/rto3) + &
          exp(-HI_alphato2*rto2)/sqrtPI*(M1_c3*rto2-M1_c4+M1_c5*rto4-M1_c6*rto2+M1_c7+M1_c8/rto2))
    m1_2=(erfc(HI_alpha*r)*(M1_c1/r-M1_c9/rto3) + &
          exp(-HI_alphato2*rto2)/sqrtPI*(M1_c10-M1_c3*rto2-M1_c5*rto4+M1_c11*rto2-M1_c12-M1_c13/rto2))

    M1(1,1)=m1_1+m1_2*rrhat(1)
    M1(1,2)=m1_2*rrhat(2);M1(2,1)=M1(1,2)
    M1(1,3)=m1_2*rrhat(3);M1(3,1)=M1(1,3)
    M1(2,2)=m1_1+m1_2*rrhat(4)
    M1(2,3)=m1_2*rrhat(5);M1(3,2)=M1(2,3)
    M1(3,3)=m1_1+m1_2*rrhat(6)
  end function M1

  subroutine M1_1d_routine(r,rto2,rrhat,M1_1d)

    use :: hi_mod, only: sqrtPI

    real(wp),intent(inout) :: M1_1d(:)
    real(wp),intent(in) :: rrhat(:)
    real(wp) :: m1_1,m1_2,r,rto2,rto3,rto4

    rto3=rto2*r; rto4=rto2*rto2
    m1_1=(erfc(HI_alpha*r)*(M1_c1/r+M1_c2/rto3) + &
          exp(-HI_alphato2*rto2)/sqrtPI*(M1_c3*rto2-M1_c4+M1_c5*rto4-M1_c6*rto2+M1_c7+M1_c8/rto2))
    m1_2=(erfc(HI_alpha*r)*(M1_c1/r-M1_c9/rto3) + &
          exp(-HI_alphato2*rto2)/sqrtPI*(M1_c10-M1_c3*rto2-M1_c5*rto4+M1_c11*rto2-M1_c12-M1_c13/rto2))

    M1_1d(1)=m1_1+m1_2*rrhat(1)
    M1_1d(2)=m1_2*rrhat(2);M1_1d(4)=M1_1d(2)
    M1_1d(3)=m1_2*rrhat(3);M1_1d(7)=M1_1d(3)
    M1_1d(5)=m1_1+m1_2*rrhat(4)
    M1_1d(6)=m1_2*rrhat(5);M1_1d(8)=M1_1d(6)
    M1_1d(9)=m1_1+m1_2*rrhat(6)

  end subroutine M1_1d_routine

  function Mstar(r,rto2,rrhat)

    real(wp) :: Mstar(3,3),rrhat(6),mstar_1,mstar_2,r,rto2,rto3

    rto3=r*rto2
    mstar_1=(1-Mstar_c1*r-M1_c1/r-M1_c2/rto3)
    mstar_2=(Mstar_c2*r-M1_c1/r+M1_c9/rto3)

    Mstar(1,1)=mstar_1+mstar_2*rrhat(1)
    Mstar(1,2)=mstar_2*rrhat(2);Mstar(2,1)=Mstar(1,2)
    Mstar(1,3)=mstar_2*rrhat(3);Mstar(3,1)=Mstar(1,3)
    Mstar(2,2)=mstar_1+mstar_2*rrhat(4)
    Mstar(2,3)=mstar_2*rrhat(5);Mstar(3,2)=Mstar(2,3)
    Mstar(3,3)=mstar_1+mstar_2*rrhat(6)

  end function Mstar

  subroutine Mstar_1d_routine(r,rto2,rrhat,Mstar_1d)

    real(wp),intent(inout) :: Mstar_1d(:)
    real(wp),intent(in) :: rrhat(:)
    real(wp) :: mstar_1,mstar_2,r,rto2,rto3

    rto3=r*rto2
    mstar_1=(1-Mstar_c1*r-M1_c1/r-M1_c2/rto3)
    mstar_2=(Mstar_c2*r-M1_c1/r+M1_c9/rto3)

    Mstar_1d(1)=mstar_1+mstar_2*rrhat(1)
    Mstar_1d(2)=mstar_2*rrhat(2);Mstar_1d(4)=Mstar_1d(2)
    Mstar_1d(3)=mstar_2*rrhat(3);Mstar_1d(7)=Mstar_1d(3)
    Mstar_1d(5)=mstar_1+mstar_2*rrhat(4)
    Mstar_1d(6)=mstar_2*rrhat(5);Mstar_1d(8)=Mstar_1d(6)
    Mstar_1d(9)=mstar_1+mstar_2*rrhat(6)

  end subroutine Mstar_1d_routine

  subroutine calcDiff_real(Rb,ntotbead,boxsize,boxorigin)

    use :: arry_mod, only: print_vector,print_matrix
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m,delrx_L

    integer,intent(in) :: ntotbead
    real(wp),intent(in) :: Rb(:)
    real(wp),intent(in) :: boxsize(3),boxorigin(3)
    real(wp),dimension(:,:),pointer :: Diff_tens_realP => null()
    integer :: ichain,jchain,ibead,ibeadub,jbead,iglob,jglob,offsetjch,offsetich
    integer :: iglobbead,jglobbead,jbeg,jend,jneig,offsetD,neig_count
    integer,parameter :: FORWARD=1,BACKWARD=-1
    real(wp),dimension(3) :: rij,rijn0,rijhat,kvec,rin0,rjn0,invboxsize,ri,rj
    integer,dimension(3) :: orient,kivec,nmax
    real(wp) :: r,r2,rhatrhat(3,3),rrhat(6),rijx_tr,invr2
    integer,target :: nvec(3)
    integer,pointer :: n1,n2,n3

    if (Dreal_sparse_mode) then

      invboxsize(1:3)=1/boxsize(1:3)
      offsetD=0
      neig_count=0
      Dreal_rowInd(1)=1
      do iglobbead=1, ntotbead-1
        iglob=(iglobbead-1)*3
        jbeg=point_D(iglobbead)
        jend=point_D(iglobbead+1)-1
        ! Check if iglobbead has neighbor:
        if (jbeg <= jend) then
          ri(1:3)=Rb(iglob+1:iglob+3)
          do jneig=jbeg, jend
            jglobbead=list_DP(jneig)
            jglob=(jglobbead-1)*3
            rj(1:3)=Rb(jglob+1:jglob+3)
            rij=ri-rj
            ! minimum image convension: (We know that rc_D < L/2)
            select case (FlowType)
              case ('Equil')
                rij(:)=rij(:)-nint(rij(:)*invboxsize(:))*boxsize(:)
              case ('PSF')
                rijx_tr=rij(1)-eps_m*rij(2)
                rijx_tr=rijx_tr-nint(rijx_tr*invboxsize(1))*boxsize(1)
                rij(2:3)=rij(2:3)-nint(rij(2:3)*invboxsize(2:3))*boxsize(2:3)
                rij(1)=rijx_tr+eps_m*rij(2)
            end select
            r2=dot_product(rij,rij);r=sqrt(r2)
            if (r <= rc_D) then
              neig_count=neig_count+1
              rijhat=rij/r
              rrhat(1)=rijhat(1)*rijhat(1)
              rrhat(2)=rijhat(1)*rijhat(2)
              rrhat(3)=rijhat(1)*rijhat(3)
              rrhat(4)=rijhat(2)*rijhat(2)
              rrhat(5)=rijhat(2)*rijhat(3)
              rrhat(6)=rijhat(3)*rijhat(3)
              Dreal_vals(offsetD+1:offsetD+9)=reshape(M1(r,r2,rrhat),shape=[9])
              ! The correction due to overlap:
              if (r < (HI_ax2)) then ! Ovelaps in Primary box (or maybe in general!)
                Dreal_vals(offsetD+1:offsetD+9)=Dreal_vals(offsetD+1:offsetD+9) + &
                                                reshape(Mstar(r,r2,rrhat),shape=[9])
              end if
              Dreal_cols(neig_count)=jglobbead
              offsetD=offsetD+9
            end if
          end do ! jneig
        end if ! jbeg.le.jend
        Dreal_rowInd(iglobbead+1)=neig_count+1
      end do ! iglobbead
      Dreal_rowInd(ntotbead+1)=neig_count+1

    else ! Dreal_sparse_mode is .false.

      ! loop over boxbeads:
!$omp parallel default(private) shared(ntotbead,Rb,Diff_tens_real,boxsize,HI_ax2,rc_Dto2) &
!$omp shared(FlowType,delrx_L,rc_D)
!$omp do schedule(auto)
      do jglobbead=1, ntotbead
        jglob=(jglobbead-1)*3
        do iglobbead=1, jglobbead
          iglob=(iglobbead-1)*3
          ! Calculation of the pair particle distance:
          rin0(1:3)=Rb(iglob+1:iglob+3)
          rjn0(1:3)=Rb(jglob+1:jglob+3)
          rijn0=rin0-rjn0
          ! Calculating Dij,PBC:
          Diff_tens_realP => Diff_tens_real(iglob+1:iglob+3,jglob+1:jglob+3)
          Diff_tens_realP=0._wp
          ! The summation over real space (loop over primary and all periodic image boxes):
          nmax=ceiling(rc_D/boxsize)
          n1 => nvec(1)
          n2 => nvec(2)
          n3 => nvec(3)
n1lp:     do n1=-nmax(1),nmax(1)
n2lp:       do n2=-nmax(2),nmax(2)
n3lp:         do n3=-nmax(3),nmax(3)
                if ((iglobbead == jglobbead).and.(all(nvec == 0))) cycle n3lp
                select case (FlowType)
                  case ('Equil')
                    rij=rijn0+nvec*boxsize
                  case ('PSF')
                    rij(1)=rijn0(1)+nvec(1)*boxsize(1)+nvec(2)*delrx_L
                    rij(2:3)=rijn0(2:3)+nvec(2:3)*boxsize(2:3)
                end select
                r2=dot_product(rij,rij)
                if (r2-rc_Dto2 <= 1.e-6) then
                  invr2=1/r2
                  rrhat(1)=rij(1)*rij(1)*invr2
                  rrhat(2)=rij(1)*rij(2)*invr2
                  rrhat(3)=rij(1)*rij(3)*invr2
                  rrhat(4)=rij(2)*rij(2)*invr2
                  rrhat(5)=rij(2)*rij(3)*invr2
                  rrhat(6)=rij(3)*rij(3)*invr2
                  r=sqrt(r2)
                  ! Adding the M1 contribution:
                  Diff_tens_realP=Diff_tens_realP+M1(r,r2,rrhat)
                  ! Overlaps in Primary box:
                  if (r < (HI_ax2)) Diff_tens_realP=Diff_tens_realP+Mstar(r,r2,rrhat)
                end if
              end do n3lp
            end do n2lp
          end do n1lp
        end do !iglobbead
      end do ! jglobbead

!$omp end do
!$omp end parallel

    end if ! Dreal_sparse_mode

  end subroutine calcDiff_real

  subroutine calcDiff_recip(Rb,ntotbead,boxsize,boxorigin)

    use :: arry_mod, only: print_spmatrix,print_vector,print_matrix
    use :: flow_mod, only: FlowType
    use :: trsfm_mod, only: eps_m

    integer,intent(in) :: ntotbead
    real(wp),intent(in) :: Rb(:)
    real(wp) :: boxsize(3),boxorigin(3),xi(3),Mx,My,Mz
    integer :: iglobbead,iglob,icoor,nearestMesh(3),igrid,jgrid,kgrid
    integer :: k_ind,k1,k2,k3,elem_count,elem_counttmp,k2Kx,k3KxKy
    integer,allocatable,dimension(:,:),save :: grid

!!$omp threadprivate (grid)

!!$omp parallel default(private) copyin(grid) &
!!$omp shared(ntotbead,Rb,boxsize,boxorigin,K_mesh,InterpMethod,p_PME,p_PMEto3,r_str,P_vals) &
!!$omp shared(P_cols,eps_m,FlowType)

    allocate(grid(p_PME,3))
!   elem_count=0 ! used for sparse P.
!!$omp do schedule(auto)
    do iglobbead=1, ntotbead
      iglob=(iglobbead-1)*3
      select case (FlowType)
        case ('Equil')
          xi(1:3)=K_mesh(1:3)/boxsize(1:3)*(Rb(iglob+1:iglob+3)-boxorigin(1:3))
        case ('PSF')
          xi(1)=K_mesh(1)/boxsize(1)*(Rb(iglob+1)-boxorigin(1)- &
                                     (Rb(iglob+2)-boxorigin(2))*eps_m)
          xi(2:3)=K_mesh(2:3)/boxsize(2:3)*(Rb(iglob+2:iglob+3)-boxorigin(2:3))
      end select

      select case (InterpMethod)
        case('BSpline')
          do icoor=1, 3 ! The sweeping direction
            ! The nearest mesh which is less than fractional coordinate:
            nearestMesh(icoor)=floor(xi(icoor))
            grid(1,icoor)=nearestMesh(icoor) ! The nearest mesh point.
            do igrid=1, p_PME-1
              if ((nearestMesh(icoor)-igrid) < 0) then
                grid(igrid+1,icoor)=nearestMesh(icoor)-igrid+K_mesh(icoor) ! Wrap around periodic box
              else
                grid(igrid+1,icoor)=nearestMesh(icoor)-igrid
              end if
            end do ! igrid
          end do ! icoor
        case('Lagrange')
          ! Not implemented yet!
      end select

      ! For Fortran-style arrangement:

!      elem_counttmp=elem_count
      elem_count=(iglobbead-1)*p_PMEto3
      elem_counttmp=elem_count
      kg: do kgrid=1, p_PME
        k3=grid(kgrid,3)
        select case(InterpMethod)
          case('BSpline')
            Mz=M_spl(xi(3)-(nearestMesh(3)-kgrid+1),p_PME)
          case('Lagrange')
            ! Not implemented yet!
        end select
        k3KxKy=k3*r_str(3)
        jg: do jgrid=1, p_PME
          k2=grid(jgrid,2)
          select case(InterpMethod)
            case('BSpline')
              My=M_spl(xi(2)-(nearestMesh(2)-jgrid+1),p_PME)
            case('Lagrange')
              ! Not implemented yet!
          end select
          k2Kx=k2*r_str(2)
          ig: do igrid=1, p_PME
            k1=grid(igrid,1)
            select case(InterpMethod)
              case('BSpline')
                Mx=M_spl(xi(1)-(nearestMesh(1)-igrid+1),p_PME)
              case('Lagrange')
                ! Not implemented yet!
            end select
              ! Calculate the scalar k-index for grid points:
            k_ind=k1+k2Kx+k3KxKy
            P_vals(elem_counttmp+p_PMEto3-mod(elem_count,p_PMEto3))=Mx*My*Mz
            P_cols(elem_counttmp+p_PMEto3-mod(elem_count,p_PMEto3))=k_ind
            elem_count=elem_count+1
          end do ig
        end do jg
      end do kg

    end do ! iglobbead

!!$omp end do
    deallocate(grid)
!!$omp end parallel

  end subroutine calcDiff_recip

! %=========================== These Routines are for the case of PME ======================%

  subroutine PME_cpu(F,ntotbead,boxsize,DF_tot)

    use :: arry_mod, only: print_vector,print_matrix,print_spmatrix

    integer,intent(in) :: ntotbead
    real(wp),dimension(:),intent(in),target :: F
    real(wp),dimension(:),intent(inout) :: DF_tot
    real(wp) :: boxsize(3)
    integer(long) :: count0,count1

    if (doTiming) then
      PMEcount=PMEcount+1
      call tick(count0)
    end if
    ! Self part:
    call calcDF_self() ! Calculation of correction of D.F due to self interaction.
    ! Real part:
    if (doTiming) call tick(count1)
    call calcDF_real() ! Calculation of real part of D.F.
    if (doTiming) et_R=et_R+tock(count1)
    ! Recip part:
    if (doTiming) call tick(count1)
    call calcDF_recip() ! Calculation of reciprocal part of D.F.
    if (doTiming) et_K=et_K+tock(count1)

    if (doTiming) et_PME=et_PME+tock(count0)

    ! Tot DF calculation:
    DF_tot=DF_self+DF_real+DF_recip

    ! call print_vector(DF_recip,'df_recip')
    ! call print_vector(DF_self,'df_self')
    ! call print_vector(DF_real,'df_real')

  contains

    subroutine calcDF_self()

      DF_self=HI_c0*F

    end subroutine calcDF_self

    subroutine calcDF_real()

      if (Dreal_sparse_mode) then
#ifdef USE_DP
          call mkl_dbsrsymv('U',ntotbead,3,Dreal_vals,Dreal_rowInd,Dreal_cols,F,DF_real)
#elif USE_SP
          call mkl_sbsrsymv('U',ntotbead,3,Dreal_vals,Dreal_rowInd,Dreal_cols,F,DF_real)
#endif
      else
        call symv(Diff_tens_real,F,DF_real)
      end if

    end subroutine calcDF_real

    subroutine calcDF_recip()

      use :: mkl_dfti

     ! complex(wp) :: Infl(3,3),C_mat(3),D_mat(3)
      integer :: i,j,k,ic,mivecx,mivecy,mivecz,miveczy,mpivecx,mpivecy,mpivecz
      integer :: k1,k2,k3,myStatus,m_ind,mtot,mtmp
      real(wp) :: mpx,mpy,mpz,m2,m2tmp,invm2,m2_vectmp,mpmphat(6),mpvec(3),mphat(3)
      real(wp) :: Infl(6)
      complex(wp) :: C_mat(3)
      real(wp),pointer,dimension(:) :: FPx,FPy,FPz,DFPx,DFPy,DFPz
      complex(wp),dimension(0:Kcto3/2-1) :: Cx_mat,Cy_mat,Cz_mat
      pointer(Cx_ptr,Cx_mat), (Cy_ptr,Cy_mat), (Cz_ptr,Cz_mat)
      complex(wp) :: c1,c2,c3
      integer(long) :: count2,count3

      ! Associate complex view with F_mesh for root thread.
      Cx_ptr=loc(F_mesh(0:Kcto3-1,1))
      Cy_ptr=loc(F_mesh(0:Kcto3-1,2))
      Cz_ptr=loc(F_mesh(0:Kcto3-1,3))
      !-------------------------------------------
      !>>> Spreading the forces onto regular mesh:
      !-------------------------------------------
      if (doTiming) call tick(count2)
      FPx => F(1:ntotbead*3-2:3)
      FPy => F(2:ntotbead*3-1:3)
      FPz => F(3:ntotbead*3:3)
#ifdef USE_DP
        call mkl_dcsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPx,0._wp,F_meshPx)
        call mkl_dcsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPy,0._wp,F_meshPy)
        call mkl_dcsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPz,0._wp,F_meshPz)
#elif USE_SP
        call mkl_scsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPx,0._wp,F_meshPx)
        call mkl_scsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPy,0._wp,F_meshPy)
        call mkl_scsrmv('T',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),FPz,0._wp,F_meshPz)
#endif
      if (doTiming) et_SPR=et_SPR+tock(count2)
      !-----------------------------------
      !>>> FWD in-place Fourier Transform:
      !-----------------------------------
      if (doTiming) call tick(count2)
      FFTStatus=DftiComputeForward(FFTfwDescHand,F_meshPx)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeForward; x-component of force.")'
        print '("  Error, status = ",i0)', FFTStatus
        print *, DftiErrorMessage(FFTStatus)
        stop
      end if
      FFTStatus=DftiComputeForward(FFTfwDescHand,F_meshPy)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeForward; y-component of force.")'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if
      FFTStatus=DftiComputeForward(FFTfwDescHand,F_meshPz)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeForward; z-component of force.")'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if
      if (doTiming) et_FFT=et_FFT+tock(count2)
      !-------------------------------------
      !>>>> Applying the Influence function:
      !-------------------------------------
      ! The prime indecis, i.e. mpivecx~z, is to consider the original frequency to be used in M2.
      if (doTiming) call tick(count2)
     ! mtot=0
      !$omp parallel default(private) shared(K_mesh,c_str,mpvecx,mpvecy,mpvecz,Cx_ptr,Cy_ptr,Cz_ptr) &
      !$omp shared(InterpMethod,m2_vec)
      !$omp do schedule(auto)
      mzy: do miveczy=0, K_mesh(3)*K_mesh(2)-1
        mivecz=miveczy/K_mesh(2)
        mivecy=mod(miveczy,K_mesh(2))
        mpz=mpvecz(mivecz)
        mpy=mpvecy(mivecy)
        m2tmp=mpz*mpz+mpy*mpy
        mtmp=miveczy*c_str(2) ! miveczy*(K1/2+1)
        mx: do mivecx=0, K_mesh(1)/2
          mpx=mpvecx(mivecx)
          m2=m2tmp+mpx*mpx
          if (m2 == 0) then
            Cx_mat(0)=0._wp;Cy_mat(0)=0._wp;Cz_mat(0)=0._wp
          else
            invm2=1._wp/m2
           ! mtot=mtot+1
            m_ind=mivecx+mivecy*c_str(2)+mivecz*c_str(3)
            C_mat=[Cx_mat(m_ind),Cy_mat(m_ind),Cz_mat(m_ind)]
            select case(InterpMethod)
              case('BSpline')
                mpmphat(1)=invm2*mpx*mpx
                mpmphat(2)=invm2*mpx*mpy
                mpmphat(3)=invm2*mpx*mpz
                mpmphat(4)=invm2*mpy*mpy
                mpmphat(5)=invm2*mpy*mpz
                mpmphat(6)=invm2*mpz*mpz
                m2_vectmp=m2_vec(mtmp+mivecx)
                Infl(1)=(1._wp-mpmphat(1))*m2_vectmp
                Infl(2)=-mpmphat(2)*m2_vectmp
                Infl(3)=-mpmphat(3)*m2_vectmp
                Infl(4)=(1._wp-mpmphat(4))*m2_vectmp
                Infl(5)=-mpmphat(5)*m2_vectmp
                Infl(6)=(1._wp-mpmphat(6))*m2_vectmp
              case('Lagrange')
                ! Not implemented yet!
            end select
            Cx_mat(m_ind)=Infl(1)*C_mat(1)+Infl(2)*C_mat(2)+Infl(3)*C_mat(3)
            Cy_mat(m_ind)=Infl(2)*C_mat(1)+Infl(4)*C_mat(2)+Infl(5)*C_mat(3)
            Cz_mat(m_ind)=Infl(3)*C_mat(1)+Infl(5)*C_mat(2)+Infl(6)*C_mat(3)
           ! print *,'id:',omp_get_thread_num()
           ! print *,'mp:',mpx,mpy,mpz
           ! call print_vector(Infl,'infl')
           ! print *,'dmat:'
           ! print *,Cx_mat(m_ind)
           ! print *,Cy_mat(m_ind)
           ! print *,Cz_mat(m_ind)
          end if ! m2.eq.0
        end do mx
      end do mzy
      !$omp end do
      !$omp end parallel
      ! Note!!: We have to symmetrize the resulting vector after applying influence function, in case of even K's.
      !  This is because the CCE symmetry breaks for (i=K/2,j,k; i=K/2,K-j,K-k) and similar cases, as we apply M2.
      !  So in order to keep haveing C->R while having the correct backward FFT we need to symmetrize these cases.
      !  This happens because we essentially neglect the Fourier vectors like i=-K/2,j,k whose imaginary parts ca-
      !  ncels that of i=K/2,j,k.
      if (any(mod(K_mesh,2) == 0)) then
        do k=0, K_mesh(3)-1
          do j=0, K_mesh(2)-1
            do i=0, K_mesh(1)/2
              if (i == mod(K_mesh(1)-i,K_mesh(1))) then
                ! X component:
                c1=Cx_mat(i+c_str(2)*j+c_str(3)*k)
                c2=Cx_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))
                if (c1 /= conjg(c2)) then
                  c3=0.5*(c1+conjg(c2))
                  Cx_mat(i+c_str(2)*j+c_str(3)*k)=c3
                  Cx_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))=conjg(c3)
                endif
                ! Y component:
                c1=Cy_mat(i+c_str(2)*j+c_str(3)*k)
                c2=Cy_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))
                if (c1 /= conjg(c2)) then
                  c3=0.5*(c1+conjg(c2))
                  Cy_mat(i+c_str(2)*j+c_str(3)*k)=c3
                  Cy_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))=conjg(c3)
                endif
                ! Z component:
                c1=Cz_mat(i+c_str(2)*j+c_str(3)*k)
                c2=Cz_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))
                if (c1 /= conjg(c2)) then
                  c3=0.5*(c1+conjg(c2))
                  Cz_mat(i+c_str(2)*j+c_str(3)*k)=c3
                  Cz_mat(i+c_str(2)*mod(K_mesh(2)-j,K_mesh(2))+c_str(3)*mod(K_mesh(3)-k,K_mesh(3)))=conjg(c3)
                endif
              endif
            enddo
          enddo
        enddo
      end if
      if (doTiming) et_INF=et_INF+tock(count2)
      !-----------------------------------
      !>>> BWD in-place Fourier Transform:
      !-----------------------------------
      if (doTiming) call tick(count2)
      FFTStatus=DftiComputeBackward(FFTbwDescHand,Cx_mat)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeBackward; x-component of force.")'
        print '(" Error with status = ",i0)', FFTStatus
        stop
      end if
      FFTStatus=DftiComputeBackward(FFTbwDescHand,Cy_mat)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeBackward; y-component of force.")'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if
      FFTStatus=DftiComputeBackward(FFTbwDescHand,Cz_mat)
      if (FFTStatus /= 0) then
        print '(" Error!!: Problem in DftiComputeBackward; z-component of force.")'
        print '("  Error, status = ",i0)', FFTStatus
        stop
      end if
      if (doTiming) et_IFFT=et_IFFT+tock(count2)
      !--------------------------------------------------------------------------
      !>>> Interpolating back the forces from regular mesh to particle positions:
      !--------------------------------------------------------------------------
      if (doTiming) call tick(count2)
      DFPx => DF_recip(1:ntotbead*3-2:3)
      DFPy => DF_recip(2:ntotbead*3-1:3)
      DFPz => DF_recip(3:ntotbead*3:3)
#ifdef USE_DP
        call mkl_dcsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPx,0._wp,DFPx)
        call mkl_dcsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPy,0._wp,DFPy)
        call mkl_dcsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPz,0._wp,DFPz)
#elif USE_SP
        call mkl_scsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPx,0._wp,DFPx)
        call mkl_scsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPy,0._wp,DFPy)
        call mkl_scsrmv('N',ntotbead,Kcto3,1._wp,'GIIC',P_vals,P_cols,P_rowInd,P_rowInd(2),F_meshPz,0._wp,DFPz)
#endif
      if (doTiming) et_INT=et_INT+tock(count2)

    end subroutine calcDF_recip

  end subroutine PME_cpu

end module Diffcalc_mod