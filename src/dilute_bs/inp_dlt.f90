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
module inp_dlt

  use :: prcn_mod

  implicit none

  save

  ! Parameters from user:
  integer,protected :: nchain
  integer,protected :: nseg
  character(len=10),protected :: tplgy
  integer,protected :: Na
  character(len=10),protected :: arm_plc
  integer,allocatable :: Ia(:)
  integer,protected :: nseg_ar
  character(len=10),protected :: LambdaMethod
  real(wp),protected :: lambda
  character(len=10),protected :: ForceLaw
  character(len=10),protected :: TruncMethod
  real(wp),protected :: qr_l
  real(wp),protected :: N_Ks,RWS_v,WLC_v
  real(wp),protected :: b
  logical,protected :: applFext
  logical,protected :: srf_tet
  real(wp),protected :: Fext0
  integer,protected :: iflow
  integer,protected :: nWi
  real(wp),protected :: Wii,Wif
  character(len=10),protected :: WiSpacing
  character(len=10),protected :: HITens
  real(wp),protected :: hstar
  character(len=10),protected :: DecompMeth
  integer,protected :: ncols
  integer :: mBlLan,mubBlLan
  logical,protected :: mset
  integer :: LCheb,Lub
  logical,protected :: Lset
  real(wp),protected :: errormin
  integer,protected :: upfactr
  logical,protected :: AveIterCalc
  character(len=10),protected :: EV_bb
  character(len=10),protected :: EV_bw
  real(wp),protected :: zstar,dstar
  real(wp),protected :: Aw
  character(len=10),protected :: dstarCalc
  real(wp),protected :: LJ_eps,LJ_sig,LJ_rtr,LJ_rc
  integer,protected :: minNonBond
  character(len=10),protected :: initmode
  real(wp),protected :: infrx,infry,infrz
  real(wp),protected :: tend,tss
  real(wp) :: trst
  integer,protected :: ndt
  real(wp),protected :: dti,dtf
  character(len=10),protected :: dtSpacing
  character(len=10),protected :: dtCalc
  logical,protected :: dtScale
  integer,protected :: nAdjSeq
  logical,protected :: Adjust_dt
  real(wp),allocatable,protected :: AdjSeq(:)
  real(wp),allocatable,protected :: AdjFact(:)
  real(wp),protected :: tol
  integer,protected :: nroots
  integer,protected :: PrScale
  logical,protected :: CoM,CoHR
  real(wp),protected :: frm_rt_rep,frm_rt_pp,frm_rt_rst,frm_rt_dmp
  logical,protected :: StrCalc
  logical,protected :: cnf_srt
  logical,protected :: indvlext
  integer,protected :: residx
  logical,protected :: TimerA
  logical,protected :: RgCalc
  logical,protected :: cosThCalc
  character(len=10),protected :: cosmode
  logical,protected :: DumpstrCalc
  integer,protected :: nstr
  real(wp),allocatable,protected :: dumpstr(:)
  ! Parameters based on input data:
  real(wp),allocatable,protected :: Wi(:),Pe(:)
  real(wp),allocatable :: dttemp1(:),dt(:,:),dt_tmp(:,:)
  integer,allocatable,protected :: ntime(:,:),itime_AdjSeq(:,:,:)
  integer,protected :: nbead,nbeadx3,nsegx3
  integer,protected :: ntotang
  real(wp),protected :: qmax
  real(wp),protected :: q_l,qr_lto2
  integer,protected :: npchain
  real(wp),protected :: RWS_C,RWS_D,WLC_A,WLC_B,WLC_C
  logical,protected :: unif_flow, sph_flow

contains

  subroutine read_dlt(id,inpFile)

    use,intrinsic :: iso_fortran_env
    use :: strg_mod, only: parse,value
    use :: arry_mod, only: sort

    integer :: il,j,ntokens,u1,stat,ios,id,k,i,n,msec,tm_inf(8),icnt
    real(wp),allocatable :: u(:)
    character(len=1024) :: line
    character(len=100) :: tokens(10)
    character(len=20) :: inpFile

    ! Default values
    tplgy='Linear'
    LambdaMethod='Rouse'
    ForceLaw='Hookean'
    TruncMethod='Linear';qr_l=1._wp
    applFext=.false.
    srf_tet=.false.;arm_plc='Random'
    unif_flow=.false.
    sph_flow=.false.
    iflow=1
    nWi=1;Wii=0._wp;Wif=0._wp;WiSpacing='Linear'
    hstar=0._wp;HITens='RPY';DecompMeth='Cholesky';ncols=1
    mBlLan=3;mubBlLan=15;mset=.false.
    LCheb=2;Lub=20;Lset=.false.
    AveIterCalc=.false.
    errormin=real(1.e-2,kind=wp)
    upfactr=50
    EV_bb='NoEV';dstar=1._wp;dstarCalc='Kumar';minNonBond=1
    EV_bw='NoEV';Aw=25._wp
    initmode='st';infrx=0.7_wp;infry=0._wp;infrz=0._wp
    tend=10._wp;tss=5._wp;trst=0._wp
    ndt=1;dti=0.01_wp;dtf=0.01_wp;dtSpacing='Linear';dtCalc='Self';dtScale=.false.
    Adjust_dt=.false.
    tol=real(1.e-4,kind=wp)
    nroots=10**6;PrScale=1
    CoM=.false.;CoHR=.false.
    frm_rt_rep=0.1_wp;frm_rt_pp=0.002_wp;frm_rt_rst=0.1;frm_rt_dmp=0.1
    StrCalc=.true.
    cnf_srt=.false.;residx=10
    indvlext=.false.
    TimerA=.false.
    RgCalc=.false.
    cosThCalc=.false.;cosmode='angle'
    DumpstrCalc=.false.

    open (newunit=u1,action='read',file=inpFile,status='old')
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
            case ('nchain')
              call value(tokens(j+1),nchain,ios)
            case ('nseg')
              call value(tokens(j+1),nseg,ios)
            case ('tplgy')
              tplgy=trim(adjustl(tokens(j+1)))
            case ('nseg_ar')
              call value(tokens(j+1),nseg_ar,ios)
            case ('Arms')
              call value(tokens(j+1),Na,ios)
              arm_plc=trim(adjustl(tokens(j+2)))
              allocate(Ia(Na+1))
              Ia(1)=1 ! The first index is always one.
              if (arm_plc == 'Fixed') then
                do k=1, Na
                  call value(tokens(j+2+k),Ia(k+1),ios)
                end do
              elseif (arm_plc == 'Stored') then
                ! Will be read in dlt_bs routine
              else
                allocate(u(Na))
                call date_and_time(values=tm_inf)
                msec=(1000*tm_inf(7)+tm_inf(8))*((id-83)*359)
                call random_seed(size=n)
                call random_seed(put=(/(i*msec,i=1,n)/))
                icnt=0
rndmlp:         do
                  call random_number(u)
                  do k=1, Na
                    Ia(k+1)=2+floor((nseg-Na*nseg_ar-1)*u(k))
                  end do
                  do i=1, Na
                    do k=1, i-1
                      if (Ia(i+1) == Ia(k+1)) then
                        icnt=icnt+1
                        if (icnt >= 100) then
                          print '(" Random placement of arms did not converge.")'
                          stop
                        end if
                        cycle rndmlp
                      end if
                    end do
                  end do
                  exit rndmlp
                end do rndmlp
                deallocate(u)
              end if
              call sort(Ia)
!              print *,'id:',id,Ia
            case ('Rel-Model')
              LambdaMethod=trim(adjustl(tokens(j+1)))
              if (LambdaMethod == 'Self') then
                call value(tokens(j+2),lambda,ios)
              end if
            case ('SPR-Force')
              ForceLaw=trim(adjustl(tokens(j+1)))
            case ('Truncation')
              TruncMethod=trim(adjustl(tokens(j+1)))
              if (TruncMethod /= 'None') then
                call value(tokens(j+2),qr_l,ios)
              end if
            case ('N_Ks')
              call value(tokens(j+1),N_Ks,ios)
              select case (ForceLaw)
                case ('RWS')
                  call value(tokens(j+1),RWS_v,ios)
                case ('WLC_UD','WLC_GEN')
                  call value(tokens(j+1),WLC_v,ios)
              end select
            case ('b')
              call value(tokens(j+1),b,ios)
            case ('EXT-Force')
              if(tokens(j+1) == 'TRUE') then
                applFext=.true.
                call value(tokens(j+2),Fext0,ios)
              elseif(tokens(j+1) == 'FALSE') then
                applFext=.false.
              end if
            case ('Surf-Tethered')
              if(tokens(j+1) == 'TRUE') then
                srf_tet=.true.
              elseif(tokens(j+1) == 'FALSE') then
                srf_tet=.false.
              end if
            case ('Unif-Flow')
              if(tokens(j+1) == 'TRUE') then
                unif_flow=.true.
              elseif(tokens(j+1) == 'FALSE') then
                unif_flow=.false.
              end if
            case ('Sph-Flow')
              if(tokens(j+1) == 'TRUE') then
                sph_flow=.true.
              elseif(tokens(j+1) == 'FALSE') then
                sph_flow=.false.
              end if
            case ('Flow-Type')
              call value(tokens(j+1),iflow,ios)
            case ('nWi')
              call value(tokens(j+1),nWi,ios)
            case ('Wi')
              call value(tokens(j+1),Wii,ios)
              call value(tokens(j+2),Wif,ios)
              WiSpacing=trim(adjustl(tokens(j+3)))
            case ('HITens')
              HITens=trim(adjustl(tokens(j+1)))
            case ('hstar')
              call value(tokens(j+1),hstar,ios)
            case ('DecompMeth')
              DecompMeth=trim(adjustl(tokens(j+1)))
            case ('ncols')
              call value(tokens(j+1),ncols,ios)
            case ('m')
              call value(tokens(j+1),mBlLan,ios)
              call value(tokens(j+2),mubBlLan,ios)
              if(tokens(j+3) == 'TRUE') then
                mset=.true.
              elseif(tokens(j+3) == 'FALSE') then
                mset=.false.
              end if
            case ('L')
              call value(tokens(j+1),LCheb,ios)
              call value(tokens(j+2),Lub,ios)
              if(tokens(j+3) == 'TRUE') then
                Lset=.true.
              elseif(tokens(j+3) == 'FALSE') then
                Lset=.false.
              end if
            case ('errormin')
              call value(tokens(j+1),errormin,ios)
            case ('upfactr')
              call value(tokens(j+1),upfactr,ios)
            case ('AveIter-rep')
              if(tokens(j+1) == 'TRUE') then
                AveIterCalc=.true.
              elseif(tokens(j+1) == 'FALSE') then
                AveIterCalc=.false.
              end if
            case ('EV-bb')
              EV_bb=trim(adjustl(tokens(j+1)))
            case ('EV-bw')
              EV_bw=trim(adjustl(tokens(j+1)))
            case ('zstar')
              call value(tokens(j+1),zstar,ios)
            case ('dstar')
              call value(tokens(j+1),dstar,ios)
              dstarCalc=trim(adjustl(tokens(j+2)))
            case ('LJ-Par')
              call value(tokens(j+1),LJ_eps,ios)
              call value(tokens(j+2),LJ_sig,ios)
              call value(tokens(j+3),LJ_rtr,ios)
              call value(tokens(j+4),LJ_rc,ios)
            case ('minNonBond')
              call value(tokens(j+1),minNonBond,ios)
            case ('initmode')
              initmode=trim(adjustl(tokens(j+1)))
            case ('init-ext')
              call value(tokens(j+1),infrx,ios)
              call value(tokens(j+2),infry,ios)
              call value(tokens(j+3),infrz,ios)
            case ('tend')
              call value(tokens(j+1),tend,ios)
            case ('tss')
              call value(tokens(j+1),tss,ios)
            case ('trst')
              call value(tokens(j+1),trst,ios)
            case ('ndt')
              call value(tokens(j+1),ndt,ios)
            case ('dt')
              call value(tokens(j+1),dti,ios)
              call value(tokens(j+2),dtf,ios)
              dtSpacing=trim(adjustl(tokens(j+3)))
              dtCalc=trim(adjustl(tokens(j+4)))
              if (dtCalc == 'Self') then
                if(tokens(j+5) == 'TRUE') then
                  dtScale=.true.
                elseif(tokens(j+5) == 'FALSE') then
                  dtScale=.false.
                end if
              end if
            case ('Adjust_dt')
              if(tokens(j+1) == 'TRUE') then
                Adjust_dt=.true.
              elseif(tokens(j+1) == 'FALSE') then
                Adjust_dt=.false.
              end if
            case ('nAdjSeq')
              call value(tokens(j+1),nAdjSeq,ios)
            case ('AdjSeq')
              allocate(AdjSeq(nAdjSeq))
              do k=1, nAdjSeq
                call value(tokens(j+k),AdjSeq(k),ios)
              end do
            case ('AdjFact')
              allocate(AdjFact(nAdjSeq))
              do k=1, nAdjSeq
                call value(tokens(j+k),AdjFact(k),ios)
              end do
            case ('tol')
              call value(tokens(j+1),tol,ios)
            case ('nroots')
              call value(tokens(j+1),nroots,ios)
            case ('PrScale')
              call value(tokens(j+1),PrScale,ios)
            case ('CoM')
              if(tokens(j+1) == 'TRUE') then
                CoM=.true.
              elseif(tokens(j+1) == 'FALSE') then
                CoM=.false.
              end if
            case ('CoHR')
              if(tokens(j+1) == 'TRUE') then
                CoHR=.true.
              elseif(tokens(j+1) == 'FALSE') then
                CoHR=.false.
              end if
            case ('frm-rt-rep')
              call value(tokens(j+1),frm_rt_rep,ios)
            case ('frm-rt-pp')
              call value(tokens(j+1),frm_rt_pp,ios)
            case ('frm-rt-rst')
              call value(tokens(j+1),frm_rt_rst,ios)
            case ('frm-rt-dmp')
              call value(tokens(j+1),frm_rt_dmp,ios)
            case ('Conf-anal')
              if(tokens(j+1) == 'TRUE') then
                StrCalc=.true.
              elseif(tokens(j+1) == 'FALSE') then
                StrCalc=.false.
              end if
            case ('Conf-sort')
              if(tokens(j+1) == 'TRUE') then
                cnf_srt=.true.
              elseif(tokens(j+1) == 'FALSE') then
                cnf_srt=.false.
              end if
            case ('Indvl-ext')
              if(tokens(j+1) == 'TRUE') then
                indvlext=.true.
              elseif(tokens(j+1) == 'FALSE') then
                indvlext=.false.
              end if
            case ('resltn-idx')
              call value(tokens(j+1),residx,ios)
            case ('Timer-rep')
              if(tokens(j+1) == 'TRUE') then
                TimerA=.true.
              elseif(tokens(j+1) == 'FALSE') then
                TimerA=.false.
              end if
            case ('Rg')
              if(tokens(j+1) == 'TRUE') then
                RgCalc=.true.
              elseif(tokens(j+1) == 'FALSE') then
                RgCalc=.false.
              end if
            case ('cosTh')
              if(tokens(j+1) == 'TRUE') then
                cosThCalc=.true.
                cosmode=trim(adjustl(tokens(j+2)))
              elseif(tokens(j+1) == 'FALSE') then
                cosThCalc=.false.
              end if
            case ('Dumpstr')
              if(tokens(j+1) == 'TRUE') then
                DumpstrCalc=.true.
                call value(tokens(j+2),nstr,ios)
                allocate(dumpstr(nstr))
                do k=1, nstr
                  call value(tokens(j+2+k),dumpstr(k),ios)
                end do
              elseif(tokens(j+1) == 'FALSE') then
                DumpstrCalc=.false.
              end if
          end select
        end do ! j
      end if ! ntokens
    end do ef

    close(u1)

  end subroutine read_dlt

  subroutine prcs_inp(id,p)

    use :: arry_mod, only: linspace,logspace

    integer :: id,p,iWi,idt,iAdjSeq
    real(wp),parameter :: PI=4*atan(1._wp)
    real(wp) :: factor,dttemp2,delAdjSeq


    allocate(Wi(nWi),Pe(nWi),dttemp1(ndt),dt(nWi,ndt),dt_tmp(nWi,ndt))
    allocate(ntime(nWi,ndt),itime_AdjSeq(nWi,ndt,nAdjSeq))
    if (WiSpacing == 'Log') then
      call logspace(Wii,Wif,Wi)
    elseif (WiSpacing == 'Linear') then
      call linspace(Wii,Wif,Wi)
    end if
    nbead=nseg+1
    nbeadx3=nbead*3;nsegx3=nseg*3
    ntotang=nchain*(nbead-2) ! used for calculation of <cosTh>
    if (dstarCalc == 'Kumar') dstar=dstar*zstar**(1.0_wp/5)

    select case (LambdaMethod)
      case ('Tanner')
        lambda=(b*nseg+7)/(15*b*nseg)*b/(b+5)*(2*(nseg+1)**2+7-&
               12*((nseg+1)**2+1)/real(nseg+1,kind=wp)/(b+7))
      case ('Zimm')
        if (hstar /= 0.0_wp) then
          ! lambda_1:
          lambda=nbead**1.5*1.22_wp/(hstar*PI**2)
          ! lambda_eta:
!         lambda=2.39_wp*(nbead**1.5*1.22_wp/(hstar*PI**2))
        end if
      case ('Rouse')
        ! lambda_1:
        lambda=0.5*1.0/(sin(PI/(2*nbead)))**2
        ! lambda_eta:
!       lambda=1.64_wp*0.5*1.0/(sin(PI/(2*nbead)))**2
      case ('Prabhakar')
        if (hstar /= 0.0_wp) then
          lambda=0.42_wp*nseg**(1.5)*sqrt(12.0_wp)/(hstar*PI**1.5)
        end if
      case ('Self')
      case default
        if (id == 0) print '(" Error: Incorrect Rel-Model.")'
        stop
    end select

    if (id == 0) then
      print *
      print '(" Longest Relaxation Time:",f10.4)',lambda
    end if

    if (lambda == 0._wp) then
      if (id == 0) print '(" Error: Incorrect LambdaMethod.")'
      stop
    end if

    Pe(:)=Wi(:)/lambda
    if (id == 0) then
      print *
      print '(7x,a)', 'Wi            Pe            dt          ntime'
      print '(7x,a)', '---------------------------------------------'
    end if
    if (dtCalc == 'Self') then
      if (dtSpacing == 'Log') then
        call logspace(dti,dtf,dttemp1)
      elseif (dtSpacing == 'Linear') then
        call linspace(dti,dtf,dttemp1)
      end if
    endif
    do iWi=1, nWi
      if (Wi(iWi) < 0.5_wp) then
        if (dtCalc == 'Hsieh') then
          dt(iWi,:)=lambda/nbead**2
        elseif (dtCalc == 'Zimm') then
          dt(iWi,:)=lambda*hstar/(3*nbead**1.5)
        elseif (dtCalc == 'Self') then
          dt(iWi,:)=dttemp1
        endif
      else
        if (dtCalc == 'Hsieh') then
          dt(iWi,:)=lambda/(nbead**2*2*Wi(iWi))
        elseif (dtCalc == 'Zimm') then
          dt(iWi,:)=lambda*hstar/(3*nbead**1.5*2*Wi(iWi))
        elseif (dtCalc == 'Self') then
          if (dtScale) then
            dt(iWi,:)=dttemp1/(2*Wi(iWi))
          else
            dt(iWi,:)=dttemp1
          end if
        end if
      end if
!     To make the dt's in appropriate format
      do idt=1, ndt
        factor=1.0_wp
        do while (floor(factor*dt(iWi,idt)) <= 100)
          factor=factor*10
        end do
        dttemp2=floor(factor*dt(iWi,idt))
        dt(iWi,idt)=dttemp2/factor
        dt_tmp(iWi,idt)=dt(iWi,idt) ! It is used to store the original value of dt.
        if (Adjust_dt) then
          ! Calculating the ntime based on AdjSeq and AdjFact:
          do iAdjSeq=1, nAdjSeq
            ! itime_AdjSeq specifies the itime at which time steps are adjusted.
            if (iAdjSeq.ne.1) then
              delAdjSeq=AdjSeq(iAdjSeq)-AdjSeq(iAdjSeq-1)
              itime_AdjSeq(iWi,idt,iAdjSeq)=itime_AdjSeq(iWi,idt,iAdjSeq-1) + &
                     ceiling( delAdjSeq*lambda/(dt(iWi,idt)*AdjFact(iAdjSeq)) )
            else
              delAdjSeq=AdjSeq(1)
              itime_AdjSeq(iWi,idt,1)=ceiling( delAdjSeq*lambda/(dt(iWi,idt)*&
                                               AdjFact(1)) )
            end if
          end do ! iAdjSeq
          ! We use dt*AdjFact(nAdjSeq) for the period after AdjSeq(nAdjSeq):
          ntime(iWi,idt)=itime_AdjSeq(iWi,idt,nAdjSeq) + &
           ceiling( (tend-AdjSeq(nAdjSeq))*lambda/(dt(iWi,idt)*AdjFact(nAdjSeq)) )
        else
          ntime(iWi,idt)=ceiling(tend*lambda/dt(iWi,idt))
        end if
        if (id == 0) then
          write(*,'(f14.7,1x,f14.7,1x,e10.2,1x,i10)') Wi(iWi),Pe(iWi),dt(iWi,idt),ntime(iWi,idt)
        end if
      end do ! idt
    end do ! iWi

    qmax=sqrt(b) ! Segmental max extension
    q_l=qr_l*qmax
    qr_lto2=qr_l*qr_l
    ! Figure out how many chains are going to each process:
    npchain = (nchain - 1) / p + 1
    ! Checking the number of processors:
    if ((mod(nchain,p) /= 0) .and. (id == 0)) then
      print '(" No. processors not evenly chosen; mod(nchain,p)/=0.")'
      print '(" chains per processor: ",i0)',npchain
    end if
    ! npchain=int(nchain/p)

    ! Specifying 'RWS','WLC-UD', or 'WLC-GEN' parameters:
    select case (ForceLaw)
      case ('RWS')
        RWS_C=3-10._wp/(3*RWS_v)+10._wp/(27*RWS_v*RWS_v)
        RWS_D=1+2._wp/(3*RWS_v)+10._wp/(27*RWS_v*RWS_v)
      case ('WLC_UD')
        WLC_A=3._wp/32-3/(8*WLC_v)-3/(2*WLC_v**2)
        WLC_B=(13._wp/32+0.4086_wp/WLC_v-14.79_wp/(4*WLC_v**2))/ &
              (1-4.225_wp/(2*WLC_v)+4.87_wp/(4*WLC_v**2))
      case ('WLC_GEN')
        WLC_A=3._wp/32-3/(8*WLC_v)-3/(2*WLC_v**2)
!        WLC_B=(13._wp/32+4.6765_wp/(2*WLC_v)+13.7004_wp/(4*WLC_v**2))/ &
!              (1+1.7228_wp/(2*WLC_v)-3.3722_wp/(4*WLC_v**2))
!        WLC_C=(1-1.4538_wp*(2*WLC_v)+0.4056_wp*(4*WLC_v**2))/(2*WLC_&
!               &v- 0.9709_wp*(4*WLC_v**2)+0.2296_wp*(8*WLC_v**3))
        WLC_B=(13._wp/32+3.4719_wp/(2*WLC_v)+2.5064_wp/(4*WLC_v**2))/ &
              (1-1.2906_wp/(2*WLC_v)+0.6482_wp/(4*WLC_v**2))
        WLC_C=(1-1.2370_wp*(2*WLC_v)+0.8105_wp*(4*WLC_v**2))/(2*WLC_&
               &v- 1.0243_wp*(4*WLC_v**2)+0.4595_wp*(8*WLC_v**3))
    end select

  end subroutine

  subroutine del_inp

    deallocate(Wi,Pe,dttemp1,dt,dt_tmp,ntime,itime_AdjSeq)

  end subroutine del_inp

end module inp_dlt
