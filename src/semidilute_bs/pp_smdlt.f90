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
module pp_smdlt
 
  use :: prcn_mod
  
  implicit none

  save ! These parameters remains in the scope of the main program
  private ! All defined prameters are private

  !> The counter for time averaging
  integer :: jcount
  !> The counter for average iteration in Lanczos
  integer :: kcount
  !> The counter for run averaging 
  integer :: lcount
  !> @name Group1
  !! File units
  !> @{
  integer :: u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20
  integer :: u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38
  integer :: u39
  !> @}
  !> @name Group1
  !! The variables regarding end to end distance
  !> @{
  real(wp) :: sqqetoeAve,qetoeAve,sdsqqetoeAve,sdqetoeAve
  real(wp) :: tAvesqqetoeAve,tAveqetoeAve,tAvesdsqqetoeAve,tAvesdqetoeAve
  real(wp) :: sdtAvesqqetoeAve,sdtAveqetoeAve,rAvetAvesqqetoe,rAvetAveqetoe
  real(wp) :: rAvesdtAvesqqetoe,rAvesdtAveqetoe,sdrAvetAvesqqetoe,sdrAvetAveqetoe
  real(wp) :: rAvetAvesqqetoeTot,rAvetAveqetoeTot,rAvesdtAvesqqetoeTot
  real(wp) :: rAvesdtAveqetoeTot,sdrAvetAvesqqetoeTot,sdrAvetAveqetoeTot
  !> @}
  !> @name Group2
  !! The variables regarding stress tensor: 
  !> @{
  real(wp) :: tauxx,tauxy,tauyy,tauzz,tauxxyy,tauyyzz
  real(wp) :: tAvetauxy,tAvetauxxyy,tAvetauyyzz,sdtAvetauxy,sdtAvetauxxyy,sdtAvetauyyzz
  real(wp) :: rAvetAvetauxy,rAvetAvetauxxyy,rAvetAvetauyyzz,rAvesdtAvetauxy,rAvesdtAvetauxxyy
  real(wp) :: rAvesdtAvetauyyzz,sdrAvetAvetauxy,sdrAvetAvetauxxyy,sdrAvetAvetauyyzz
  real(wp) :: rAvetAvetauxyTot,rAvetAvetauxxyyTot,rAvetAvetauyyzzTot,rAvesdtAvetauxyTot
  real(wp) :: rAvesdtAvetauxxyyTot,rAvesdtAvetauyyzzTot,sdrAvetAvetauxyTot
  real(wp) :: sdrAvetAvetauxxyyTot,sdrAvetAvetauyyzzTot
  !> @}
  !> @name Group3
  !! The variables regarding diffusion coefficient:
  !> @{
  real(wp) :: DcmAve,sdDcmAve,MSDAve,sdMSDAve
  real(wp) :: tAveDcmAve,tAvesdDcmAve,sdtAveDcmAve
  real(wp) :: rAvetAveDcm,rAvesdtAveDcm,sdrAvetAveDcm
  real(wp) :: rAvetAveDcmTot,rAvesdtAveDcmTot,sdrAvetAveDcmTot
  !> @}
  !> @name Group4
  !! The variables ragarding radius of gyration:
  !> @{
  real(wp) :: RgSqAve,sdRgSqAve,AspherAve,sdAspherAve,Aspher,RgSqTensEVbar,traceRgSqTens
  real(wp) :: tAveRgSqAve,tAvesdRgSqAve,sdtAveRgSqAve,tAveAspherAve,tAvesdAspherAve
  real(wp) :: sdtAveAspherAve
  real(wp) :: rAvetAveRgSq,rAvesdtAveRgSq,sdrAvetAveRgSq,rAvetAveAspher,rAvesdtAveAspher
  real(wp) :: sdrAvetAveAspher,rAvetAveRgSqTot,rAvesdtAveRgSqTot,sdrAvetAveRgSqTot
  real(wp) :: rAvetAveAspherTot,rAvesdtAveAspherTot,sdrAvetAveAspherTot
  !> @}
  !> @name Group5
  !! The variables regarding iteration number of Lanczos method:
  !> @{
  real(wp) :: mAve,sdmAve,LAve,sdLAve
  real(wp) :: tAvem,sdtAvem
  !> @}
  !> @name Group6
  !! Arrays for averaging the results of all runs
  !> @{
  real(wp),allocatable,dimension(:) :: rAvesqqetoe,rAveqetoe,rAvetauxy,rAvetauxxyy
  real(wp),allocatable,dimension(:) :: rAvetauyyzz,rAveRgSq,sdrAvesqqetoe,sdrAveqetoe
  real(wp),allocatable,dimension(:) :: sdrAvetauxy,sdrAvetauxxyy,sdrAvetauyyzz,sdrAveRgSq
  real(wp),allocatable,dimension(:) :: rAveMSD,sdrAveMSD,rAvesqqetoeTot,rAveqetoeTot
  real(wp),allocatable,dimension(:) :: rAvetauxyTot,rAvetauxxyyTot,rAvetauyyzzTot
  real(wp),allocatable,dimension(:) :: rAveRgSqTot,sdrAvesqqetoeTot,sdrAveqetoeTot
  real(wp),allocatable,dimension(:) :: sdrAvetauxyTot,sdrAvetauxxyyTot,sdrAvetauyyzzTot
  real(wp),allocatable,dimension(:) :: sdrAveRgSqTot,rAveMSDTot,sdrAveMSDTot
  !> @}
  !> A counter used in correlation function
  integer :: kchk
  !> The array to store the Ree correlation function
  real(wp),allocatable,protected :: cf_ree(:,:,:,:)
  !> The array to store the stress correlation function
  real(wp),allocatable,target,protected :: cf_rg(:,:,:)
  

  public :: init_pp        ,&
            data_time_init ,&
            data_run_init  ,&
            material_func  ,&
            CorrFcn        ,&
            del_pp

  public :: StrCalc,StrPr_mode,RgCalc,CoMDiff,AveIterCalc,&
               CorFun,CF_mode,cf_ree,cf_rg,kchk,add_cmb
!  protected :: StrCalc,StrPr_mode,RgCalc,CoMDiff,AveIterCalc,&
!               CorFun,CF_mode

  !> If the Comb chains exist
  logical,protected :: add_cmb
  !> If the stress calculation is desired
  logical,protected :: StrCalc
  !> If the timing report is intended
  logical,protected :: doTiming
  !> If the center of mass diffusivity calculation is intended
  logical,protected :: CoMDiff
  !> If the calculation of radius of gyration is intended
  logical,protected :: RgCalc
  !> If the calculation of average iteration in calculation of D is intended
  logical,protected :: AveIterCalc
  !> If the calculation of correlation function is intended
  logical,protected :: CorFun
  !> The type of correlation function
  character(len=10),protected :: CF_mode
  !> The printing type for stress, stress or viscosity
  character(len=10),protected :: StrPr_mode
  

  
contains

  !> Initializes the processing of the data in pp_smdlt
  !! \myrank the rank of the process
  subroutine init_pp(myrank,nrun)

    use :: strg_mod
    use,intrinsic :: iso_fortran_env
    use :: hi_mod, only: DecompMeth
    use :: flow_mod, only: FlowType
    use :: mpi
    !include 'mpif.h'
  
    integer,intent(in) :: myrank,nrun
    integer :: stat1,stat2,stat3,stat4,stat5,stat6,stat7,stat8,stat9,stat10
    integer :: stat11,stat12,stat13,stat14,stat15,stat16,stat17,stat18
    character(len=1024) :: file1,file2,file3,file4,file5,file6,file7,file8,file9
    character(len=1024) :: file10,file11,file12,file13,file14,file15,file16,file17,file18
    character(len=1024) :: format_str,line
    character(len=100) :: tokens(50)
    integer :: i,j,ntokens,u0,il,stat,ierr
	write(*,*) "module:pp_smdlt:init_pp"
    ! Default values:
    StrCalc=.false.;StrPr_mode='Visc'
    doTiming=.false.
    CoMDiff=.false.
    RgCalc=.false.
    CorFun=.false.;CF_mode='Rg'

    open (newunit=u0,action='read',file='input.dat',status='old')
    il=1
ef: do
      read(u0,'(A)',iostat=stat) line
      if (stat == iostat_end) then
        exit ef ! end of file
      elseif (stat > 0) then
        print '(" io_mod: Error reading line ",i0, " Process ID ",i0," in pp")',il,myrank
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
            case ('Conf-anal')
              if(tokens(j+1) == 'TRUE') then
                StrCalc=.true.
                StrPr_mode=trim(adjustl(tokens(j+2)))
              elseif(tokens(j+1) == 'FALSE') then
                StrCalc=.false.
              else
                print '(" Inconsistent Conf-anal.")'
                stop
              end if
            case ('Timer-rep')
              if(tokens(j+1) == 'TRUE') then
                doTiming=.true.
              elseif(tokens(j+1) == 'FALSE') then
                doTiming=.false.
              else
                print '(" Inconsistent Timing-rep.")'
                stop
              end if
            case ('CoM-diff')
              if(tokens(j+1) == 'TRUE') then
                CoMDiff=.true.
              elseif(tokens(j+1) == 'FALSE') then
                CoMDiff=.false.
              else
                print '(" Inconsistent CoM-diff.")'
                stop
              end if
            case ('Rg-calc')
              if(tokens(j+1) == 'TRUE') then
                RgCalc=.true.
              elseif(tokens(j+1) == 'FALSE') then
                RgCalc=.false.
              else
                print '(" Inconsistent Rg-calc.")'
                stop
              end if
            case ('AveIter-rep')
              if(tokens(j+1) == 'TRUE') then
                AveIterCalc=.true.
              elseif(tokens(j+1) == 'FALSE') then
                AveIterCalc=.false.
              else
                print '(" Inconsistent AveIter-rep.")'
                stop
              end if
            case ('CorrFcn-calc')
              if(tokens(j+1) == 'TRUE') then
                CorFun=.true.
                CF_mode=trim(adjustl(tokens(j+2)))
              elseif(tokens(j+1) == 'FALSE') then
                CorFun=.false.
              else
                print '(" Inconsistent CorrFcn-calc.")'
                stop
              end if
          end select
        end do ! j
      end if ! ntokens
    end do ef
    close(u0)
    if (CorFun) then
      kchk=0
      ! This file is used to save the correlation function information:
      select case (CF_mode)
        case ('Ree')
          open (newunit=u39,file='data/Ct-Ree.dat',status='replace',position='append')
        case ('Rg')
          open (newunit=u39,file='data/Ct-Rg.dat',status='replace',position='append')
      end select
      StrCalc=.true. 
    end if

    format_str="(A,i0.3,'.dat')"
    ! Deleting all files in "Ree,Rg,CoM,rheol" directories.
    if (myrank == 0) then
      i=0
      do
        if (StrCalc) then
          write(file1,format_str)'data/Ree/ReeSqvsTimeID',i
          write(file2,format_str)'data/Ree/ReeRelvsTimeID',i
          write(file8,format_str)'data/Ree/tAveReeSqID',i
          write(file9,format_str)'data/Ree/tAveReeRelID',i
          select case (StrPr_mode)
            case ('Visc')
              write(file3,format_str)'data/rheol/EtavsTimeID',i
              write(file4,format_str)'data/rheol/Psi1vsTimeID',i
              write(file5,format_str)'data/rheol/Psi2vsTimeID',i
              write(file10,format_str)'data/rheol/tAveEtavsWiID',i
              write(file11,format_str)'data/rheol/tAvePsi1vsTimeID',i
              write(file12,format_str)'data/rheol/tAvePsi2vsTimeID',i
            case ('Stress')
              write(file3,format_str)'data/rheol/TauxyvsTimeID',i
              write(file4,format_str)'data/rheol/N1vsTimeID',i
              write(file5,format_str)'data/rheol/N2vsTimeID',i
              write(file10,format_str)'data/rheol/tAveTauxyvsWiID',i
              write(file11,format_str)'data/rheol/tAveN1vsWiID',i
              write(file12,format_str)'data/rheol/tAveN2vsWiID',i
          end select
          
        end if
        if (RgCalc) then
          write(file6,format_str)'data/Rg/RgSqvsTimeID',i
          write(file13,format_str)'data/Rg/tAveRgSqID',i
          write(file14,format_str)'data/Rg/tAveAsphericityID',i
        end if 
        if ((FlowType == 'Equil').and.CoMDiff) then
          write(file7,format_str)'data/CoM/DcmvsTimeID',i
          write(file18,format_str)'data/CoM/MSDvsTimeID',i
          write(file15,format_str)'data/CoM/tAveDcmID',i
        end if
        if (AveIterCalc) then
          write(file16,format_str)'data/m/mvsTimeID',i
          write(file17,format_str)'data/m/tAvemID',i
        end if

        open(newunit=u1,iostat=stat1,file=file1,status='old')
        open(newunit=u2,iostat=stat2,file=file2,status='old')
        open(newunit=u3,iostat=stat3,file=file3,status='old')
        open(newunit=u25,iostat=stat4,file=file4,status='old')
        open(newunit=u26,iostat=stat5,file=file5,status='old')
        open(newunit=u18,iostat=stat6,file=file6,status='old')
        open(newunit=u19,iostat=stat7,file=file7,status='old')

        open(newunit=u16,iostat=stat8,file=file8,status='old')
        open(newunit=u4,iostat=stat9,file=file9,status='old')
        open(newunit=u5,iostat=stat10,file=file10,status='old')
        open(newunit=u6,iostat=stat11,file=file11,status='old')
        open(newunit=u7,iostat=stat12,file=file12,status='old')
        open(newunit=u12,iostat=stat13,file=file13,status='old')
        open(newunit=u13,iostat=stat14,file=file14,status='old')
        open(newunit=u10,iostat=stat15,file=file15,status='old')
        open(newunit=u14,iostat=stat16,file=file16,status='old')
        open(newunit=u15,iostat=stat17,file=file17,status='old')
        open(newunit=u37,iostat=stat18,file=file18,status='old')

        if (stat1 == 0) close((u1),status='delete')
        if (stat2 == 0) close((u2),status='delete')
        if (stat3 == 0) close((u3),status='delete')
        if (stat4 == 0) close((u25),status='delete')
        if (stat5 == 0) close((u26),status='delete')
        if (stat6 == 0) close((u18),status='delete')
        if (stat7 == 0) close((u19),status='delete')

        if (stat8 == 0) close((u16),status='delete')
        if (stat9 == 0) close((u4),status='delete')
        if (stat10 == 0) close((u5),status='delete')
        if (stat11 == 0) close((u6),status='delete')
        if (stat12 == 0) close((u7),status='delete')
        if (stat13 == 0) close((u12),status='delete')
        if (stat14 == 0) close((u13),status='delete')
        if (stat15 == 0) close((u10),status='delete')
        if (stat16 == 0) close((u14),status='delete')
        if (stat17 == 0) close((u15),status='delete')
        if (stat18 == 0) close((u37),status='delete')

        if ((stat1 /= 0) .and. (stat2 /= 0) .and. (stat3 /= 0) .and. &
            (stat4 /= 0) .and. (stat5 /= 0) .and. (stat6 /= 0) .and. &
            (stat7 /= 0) .and. (stat8 /= 0) .and. (stat9 /= 0) .and. &
            (stat10 /= 0) .and. (stat11 /= 0) .and. (stat12 /= 0) .and. &
            (stat13 /= 0) .and. (stat14 /= 0) .and. (stat15 /= 0) .and. &
            (stat16 /= 0) .and. (stat17 /= 0) .and. (stat18 /= 0) ) exit
        i=i+1
      end do
    end if ! myrank
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if (StrCalc) then
      write(file1,format_str)'data/Ree/ReeSqvsTimeID',myrank
      write(file2,format_str)'data/Ree/ReeRelvsTimeID',myrank
      open (newunit=u1,file=trim(adjustl(file1)),status='replace',position='append')
      write(u1,*) "# Wi, dt, Time, <Ree.Ree>, sd<Ree.Ree> #"
      write(u1,*) "# ------------------------------------ #"
      open (newunit=u2,file=trim(adjustl(file2)),status='replace',position='append')
      write(u2,*) "# Wi, dt, Strain, <|Ree|>/|Ree|_max, sd<Ree^2>/|Ree^2|_max #" 
      write(u2,*) "# -------------------------------------------------------- #"

      select case (StrPr_mode)
        case ('Visc')
          write(file3,format_str)'data/rheol/EtavsTimeID',myrank
          write(file4,format_str)'data/rheol/Psi1vsTimeID',myrank
          write(file5,format_str)'data/rheol/Psi2vsTimeID',myrank
          open (newunit=u3,file=trim(adjustl(file3)),status='replace',position='append')
          write(u3,*) "# Wi, dt, Time, Etap=(Taupxy)/Pe #"
          write(u3,*) "# ------------------------------ #"
          open (newunit=u25,file=trim(adjustl(file4)),status='replace',position='append')
          write(u25,*) "# Wi, dt, Time, Psi1=-(Taupxx-Taupyy)/Pe^2 #"
          write(u25,*) "# ---------------------------------------- #"
          open (newunit=u26,file=trim(adjustl(file5)),status='replace',position='append')
          write(u26,*) "# Wi, dt, Time, Psi2=-(Taupyy-Taupzz)/Pe^2 #"
          write(u26,*) "# ---------------------------------------- #"
        case ('Stress')
          write(file3,format_str)'data/rheol/TauxyvsTimeID',myrank
          write(file4,format_str)'data/rheol/N1vsTimeID',myrank
          write(file5,format_str)'data/rheol/N2vsTimeID',myrank
          open (newunit=u3,file=trim(adjustl(file3)),status='replace',position='append')
          write(u3,*) "# Wi, dt, Time, Taupxy #"
          write(u3,*) "# -------------------- #"
          open (newunit=u25,file=trim(adjustl(file4)),status='replace',position='append')
          write(u25,*) "# Wi, dt, Time, Taupxx-Taupyy #"
          write(u25,*) "# --------------------------- #"
          open (newunit=u26,file=trim(adjustl(file5)),status='replace',position='append')
          write(u26,*) "# Wi, dt, Taupyy-Taupzz #"
          write(u26,*) "# --------------------- #"
      end select
    end if ! StrCalc
    if (RgCalc) then
      write(file6,format_str)'data/Rg/RgSqvsTimeID',myrank
      open (newunit=u18,file=trim(adjustl(file6)),status='replace',position='append')
      write(u18,*) "# Wi, dt, Time, <Rg.Rg>, sd<Rg.Rg> #"
      write(u18,*) "# -------------------------------- #"
    end if
    if ((FlowType == 'Equil').and.CoMDiff) then
      write(file7,format_str)'data/CoM/DcmvsTimeID',myrank
      write(file18,format_str)'data/CoM/MSDvsTimeID',myrank
      open (newunit=u19,file=trim(adjustl(file7)),status='replace',position='append')
      write(u19,*) "# Wi, dt, Time, <Dcm>, sd<Dcm> #"
      write(u19,*) "# ---------------------------- #"
      open (newunit=u37,file=trim(adjustl(file18)),status='replace',position='append')
      write(u37,*) "# Wi, dt, Time, <MSD>, sd<MSD> #"
      write(u37,*) "# ---------------------------- #"
    end if

    if (StrCalc) then
      if (nrun > 1) then
        if (myrank == 0) then
          open (newunit=u20,file='data/rAveReeSqvsTime.dat',status='replace',position='append')
          write(u20,*) "# nbead, Wi, dt, Time, <Ree.Ree>r, sd<Ree.Ree>r #"
          write(u20,*) "# --------------------------------------------- #"
          open (newunit=u21,file='data/rAveRgSqvsTime.dat',status='replace',position='append')
          write(u21,*) "# nbead, Wi, dt, Time, <Rg.Rg>r, sd<Rg.Rg>r #"
          write(u21,*) "# ----------------------------------------- #"
          open (newunit=u23,file='data/rAveReeRelvsTime.dat',status='replace',position='append')
          write(u23,*) "# Wi, dt, Strain, <|Ree|>r/|Ree|_mx, sd<|Ree|>r/|Ree|_mx #" 
          write(u23,*) "# ------------------------------------------------------ #"
          select case (StrPr_mode)
            case ('Visc')
              open (newunit=u24,file='data/rAveEtavsTime.dat',status='replace',position='append')
              write(u24,*) "# Wi, dt, Time, Etap=<Taupxy>r/Pe, sd<Taupxy>r/Pe #"
              write(u24,*) "# ----------------------------------------------- #"
              open (newunit=u27,file='data/rAvePsi1vsTime.dat',status='replace',position='append') 
              write(u27,*) "# Wi, dt, <Psi1=-(Taupxx-Taupyy)/Pe^2>r,  sd<Psi1>r #"
              write(u27,*) "# ------------------------------------------------- #"
              open (newunit=u28,file='data/rAvePsi2vsTime.dat',status='replace',position='append') 
              write(u28,*) "# Wi, dt, <Psi2=-(Taupyy-Taupzz)/Pe^2>r, sd<Psi2>r #"
              write(u28,*) "# ------------------------------------------------ #"
            case ('Stress')
              open (newunit=u24,file='data/rAveTauxyvsTime.dat',status='replace',position='append')
              write(u24,*) "# Wi, dt, Time, <Taupxy>r, sd<Taupxy>r #"
              write(u24,*) "# ------------------------------------ #"
              open (newunit=u27,file='data/rAveN1vsTime.dat',status='replace',position='append') 
              write(u27,*) "# Wi, dt, <N1=(Taupxx-Taupyy)>r, sd<N1>r #"
              write(u27,*) "# -------------------------------------- #"
              open (newunit=u28,file='data/rAveN2vsTime.dat',status='replace',position='append') 
              write(u28,*) "# Wi, dt, <N2=(Taupyy-Taupzz)>r, sd<N2>r #"
              write(u28,*) "# -------------------------------------- #"
          end select  
        end if ! myrank
      end if ! nrun

      write(file8,format_str)'data/Ree/tAveReeSqID',myrank
      write(file9,format_str)'data/Ree/tAveReeRelID',myrank
      open (newunit=u16,file=trim(adjustl(file8)),status='replace',position='append')
      write(u16,*) "# nbead, Wi, dt, <Ree2>t, <sd(Ree2)>t, sd<Ree2>t #"
      write(u16,*) "# ---------------------------------------------- #"            
      open (newunit=u4,file=trim(adjustl(file9)),status='replace',position='append')
      write(u4,*) "# Wi, dt, <<|Ree|>/|Ree|_mx>t, <sd<|Ree|>/|Ree|_mx>t, sd<<|Ree|>/|Ree|_mx>t #"
      write(u4,*) "# ------------------------------------------------------------------------- #"
      select case (StrPr_mode)
        case ('Visc')
          write(file10,format_str)'data/rheol/tAveEtavsWiID',myrank
          write(file11,format_str)'data/rheol/tAvePsi1vsWiID',myrank
          write(file12,format_str)'data/rheol/tAvePsi2vsWiID',myrank
          open (newunit=u5,file=trim(adjustl(file10)),status='replace',position='append')
          write(u5,*) "# Wi, dt, <Etap=-(Taupxy)/Pe>t, sd<Etap>t #"
          write(u5,*) "# --------------------------------------- #"
          open (newunit=u6,file=trim(adjustl(file11)),status='replace',position='append')
          write(u6,*) "# Wi, dt, <Psi1=-(Taupxx-Taupyy)/Pe^2>t, sd<Psi1>t #"
          write(u6,*) "# ------------------------------------------------ #"
          open (newunit=u7,file=trim(adjustl(file12)),status='replace',position='append')
          write(u7,*) "# Wi, dt, <Psi2=-(Taupyy-Taupzz)/Pe^2>t, sd<Psi2>t #"
          write(u7,*) "# ------------------------------------------------ #"
        case ('Stress')
          write(file10,format_str)'data/rheol/tAveTauxyvsWiID',myrank
          write(file11,format_str)'data/rheol/tAveN1vsWiID',myrank
          write(file12,format_str)'data/rheol/tAveN2vsWiID',myrank
          open (newunit=u5,file=trim(adjustl(file10)),status='replace',position='append')
          write(u5,*) "# Wi, dt, <Taupxy>t, sd<Etap>t #"
          write(u5,*) "# ---------------------------- #"
          open (newunit=u6,file=trim(adjustl(file11)),status='replace',position='append')
          write(u6,*) "# Wi, dt, <N1=(Taupxx-Taupyy)>t, sd<N1>t #"
          write(u6,*) "# -------------------------------------- #"
          open (newunit=u7,file=trim(adjustl(file12)),status='replace',position='append')
          write(u7,*) "# Wi, dt, <N2=(Taupyy-Taupzz)>t, sd<N2>t #"
          write(u7,*) "# -------------------------------------- #"
      end select
      if (nrun > 1) then
        if (myrank == 0) then
          open (newunit=u36,file='data/rAveReeSq.dat',status='unknown',position='append')
          write(u36,*) "# nbead, Wi, dt, <<Ree2>t>r, <sd(<Ree2>t)>r, sd<Ree2>r #"
          write(u36,*) "# ---------------------------------------------------- #"            
          open (newunit=u29,file='data/rAveRee.dat',status='unknown',position='append')
          write(u29,*) "# Wi, dt, <<|Qee|>/|Qee|_mx>r, <sd<|Qee|>/|Qee|_mx>r, sd<<|Qee|>/|Qee|_mx>r #"
          write(u29,*) "# ------------------------------------------------------------------------- #"
          select case (StrPr_mode)
            case ('Visc')
              open (newunit=u30,file='data/rAveEtavsWi.dat',status='unknown',position='append')
              write(u30,*) "# Wi, dt, <Etap=-<Taupxy>t/Pe>r, <sd<Taupxy>t/Pe>r, sd<<Etap>t>r #"
              write(u30,*) "# -------------------------------------------------------------- #"
              open (newunit=u31,file='data/rAvePsi1vsWi.dat',status='unknown',position='append') 
              write(u31,*) "# Wi, dt, <Psi1=-<Taupxx-Taupyy>t/Pe^2>r, <sd<Taupxx-Taupyy>t/Pe^2>r, sd<<Psi1>t>r #"
              write(u31,*) "# -------------------------------------------------------------------------------- #"
              open (newunit=u32,file='data/rAvePsi2vsWi.dat',status='unknown',position='append') 
              write(u32,*) "# Wi, dt, <Psi2=-<Taupyy-Taupzz>t/Pe^2>r, <sd<Taupyy-Taupzz>t/Pe^2>r, sd<<Psi2>t>r #"
              write(u32,*) "# -------------------------------------------------------------------------------- #"
            case ('Stress')
              open (newunit=u30,file='data/rAveTauxyvsWi.dat',status='unknown',position='append')
              write(u30,*) "# Wi, dt, <<Taupxy>t>r, <sd<Taupxy>t>r, sd<<Taupxy>t>r #"
              write(u30,*) "# ---------------------------------------------------- #"
              open (newunit=u31,file='data/rAveN1vsWi.dat',status='unknown',position='append') 
              write(u31,*) "# Wi, dt, <N1=<Taupxx-Taupyy>t>r, <sd<Taupxx-Taupyy>t>r, sd<<N1>t>r #"
              write(u31,*) "# ----------------------------------------------------------------- #"
              open (newunit=u32,file='data/rAveN2vsWi.dat',status='unknown',position='append') 
              write(u32,*) "# Wi, dt, <N2=<Taupyy-Taupzz>t>r, <sd<Taupyy-Taupzz>t>r, sd<<N2>t>r #"
              write(u32,*) "# ----------------------------------------------------------------- #"
          end select
        end if ! myrank
      end if ! nrun
    end if ! StrCalc
    if (RgCalc) then
      write(file13,format_str)'data/Rg/tAveRgSqID',myrank
      write(file14,format_str)'data/Rg/tAveAsphericityID',myrank
      open (newunit=u12,file=trim(adjustl(file13)),status='replace',position='append')
      write(u12,*) "# nbead, Wi, dt, <Rg2>t, <sd(Rg2)>t, sd<Rg2>t #"
      write(u12,*) "# ------------------------------------------- #"            
      open (newunit=u13,file=trim(adjustl(file14)),status='replace',position='append')
      write(u13,*) "# nbead, Wi, dt, <Asphericity>t, <sd(Asphericity)>t, sd<Asphericity>t #"
      write(u13,*) "# ------------------------------------------------------------------- #"
      if (myrank == 0) then
        open (newunit=u34,file='data/rAveRgSq.dat',status='unknown',position='append')
        write(u34,*) "# nbead, Wi, dt, <Rg2>r, <sd(Rg2)>r, sd<Rg2>r #"
        write(u34,*) "# ------------------------------------------- #"            
        open (newunit=u35,file='data/rAveAsphericity.dat',status='unknown',position='append')
        write(u35,*) "# nbead, Wi, dt, <Asphericity>r, <sd(Asphericity)>r, sd<Asphericity>r #"
        write(u35,*) "# ------------------------------------------------------------------- #"            
      end if
    end if
    if ((FlowType == 'Equil').and.CoMDiff) then
      write(file15,format_str)'data/CoM/tAveDcmID',myrank
      open (newunit=u10,file=trim(adjustl(file15)),status='replace',position='append')
      write(u10,*) "# nbead, dt, <Dcm>t, <sd(Dcm)>t, sd<Dcm>t #"
      write(u10,*) "# --------------------------------------- #"
      if (myrank == 0) then
        open (newunit=u33,file='data/rAveDcm.dat',status='unknown',position='append')
        write(u33,*) "# nbead, dt, <Dcm>r, <sd(Dcm)>r, sd<Dcm>r #"
        write(u33,*) "# --------------------------------------- #"            
        open (newunit=u38,file='data/rAveMSD.dat',status='replace',position='append')
        write(u38,*) "# nbead, dt, Time, <MSD>r, sd<MSD>r #"
        write(u38,*) "# --------------------------------- #"            
      end if
    end if
    if (AveIterCalc) then
      if (DecompMeth == 'Lanczos') then
        write(file16,format_str)'data/m/tAvemID',myrank
        write(file17,format_str)'data/m/mvsTimeID',myrank
        open (newunit=u14,file=trim(adjustl(file16)),status='replace',position='append')
        write(u14,*) "# nbead, Wi, dt, <m>t, sd<m>t #"
        write(u14,*) "# ----------------------------#"
        open (newunit=u15,file=trim(adjustl(file17)),status='replace',position='append')
        write(u15,*) "# FlowType, Wi, dt, Time, m #"
        write(u15,*) "# --------------------------#"
      end if ! DecompMeth
    end if

  end subroutine init_pp

  !> Initializes the processing of the data for material functions at initial time step
  !! \myrank the rank of the process
  subroutine data_time_init(myrank)
 
    use :: hi_mod, only: DecompMeth
    use :: flow_mod, only: FlowType

    integer,intent(in) :: myrank
	write(*,*) "module:pp_smdlt:data_time_init"
    jcount=0;kcount=0;lcount=0
    ! Variables for time averaging:
    if (StrCalc) then
      ! related to end to end distance qetoe.
      tAvesqqetoeAve=0._wp;tAveqetoeAve=0._wp
      tAvesdsqqetoeAve=0._wp;tAvesdqetoeAve=0._wp
      sdtAvesqqetoeAve=0._wp;sdtAveqetoeAve=0._wp
      ! related to stress tensor.
      tAvetauxy=0._wp;tAvetauxxyy=0._wp;tAvetauyyzz=0._wp
      sdtAvetauxy=0._wp;sdtAvetauxxyy=0._wp;sdtAvetauyyzz=0._wp
    end if
    ! Diffusion coefficient of center of mass.
    if ((FlowType == 'Equil').and.CoMDiff) then
      tAveDcmAve=0._wp;tAvesdDcmAve=0._wp;sdtAveDcmAve=0._wp
    end if
    ! Radius of gyration.
    if (RgCalc) then
      tAveRgSqAve=0._wp;tAvesdRgSqAve=0._wp;sdtAveRgSqAve=0._wp
      tAveAspherAve=0._wp;tAvesdAspherAve=0._wp;sdtAveAspherAve=0._wp
    end if           
    if (AveIterCalc) then
      if (DecompMeth == 'Lanczos') then
        tAvem=0._wp;sdtAvem=0._wp
      end if
    end if

  end subroutine data_time_init

  !> Initializes the processing of the data for material functions at initial run
  !! \myrank the rank of the process
  subroutine data_run_init(myrank,ntime,tgap,ndmp,nprun,ntotchain,ntotbeadx3)
 
    use :: hi_mod, only: DecompMeth
    use :: flow_mod, only: FlowType

    integer,intent(in) :: myrank,ntime,tgap,ndmp,nprun,ntotchain,ntotbeadx3
!    integer :: ntime,tgap
	write(*,*) "module:pp_smdlt:data_run_init"
    ! Variables for run averaging:
    if (StrCalc) then
      ! related to end to end distance qetoe.
      rAvesqqetoe=0._wp;rAveqetoe=0._wp
      sdrAvesqqetoe=0._wp;sdrAveqetoe=0._wp
      rAvetAvesqqetoe=0._wp;rAvetAveqetoe=0._wp
      rAvesdtAvesqqetoe=0._wp;rAvesdtAveqetoe=0._wp
      sdrAvetAvesqqetoe=0._wp;sdrAvetAveqetoe=0._wp
      ! related to stress tensor.
      rAvetauxy=0._wp;rAvetauxxyy=0._wp;rAvetauyyzz=0._wp
      sdrAvetauxy=0._wp;sdrAvetauxxyy=0._wp;sdrAvetauyyzz=0._wp
      rAvetAvetauxy=0._wp;rAvetAvetauxxyy=0._wp;rAvetAvetauyyzz=0._wp
      rAvesdtAvetauxy=0._wp;rAvesdtAvetauxxyy=0._wp;rAvesdtAvetauyyzz=0._wp
      sdrAvetAvetauxy=0._wp;sdrAvetAvetauxxyy=0._wp;sdrAvetAvetauyyzz=0._wp
      if (allocated(rAvesqqetoe))   deallocate(rAvesqqetoe)
      if (allocated(rAveqetoe))     deallocate(rAveqetoe)
      if (allocated(rAvetauxy))     deallocate(rAvetauxy)
      if (allocated(rAvetauxxyy))   deallocate(rAvetauxxyy)
      if (allocated(rAvetauyyzz))   deallocate(rAvetauyyzz)
      if (allocated(sdrAvesqqetoe)) deallocate(sdrAvesqqetoe)
      if (allocated(sdrAveqetoe))   deallocate(sdrAveqetoe)
      if (allocated(sdrAvetauxy))   deallocate(sdrAvetauxy)
      if (allocated(sdrAvetauxxyy)) deallocate(sdrAvetauxxyy)
      if (allocated(sdrAvetauyyzz)) deallocate(sdrAvetauyyzz)
      allocate(rAvesqqetoe(ndmp),rAveqetoe(ndmp),rAvetauxy(ndmp),rAvetauxxyy(ndmp))
      allocate(rAvetauyyzz(ndmp),sdrAvesqqetoe(ndmp),sdrAveqetoe(ndmp))
      allocate(sdrAvetauxy(ndmp),sdrAvetauxxyy(ndmp),sdrAvetauyyzz(ndmp))
      if (myrank == 0) then
        if (allocated(rAvesqqetoeTot))   deallocate(rAvesqqetoeTot)
        if (allocated(rAveqetoeTot))     deallocate(rAveqetoeTot)
        if (allocated(rAvetauxyTot))     deallocate(rAvetauxyTot)
        if (allocated(rAvetauxxyyTot))   deallocate(rAvetauxxyyTot)
        if (allocated(rAvetauyyzzTot))   deallocate(rAvetauyyzzTot)
        if (allocated(sdrAvesqqetoeTot)) deallocate(sdrAvesqqetoeTot)
        if (allocated(sdrAveqetoeTot))   deallocate(sdrAveqetoeTot)
        if (allocated(sdrAvetauxyTot))   deallocate(sdrAvetauxyTot)
        if (allocated(sdrAvetauxxyyTot)) deallocate(sdrAvetauxxyyTot)
        if (allocated(sdrAvetauyyzzTot)) deallocate(sdrAvetauyyzzTot)
        allocate(rAvesqqetoeTot(ndmp),rAveqetoeTot(ndmp),rAvetauxyTot(ndmp))
        allocate(rAvetauxxyyTot(ndmp),rAvetauyyzzTot(ndmp),sdrAvesqqetoeTot(ndmp))
        allocate(sdrAveqetoeTot(ndmp),sdrAvetauxyTot(ndmp),sdrAvetauxxyyTot(ndmp))
        allocate(sdrAvetauyyzzTot(ndmp))
      end if ! myrank
    end if
    ! Diffusion coefficient of center of mass.
    if ((FlowType == 'Equil').and.CoMDiff) then
      rAveMSD=0._wp;sdrAveMSD=0._wp
      rAvetAveDcm=0._wp;rAvesdtAveDcm=0._wp;sdrAvetAveDcm=0._wp        
      if (allocated(rAveMSD))   deallocate(rAveMSD)
      if (allocated(sdrAveMSD)) deallocate(sdrAveMSD)
      allocate(rAveMSD(ndmp),sdrAveMSD(ndmp))
      if (myrank == 0) then
        if (allocated(rAveMSDTot))   deallocate(rAveMSDTot)
        if (allocated(sdrAveMSDTot)) deallocate(sdrAveMSDTot)
        allocate(rAveMSDTot(ndmp),sdrAveMSDTot(ndmp))
      end if
    end if
    ! Radius of gyration.
    if (RgCalc) then
      rAveRgSq=0._wp;sdrAveRgSq=0._wp
      rAvetAveRgSq=0._wp;rAvesdtAveRgSq=0._wp;sdrAvetAveRgSq=0._wp
      rAvetAveAspher=0._wp;rAvesdtAveAspher=0._wp;sdrAvetAveAspher=0._wp
      if (allocated(rAveRgSq))   deallocate(rAveRgSq)
      if (allocated(sdrAveRgSq)) deallocate(sdrAveRgSq)
      allocate(rAveRgSq(ndmp),sdrAveRgSq(ndmp))
      if (myrank == 0) then
        if (allocated(rAveRgSqTot))   deallocate(rAveRgSqTot)
        if (allocated(sdrAveRgSqTot)) deallocate(sdrAveRgSqTot)
        allocate(rAveRgSqTot(ndmp),sdrAveRgSqTot(ndmp))
      end if ! myrank
    end if
    if (CorFun) then
      select case (CF_mode)
        case ('Ree')
          if (allocated(cf_ree)) deallocate(cf_ree)
          allocate(cf_ree(3,ntime/tgap,ntotchain,nprun))
        case ('Rg')
          if (allocated(cf_rg)) deallocate(cf_rg)
          allocate(cf_rg(ntotbeadx3,ntime/tgap,nprun))
        end select
    end if

  end subroutine data_run_init



  !!> material_func 
  subroutine material_func(myrank,irun,itime,time,Wi,Pe,dt,tgap,ntime,nchain,nseg,nbead,nsegx3,nbeadx3,nchain_cmb,nseg_cmb,nseg_cmbbb,&
               add_cmb,ntotchain,lambda,tss,trst,MPI_REAL_WP,nrun,nprun,tend,bs       ,chains        ,R)
 !      call material_func(id,irun,itime,time,Wi,Pe,dt,tgap,ntime,nchain,nseg,nbead,nsegx3,nbeadx3,nchain_cmb,nseg_cmb,&
 !             add_cmb,ntotchain,lambda,tss,trst,MPI_REAL_WP,nrun,nprun,tend,this%size,this%BoxChains,this%R_tilde)
    use :: hi_mod, only: DecompMeth,mst
    use :: chain_mod, only: chain
    use :: arry_mod, only: print_vector,print_matrix
    use :: flow_mod, only: FlowType
    use :: sprforce_mod, only: qmx
    use :: io_mod, only: rcmst
    use :: mpi
    !include 'mpif.h'

    integer,intent(in) :: myrank,irun,itime,ntime,tgap,nchain,nseg,nbead,nsegx3,nbeadx3
	integer,intent(in) :: nchain_cmb,nseg_cmb,ntotchain,nseg_cmbbb
	integer            :: nbead_cmb,nbead_cmbx3
	logical,intent(in) :: add_cmb
    integer,intent(in) :: MPI_REAL_WP,nrun,nprun
    type(chain),intent(in) :: chains(:) ! this%BoxChains
    real(wp),intent(in) :: Wi,Pe,dt,time,R(:),lambda,tss,trst,tend,bs(3)
    integer :: jchain,ierr,i,j,info,icount
    real(wp) :: rcmtmp(3),rcmExcess(3),rcmTot(3),RgSqTensEV(3),RgSqEVdiff(3)
    real(wp) :: RgSqTens(3,3),count_time
    real(wp),dimension(:),pointer  :: qP,RPi,RPj
    character(len=99),parameter :: fmtfef="(f8.2,1x,e11.3,1x,f14.7)"
    character(len=99),parameter :: fmtfeff="(f8.2,1x,e11.3,1x,f14.7,2x,f18.10)"
    character(len=99),parameter :: fmtfefff="(f8.2,1x,e11.3,1x,f14.7,2x,f18.10,2x,f14.7)"
    character(len=99),parameter :: fmtfe2f="(f8.2,1x,e11.3,1x,2(f20.7,2x))"
    character(len=99),parameter :: fmtfe3f="(f8.2,1x,e11.3,1x,3(f14.7,2x))"
    character(len=99),parameter :: fmtie3f="(i4,1x,e11.3,1x,3(f14.7,2x))"
    character(len=99),parameter :: fmtife2f="(i4,1x,f14.7,1x,e11.3,1x,2(f14.7,2x))"
    character(len=99),parameter :: fmtife3f="(i4,1x,f14.7,1x,e11.3,1x,3(f14.7,2x))"
    character(len=99),parameter :: fmtae3f="(a,1x,e11.3,1x,3(f14.7,2x))"
    character(len=99),parameter :: fmtie2f="(i4,1x,e11.3,1x,2(f14.7,2x))"
	write(*,*) "module:pp_smdlt:material_func"
	if (add_cmb) then
    nbead_cmb=nseg_cmb+1
	nbead_cmbx3=nbead_cmb*3
	end if
	
    lcount=lcount+1
    if (StrCalc) then
      call StressCalc(itime,time,ntime,irun,myrank,nchain,nseg,nsegx3,nbeadx3,nchain_cmb,nseg_cmb,nseg_cmbbb,ntotchain,add_cmb,tss,lambda,&
                        nrun,tend,chains,R)
    end if ! StrCalc

    ! Center of mass diffusion calculation:
    if ((FlowType == 'Equil').and.CoMDiff) then
      DcmAve=0._wp;sdDcmAve=0._wp
      do jchain=1, nchain
        rcmtmp=chains(jchain)%chain_rcm(:)!rcm(:,jchain)
        rcmExcess(:)=bs(:)*chains(jchain)%chain_cmif(:)
        rcmtot=rcmtmp+rcmExcess-rcmst(jchain,:,irun)
        DcmAve=DcmAve+dot_product(rcmtot,rcmtot)             
        sdDcmAve=sdDcmAve+dot_product(rcmtot,rcmtot)*dot_product(rcmtot,rcmtot)
      end do
	  if (add_cmb) then
        do jchain=1, nchain_cmb
            rcmtmp=chains(jchain+nchain)%chain_rcm(:)!rcm(:,jchain)
            rcmExcess(:)=bs(:)*chains(jchain+nchain)%chain_cmif(:)
            rcmtot=rcmtmp+rcmExcess-rcmst(jchain+nchain,:,irun)
            DcmAve=DcmAve+dot_product(rcmtot,rcmtot)             
            sdDcmAve=sdDcmAve+dot_product(rcmtot,rcmtot)*dot_product(rcmtot,rcmtot)
        end do
	  end if !add_cmb
    end if ! CoMDiff
	
    ! Radius of gyration:
    if (RgCalc) then
      ! R=rv-rcm
      RgSqAve=0._wp;sdRgSqAve=0._wp;AspherAve=0._wp;sdAspherAve=0._wp
      do jchain=1, nchain
        do j=1, 3
          RPj => chains(jchain)%chain_R(j:nbeadx3-(3-j):3)
          do i=1, 3
            RPi => chains(jchain)%chain_R(i:nbeadx3-(3-i):3)
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
        Aspher=1._wp/6/RgSqTensEVbar**2*dot_product(RgSqEVdiff,RgSqEVdiff)
        RgSqAve=RgSqAve+traceRgSqTens
        sdRgSqAve=sdRgSqAve+traceRgSqTens*traceRgSqTens
        AspherAve=AspherAve+Aspher
        sdAspherAve=sdAspherAve+Aspher*Aspher
      end do
	  
	  if (add_cmb) then
	  do jchain=1, nchain_cmb
        do j=1, 3
          RPj => chains(jchain+nchain)%chain_R(j:nbead_cmbx3-(3-j):3)
          do i=1, 3
            RPi => chains(jchain+nchain)%chain_R(i:nbead_cmbx3-(3-i):3)
            RgSqTens(i,j)=1._wp/nbead_cmb*dot(RPi,RPj)
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
        Aspher=1._wp/6/RgSqTensEVbar**2*dot_product(RgSqEVdiff,RgSqEVdiff)
        RgSqAve=RgSqAve+traceRgSqTens
        sdRgSqAve=sdRgSqAve+traceRgSqTens*traceRgSqTens
        AspherAve=AspherAve+Aspher
        sdAspherAve=sdAspherAve+Aspher*Aspher
      end do !
	  end if !add_cmb
	  
    end if ! RgCalc

    ! Evaluations as a function of time:
    if (StrCalc) then
      sdsqqetoeAve=sqrt(abs(sdsqqetoeAve-sqqetoeAve*sqqetoeAve)/(ntotchain-1))
      sdqetoeAve=sqrt(abs(sdqetoeAve-qetoeAve*qetoeAve)/(ntotchain-1))
      ! For run average:
      if (nrun > 1) then
        rAvesqqetoe(lcount)=rAvesqqetoe(lcount)+sqqetoeAve
        rAveqetoe(lcount)=rAveqetoe(lcount)+qetoeAve
        rAvetauxy(lcount)=rAvetauxy(lcount)+tauxy
        rAvetauxxyy(lcount)=rAvetauxxyy(lcount)+tauxxyy
        rAvetauyyzz(lcount)=rAvetauyyzz(lcount)+tauyyzz
        sdrAvesqqetoe(lcount)=sdrAvesqqetoe(lcount)+sqqetoeAve*sqqetoeAve
        sdrAveqetoe(lcount)=sdrAveqetoe(lcount)+qetoeAve*qetoeAve
        sdrAvetauxy(lcount)=sdrAvetauxy(lcount)+tauxy*tauxy
        sdrAvetauxxyy(lcount)=sdrAvetauxxyy(lcount)+tauxxyy*tauxxyy
        sdrAvetauyyzz(lcount)=sdrAvetauyyzz(lcount)+tauyyzz*tauyyzz
        ! Comment out if you want in detail information of each process.
        write(u1,fmtfefff) Wi,dt,time,sqqetoeAve,sdsqqetoeAve
        write(u2,fmtfefff) Wi,dt,Pe*time,qetoeAve/(qmx*nseg),sdqetoeAve/(qmx*nseg)
      else ! nrun=1
        write(u1,fmtfefff) Wi,dt,time,sqqetoeAve,sdsqqetoeAve
        write(u2,fmtfefff) Wi,dt,Pe*time,qetoeAve/(qmx*nseg),sdqetoeAve/(qmx*nseg)
      end if
      ! To investigate the transient behavior of Ree in rank 0:
!      if ((nrun > 1) .and. (myrank == 0)) write(u1,fmtfefff) Wi,dt,time,sqqetoeAve,sdsqqetoeAve
      select case (StrPr_mode)
        case ('Visc')
          if (FlowType /= 'Equil') then
            if (FlowType == 'PSF') then
              if (nrun > 1) then
!               Comment out if you want in detail information of each process.
!                write(u3,fmtfeff) Wi,dt,time,-tauxy/Pe
!                write(u25,fmtfeff) Wi,dt,time,-tauxxyy/Pe**2
!                write(u26,fmtfeff) Wi,dt,time,-tauyyzz/Pe**2
              else ! nrun=1
                write(u3,fmtfeff) Wi,dt,time,-tauxy/Pe
                write(u25,fmtfeff) Wi,dt,time,-tauxxyy/Pe**2
                write(u26,fmtfeff) Wi,dt,time,-tauyyzz/Pe**2
              end if ! nrun
            !else
			end if
			if (FlowType == 'PEF') then
!             Not implemented yet!
              write(u9,fmtfef) Pe*time,dt,-tauxxyy/Pe
            end if
          end if
        case ('Stress')
          if (nrun > 1) then
!           Comment out if you want in detail information of each process.
!            write(u3,fmtfeff) Wi,dt,time,tauxy
!            write(u25,fmtfeff) Wi,dt,time,tauxxyy
!            write(u26,fmtfeff) Wi,dt,time,tauyyzz
          else ! nrun=1
            write(u3,fmtfeff) Wi,dt,time,tauxy
            write(u25,fmtfeff) Wi,dt,time,tauxxyy
            write(u26,fmtfeff) Wi,dt,time,tauyyzz
          end if ! nrun
      end select
    end if ! StrCalc
    if (AveIterCalc) then
      if (DecompMeth == 'Lanczos') then
        if (nrun > 1) then
!         Comment out if you want in detail information of each process.
!          write(u15,fmtae3f) FlowType,Wi,dt,time,real(mst)
        else ! nrun=1
          write(u15,fmtae3f) FlowType,Wi,dt,time,real(mst)
        end if
        kcount=kcount+1
        tAvem=tAvem+mst
        sdtAvem=sdtAvem+mst*mst
      end if
    end if ! AveIterCalc
    if ((FlowType == 'Equil').and.CoMDiff) then
      MSDAve=DcmAve/ntotchain
      sdMSDAve=sdDcmAve/ntotchain
      sdMSDAve=sqrt(abs(sdMSDAve-MSDAve**2)/(ntotchain-1))
      DcmAve=DcmAve/(6*time*ntotchain)
      sdDcmAve=sdDcmAve/(36*time*time*ntotchain)
      sdDcmAve=sqrt(abs(sdDcmAve-DcmAve**2)/(ntotchain-1))
      ! For run averages:
      if (nrun > 1) then
        rAveMSD(lcount)=rAveMSD(lcount)+MSDAve
        sdrAveMSD(lcount)=sdrAveMSD(lcount)+MSDAve*MSDAve
!       Comment out if you want in detail information of each process.
!        write(u37,fmtfefff) Wi,dt,time,MSDAve,sdMSDAve
!        write(u19,fmtfefff) Wi,dt,time,DcmAve,sdDcmAve
      else ! nrun=1
        write(u37,fmtfefff) Wi,dt,time,MSDAve,sdMSDAve
        write(u19,fmtfefff) Wi,dt,time,DcmAve,sdDcmAve
      end if ! nrun
    end if
	
    if (RgCalc) then
      RgSqAve=RgSqAve/ntotchain
      AspherAve=AspherAve/ntotchain
      sdRgSqAve=sdRgSqAve/ntotchain
      sdAspherAve=sdAspherAve/ntotchain
      sdRgSqAve=sqrt(abs(sdRgSqAve-RgSqAve**2)/(ntotchain-1))
      sdAspherAve=sqrt(abs(sdAspherAve-AspherAve**2)/(ntotchain-1))
!     For run average:
      if (nrun > 1) then
        rAveRgSq(lcount)=rAveRgSq(lcount)+RgSqAve
        sdrAveRgSq(lcount)=sdrAveRgSq(lcount)+RgSqAve*RgSqAve
!       Comment out if you want in detail information of each process.
!        write(u18,fmtfefff) Wi,dt,time,RgSqAve,sdRgSqAve
      else ! nrun=1
        write(u18,fmtfefff) Wi,dt,time,RgSqAve,sdRgSqAve  !Rg
      end if
    end if

    ! After ?*'chain-largest (or) characteristic-relaxation-time' average over time:
    if (time >= tss*lambda) then
      jcount=jcount+1
      if (StrCalc) then
        ! Terms for stress and relative extension:
        tAvesqqetoeAve=tAvesqqetoeAve+sqqetoeAve
        tAveqetoeAve=tAveqetoeAve+qetoeAve
        tAvetauxy=tAvetauxy+tauxy
        tAvetauxxyy=tAvetauxxyy+tauxxyy
        tAvetauyyzz=tAvetauyyzz+tauyyzz
        ! Terms for standard deviation of stress and relative extension
        tAvesdsqqetoeAve=tAvesdsqqetoeAve+sdsqqetoeAve
        tAvesdqetoeAve=tAvesdqetoeAve+sdqetoeAve
        ! Standard deviation of terms for stress and relative extension
        sdtAvesqqetoeAve=sdtAvesqqetoeAve+sqqetoeAve*sqqetoeAve
        sdtAveqetoeAve=sdtAveqetoeAve+qetoeAve*qetoeAve
        sdtAvetauxy=sdtAvetauxy+tauxy*tauxy
        sdtAvetauxxyy=sdtAvetauxxyy+tauxxyy*tauxxyy
        sdtAvetauyyzz=sdtAvetauyyzz+tauyyzz*tauyyzz
      end if ! StrCalc
      if ((FlowType == 'Equil').and.CoMDiff) then
        tAveDcmAve=tAveDcmAve+DcmAve
        tAvesdDcmAve=tAvesdDcmAve+sdDcmAve
        sdtAveDcmAve=sdtAveDcmAve+DcmAve*DcmAve
      end if ! CoMDiff
      if (RgCalc) then
        tAveRgSqAve=tAveRgSqAve+RgSqAve
        tAveAspherAve=tAveAspherAve+AspherAve
        tAvesdRgSqAve=tAvesdRgSqAve+sdRgSqAve
        tAvesdAspherAve=tAvesdAspherAve+sdAspherAve
        sdtAveRgSqAve=sdtAveRgSqAve+RgSqAve*RgSqAve
        sdtAveAspherAve=sdtAveAspherAve+AspherAve*AspherAve
      end if ! RgCalc
    end if ! time >= tss*lambda
    ! Finalize the time averages in the last iteration: 
    if (itime == ntime) then
      if (jcount /= 0) then
        if (StrCalc) then
          ! Terms for stress and relative extension:
          tAvesqqetoeAve=tAvesqqetoeAve/jcount
          tAveqetoeAve=tAveqetoeAve/jcount
          tAvetauxy=tAvetauxy/jcount
          tAvetauxxyy=tAvetauxxyy/jcount
          tAvetauyyzz=tAvetauyyzz/jcount
          ! Terms for standard deviation of stress and relative extension:
          tAvesdsqqetoeAve=tAvesdsqqetoeAve/jcount
          tAvesdqetoeAve=tAvesdqetoeAve/jcount
          ! Standard deviation of terms for stress and relative extension:
          sdtAvesqqetoeAve=sdtAvesqqetoeAve/jcount
          sdtAveqetoeAve=sdtAveqetoeAve/jcount
          sdtAvetauxy=sdtAvetauxy/jcount
          sdtAvetauxxyy=sdtAvetauxxyy/jcount
          sdtAvetauyyzz=sdtAvetauyyzz/jcount
          sdtAvesqqetoeAve=sqrt(abs(sdtAvesqqetoeAve-tAvesqqetoeAve**2)/(jcount-1))
          sdtAveqetoeAve=sqrt(abs(sdtAveqetoeAve-tAveqetoeAve**2)/(jcount-1))
          sdtAvetauxy=sqrt(abs(sdtAvetauxy-tAvetauxy**2)/(jcount-1))
          sdtAvetauxxyy=sqrt(abs(sdtAvetauxxyy-tAvetauxxyy**2)/(jcount-1))
          sdtAvetauyyzz=sqrt(abs(sdtAvetauyyzz-tAvetauyyzz**2)/(jcount-1))
          ! For run averaging:
          if (nrun > 1) then
            rAvetAvesqqetoe=rAvetAvesqqetoe+tAvesqqetoeAve
            rAvetAveqetoe=rAvetAveqetoe+tAveqetoeAve
            rAvetAvetauxy=rAvetAvetauxy+tAvetauxy
            rAvetAvetauxxyy=rAvetAvetauxxyy+tAvetauxxyy
            rAvetAvetauyyzz=rAvetAvetauyyzz+tAvetauyyzz
            rAvesdtAvesqqetoe=rAvesdtAvesqqetoe+sdtAvesqqetoeAve
            rAvesdtAveqetoe=rAvesdtAveqetoe+sdtAveqetoeAve
            rAvesdtAvetauxy=rAvesdtAvetauxy+sdtAvetauxy
            rAvesdtAvetauxxyy=rAvesdtAvetauxxyy+sdtAvetauxxyy
            rAvesdtAvetauyyzz=rAvesdtAvetauyyzz+sdtAvetauyyzz            
            sdrAvetAvesqqetoe=sdrAvetAvesqqetoe+tAvesqqetoeAve*tAvesqqetoeAve
            sdrAvetAveqetoe=sdrAvetAveqetoe+tAveqetoeAve*tAveqetoeAve
            sdrAvetAvetauxy=sdrAvetAvetauxy+tAvetauxy*tAvetauxy
            sdrAvetAvetauxxyy=sdrAvetAvetauxxyy+tAvetauxxyy*tAvetauxxyy
            sdrAvetAvetauyyzz=sdrAvetAvetauyyzz+tAvetauyyzz*tAvetauyyzz
          end if
          write(u4,fmtfe3f) Wi,dt,tAveqetoeAve/(qmx*nseg),tAvesdqetoeAve/(qmx*nseg),sdtAveqetoeAve/(qmx*nseg)
          write(u16,fmtife3f) nbead,Wi,dt,tAvesqqetoeAve,tAvesdsqqetoeAve,sdtAvesqqetoeAve
          select case (StrPr_mode)
            case ('Visc')
              if (FlowType /= 'Equil') then
                write(u5,fmtfe2f) Wi,dt,-tAvetauxy/Pe,sdtAvetauxy/Pe
                write(u6,fmtfe2f) Wi,dt,-tAvetauxxyy/Pe**2,sdtAvetauxxyy/Pe**2
                write(u7,fmtfe2f) Wi,dt,-tAvetauyyzz/Pe**2,sdtAvetauyyzz/Pe**2
              end if ! FlowType /= 'Equil'
            case ('Stress')
              write(u5,fmtfe2f) Wi,dt,tAvetauxy,sdtAvetauxy
              write(u6,fmtfe2f) Wi,dt,tAvetauxxyy,sdtAvetauxxyy
              write(u7,fmtfe2f) Wi,dt,tAvetauyyzz,sdtAvetauyyzz
          end select
        end if ! StrCalc
        if ((FlowType == 'Equil').and.CoMDiff) then
          tAveDcmAve=tAveDcmAve/jcount
          tAvesdDcmAve=tAvesdDcmAve/jcount
          sdtAveDcmAve=sdtAveDcmAve/jcount
          sdtAveDcmAve=sqrt(abs(sdtAveDcmAve-tAveDcmAve**2)/(jcount-1))
          ! For run averaging:
          if (nrun > 1) then
            rAvetAveDcm=rAvetAveDcm+tAveDcmAve
            rAvesdtAveDcm=rAvesdtAveDcm+sdtAveDcmAve
            sdrAvetAveDcm=sdrAvetAveDcm+tAveDcmAve*tAveDcmAve
          end if
          write(u10,fmtie3f) nbead,dt,tAveDcmAve,tAvesdDcmAve,sdtAveDcmAve
        end if ! CoMDiff
        if (RgCalc) then
          tAveRgSqAve=tAveRgSqAve/jcount
          tAveAspherAve=tAveAspherAve/jcount
          tAvesdRgSqAve=tAvesdRgSqAve/jcount
          tAvesdAspherAve=tAvesdAspherAve/jcount
          sdtAveRgSqAve=sdtAveRgSqAve/jcount
          sdtAveAspherAve=sdtAveAspherAve/jcount
          sdtAveRgSqAve=sqrt(abs(sdtAveRgSqAve-tAveRgSqAve**2)/(jcount-1))
          sdtAveAspherAve=sqrt(abs(sdtAveAspherAve-tAveAspherAve**2)/(jcount-1))
          ! For run averaging:
          if (nrun > 1) then
            rAvetAveRgSq=rAvetAveRgSq+tAveRgSqAve
            rAvetAveAspher=rAvetAveAspher+tAveAspherAve
            rAvesdtAveRgSq=rAvesdtAveRgSq+sdtAveRgSqAve
            rAvesdtAveAspher=rAvesdtAveAspher+sdtAveAspherAve
            sdrAvetAveRgSq=sdrAvetAveRgSq+tAveRgSqAve*tAveRgSqAve
            sdrAvetAveAspher=sdrAvetAveAspher+tAveAspherAve*tAveAspherAve
          end if
          write(u12,fmtife3f) nbead,Wi,dt,tAveRgSqAve,tAvesdRgSqAve,sdtAveRgSqAve
          write(u13,fmtife3f) nbead,Wi,dt,tAveAspherAve,tAvesdAspherAve,sdtAveAspherAve
        end if ! RgCalc
      else ! jcount=0
        if (StrCalc) then
          ! For run averaging:
          if (nrun > 1) then
            rAvetAvesqqetoe=rAvetAvesqqetoe+sqqetoeAve
            rAvetAveqetoe=rAvetAveqetoe+qetoeAve
            rAvetAvetauxy=rAvetAvetauxy+tauxy
            rAvetAvetauxxyy=rAvetAvetauxxyy+tauxxyy
            rAvetAvetauyyzz=rAvetAvetauyyzz+tauyyzz
            rAvesdtAvesqqetoe=rAvesdtAvesqqetoe+sdsqqetoeAve
            rAvesdtAveqetoe=rAvesdtAveqetoe+sdqetoeAve
            sdrAvetAvesqqetoe=sdrAvetAvesqqetoe+sqqetoeAve*sqqetoeAve
            sdrAvetAveqetoe=sdrAvetAveqetoe+qetoeAve*qetoeAve
            sdrAvetAvetauxy=sdrAvetAvetauxy+tauxy*tauxy
            sdrAvetAvetauxxyy=sdrAvetAvetauxxyy+tauxxyy*tauxxyy
            sdrAvetAvetauyyzz=sdrAvetAvetauyyzz+tauyyzz*tauyyzz
          end if
          write(u4,fmtfe2f) Wi,dt,qetoeAve/(qmx*nseg),sdqetoeAve/(qmx*nseg)  !normalised based on linear chainsegment
          select case (StrPr_mode)
            case ('Visc')
              if (FlowType /= 'Equil') then
                write(u5,fmtfef) Wi,dt,-tauxy/Pe
                write(u6,fmtfef) Wi,dt,-tauxxyy/(Pe**2)
                write(u7,fmtfef) Wi,dt,-tauyyzz/(Pe**2)
              end if ! FlowType.ne.'Equil'
            case ('Stress')
              write(u5,fmtfef) Wi,dt,tauxy
              write(u6,fmtfef) Wi,dt,tauxxyy
              write(u7,fmtfef) Wi,dt,tauyyzz
          end select
        end if ! StrCalc
        if ((FlowType == 'Equil').and.CoMDiff) then
          ! For run averaging:
          if (nrun > 1) then
            rAvetAveDcm=rAvetAveDcm+DcmAve
            rAvesdtAveDcm=rAvesdtAveDcm+sdDcmAve
            sdrAvetAveDcm=sdrAvetAveDcm+DcmAve*DcmAve
          end if
          write(u10,fmtie2f) nbead,dt,DcmAve,sdDcmAve
        end if
        if (RgCalc) then
          ! For run averaging:
          if (nrun > 1) then
            rAvetAveRgSq=rAvetAveRgSq+RgSqAve
            rAvetAveAspher=rAvetAveAspher+AspherAve
            rAvesdtAveRgSq=rAvesdtAveRgSq+sdRgSqAve
            rAvesdtAveAspher=rAvesdtAveAspher+sdAspherAve
            sdrAvetAveRgSq=sdrAvetAveRgSq+RgSqAve*RgSqAve
            sdrAvetAveAspher=sdrAvetAveAspher+AspherAve*AspherAve
          end if
          write(u12,fmtife2f) nbead,Wi,dt,RgSqAve,sdRgSqAve
          write(u13,fmtife2f) nbead,Wi,dt,AspherAve,sdAspherAve
        end if ! RgCalc
      end if ! jcount /= 0
      if (AveIterCalc) then
        if (kcount /= 0) then
          if (DecompMeth == 'Lanczos') then
            tAvem=tAvem/kcount
            sdtAvem=sdtAvem/kcount
            sdtAvem=sqrt(abs(sdtAvem-tAvem*tAvem)/(kcount-1))
            write(u14,fmtife2f) nbead,Wi,dt,tAvem,sdtAvem
          end if
        end if ! kcount /= 0
      end if ! AveIterCalc
    end if ! itime == ntime
    ! Finalize the run averages in the last iteration:
    if ((nrun > 1) .and. (irun == nprun) .and. (itime == ntime)) then

      ! Terms for stress and relative extension:
      if (StrCalc) then

        do icount=1, lcount

          rAvesqqetoe(icount)=rAvesqqetoe(icount)/nrun
          rAveqetoe(icount)=rAveqetoe(icount)/nrun
          rAvetauxy(icount)=rAvetauxy(icount)/nrun
          rAvetauxxyy(icount)=rAvetauxxyy(icount)/nrun
          rAvetauyyzz(icount)=rAvetauyyzz(icount)/nrun
          sdrAvesqqetoe(icount)=sdrAvesqqetoe(icount)/nrun
          sdrAveqetoe(icount)=sdrAveqetoe(icount)/nrun
          sdrAvetauxy(icount)=sdrAvetauxy(icount)/nrun
          sdrAvetauxxyy(icount)=sdrAvetauxxyy(icount)/nrun
          sdrAvetauyyzz(icount)=sdrAvetauyyzz(icount)/nrun
       
        end do ! icount
!
        rAvetAvesqqetoe=rAvetAvesqqetoe/nrun
        rAvetAveqetoe=rAvetAveqetoe/nrun
        rAvetAvetauxy=rAvetAvetauxy/nrun
        rAvetAvetauxxyy=rAvetAvetauxxyy/nrun
        rAvetAvetauyyzz=rAvetAvetauyyzz/nrun
        rAvesdtAvesqqetoe=rAvesdtAvesqqetoe/nrun
        rAvesdtAveqetoe=rAvesdtAveqetoe/nrun
        rAvesdtAvetauxy=rAvesdtAvetauxy/nrun
        rAvesdtAvetauxxyy=rAvesdtAvetauxxyy/nrun
        rAvesdtAvetauyyzz=rAvesdtAvetauyyzz/nrun
        sdrAvetAvesqqetoe=sdrAvetAvesqqetoe/nrun
        sdrAvetAveqetoe=sdrAvetAveqetoe/nrun
        sdrAvetAvetauxy=sdrAvetAvetauxy/nrun
        sdrAvetAvetauxxyy=sdrAvetAvetauxxyy/nrun
        sdrAvetAvetauyyzz=sdrAvetAvetauyyzz/nrun

        ! Reduction of terms for transient stress and relative extension:
        call MPI_Reduce(rAvesqqetoe,rAvesqqetoeTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAveqetoe,rAveqetoeTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetauxy,rAvetauxyTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetauxxyy,rAvetauxxyyTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetauyyzz,rAvetauyyzzTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvesqqetoe,sdrAvesqqetoeTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAveqetoe,sdrAveqetoeTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetauxy,sdrAvetauxyTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetauxxyy,sdrAvetauxxyyTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetauyyzz,sdrAvetauyyzzTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(rAvetAvesqqetoe,rAvetAvesqqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetAveqetoe,rAvetAveqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetAvetauxy,rAvetAvetauxyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetAvetauxxyy,rAvetAvetauxxyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetAvetauyyzz,rAvetAvetauyyzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(rAvesdtAvesqqetoe,rAvesdtAvesqqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAveqetoe,rAvesdtAveqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAvetauxy,rAvesdtAvetauxyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAvetauxxyy,rAvesdtAvetauxxyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAvetauyyzz,rAvesdtAvetauyyzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(sdrAvetAvesqqetoe,sdrAvetAvesqqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAveqetoe,sdrAvetAveqetoeTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAvetauxy,sdrAvetAvetauxyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAvetauxxyy,sdrAvetAvetauxxyyTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAvetauyyzz,sdrAvetAvetauyyzzTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (myrank == 0) then

          do icount=1, lcount 
            if (icount /= lcount) then
              count_time=tgap*dt*icount
            else
              count_time=ntime*dt
            end if
            sdrAvesqqetoeTot(icount)=sqrt(abs(sdrAvesqqetoeTot(icount)-rAvesqqetoeTot(icount)**2)/(nrun-1))
            sdrAveqetoeTot(icount)=sqrt(abs(sdrAveqetoeTot(icount)-rAveqetoeTot(icount)**2)/(nrun-1))
            sdrAvetauxyTot(icount)=sqrt(abs(sdrAvetauxyTot(icount)-rAvetauxyTot(icount)**2)/(nrun-1))
            sdrAvetauxxyyTot(icount)=sqrt(abs(sdrAvetauxxyyTot(icount)-rAvetauxxyyTot(icount)**2)/(nrun-1))
            sdrAvetauyyzzTot(icount)=sqrt(abs(sdrAvetauyyzzTot(icount)-rAvetauyyzzTot(icount)**2)/(nrun-1)) 
            write(u23,fmtfe3f) Wi,dt,count_time,rAveqetoeTot(icount)/(qmx*nseg),sdrAveqetoeTot(icount)/(qmx*nseg)
            write(u20,fmtife3f) nbead,Wi,dt,count_time,rAvesqqetoeTot(icount),sdrAvesqqetoeTot(icount)
            select case (StrPr_mode)
              case ('Visc')
                if (FlowType /= 'Equil') then
                  write(u24,fmtfe3f) Wi,dt,count_time,-rAvetauxyTot(icount)/Pe,sdrAvetauxyTot(icount)/Pe
                  write(u27,fmtfe3f) Wi,dt,count_time,-rAvetauxxyyTot(icount)/Pe**2,sdrAvetauxxyyTot(icount)/Pe**2
                  write(u28,fmtfe3f) Wi,dt,count_time,-rAvetauyyzzTot(icount)/Pe**2,sdrAvetauyyzzTot(icount)/Pe**2
                end if ! FlowType /= 'Equil'
              case ('Stress')
                write(u24,fmtfe3f) Wi,dt,count_time,rAvetauxyTot(icount),sdrAvetauxyTot(icount)
                write(u27,fmtfe3f) Wi,dt,count_time,rAvetauxxyyTot(icount),sdrAvetauxxyyTot(icount)
                write(u28,fmtfe3f) Wi,dt,count_time,rAvetauyyzzTot(icount),sdrAvetauyyzzTot(icount)
            end select

          end do ! icount

          sdrAvetAvesqqetoeTot=sqrt(abs(sdrAvetAvesqqetoeTot-rAvetAvesqqetoeTot**2)/(nrun-1))
          sdrAvetAveqetoeTot=sqrt(abs(sdrAvetAveqetoeTot-rAvetAveqetoeTot**2)/(nrun-1))
          sdrAvetAvetauxyTot=sqrt(abs(sdrAvetAvetauxyTot-rAvetAvetauxyTot**2)/(nrun-1))
          sdrAvetAvetauxxyyTot=sqrt(abs(sdrAvetAvetauxxyyTot-rAvetAvetauxxyyTot**2)/(nrun-1))
          sdrAvetAvetauyyzzTot=sqrt(abs(sdrAvetAvetauyyzzTot-rAvetAvetauyyzzTot**2)/(nrun-1))
          write(u29,fmtfe3f) Wi,dt,rAvetAveqetoeTot/(qmx*nseg),rAvesdtAveqetoeTot/(qmx*nseg),&
                            sdrAvetAveqetoeTot/(qmx*nseg)
          write(u36,fmtife3f) nbead,Wi,dt,rAvetAvesqqetoeTot,rAvesdtAvesqqetoeTot,sdrAvetAvesqqetoeTot
          select case (StrPr_mode)
            case ('Visc')
              if (FlowType /= 'Equil') then
                write(u30,fmtfe3f) Wi,dt,-rAvetAvetauxyTot/Pe,rAvesdtAvetauxyTot/Pe,sdrAvetAvetauxyTot/Pe
                write(u31,fmtfe3f) Wi,dt,-rAvetAvetauxxyyTot/Pe**2,rAvesdtAvetauxxyyTot/Pe**2,&
                                   sdrAvetAvetauxxyyTot/Pe**2
                write(u32,fmtfe3f) Wi,dt,-rAvetAvetauyyzzTot/Pe**2,rAvesdtAvetauyyzzTot/Pe**2,&
                                   sdrAvetAvetauyyzzTot/Pe**2
              end if ! FlowType /= 'Equil'
            case ('Stress')
              write(u30,fmtfe3f) Wi,dt,rAvetAvetauxyTot,rAvesdtAvetauxyTot,sdrAvetAvetauxyTot
              write(u31,fmtfe3f) Wi,dt,rAvetAvetauxxyyTot,rAvesdtAvetauxxyyTot,sdrAvetAvetauxxyyTot
              write(u32,fmtfe3f) Wi,dt,rAvetAvetauyyzzTot,rAvesdtAvetauyyzzTot,sdrAvetAvetauyyzzTot
          end select

        end if ! myrank

      end if ! StrCalc

      if (RgCalc) then

        do icount=1, lcount
          rAveRgSq(icount)=rAveRgSq(icount)/nrun
          sdrAveRgSq(icount)=sdrAveRgSq(icount)/nrun
        end do

        rAvetAveRgSq=rAvetAveRgSq/nrun
        rAvetAveAspher=rAvetAveAspher/nrun
        rAvesdtAveRgSq=rAvesdtAveRgSq/nrun
        rAvesdtAveAspher=rAvesdtAveAspher/nrun
        sdrAvetAveRgSq=sdrAvetAveRgSq/nrun
        sdrAvetAveAspher=sdrAvetAveAspher/nrun


        call MPI_Reduce(rAveRgSq,rAveRgSqTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAveRgSq,sdrAveRgSqTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(rAvetAveRgSq,rAvetAveRgSqTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvetAveAspher,rAvetAveAspherTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(rAvesdtAveRgSq,rAvesdtAveRgSqTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAveAspher,rAvesdtAveAspherTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(sdrAvetAveRgSq,sdrAvetAveRgSqTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAveAspher,sdrAvetAveAspherTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (myrank == 0) then

          do icount =1, lcount
            if (icount /= lcount) then
              count_time=tgap*dt*icount
            else
              count_time=ntime*dt
            end if
            sdrAveRgSqTot(icount)=sqrt(abs(sdrAveRgSqTot(icount)-rAveRgSqTot(icount)**2)/(nrun-1))
            write(u21,fmtife3f) nbead,Wi,dt,count_time,rAveRgSqTot(icount),sdrAveRgSqTot(icount)
          end do ! icount

          sdrAvetAveRgSqTot=sqrt(abs(sdrAvetAveRgSqTot-rAvetAveRgSqTot**2)/(nrun-1))
          sdrAvetAveAspherTot=sqrt(abs(sdrAvetAveAspherTot-rAvetAveAspherTot**2)/(nrun-1))
          write(u34,fmtife3f) nbead,Wi,dt,rAvetAveRgSqTot,rAvesdtAveRgSqTot,sdrAvetAveRgSqTot
          write(u35,fmtife3f) nbead,Wi,dt,rAvetAveAspherTot,rAvesdtAveAspherTot,sdrAvetAveAspherTot

        end if ! myrank

      end if ! RgCalc

      if ((FlowType == 'Equil').and.CoMDiff) then

        do icount=1, lcount
          rAveMSD(icount)=rAveMSD(icount)/nrun
          sdrAveMSD(icount)=sdrAveMSD(icount)/nrun
        end do

        rAvetAveDcm=rAvetAveDcm/nrun
        rAvesdtAveDcm=rAvesdtAveDcm/nrun
        sdrAvetAveDcm=sdrAvetAveDcm/nrun

        call MPI_Reduce(rAveMSD,rAveMSDTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAveMSD,sdrAveMSDTot,lcount,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_Reduce(rAvetAveDcm,rAvetAveDcmTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(rAvesdtAveDcm,rAvesdtAveDcmTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_Reduce(sdrAvetAveDcm,sdrAvetAveDcmTot,1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (myrank == 0) then

          do icount=1, lcount

            if (icount /= lcount) then
              count_time=tgap*dt*icount
            else
              count_time=ntime*dt
            end if
            sdrAveMSDTot(icount)=sqrt(abs(sdrAveMSDTot(icount)-rAveMSDTot(icount)**2)/(nrun-1))
            write(u38,fmtife3f) nbead,Wi,dt,count_time,rAveMSDTot(icount),sdrAveMSDTot(icount)

          end do

          sdrAvetAveDcmTot=sqrt(abs(sdrAvetAveDcmTot-rAvetAveDcmTot**2)/(nrun-1))
          write(u33,fmtife3f) nbead,Wi,dt,rAvetAveDcmTot,rAvesdtAveDcmTot,sdrAvetAveDcmTot
        end if ! myrank
      end if ! CoMDiff


    end if ! irun > 1 & nrun == nprun & itime == ntime

  end subroutine material_func


  !!> StressCalc
  subroutine StressCalc(itime,time,ntime,irun,myrank,nchain,nseg,nsegx3,nbeadx3,nchain_cmb,nseg_cmb,nseg_cmbbb,ntotchain,add_cmb,tss,lambda,&
                        nrun,tend,chains,R)

    use :: chain_mod, only: chain
    use :: arry_mod, only: print_matrix,print_vector
    use :: force_smdlt, only: Fphi,rFphi
    use :: conv_mod, only: QtoR

    integer,intent(in) :: itime,ntime,irun,myrank,nchain,nseg,nsegx3,nbeadx3,nrun
    type(chain),intent(in) :: chains(:)
    real(wp),intent(in) :: time,tss,lambda,tend,R(:)
    real(wp),dimension(:),pointer :: qP,qPx,qPy,qPz,RPx,RPy,RPz
    real(wp),dimension(:),pointer :: FPx,FPy,FPz
    real(wp),dimension(:),pointer :: FphiPx => null(),FphiPy => null(),FphiPz => null()
    real(wp),dimension(:),pointer :: FsegPx => null(),FsegPy => null(),FsegPz => null()
    real(wp) :: txx,txy,tyy,tzz,txxyy,tyyzz,qetoex,qetoey,qetoez,sqqetoe,qetoe
    real(wp) :: txz,tyx,tyz,tzx,tzy,tauxz,tauyx,tauyz,tauzx,tauzy,rFev(4)
    integer :: ichain,offsetb,offsets
    integer,intent(in) :: nchain_cmb,nseg_cmb,ntotchain,nseg_cmbbb  !,add_cmb
	integer            :: nbead_cmb,nbead_cmbx3
	logical,intent(in) :: add_cmb
	write(*,*) "module:pp_smdlt:StressCalc"
	if (add_cmb) then
    nbead_cmb=nseg_cmb+1
	nbead_cmbx3=nbead_cmb*3
	end if
	
!   !-------------------------------------------------------------------------!
!   ! Stress is calculated for each chain by summing up the contribution      !
!   !   of all segments. Then all chain stresses in each process will be      !
!   !   added up.                                                             !
!   ! Note that standard deviation is calculated based on the formula for     !
!   !   sample standard deviation of the sampled mean:                        !
!   !   sigma_mean = 1/sqrt(N-1) * sigma ;                                    !
!   !   where sigma is the standard deviation of the distribution:            !
!   !   sigma = sqrt(ave(X^2)-(ave(x))^2)                                     !
!   ! Therefore, we will calculate tauxx and sdxx in each process and finally !
!   !   they will be summed up to give total stress and also appropriate term !
!   !   for getting the deviation.                                            !
!   !-------------------------------------------------------------------------!
    sqqetoeAve=0._wp;qetoeAve=0._wp;sdsqqetoeAve=0._wp;sdqetoeAve=0._wp
    tauxx=0._wp;tauxy=0._wp;tauyy=0._wp;tauzz=0._wp;tauxxyy=0._wp;tauyyzz=0._wp
    ! The data at ntime is excluded since a consistent interval isn't guaranteed.
    if ((CorFun) .and. (time < tend*lambda) .and. (itime /= ntime)) then
      if (CF_mode == 'Rg') then
        tauxz=0._wp;tauyx=0._wp;tauyz=0._wp;tauzx=0._wp;tauzy=0._wp
      end if
      kchk=kchk+1
    end if
	!linear chians firts
    do ichain=1, nchain
      offsetb=(chains(ichain)%chain_ID-1)*nbeadx3
      offsets=(chains(ichain)%chain_ID-1)*nsegx3
      qP => chains(ichain)%chain_Q(:)

      qPx => chains(ichain)%chain_Q(1:nsegx3-2:3)
      qPy => chains(ichain)%chain_Q(2:nsegx3-1:3) 
      qPz => chains(ichain)%chain_Q(3:nsegx3:3)

      qetoex=sum(qPx);qetoey=sum(qPy);qetoez=sum(qPz)
      sqqetoe=qetoex*qetoex+qetoey*qetoey+qetoez*qetoez
      qetoe=sqrt(sqqetoe)
      ! This data is stored for correlation function calculation.
      if ((CorFun) .and. (CF_mode == 'Ree') .and. (itime /= ntime)) then
        cf_ree(1:3,kchk,ichain,irun)=[qetoex,qetoey,qetoez]
      end if
      qetoeAve=qetoeAve+qetoe
      sqqetoeAve=sqqetoeAve+sqqetoe
      sdsqqetoeAve=sdsqqetoeAve+sqqetoe*sqqetoe
      sdqetoeAve=sdqetoeAve+qetoe*qetoe
    end do! chain loop
    if (add_cmb) then
      do ichain=1, nchain_cmb
        offsetb=(chains(nchain)%chain_ID-1)*nbeadx3+(chains(ichain+nchain)%chain_ID-chains(nchain)%chain_ID-1)*nbead_cmbx3
        offsets=(chains(ichain)%chain_ID-1)*nsegx3+(chains(ichain+nchain)%chain_ID-chains(ichain)%chain_ID-1)*(nseg_cmb*3)
        qP => chains(ichain+nchain)%chain_Q(:)
		
        qPx => chains(ichain+nchain)%chain_Q(1:(nseg_cmbbb*3)-2:3)
        qPy => chains(ichain+nchain)%chain_Q(2:(nseg_cmbbb*3)-1:3) 
        qPz => chains(ichain+nchain)%chain_Q(3:(nseg_cmbbb*3):3)

        qetoex=sum(qPx);qetoey=sum(qPy);qetoez=sum(qPz)
        sqqetoe=qetoex*qetoex+qetoey*qetoey+qetoez*qetoez
        qetoe=sqrt(sqqetoe)
      ! This data is stored for correlation function calculation.
        if ((CorFun) .and. (CF_mode == 'Ree') .and. (itime /= ntime)) then
          cf_ree(1:3,kchk,ichain+nchain,irun)=[qetoex,qetoey,qetoez]
        end if
        qetoeAve=qetoeAve+qetoe
        sqqetoeAve=sqqetoeAve+sqqetoe
        sdsqqetoeAve=sdsqqetoeAve+sqqetoe*sqqetoe
        sdqetoeAve=sdqetoeAve+qetoe*qetoe
      end do! chain_cmb loop
	end if !add_cmb
	
    sqqetoeAve=sqqetoeAve/ntotchain
    qetoeAve=qetoeAve/ntotchain
    sdsqqetoeAve=sdsqqetoeAve/ntotchain
    sdqetoeAve=sdqetoeAve/ntotchain

    tauxx = rFphi(1)/ntotchain
    tauxy = rFphi(2)/ntotchain 
    tauyy = rFphi(3)/ntotchain
    tauzz = rFphi(4)/ntotchain

    if ((CorFun) .and. (CF_mode == 'Rg') .and. (itime /= ntime)) then
      cf_rg(:,kchk,irun)=R
    end if

    tauxxyy = (rFphi(1)-rFphi(3))/ntotchain
    tauyyzz = (rFphi(3)-rFphi(4))/ntotchain
  
  end subroutine StressCalc


  !!> CorrFcn
  subroutine CorrFcn(dt,Wi,tgap,myrank,p,nchain,nchain_cmb,nbead,tend,lambda,nprun,MPI_REAL_WP)

    use :: arry_mod, only: print_vector,print_matrix
    use :: mpi
    !include 'mpif.h'

    integer,intent(in) :: tgap,myrank,p,nchain,nchain_cmb,nbead,nprun,MPI_REAL_WP
    real(wp),intent(in) :: dt,Wi,tend,lambda
    real(wp) :: CF0(3),term,term1,term2
	real(wp) :: CF0_cmb(3),term_cmb,term1_cmb,term2_cmb
    integer :: irun,ichain,t0,t0tmp,tt0tmp,ttmp,tcor,trun,tt0max,tt0,t,ierr,ttmpmx
    real(wp),allocatable :: CF(:),CF0t(:),CF00(:),CF0tTot(:),CF00Tot(:),sdCF(:)
    real(wp),allocatable :: CF_cmb(:),CF0t_cmb(:),CF00_cmb(:),CF0tTot_cmb(:),CF00Tot_cmb(:),sdCF_cmb(:)
    integer,allocatable :: norm(:), norm_cmb(:)
    real(wp),pointer :: rCFT(:),rCF0(:) ,rCFT_cmb(:),rCF0_cmb(:)
    character(len=99),parameter :: fmtfe2f="(f8.2,1x,e11.3,1x,2(f20.7,2x))"
    character(len=99),parameter :: fmtfe3f="(f8.2,1x,e11.3,1x,3(f14.7,2x))"
	write(*,*) "module:pp_smdlt:CorrFcn"
    select case (CF_mode)
      case ('Ree')
        write(u39,*) "# Wi, dt, time, <Qee(t).Qee(0)>.<Qee(0).Qee(0)>"
        write(u39,*) "#----------------------------------------------#"
      case ('Rg')
        write(u39,*) "# Wi, dt, time, <[r(t)-rcm(t)].[r(0)-rcm(0)]>/Rg^2 #"
        write(u39,*) "#--------------------------------------------------#"
    end select

    ! Total number of data points recorded for calculating CorFun.
    trun=kchk*tgap
    ! The number of points for representation of CorFun (Here it is practically ntime).
    tcor=ceiling(tend*lambda/dt)
    ! The maximum distance between two correlated values.
    ttmpmx=min(trun,tgap+tcor)-tgap
    allocate(CF(0:ttmpmx/tgap),sdCF(0:ttmpmx/tgap),norm(0:ttmpmx/tgap))
    allocate(CF0t(0:ttmpmx/tgap),CF00(0:ttmpmx/tgap))
    if (myrank == 0) allocate(CF0tTot(0:ttmpmx/tgap),CF00Tot(0:ttmpmx/tgap))
    CF=0._wp;sdCF=0._wp;norm=0
    CF0t=0._wp;CF00=0._wp
	
    if (add_cmb) then
	    allocate(CF_cmb(0:ttmpmx/tgap),sdCF_cmb(0:ttmpmx/tgap),norm_cmb(0:ttmpmx/tgap))
		allocate(CF0t_cmb(0:ttmpmx/tgap),CF00_cmb(0:ttmpmx/tgap))
		if (myrank == 0) allocate(CF0tTot_cmb(0:ttmpmx/tgap),CF00Tot_cmb(0:ttmpmx/tgap))
		CF_cmb=0._wp;sdCF_cmb=0._wp;norm_cmb=0
		CF0t_cmb=0._wp;CF00_cmb=0._wp
	end if
	
    select case (CF_mode)

      case ('Ree')

        do irun=1, nprun
          do ichain=1, nchain

            do t0tmp=1, kchk
              t0=t0tmp*tgap ! TAU0
              CF0(:)=cf_ree(:,t0tmp,ichain,irun) ! cf(TAU0)
              tt0max=min(trun,t0+tcor) ! ( TAU+TAU0 )max

              do tt0=t0, tt0max, tgap ! TAU+TAU0
                tt0tmp=tt0/tgap
                t=tt0-t0
                ttmp=t/tgap
                term1=dot_product(CF0(:),cf_ree(:,tt0tmp,ichain,irun))
                term2=dot_product(CF0(:),CF0(:))
                CF0t(ttmp)=CF0t(ttmp)+term1
                CF00(ttmp)=CF00(ttmp)+term2
                norm(ttmp)=norm(ttmp)+1
              end do ! tt0

            end do ! t0tmp

          end do ! ichain
		  
		  if (add_cmb) then
         do ichain=1, nchain_cmb

            do t0tmp=1, kchk
              t0=t0tmp*tgap ! TAU0
              CF0_cmb(:)=cf_ree(:,t0tmp,ichain+nchain,irun) ! cf(TAU0)
              tt0max=min(trun,t0+tcor) ! ( TAU+TAU0 )max

              do tt0=t0, tt0max, tgap ! TAU+TAU0
                tt0tmp=tt0/tgap
                t=tt0-t0
                ttmp=t/tgap
                term1_cmb=dot_product(CF0_cmb(:),cf_ree(:,tt0tmp,ichain+nchain,irun))
                term2_cmb=dot_product(CF0_cmb(:),CF0_cmb(:))
                CF0t_cmb(ttmp)=CF0t_cmb(ttmp)+term1_cmb
                CF00_cmb(ttmp)=CF00_cmb(ttmp)+term2_cmb
                norm_cmb(ttmp)=norm_cmb(ttmp)+1
              end do ! tt0

            end do ! t0tmp

          end do ! ichain_cmb
		  
		  end if !add_comb
        end do ! irun

      case ('Rg')

        do irun=1, nprun
          do t0tmp=1, kchk
            t0=t0tmp*tgap ! TAU0
            rCF0 => cf_rg(:,t0tmp,irun) ! cf(TAU0)
            tt0max=min(trun,t0+tcor) ! ( TAU+TAU0 )max
            do tt0=t0, tt0max, tgap ! TAU+TAU0
              tt0tmp=tt0/tgap
              t=tt0-t0
              ttmp=t/tgap
              rCFT => cf_rg(:,tt0tmp,irun) ! cf(TAU+TAU0)
              term1=dot(rCF0,rCFT)
              CF0t(ttmp)=CF0t(ttmp)+term1 ! Correlation function
              norm(ttmp)=norm(ttmp)+1
            end do ! tt0
          end do ! t0tmp
        end do ! irun

    end select

    call MPI_Reduce(CF0t,CF0tTot,ttmpmx/tgap+1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(CF00,CF00Tot,ttmpmx/tgap+1,MPI_REAL_WP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (myrank == 0) then

      do ttmp=0, ttmpmx/tgap ! Since ntime is not considered, kchk may be 1 unit smaller than tcor/tgap.
        select case (CF_mode)
          case ('Ree')
            CF0tTot(ttmp)=CF0tTot(ttmp)/(p*norm(ttmp))
            CF00Tot(ttmp)=CF00Tot(ttmp)/(p*norm(ttmp))
			if (add_cmb) then
            CF0tTot_cmb(ttmp)=CF0tTot_cmb(ttmp)/(p*norm_cmb(ttmp))
            CF00Tot_cmb(ttmp)=CF00Tot_cmb(ttmp)/(p*norm_cmb(ttmp))
			write(u39,fmtfe2f) Wi,dt,ttmp*tgap*dt,CF0tTot(ttmp)/CF00Tot(ttmp),CF0tTot_cmb(ttmp)/CF00Tot_cmb(ttmp)
			else
			write(u39,fmtfe2f) Wi,dt,ttmp*tgap*dt,CF0tTot(ttmp)/CF00Tot(ttmp)
			end if

          case ('Rg')
            CF0tTot(ttmp)=CF0tTot(ttmp)/(p*norm(ttmp))
			if (add_cmb) then
			CF0tTot_cmb(ttmp)=CF0tTot_cmb(ttmp)/(p*norm_cmb(ttmp))
			write(u39,fmtfe2f) Wi,dt,ttmp*tgap*dt,CF0tTot(ttmp)/CF0tTot(0),CF0tTot_cmb(ttmp)/CF0tTot_cmb(0)
			else
            write(u39,fmtfe2f) Wi,dt,ttmp*tgap*dt,CF0tTot(ttmp)/CF0tTot(0)
			end if
        end select
      end do

!      print *,'Rg:',CF0tTot(0)/(nbead*nchain)

    end if ! myrank
    ! reset kchk:
    kchk=0
    deallocate(CF,sdCF,norm)
    deallocate(CF0t,CF00)
	if (add_cmb) then
	deallocate(CF_cmb,sdCF_cmb,norm_cmb)
    deallocate(CF0t_cmb,CF00_cmb)
	end if

  end subroutine CorrFcn



  subroutine del_pp(myrank)

    use :: flow_mod, only: FlowType

    integer,intent(in) :: myrank
	write(*,*) "module:pp_smdlt:del_pp"
    close(u1);close(u2);close(u3);close(u4);close(u5);close(u6);close(u7);close(u8);close(u9)
    close(u10);close(u11);close(u12);close(u13);close(u14);close(u15);close(u16);close(u17)
    close(u18);close(u19);close(u20);close(u21);close(u22);close(u23);close(u24);close(u25)
    close(u26);close(u27);close(u28);close(u29);close(u30);close(u31);close(u32);close(u33)
    close(u34);close(u35);close(u36);close(u37);close(u38);close(u39)

    if (StrCalc) then
      deallocate(rAvesqqetoe,rAveqetoe,rAvetauxy,rAvetauxxyy,rAvetauyyzz,sdrAvesqqetoe)
      deallocate(sdrAveqetoe,sdrAvetauxy,sdrAvetauxxyy,sdrAvetauyyzz)
      if (myrank == 0) then
        deallocate(rAvesqqetoeTot,rAveqetoeTot,rAvetauxyTot,rAvetauxxyyTot,rAvetauyyzzTot)
        deallocate(sdrAvesqqetoeTot,sdrAveqetoeTot,sdrAvetauxyTot,sdrAvetauxxyyTot,sdrAvetauyyzzTot)
      end if ! myrank
    end if ! StrCalc
    if (RgCalc) then
      deallocate(rAveRgSq,sdrAveRgSq)
      if (myrank == 0) deallocate(rAveRgSqTot,sdrAveRgSqTot)
    end if
    if ((FlowType == 'Equil').and.CoMDiff) then
      deallocate(rAveMSD,sdrAveMSD)
      if (myrank == 0) deallocate(rAveMSDTot,sdrAveMSDTot)
    end if
    if (CorFun) then
      select case (CF_mode)
        case ('Ree')
          deallocate(cf_ree)
        case ('Rg')
          deallocate(cf_rg)
      end select
    end if

  end subroutine del_pp

end module pp_smdlt
