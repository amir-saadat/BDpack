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
module force_mod

  use :: prcn_mod

  implicit none

contains

  subroutine sprforce(id,q,nseg,ForceLaw,TruncMethod,Fseg)

    use :: mpi_mod, only: forced_exit
    use :: inp_dlt, only: b,qr_l,RWS_v,WLC_v,qmax,qr_lto2,q_l,RWS_C,RWS_D,&
                          WLC_A,WLC_B

    integer,intent(in) :: id,nseg
    character(len=10) :: ForceLaw,TruncMethod
    real(wp),dimension(3*nseg),intent(in),target :: q
    real(wp),dimension(3*nseg),intent(inout),target :: Fseg
    real(wp),dimension(:), pointer    :: qP,FsegP
    real(wp) :: qmag,F,qmagsq,qr,qrsq
    integer :: iseg,offset

    do iseg=1,nseg

      offset=3*(iseg-1)
      qP => q(offset+1:offset+3)
      qmagsq=dot(qP,qP)!;qmag=sqrt(qmag)
      FsegP => Fseg(offset+1:offset+3)
      qrsq=qmagsq/b

      select case (ForceLaw)
        case ('FENE')
          if (qmagsq.lt.q_l**2) then
            F=1/(1-qrsq)
          else
            qr=sqrt(qrsq)
            if (TruncMethod.eq.'Linear') then
              F=(1/qr)*( qr_l/(1-qr_lto2) + (1+qr_lto2)/(1-qr_lto2)**2*&
                (qr-qr_l) )
            elseif (TruncMethod.eq.'Cnst') then
              F=(1/qr)*( qr_l/(1-qr_lto2) )
            else
               print '(" Warning!!: Wrong Truncation method for spring F-E.")'
            end if
          end if
        case ('WLC_MS')
          qr=sqrt(qrsq)
          F=2/(3*qr)*(0.25/((1-qr)**2)-0.25+qr)
        case ('WLC_UD','WLC_GEN')
          F=2._wp/3*(1/(1-qrsq)**2-7/(2*WLC_v*(1-qrsq))+WLC_A+WLC_B*(1-qrsq))
        case ('ILCCP')
          if (qmagsq < q_l**2) then
            F=(1-qrsq/3)/(1-qrsq)
          else
            qr=sqrt(qrsq)
            if (TruncMethod == 'Linear') then
              F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) + (3+qr_lto2*&
                              qr_lto2)/(3*(1-qr_lto2)**2)*(qr-qr_l) )
            elseif (TruncMethod == 'Cnst') then
              F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) )
            else
               print '(" Warning!!: Wrong Truncation method for spring F-E.")'
            end if
          end if
        case ('RWS')
          F=RWS_C/3*(1-RWS_D/RWS_C*qrsq)/(1-qrsq)
        case ('Hookean')
          F = 1._wp
        case default
          print '(" Not a valid Force Law.")'
          ! stop
          call forced_exit(id)
      end select

      call copy(qP,FsegP)! FsegP=qP
      call scal(FsegP,F) ! FsegP:=F*FsegP

    end do

  end subroutine sprforce

  subroutine sprupdate(id,root_f,PrScale,nroots,dt,RHSP,qstar,iseg,nseg,&
                       ForceLaw,TruncMethod,qbar,Fbarseg,Fbarbead,tplgy,&
                       Amat,nseg_bb,nseg_ar,Ia,Na,itime)

    use :: mpi_mod, only: forced_exit
    use :: arry_mod, only: print_vector,print_matrix
    use :: root_mod, only: CubeRoot,root_fndr
    use :: inp_dlt, only: b,qr_l,RWS_v,WLC_v,qmax,qr_lto2,q_l,RWS_C,RWS_D,&
                          WLC_A,WLC_B

    character(len=10),intent(in) :: ForceLaw,TruncMethod,tplgy
    real(wp),intent(in) :: dt
    integer,intent(in) :: id,PrScale,nroots,iseg,nseg,nseg_bb,nseg_ar,Ia(:),Na,itime
    real(wp),dimension(3),intent(in) :: RHSP
    real(wp),dimension(3*(nseg+1)),intent(inout)    :: Fbarbead
    real(wp),dimension(3*(nseg)),intent(in),target  :: qstar
    real(wp),dimension(:),pointer                   :: qcheck
    real(wp),dimension(3*nseg),intent(inout),target :: qbar,Fbarseg
    real(wp),dimension(:), pointer                  :: qbarP,FbarsegP
    real(wp),dimension(PrScale*nroots) :: root_f
    real(wp),dimension(3*nseg,3*(nseg+1)),intent(in) :: Amat
    real(wp) :: denom
    real(wp) :: RHSmag,qcheckmag,a1,a2,a3,qmag,F,coeffs(8),qr,slope,qden
    integer :: iRHS,offset,iarm,iseg_ar,offset_ar

    RHSmag=dot(RHSP,RHSP);RHSmag=sqrt(RHSmag)
    iRHS=int(RHSmag/(0.01_wp/PrScale)) ! Based on LookUp table definitions
    offset=3*(iseg-1)
    qbarP => qbar(offset+1:offset+3)
    FbarsegP => Fbarseg(offset+1:offset+3)
    !-------------------------------------------------------------------!
    ! qcheck is introduced to check the extent of extension of segments.!
    ! So in high extensions, the coefficients a1, a2, a3 are going to   !
    ! be evaluated exactly.                                             !
    !-------------------------------------------------------------------!
    qcheck => qstar(offset+1:offset+3)
    qcheckmag=dot(qcheck,qcheck);qcheckmag=sqrt(qcheckmag)

    if (ForceLaw /= 'Hookean') then
      if ((iRHS > (nroots*PrScale)) .or. (qcheckmag > (0.9_wp*qmax)) .or. &
          (PrScale > 100)) then
        select case (ForceLaw)
          case ('FENE')
            ! FENE
            a1=-RHSmag
            a2=-b*(1+dt/2)
            a3=b*RHSmag
            call CubeRoot(real(a1,kind=double),real(a2,kind=double),&
                          real(a3,kind=double),real(qmax,kind=double),qmag)
            F=1/(1-(qmag/qmax)**2)
          case ('WLC_MS')
            ! WLC, Worm Like Chain proposed by Marko and Siggia
            denom=1+dt/3
            a1=-(((2+0.75*dt)*qmax)+RHSmag)/denom
            a2=((1+0.5*dt)*b+2*qmax*RHSmag)/denom
            a3=-RHSmag*b/denom
            call CubeRoot(real(a1,kind=double),real(a2,kind=double),&
                          real(a3,kind=double),real(qmax,kind=double),qmag)
            F=2*qmax/(3*qmag)*(0.25/((1-qmag/qmax)**2)-0.25+qmag/qmax)
          case ('WLC_UD','WLC_GEN')
            ! WLC, Worm Like Chain modified by Underhill and Doyle
            coeffs(1)=RHSmag
            coeffs(2)=-1-dt/3+7*dt/(6*WLC_v)-dt/3*(WLC_A+WLC_B)
            coeffs(3)=-2*RHSmag/b
            coeffs(4)=2/b-7*dt/(6*WLC_v*b)+dt/b*(2*WLC_A/3+WLC_B)
            coeffs(5)=RHSmag/b**2
            coeffs(6)=-1/b**2-dt/b**2*(WLC_A/3+WLC_B)
            coeffs(7)=0._wp
            coeffs(8)=WLC_B*dt/(3*b**3)
            call root_fndr(real(coeffs,kind=double),real(qmax,kind=double),qmag)
            qr=qmag/qmax
            F=2._wp/3*(1/(1-qr**2)**2-7/(2*WLC_v*(1-qr**2))+WLC_A+WLC_B*(1-qr**2))
          case ('ILCCP')
            ! ILCCP, Inverse Langevin Chain (Cohen-Pade approximation)
            denom=1+dt/6
            a1=-RHSmag/denom
            a2=-(1+0.5*dt)*b/denom
            a3=RHSmag*b/denom
            call CubeRoot(real(a1,kind=double),real(a2,kind=double),&
                          real(a3,kind=double),real(qmax,kind=double),qmag)
            F=(1-(qmag/qmax)**2/3)/(1-(qmag/qmax)**2)
          case ('RWS')
            ! RWS, From Underhill and Doyle
            denom=1+RWS_D*dt/6
            a1=-RHSmag/denom
            a2=-(1+RWS_C/6*dt)*b/denom
            a3=RHSmag*b/denom
            call CubeRoot(real(a1,kind=double),real(a2,kind=double),&
                          real(a3,kind=double),real(qmax,kind=double),qmag)
            F = RWS_C/3*(1-RWS_D/RWS_C*(qmag/qmax)**2)/(1-(qmag/qmax)**2)
        end select
        if ((qmag<=0._wp).or.(qmag>=qmax)) then
          print '(" Oops! wrong root in SegForce Routine.")'
          print '(" Root No: ",i20," |q|: ",f24.7)',iRHS,qmag
          ! stop
          call forced_exit(id)
        end if
      else
        slope=root_f(iRHS+2) - root_f(iRHS+1)
        qmag=root_f(iRHS+1)+(slope*((RHSmag/(0.01_wp/PrScale))-iRHS))
        qr=qmag/qmax
        select case (ForceLaw)
          case ('FENE')
            if (qmag < q_l) then
              F = 1/(1-qr**2)
            else
              if (TruncMethod == 'Linear') then
                F=(1/qr)*( qr_l/(1-qr_lto2) + &
                  (1+qr_lto2)/(1-qr_lto2)**2*(qr-qr_l) )
              elseif (TruncMethod == 'Cnst') then
                F=(1/qr)*( qr_l/(1-qr_lto2) )
              else
                 print '(" Wrong Truncation method for spring F-E.")'
              end if
            end if
          case ('WLC_MS')
            F=2/(3*qr)*(0.25/((1-qr)**2)-0.25+qr)
          case ('WLC_UD','WLC_GEN')
            F=2._wp/3*(1/(1-qr**2)**2-7/(2*WLC_v*(1-qr**2))+WLC_A+WLC_B*(1-qr**2))
          case ('ILCCP')
            if (qmag < q_l) then
              F=(1-qr**2/3)/(1-qr**2)
            else
              if (TruncMethod == 'Linear') then
                F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) + &
                  (3+qr_lto2*qr_lto2)/(3*(1-qr_lto2)**2)*(qr-qr_l) )
              elseif (TruncMethod == 'Cnst') then
                F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) )
              else
                 print '(" Wrong Truncation method for spring F-E.")'
              end if
            end if
          case ('RWS')
            F=RWS_C/3*(1-RWS_D/RWS_C*qr**2)/(1-qr**2)
        end select
      end if
    else
      ! Note that for the Hookean case the force is always 1.
      F=1._wp
    end if

    ! Extra caution, Can be waived.
    qden=1+dt*F/2
    call copy(RHSP,qbarP)! qbarP=RHSP
    qbarP=1/qden*qbarP
    FbarsegP=qbarP
    qmag=dot(qbarP,qbarP);qmag=sqrt(qmag)
    qr=qmag/qmax
    select case (ForceLaw)
      case ('FENE')
        if (qmag < q_l) then
          F=1/(1-qr**2)
        else
          if (TruncMethod == 'Linear') then
            F=(1/qr)*( qr_l/(1-qr_lto2) + (1+qr_lto2)/(1-qr_lto2)**2*(qr-qr_l) )
          elseif (TruncMethod == 'Cnst') then
            F=(1/qr)*( qr_l/(1-qr_lto2) )
          else
             print '(" Wrong Truncation method for spring F-E.")'
          end if
        end if
      case ('WLC_MS')
        F=2/(3*qr)*(0.25/((1-qr)**2)-0.25+qr)
      case ('WLC_UD','WLC_GEN')
        F=2._wp/3*(1/(1-qr**2)**2-7/(2*WLC_v*(1-qr**2))+WLC_A+WLC_B*(1-qr**2))
      case ('ILCCP')
        if (qmag < q_l) then
          F=(1-qr**2/3)/(1-qr**2)
        else
          if (TruncMethod == 'Linear') then
            F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) + &
              (3+qr_lto2*qr_lto2)/(3*(1-qr_lto2)**2)*(qr-qr_l) )
          elseif (TruncMethod == 'Cnst') then
            F=(1/qr)*( qr_l*(3-qr_lto2)/(3*(1-qr_lto2)) )
          else
             print '(" Wrong Truncation method for spring F-E.")'
          end if
        end if
      case ('RWS')
        F=RWS_C/3*(1-RWS_D/RWS_C*qr**2)/(1-qr**2)
    end select
    call scal(FbarsegP,F)       ! FbarsegP:=F*qbarP  !

    select case (tplgy)

      case ('Linear')
        ! As the change in Fseg(i) only affects Fbead(i:i+1):
        if (nseg == 1) then
          Fbarbead(1:3)=Fbarseg(1:3)
          Fbarbead(4:6)=-Fbarseg(1:3)
        else
          if (iseg == 1) then
            Fbarbead(1:3)=Fbarseg(1:3)
            Fbarbead(4:6)=Fbarseg(4:6)-Fbarseg(1:3)
          elseif (iseg == nseg) then
            Fbarbead(offset+1:offset+3)=Fbarseg(offset+1:offset+3)-&
                                        Fbarseg(offset-2:offset)
            Fbarbead(offset+4:offset+6)=-Fbarseg(offset+1:offset+3)
          else
            Fbarbead(offset+1:offset+3)=Fbarseg(offset+1:offset+3)-&
                                        Fbarseg(offset-2:offset)
            Fbarbead(offset+4:offset+6)=Fbarseg(offset+4:offset+6)-&
                                        Fbarseg(offset+1:offset+3)
          end if
        end if

      case ('Comb')
!        call gemv(Amat,Fbarseg,Fbarbead,alpha=-1.d0,trans='T')
        if (iseg <= nseg_bb) then
          if (iseg == 1) then
            Fbarbead(1:3)=Fbarseg(1:3)
            Fbarbead(4:6)=Fbarseg(4:6)-Fbarseg(1:3)
          elseif (iseg == nseg_bb) then
            Fbarbead(offset+1:offset+3)=Fbarseg(offset+1:offset+3)-&
                                        Fbarseg(offset-2:offset)
            Fbarbead(offset+4:offset+6)=-Fbarseg(offset+1:offset+3)
          else
            Fbarbead(offset+1:offset+3)=Fbarseg(offset+1:offset+3)-&
                                        Fbarseg(offset-2:offset)
            Fbarbead(offset+4:offset+6)=Fbarseg(offset+4:offset+6)-&
                                        Fbarseg(offset+1:offset+3)
          end if
          do iarm=1, Na
            offset_ar=(nseg_bb+(iarm-1)*nseg_ar)*3
            if (Ia(iarm+1) == iseg) then
              Fbarbead(offset+1:offset+3)=Fbarbead(offset+1:offset+3)+Fbarseg(offset_ar+1:offset_ar+3)
            elseif (Ia(iarm+1) == iseg+1) then
              Fbarbead(offset+4:offset+6)=Fbarbead(offset+4:offset+6)+Fbarseg(offset_ar+1:offset_ar+3)
            end if
          end do
        else ! iseg > nseg_bb
          iarm=(iseg-nseg_bb-1)/nseg_ar+1
          offset_ar=(Ia(iarm+1)-1)*3
          if (nseg_ar == 1) then
            Fbarbead(offset_ar+1:offset_ar+3)=Fbarseg(offset   +1:offset   +3)+&
                                              Fbarseg(offset_ar+1:offset_ar+3)-&
                                              Fbarseg(offset_ar-2:offset_ar)
            Fbarbead(offset+4:offset+6)=-Fbarseg(offset+1:offset+3)
          else
            iseg_ar=iseg-nseg_bb-(iarm-1)*nseg_ar
            if (iseg_ar == 1) then
              Fbarbead(offset_ar+1:offset_ar+3)=Fbarseg(offset   +1:offset   +3)+&
                                                Fbarseg(offset_ar+1:offset_ar+3)-&
                                                Fbarseg(offset_ar-2:offset_ar)
              Fbarbead(offset+4:offset+6)=Fbarseg(offset+4:offset+6)-Fbarseg(offset+1:offset+3)
            elseif (iseg_ar == nseg_ar) then
              Fbarbead(offset+1:offset+3)= Fbarseg(offset+1:offset+3)-Fbarseg(offset-2:offset)
              Fbarbead(offset+4:offset+6)=-Fbarseg(offset+1:offset+3)
            else
              Fbarbead(offset+1:offset+3)=Fbarseg(offset+1:offset+3)-Fbarseg(offset-2:offset)
              Fbarbead(offset+4:offset+6)=Fbarseg(offset+4:offset+6)-Fbarseg(offset+1:offset+3)
            end if
          end if
        end if

      case default
        call gemv(Amat,Fbarseg,Fbarbead,alpha=-1._wp,trans='T')
    end select

  end subroutine sprupdate

  subroutine bndforce(nbead_bb,q,Fbnd,itime)

    use :: arry_mod, only: print_vector
    use :: inp_dlt, only: WLC_v,WLC_C,tplgy,nseg_ar,Na,Ia,srf_tet

    integer,intent(in) :: itime
    real(wp),intent(in) :: q(:)
    real(wp),intent(inout) :: Fbnd(:)
    real(wp) :: thta(-1:1),cost(-1:1),thtal,thtar,costl,costr
    real(wp) :: thta_s,cost_s
    real(wp) :: qtmp(3,-2:1),qmg(-2:1),ehat(3,-2:1)
    real(wp) :: qtmpl(3),qtmpr(3),qmgl,qmgr,ehatl(3),ehatr(3)
    integer :: nbead,ib,os,nbead_bb,osl,iarm

!call print_vector(q,'q')

    do ib=1, nbead_bb

!      if (((ib == 1).and.(nbead_bb > 2)) .or.
!          ((ib == 1).and.srf_tet)) then
      if ((ib == 1).and.(nbead_bb > 2)) then
        qtmp(:,0)=q(1:3)
        qtmp(:,1)=q(4:6)
        qmg(0)=sqrt(dot(qtmp(:,0),qtmp(:,0)))
        qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
        ehat(:,0)=qtmp(:,0)/qmg(0)
        ehat(:,1)=qtmp(:,1)/qmg(1)
        thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
        cost(1)=cos(thta(1))
        Fbnd(1:3)=WLC_C/qmg(0)*(cost(1)*ehat(:,0)-ehat(:,1))
      elseif (ib == 2) then
        qtmp(:,-1)=qtmp(:,0)
        qtmp(:, 0)=qtmp(:,1)
        qmg(-1)=qmg(0)
        qmg( 0)=qmg(1)
        ehat(:,-1)=ehat(:,0)
        ehat(:, 0)=ehat(:,1)
        cost(0)=cost(1)
        if (nbead_bb > 3) then
          qtmp(:,1)=q(7:9)
          qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
          ehat(:,1)=qtmp(:,1)/qmg(1)
          thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
          cost(1)=cos(thta(1))
          Fbnd(4:6)=WLC_C*( (ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg( 0))- &
                             ehat(:,-1)*(1/qmg( 0)+cost(0)/qmg(-1)))+&
                             1/qmg(0)*(cost(1)*ehat(:,0)-ehat(:,1)) )
        else
          Fbnd(4:6)=WLC_C*( (ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg( 0))- &
                             ehat(:,-1)*(1/qmg( 0)+cost(0)/qmg(-1))) )
        end if
      elseif ((ib >= 3).and.(ib <= nbead_bb-2)) then
        qtmp(:,-2)=qtmp(:,-1)
        qtmp(:,-1)=qtmp(:, 0)
        qtmp(:, 0)=qtmp(:, 1)
        qmg(-2)=qmg(-1)
        qmg(-1)=qmg( 0)
        qmg( 0)=qmg( 1)
        ehat(:,-2)=ehat(:,-1)
        ehat(:,-1)=ehat(:, 0)
        ehat(:, 0)=ehat(:, 1)
        cost(-1)=cost(0)
        cost( 0)=cost(1)
        qtmp(:,1)=q(ib*3+1:ib*3+3)
        qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
        ehat(:,1)=qtmp(:,1)/qmg(1)
        thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
        cost(1)=cos(thta(1))
        os=(ib-1)*3
        Fbnd(os+1:os+3)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg(0)) - &
                                ehat(:,-1)*(1/qmg(0)+cost(0)/qmg(-1)) + &
                                1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2))+&
                                1/qmg( 0)*( cost( 1)*ehat(:, 0)-ehat(:, 1)) )
      elseif ((ib == nbead_bb-1).and.(nbead_bb > 3)) then
        qtmp(:,-2)=qtmp(:,-1)
        qtmp(:,-1)=qtmp(:, 0)
        qtmp(:, 0)=qtmp(:, 1)
        qmg(-2)=qmg(-1)
        qmg(-1)=qmg( 0)
        qmg( 0)=qmg( 1)
        ehat(:,-2)=ehat(:,-1)
        ehat(:,-1)=ehat(:, 0)
        ehat(:, 0)=ehat(:, 1)
        cost(-1)=cost(0)
        cost( 0)=cost(1)
        os=(ib-1)*3
        Fbnd(os+1:os+3)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg(0)) - &
                                ehat(:,-1)*(1/qmg(0)+cost(0)/qmg(-1)) + &
                                1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2)))
      elseif (ib == nbead_bb) then
        qtmp(:,-2)=qtmp(:,-1)
        qtmp(:,-1)=qtmp(:, 0)
        qmg(-2)=qmg(-1)
        qmg(-1)=qmg( 0)
        ehat(:,-2)=ehat(:,-1)
        ehat(:,-1)=ehat(:, 0)
        cost(-1)=cost(0)
        os=(ib-1)*3
        Fbnd(os+1:os+3)=WLC_C*(1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2)))
      end if

    end do

!print *,'1:'
!call print_vector(Fbnd,'fbnd')

    if (tplgy == 'Comb') then
      do iarm=1, Na

        do ib=1, nseg_ar+2

          ! Before and after grafting point:
          if (ib == 1) then
            ! Left and right
            osl=(Ia(iarm+1)-2)*3
            qtmpl(:)= q(osl+1:osl+3)
            qtmpr(:)=-q(osl+4:osl+6)
            os=(nbead_bb-1+(iarm-1)*nseg_ar)*3
            qtmp(:,1)=q(os+1:os+3)
            qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
            qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
            qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
            ehatl(:)=qtmpl(:)/qmgl
            ehatr(:)=qtmpr(:)/qmgr
            ehat(:,1)=qtmp(:,1)/qmg(1)
            thtal=acos(dot(qtmpl,qtmp(:,1))/(qmgl*qmg(1)))
            thtar=acos(dot(qtmpr,qtmp(:,1))/(qmgr*qmg(1)))
!print *,'theta',thtal,thtar
            costl=cos(thtal)
            costr=cos(thtar)
            Fbnd(osl+1:osl+3)=Fbnd(osl+1:osl+3)+WLC_C/qmgl*&
                              (costl*ehatl(:)-ehat(:,1))
            Fbnd(osl+7:osl+9)=Fbnd(osl+7:osl+9)+WLC_C/qmgr*&
                              (costr*ehatr(:)-ehat(:,1))
!print *,'fl',WLC_C/qmgl*(costl*ehatl(:)-ehat(:,1))
!print *,'fr',WLC_C/qmgr*(costr*ehatr(:)-ehat(:,1))
          ! Grafting point:
          elseif (ib == 2) then
            qtmp(:,0)=qtmp(:,1)
            qmg(0)=qmg(1)
            ehat(:, 0)=ehat(:,1)
            if (nseg_ar > 1) then
              qtmp(:,1)=q(os+4:os+6)
              qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
              ehat(:,1)=qtmp(:,1)/qmg(1)
              thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
              cost(1)=cos(thta(1))
              Fbnd(osl+4:osl+6)=Fbnd(osl+4:osl+6)+WLC_C*&
                                ( (ehat(:,0)*(1/qmgl  +costl/qmg(0)) -&
                                   ehatl(: )*(1/qmg(0)+costl/qmgl  ))+&
                                  (ehat(:,0)*(1/qmgr  +costr/qmg(0)) -&
                                   ehatr(: )*(1/qmg(0)+costr/qmgr  ))+&
                                   1/qmg(0)*(cost(1)*ehat(:,0)-ehat(:,1)) )
            else
              Fbnd(osl+4:osl+6)=Fbnd(osl+4:osl+6)+WLC_C*&
                                ( (ehat(:,0)*(1/qmgl  +costl/qmg(0)) -&
                                   ehatl(: )*(1/qmg(0)+costl/qmgl  ))+&
                                  (ehat(:,0)*(1/qmgr  +costr/qmg(0)) -&
                                   ehatr(: )*(1/qmg(0)+costr/qmgr  )) )
            end if
          ! First bead on the arm, next to the grafting point:
          elseif (ib == 3) then
            if (nseg_ar == 1) then
              qtmp(:,-1)=qtmp(:,0)
              qmg(-1)=qmg(0)
              ehat(:,-1)=ehat(:,0)
              Fbnd(os+4:os+6)=WLC_C*(1/qmg(-1)*(-costl*ehat(:,-1)+ehatl(:))+&
                                     1/qmg(-1)*(-costr*ehat(:,-1)+ehatr(:)) )
            elseif (nseg_ar == 2) then
              qtmp(:,-1)=qtmp(:,0)
              qtmp(:, 0)=qtmp(:,1)
              qmg(-1)=qmg(0)
              qmg( 0)=qmg(1)
              ehat(:,-1)=ehat(:,0)
              ehat(:, 0)=ehat(:,1)
              cost( 0)=cost(1)
              Fbnd(os+4:os+6)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg(0)) - &
                                      ehat(:,-1)*(1/qmg(0)+cost(0)/qmg(-1)) + &
                                      1/qmg(-1)*(-costl*ehat(:,-1)+ehatl(:))+ &
                                      1/qmg(-1)*(-costr*ehat(:,-1)+ehatr(:)) )
            else ! nseg_ar > 2
              qtmp(:,-1)=qtmp(:,0)
              qtmp(:, 0)=qtmp(:,1)
              qmg(-1)=qmg(0)
              qmg( 0)=qmg(1)
              ehat(:,-1)=ehat(:,0)
              ehat(:, 0)=ehat(:,1)
              cost( 0)=cost(1)
              qtmp(:,1)=q(os+7:os+9)
              qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
              ehat(:,1)=qtmp(:,1)/qmg(1)
              thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
              cost(1)=cos(thta(1))
              Fbnd(os+4:os+6)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg(0)) - &
                                      ehat(:,-1)*(1/qmg(0)+cost(0)/qmg(-1)) + &
                                      1/qmg(-1)*(-costl*ehat(:,-1)+ehatl(:))+ &
                                      1/qmg(-1)*(-costr*ehat(:,-1)+ehatr(:))+ &
                                      1/qmg( 0)*( cost( 1)*ehat(:, 0)-ehat(:,1)) )
            end if
          elseif ((ib > 3).and.(ib <= nseg_ar)) then
            qtmp(:,-2)=qtmp(:,-1)
            qtmp(:,-1)=qtmp(:, 0)
            qtmp(:, 0)=qtmp(:, 1)
            qmg(-2)=qmg(-1)
            qmg(-1)=qmg( 0)
            qmg( 0)=qmg( 1)
            ehat(:,-2)=ehat(:,-1)
            ehat(:,-1)=ehat(:, 0)
            ehat(:, 0)=ehat(:, 1)
            cost(-1)=cost(0)
            cost( 0)=cost(1)
            os=(nbead_bb-1+(iarm-1)*nseg_ar+ib-1)*3
            qtmp(:,1)=q(os+1:os+3)
            qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1)))
            ehat(:,1)=qtmp(:,1)/qmg(1)
            thta(1)=acos(dot(qtmp(:,0),qtmp(:,1))/(qmg(0)*qmg(1)))
            cost(1)=cos(thta(1))
            os=(nbead_bb+(iarm-1)*nseg_ar+ib-3)*3
            Fbnd(os+1:os+3)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg( 0)) - &
                                    ehat(:,-1)*(1/qmg( 0)+cost(0)/qmg(-1)) + &
                                    1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2))+&
                                    1/qmg( 0)*( cost( 1)*ehat(:, 0)-ehat(:, 1)) )
          elseif (ib == nseg_ar+1) then
            qtmp(:,-2)=qtmp(:,-1)
            qtmp(:,-1)=qtmp(:, 0)
            qtmp(:, 0)=qtmp(:, 1)
            qmg(-2)=qmg(-1)
            qmg(-1)=qmg( 0)
            qmg( 0)=qmg( 1)
            ehat(:,-2)=ehat(:,-1)
            ehat(:,-1)=ehat(:, 0)
            ehat(:, 0)=ehat(:, 1)
            cost(-1)=cost(0)
            cost( 0)=cost(1)
            os=(nbead_bb+(iarm-1)*nseg_ar+ib-3)*3
            Fbnd(os+1:os+3)=WLC_C*( ehat(:, 0)*(1/qmg(-1)+cost(0)/qmg( 0)) - &
                                    ehat(:,-1)*(1/qmg( 0)+cost(0)/qmg(-1)) + &
                                    1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2)))
          elseif (ib == nseg_ar+2) then
            qtmp(:,-2)=qtmp(:,-1)
            qtmp(:,-1)=qtmp(:, 0)
            qmg(-2)=qmg(-1)
            qmg(-1)=qmg( 0)
            ehat(:,-2)=ehat(:,-1)
            ehat(:,-1)=ehat(:, 0)
            cost(-1)=cost(0)
            os=(nbead_bb+(iarm-1)*nseg_ar+ib-3)*3
            Fbnd(os+1:os+3)=WLC_C*(1/qmg(-1)*(-cost(-1)*ehat(:,-1)+ehat(:,-2)))
          end if

        end do ! ib

      end do ! iarm
    end if ! tplgy == 'Comb'

  end subroutine bndforce

  subroutine bndupdate(nbead_bb,Fbnd,qst,Fbarbnd,itime)

    use :: arry_mod, only: print_vector,print_matrix

    integer,intent(in) :: itime
    real(wp),intent(in) :: qst(:),Fbnd(:)
    real(wp),intent(out) :: Fbarbnd(:)
    real(wp),allocatable :: Fstbnd(:)
    integer :: nbeadx3,nbead_bb

    nbeadx3=size(Fbnd)

    allocate(Fstbnd(nbeadx3))
    call bndforce(nbead_bb,qst,Fstbnd,itime)
    Fbarbnd=0.5_wp*(Fbnd+Fstbnd)
    deallocate(Fstbnd)

  end subroutine bndupdate

  subroutine tetforce(rvmrc,rc,D,dt,Ftet,rf0,itime)

    use :: arry_mod, only: print_vector,print_matrix

    real(wp),intent(in) :: rvmrc(:),rc(:),D(:,:),dt,rf0(3)
    integer,intent(in) :: itime
    real(wp),intent(inout) :: Ftet(3)
    real(wp) :: delr(3),D1(3,3)
    integer :: ipiv(6),info

    delr(1:3)=(rvmrc(1:3)+rc(1:3)-rf0(1:3))/(0.25_wp*dt)

    D1=D(1:3,1:3)
    call sysv(D1,delr,'U',ipiv,info)
    if (info == 0) then
      Ftet=-delr
    else
      print '(" Incorrect Solution for D.delr=F by sysv in force module.")'
      print '(" sysv info: ",i0)',info
    end if

  end subroutine tetforce

  subroutine tetupdate(Ftet,rvmrc,rc,D,dt,Fbartet,rf0,itime)

    real(wp),intent(in) :: Ftet(3),rvmrc(:),rc(:),D(:,:),dt,rf0(3)
    integer,intent(in) :: itime
    real(wp),intent(out) :: Fbartet(3)
    real(wp) :: Fsttet(3)

    call tetforce(rvmrc,rc,D,dt,Fsttet,rf0,itime)
    Fbartet=0.5_wp*(Ftet+Fsttet)

  end subroutine tetupdate

end module force_mod
