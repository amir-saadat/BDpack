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

  subroutine sprforce(q,nseg,ForceLaw,TruncMethod,Fseg)

    use :: inp_mod, only: b,qr_l,RWS_v,WLC_v,qmax,qr_lto2,q_l,RWS_C,RWS_D,&
                          WLC_A,WLC_B
    
    integer,intent(in) :: nseg
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
          stop
      end select
  
      call copy(qP,FsegP)! FsegP=qP
      call scal(FsegP,F) ! FsegP:=F*FsegP

    end do
  
  end subroutine sprforce

  subroutine sprupdate(root_f,PrScale,nroots,dt,RHSP,qstar,iseg,nseg,&
                       ForceLaw,TruncMethod,qbar,Fbarseg,Fbarbead,tpl&
                       &gy,Amat)

    use :: arry_mod, only: print_vector,print_matrix
    use :: root_mod, only: CubeRoot,root_fndr
    use :: inp_mod, only: b,qr_l,RWS_v,WLC_v,qmax,qr_lto2,q_l,RWS_C,RWS_D,&
                          WLC_A,WLC_B
  
    character(len=10),intent(in) :: ForceLaw,TruncMethod,tplgy
    real(wp),intent(in) :: dt
    integer,intent(in) :: PrScale,nroots,iseg,nseg
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
    integer :: iRHS,offset
  
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
            coeffs(7)=0.d0
            coeffs(8)=WLC_B*dt/(3*b**3)
            call root_fndr(coeffs,qmax,qmag)
            qr=qmag/qmax
            F=2.d0/3*(1/(1-qr**2)**2-7/(2*WLC_v*(1-qr**2))+WLC_A+WLC_B*(1-qr**2))
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
        if ((qmag<=0.0_wp).or.(qmag>=qmax)) then
          print '(" Oops! wrong root in SegForce Routine.")'
          print '(" Root No: ",i10," |q|: ",f14.7)',iRHS,qmag
          stop
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
    
    if (tplgy == 'Linear') then
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
    else
      call gemv(Amat,Fbarseg,Fbarbead,alpha=-1.d0,trans='T')
    end if

  end subroutine sprupdate

  subroutine bndforce(q,Fbnd,itime)

    use :: arry_mod, only: print_vector
    use :: inp_mod, only: WLC_v,WLC_C

    integer,intent(in) :: itime
    real(double),intent(in) :: q(:)
    real(double),intent(inout) :: Fbnd(:)
    real(double) :: thta(-1:1),cost(-1:1)
    real(double) :: qtmp(3,-2:1),qmg(-2:1),ehat(3,-2:1)
    integer :: nbead,ib,os

    nbead=size(Fbnd)/3

    do ib=1, nbead

      if (ib == 1) then
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
        if (nbead > 3) then
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
      elseif ((ib >= 3).and.(ib <= nbead-2)) then
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
                                1/qmg( 0)*( cost( 1)*ehat(:, 0)-ehat(:,1)) )
      elseif ((ib == nbead-1).and.(nbead > 3)) then
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
      elseif (ib == nbead) then
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

  end subroutine bndforce

  subroutine bndupdate(Fbnd,qst,Fbarbnd,itime)

    integer,intent(in) :: itime
    real(double),intent(in) :: qst(:),Fbnd(:)
    real(double),intent(out) :: Fbarbnd(:)
    real(double),allocatable :: Fstbnd(:)
    integer :: nbeadx3

    nbeadx3=size(Fbnd)
    allocate(Fstbnd(nbeadx3))

    call bndforce(qst,Fstbnd,itime)
    Fbarbnd=0.5d0*(Fbnd+Fstbnd)

    deallocate(Fstbnd)

  end subroutine bndupdate


end module force_mod
