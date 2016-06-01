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
subroutine HICalc(rvmrc,nseg,HITens,DiffTens,EVForceLaw,Fev)

  use :: prcn_mod
  use :: arry_mod, only: print_vector,print_matrix
  use :: inp_mod, only: hstar,zstar,dstar,LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond

  implicit real(wp) (a-h, o-z)
  real(wp),parameter :: PI=3.1415926535897958648_wp!4*atan(1.0_wp)
  real(wp),parameter :: sqrtPI=sqrt(PI)
  integer :: offset,offseti,offsetj
  character(len=10) :: HITens,EVForceLaw
  real(wp),dimension(3*(nseg+1),3*(nseg+1)),intent(out) :: DiffTens
!  real(wp),dimension(3*nseg),intent(in) :: q
  real(wp),dimension(3*(nseg+1)),intent(in) :: rvmrc
  real(wp),dimension(3*(nseg+1)),intent(out) :: Fev
!  real(wp),allocatable,dimension(:) :: qxtemp,qytemp,qztemp
  real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
  real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
  real(wp) :: rimrc(3),rjmrc(3),rijmag

  if (HITens == 'RPY') then
    ! For Rotne-Prager-Yamakawa Tensor:
    A=0.75*hstar*sqrtPI
    B=hstar**3*PI*sqrtPI/2
    C=(3.0_wp/2)*hstar**3*PI*sqrtPI
    D=2*sqrtPI*hstar
    E=9/(32*sqrtPI*hstar)
    F=3/(32*sqrtPI*hstar)
  elseif (HITens == 'OB') then
    ! For Oseen-Burgers Tensor:
    G=0.75*hstar*sqrtPI
  elseif (HITens == 'RegOB') then
    ! For Reguralized Oseen-Burgers Tensor (introduced in HCO book):
    G=0.75*hstar*sqrtPI
    O=4*PI*hstar**2/3
    P=14*PI*hstar**2/3
    R=8*PI**2*hstar**4
    S=2*PI*hstar**2
    T=R/3
  end if
  rmagmin=1.0e-7_wp ! The Minimum value accepted as the |rij|
!  qxtemp=0._wp;qytemp=0._wp;qztemp=0._wp
  if (EVForceLaw /= 'NoEV') then
    Fev=0._wp
    if (EVForceLaw == 'Gauss') then
      prefactor=zstar/dstar**5;denom=2*dstar**2
    elseif (EVForceLaw == 'LJ') then
      epsOVsig=LJ_eps/LJ_sig
      sigOVrtr=LJ_sig/LJ_rtr
      sigOVrtrto6=sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr
      LJ_prefactor_tr=4*epsOVsig*( 12*(sigOVrtrto6*sigOVrtrto6*sigOVrtr) - &
                                    6*(sigOVrtrto6*sigOVrtr))
    end if
  end if
  do jbead=1,nseg+1
    offset=3*(jbead-2)
    offsetj=3*(jbead-1)
    rjmrc=rvmrc(offsetj+1:offsetj+3)
    do ibead=1, jbead
      offseti=3*(ibead-1)
      ! Calculating the global location
      iglob=3*(ibead-1);jglob=3*(jbead-1)
      ! Calculation of rij and |rij|
      if (ibead == jbead) then
        DiffTens(iglob+1,jglob+1)=1.0_wp
        DiffTens(iglob+1,jglob+2)=0.0_wp
        DiffTens(iglob+1,jglob+3)=0.0_wp
        DiffTens(iglob+2,jglob+2)=1.0_wp
        DiffTens(iglob+2,jglob+3)=0.0_wp
        DiffTens(iglob+3,jglob+3)=1.0_wp
      else
        !----------------------------------------------------------------------------------!
        !Now to efficiently calculate the rij based on q, we construct a vector qtemp(i) to!
        !  save the conn. vectors from i. As the loop goes over i, the magnitudes of prior !
        !  indices are not needed to be calculated again:                                  !
        !----------------------------------------------------------------------------------!
        rimrc=rvmrc(offseti+1:offseti+3)
        rxij=rjmrc(1)-rimrc(1)
        ryij=rjmrc(2)-rimrc(2)
        rzij=rjmrc(3)-rimrc(3)
        rijmag=sqrt(rxij*rxij+ryij*ryij+rzij*rzij)
        rijmagto3=rijmag*rijmag*rijmag
        rijmagto5=rijmagto3*rijmag*rijmag
        if (rijmag<=rmagmin) then
          write(*,*) 'Warning: The |rij| is lower than the accepted value in &
                      routine HIEVCalc'
          write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rijmag,'|rij|min:',rmagmin
          write(*,*) ibead,jbead
          stop
        end if
        if (EVForceLaw /= 'NoEV') then
          ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
          ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.
          if (EVForceLaw == 'Gauss') then
            Fev(offseti+1)=Fev(offseti+1)-prefactor*rxij*exp(-rijmag*rijmag/denom)
            Fev(offsetj+1)=Fev(offsetj+1)+prefactor*rxij*exp(-rijmag*rijmag/denom)
            Fev(offseti+2)=Fev(offseti+2)-prefactor*ryij*exp(-rijmag*rijmag/denom)
            Fev(offsetj+2)=Fev(offsetj+2)+prefactor*ryij*exp(-rijmag*rijmag/denom)
            Fev(offseti+3)=Fev(offseti+3)-prefactor*rzij*exp(-rijmag*rijmag/denom)
            Fev(offsetj+3)=Fev(offsetj+3)+prefactor*rzij*exp(-rijmag*rijmag/denom)
          elseif (EVForceLaw == 'LJ') then
            if (abs(ibead-jbead).ge.minNonBond) then ! Only for Non-neighboring beads. Importrant!!
              if ((rijmag.gt.LJ_rtr).and.(rijmag.le.LJ_rc)) then
                sigOVrmag=LJ_sig/rijmag
                sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
                LJ_prefactor=4*epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
                                           6*(sigOVrmagto6*sigOVrmag))
                Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor*rxij/rijmag
                Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor*rxij/rijmag
                Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor*ryij/rijmag
                Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor*ryij/rijmag
                Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor*rzij/rijmag
                Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor*rzij/rijmag
              elseif (rijmag.le.LJ_rtr) then
                Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor_tr*rxij/rijmag
                Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor_tr*rxij/rijmag
                Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor_tr*ryij/rijmag
                Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor_tr*ryij/rijmag
                Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor_tr*rzij/rijmag
                Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor_tr*rzij/rijmag
              end if ! As else, nothing should be added to the corresponding force elements
            end if ! abs(ibead-jbead)
          end if ! EVForceLaw
        end if ! EVForceLaw
        if (HITens == 'RPY') then
          if (rijmag >= D) then
            Alpha=A/rijmag+B/rijmagto3;Beta=A/rijmagto3;Gamm=C/rijmagto5
            Zeta=Beta-Gamm
            Zeta12=Zeta*rxij*ryij;Zeta13=Zeta*rxij*rzij;Zeta23=Zeta*ryij*rzij
            DiffTens(iglob+1,jglob+1)=Alpha+Zeta*rxij*rxij
            DiffTens(iglob+1,jglob+2)=Zeta12;DiffTens(iglob+2,jglob+1)=Zeta12
            DiffTens(iglob+1,jglob+3)=Zeta13;DiffTens(iglob+3,jglob+1)=Zeta13        
            DiffTens(iglob+2,jglob+2)=Alpha+Zeta*ryij*ryij
            DiffTens(iglob+2,jglob+3)=Zeta23;DiffTens(iglob+3,jglob+2)=Zeta23
            DiffTens(iglob+3,jglob+3)=Alpha+Zeta*rzij*rzij
          else    
            Theta=1-E*rijmag;Xi=F/rijmag
            DiffTens(iglob+1,jglob+1)=1-E*rijmag+F*rxij*rxij/rijmag
            Xi12=F*rxij*ryij/rijmag;Xi13=F*rxij*rzij/rijmag;Xi23=F*ryij*rzij/rijmag
            DiffTens(iglob+1,jglob+2)=Xi12;DiffTens(iglob+2,jglob+1)=Xi12
            DiffTens(iglob+1,jglob+3)=Xi13;DiffTens(iglob+3,jglob+1)=Xi13
            DiffTens(iglob+2,jglob+2)=1-E*rijmag+F*ryij*ryij/rijmag
            DiffTens(iglob+2,jglob+3)=Xi23;DiffTens(iglob+3,jglob+2)=Xi23
            DiffTens(iglob+3,jglob+3)=1-E*rijmag+F*rzij*rzij/rijmag
          end if  
        elseif (HITens == 'Zimm') then
          Rho=sqrt(2.0_wp)*hstar*sqrt(1/abs(real(ibead-jbead,kind=wp)))
          DiffTens(iglob+1,jglob+1)=Rho
          DiffTens(iglob+1,jglob+2)=0.0_wp;DiffTens(iglob+2,jglob+1)=0.0_wp
          DiffTens(iglob+1,jglob+3)=0.0_wp;DiffTens(iglob+3,jglob+1)=0.0_wp
          DiffTens(iglob+2,jglob+2)=Rho
          DiffTens(iglob+2,jglob+3)=0.0_wp;DiffTens(iglob+3,jglob+2)=0.0_wp
          DiffTens(iglob+3,jglob+3)=Rho
        elseif (HITens == 'OB') then
          Psi=G/rijmag;Chi=G/rijmagto3
          Chi12=Chi*rxij*ryij;Chi13=Chi*rxij*rzij;Chi23=Chi*ryij*rzij
          DiffTens(iglob+1,jglob+1)=Psi+Chi*rxij*rxij
          DiffTens(iglob+1,jglob+2)=Chi12;DiffTens(iglob+2,jglob+1)=Chi12
          DiffTens(iglob+1,jglob+3)=Chi13;DiffTens(iglob+3,jglob+1)=Chi13
          DiffTens(iglob+2,jglob+2)=Psi+Chi*ryij*ryij
          DiffTens(iglob+2,jglob+3)=Chi23;DiffTens(iglob+3,jglob+2)=Chi23
          DiffTens(iglob+3,jglob+3)=Psi+Chi*rzij*rzij
        elseif (HITens == 'RegOB') then
          Omicron=G/(rijmag*(rijmag**2+O)**3)
          Upsilon=Omicron*(rijmag**6+P*rijmag**4+R*rijmag**2)
          Omega=Omicron*(rijmag**6+S*rijmag**4-T*rijmag**2)/(rijmag**2)
          Omega12=Omega*rxij*ryij;Omega13=Omega*rxij*rzij;Omega23=Omega*ryij*rzij
          DiffTens(iglob+1,jglob+1)=Upsilon+Omega*rxij*rxij
          DiffTens(iglob+1,jglob+2)=Omega12;DiffTens(iglob+2,jglob+1)=Omega12
          DiffTens(iglob+1,jglob+3)=Omega13;DiffTens(iglob+3,jglob+1)=Omega13
          DiffTens(iglob+2,jglob+2)=Upsilon+Omega*ryij*ryij
          DiffTens(iglob+2,jglob+3)=Omega23;DiffTens(iglob+3,jglob+2)=Omega23
          DiffTens(iglob+3,jglob+3)=Upsilon+Omega*rzij*rzij
        end if ! HITens
      end if ! ibead == jbead
    end do
  end do

end subroutine HICalc

subroutine EVCalc(rvmrc,nseg,EVForceLaw,Fev)

  use :: prcn_mod
  use :: arry_mod, only: print_vector,print_matrix
  use :: inp_mod, only: zstar,dstar,LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond

  implicit real(wp) (a-h, o-z)
  integer :: offset,offseti,offsetj 
  character(len=10) :: EVForceLaw 
!  real(wp),dimension(3*nseg),intent(in) :: q
  real(wp),dimension(3*(nseg+1)),intent(in) :: rvmrc
!  real(wp),allocatable,dimension(:) :: qxtemp,qytemp,qztemp
  real(wp),dimension(3*(nseg+1)),intent(out) :: Fev
  real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
  real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
  real(wp) :: rimrc(3),rjmrc(3),rijmag,rijmag2

  rmagmin=1.0e-7_wp ! The Minimum value accepted as the |rij|
  Fev=0._wp
  prefactor=zstar/dstar**5;denom=2*dstar**2
  if (EVForceLaw.eq.'LJ') then
    epsOVsig=LJ_eps/LJ_sig
    sigOVrtr=LJ_sig/LJ_rtr
    sigOVrtrto6=sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr
    LJ_prefactor_tr=4*epsOVsig*( 12*(sigOVrtrto6*sigOVrtrto6*sigOVrtr) - &
                                  6*(sigOVrtrto6*sigOVrtr))
  end if
  do jbead=1,nseg+1
    offset=3*(jbead-2)
    offsetj=3*(jbead-1)
    rjmrc=rvmrc(offsetj+1:offsetj+3)
    do ibead=1, jbead
      offseti=3*(ibead-1)
      if (ibead /= jbead) then
        !----------------------------------------------------------------------------------!
        !Now to efficiently calculate the rij based on q, we construct a vector qtemp(i) to!
        !  save the conn. vectors from i. As the loop goes over i, the magnitudes of prior !
        !  indecies are not needed to be calculated again:                                 !
        !----------------------------------------------------------------------------------!
        rimrc=rvmrc(offseti+1:offseti+3)
        rxij=rjmrc(1)-rimrc(1)
        ryij=rjmrc(2)-rimrc(2)
        rzij=rjmrc(3)-rimrc(3)
        rijmag2=rxij*rxij+ryij*ryij+rzij*rzij
        if (rijmag2<=rmagmin**2) then
          write(*,*) 'Warning: The |rij| is lower than the accepted value in &
                      routine EVCalc'
          write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rijmag2,'|rij|min:',rmagmin**2
          write(*,*) ibead,jbead
          stop
        end if
        ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
        ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.        
        if (EVForceLaw == 'Gauss') then
          Fev(offseti+1)=Fev(offseti+1)-prefactor*rxij*exp(-rijmag2/denom)
          Fev(offsetj+1)=Fev(offsetj+1)+prefactor*rxij*exp(-rijmag2/denom)
          Fev(offseti+2)=Fev(offseti+2)-prefactor*ryij*exp(-rijmag2/denom)
          Fev(offsetj+2)=Fev(offsetj+2)+prefactor*ryij*exp(-rijmag2/denom)
          Fev(offseti+3)=Fev(offseti+3)-prefactor*rzij*exp(-rijmag2/denom)
          Fev(offsetj+3)=Fev(offsetj+3)+prefactor*rzij*exp(-rijmag2/denom)
        elseif (EVForceLaw == 'LJ') then
          rijmag=sqrt(rijmag2)
          if (abs(ibead-jbead).ge.minNonBond) then ! Only for Non-neighboring beads. Importrant!!
            if ((rijmag.gt.LJ_rtr).and.(rijmag.le.LJ_rc)) then
              sigOVrmag=LJ_sig/rijmag
              sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
              LJ_prefactor=4*epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
                                         6*(sigOVrmagto6*sigOVrmag))
              Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor*rxij/rijmag
              Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor*rxij/rijmag
              Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor*ryij/rijmag
              Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor*ryij/rijmag
              Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor*rzij/rijmag
              Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor*rzij/rijmag
            elseif (rijmag.le.LJ_rtr) then
              Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor_tr*rxij/rijmag
              Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor_tr*rxij/rijmag
              Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor_tr*ryij/rijmag
              Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor_tr*ryij/rijmag
              Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor_tr*rzij/rijmag
              Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor_tr*rzij/rijmag
            end if ! As else, nothing should be added to the corresponding force elements
          end if ! abs(ibead-jbead)
        end if ! EVForceLaw
      end if ! ibead /= jbead
    end do
  end do

end subroutine EVCalc

subroutine EVUpdate(Fev,qstar,Bmat,rvmrc,EVForceLaw,nseg,Fbarev)

  use :: prcn_mod
  use :: inp_mod, only: zstar,dstar,LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond

  implicit real(wp) (a-h, o-z)
  integer :: offset,offseti,offsetj
  character(len=10) :: EVForceLaw
  real(wp),dimension(3*nseg),intent(in) :: qstar
  real(wp),dimension(3*(nseg+1),3*nseg),intent(in) :: Bmat
  real(wp),dimension(3*(nseg+1)) :: rvmrc
!  real(wp),allocatable,dimension(:) :: qxtemp,qytemp,qztemp
  real(wp),allocatable,dimension(:) :: Fstarev
  real(wp),dimension(3*(nseg+1)),intent(in) :: Fev
  real(wp),dimension(3*(nseg+1)),intent(out) :: Fbarev
  real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
  real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
  real(wp) :: rimrc(3),rjmrc(3),rijmag,rijmag2

  allocate(Fstarev(3*(nseg+1)))   
  rmagmin=1.0e-6_wp ! The Minimum value accepted as the |rij|
  Fstarev=0._wp
  prefactor=zstar/dstar**5;denom=2*dstar**2
  if (EVForceLaw.eq.'LJ') then
    epsOVsig=LJ_eps/LJ_sig
    sigOVrtr=LJ_sig/LJ_rtr
    sigOVrtrto6=sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr*sigOVrtr
    LJ_prefactor_tr=4*epsOVsig*( 12*(sigOVrtrto6*sigOVrtrto6*sigOVrtr) - &
                                  6*(sigOVrtrto6*sigOVrtr))
  end if
  call gemv(Bmat,qstar,rvmrc)
  do jbead=1,nseg+1
    offset=3*(jbead-2)
    offsetj=3*(jbead-1)
    rjmrc=rvmrc(offsetj+1:offsetj+3)
    do ibead=1, jbead
      offseti=3*(ibead-1)
      if (ibead /= jbead) then
        rimrc=rvmrc(offseti+1:offseti+3)
        rxij=rjmrc(1)-rimrc(1)
        ryij=rjmrc(2)-rimrc(2)
        rzij=rjmrc(3)-rimrc(3)
        rijmag2=rxij*rxij+ryij*ryij+rzij*rzij
        if (rijmag2<=rmagmin**2) then
          write(*,*) 'Warning: The |rij| is lower than the accepted value in &
                      routine EVUpdate'
          write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rijmag2,'|rij|min:',rmagmin**2
          write(*,*) ibead,jbead
          stop
        end if
        if (EVForceLaw == 'Gauss') then
          Fstarev(offseti+1)=Fstarev(offseti+1)-prefactor*rxij*exp(-rijmag2/denom)
          Fstarev(offsetj+1)=Fstarev(offsetj+1)+prefactor*rxij*exp(-rijmag2/denom)
          Fstarev(offseti+2)=Fstarev(offseti+2)-prefactor*ryij*exp(-rijmag2/denom)
          Fstarev(offsetj+2)=Fstarev(offsetj+2)+prefactor*ryij*exp(-rijmag2/denom)
          Fstarev(offseti+3)=Fstarev(offseti+3)-prefactor*rzij*exp(-rijmag2/denom)
          Fstarev(offsetj+3)=Fstarev(offsetj+3)+prefactor*rzij*exp(-rijmag2/denom)
        elseif (EVForceLaw == 'LJ') then
          rijmag=sqrt(rijmag2)
          if (abs(ibead-jbead).ge.minNonBond) then
            if ((rijmag.gt.LJ_rtr).and.(rijmag.le.LJ_rc)) then
              sigOVrmag=LJ_sig/rijmag
              sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
              LJ_prefactor=4*epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
                                         6*(sigOVrmagto6*sigOVrmag))
              Fstarev(offseti+1)=Fstarev(offseti+1)-LJ_prefactor*rxij/rijmag
              Fstarev(offsetj+1)=Fstarev(offsetj+1)+LJ_prefactor*rxij/rijmag
              Fstarev(offseti+2)=Fstarev(offseti+2)-LJ_prefactor*ryij/rijmag
              Fstarev(offsetj+2)=Fstarev(offsetj+2)+LJ_prefactor*ryij/rijmag
              Fstarev(offseti+3)=Fstarev(offseti+3)-LJ_prefactor*rzij/rijmag
              Fstarev(offsetj+3)=Fstarev(offsetj+3)+LJ_prefactor*rzij/rijmag
            elseif (rijmag.le.LJ_rtr) then
              Fstarev(offseti+1)=Fstarev(offseti+1)-LJ_prefactor_tr*rxij/rijmag
              Fstarev(offsetj+1)=Fstarev(offsetj+1)+LJ_prefactor_tr*rxij/rijmag
              Fstarev(offseti+2)=Fstarev(offseti+2)-LJ_prefactor_tr*ryij/rijmag
              Fstarev(offsetj+2)=Fstarev(offsetj+2)+LJ_prefactor_tr*ryij/rijmag
              Fstarev(offseti+3)=Fstarev(offseti+3)-LJ_prefactor_tr*rzij/rijmag
              Fstarev(offsetj+3)=Fstarev(offsetj+3)+LJ_prefactor_tr*rzij/rijmag
            end if ! As else, nothing should be added to the corresponding force elements
          end if ! abs(ibead - jbead)
        end if ! EVForceLaw 
      end if ! ibead /= jbead
    end do
  end do
  Fbarev=0.5*(Fev+Fstarev)
  
  deallocate(Fstarev)

end subroutine EVUpdate
