submodule (hiev_mod) hi_smod

  implicit none

  !--------------------------------------------
  !>>>> Interface to routines for hibb and hibw
  !--------------------------------------------
  interface
    !> initializes HI between beads
    !! \param this hibb object
    module subroutine init_hibb(this)
      class(hibb_t),intent(inout) :: this
    end subroutine init_hibb
    !> initializes HI between beads and the wall
    !! \param this hibw object
    module subroutine init_hibw(this)
      class(hibw_t),intent(inout) :: this
    end subroutine init_hibw
    !> calculates HI between beads
    !! \param this hibb object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    module subroutine calc_hibb(this,i,j,rij,DiffTens)
      class(hibb_t),intent(inout) :: this
      integer,intent(in) :: i,j 
      type(dis),intent(in) :: rij
      real(wp),intent(out) :: DiffTens(:,:)
    end subroutine calc_hibb
    !> calculates HI between beads and the wall
    !! \param this hibw object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    module subroutine calc_hibw(this,i,j,rij,DiffTens)
      class(hibw_t),intent(inout) :: this
      integer,intent(in) :: i,j 
      type(dis),intent(in) :: rij
      real(wp),intent(out) :: DiffTens(:,:)
    end subroutine calc_hibw
  end interface

contains

!  module procedure HICalc2
!
!    use :: arry_mod, only: print_vector,print_matrix
!    use :: inp_dlt, only: HITens,hstar
!    use :: inp_dlt, only: EVForceLaw,zstar,dstar
!    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond
!
!    integer :: offset,offseti,offsetj
!    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
!    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
!    real(wp) :: rimrc(3),rjmrc(3)
!    integer :: ibead,jbead,iglob,jglob
!    real(wp) :: rxij,ryij,rzij,rijmag,rijmagto3,rijmagto5
!    real(wp) :: Alpha,Beta,Gamm,Zeta,Zeta12,Zeta13,Zeta23
!    real(wp) :: Theta,Xi,Xi12,Xi13,Xi23,Rho,Psi,Chi,Chi12,Chi13,Chi23
!    real(wp) :: Omicron,Upsilon,Omega,Omega12,Omega13,Omega23
!
!    if (EVForceLaw /= 'NoEV') then
!      Fev=0._wp
!    end if
!    do jbead=1,nseg+1
!      offset=3*(jbead-2)
!      offsetj=3*(jbead-1)
!      rjmrc=rvmrc(offsetj+1:offsetj+3)
!      do ibead=1, jbead
!        offseti=3*(ibead-1)
!        ! Calculating the global location
!        iglob=3*(ibead-1);jglob=3*(jbead-1)
!        ! Calculation of rij and |rij|
!        if (ibead == jbead) then
!          DiffTens(iglob+1,jglob+1)=1._wp
!          DiffTens(iglob+1,jglob+2)=0._wp
!          DiffTens(iglob+1,jglob+3)=0._wp
!          DiffTens(iglob+2,jglob+2)=1._wp
!          DiffTens(iglob+2,jglob+3)=0._wp
!          DiffTens(iglob+3,jglob+3)=1._wp
!        else
!          rimrc=rvmrc(offseti+1:offseti+3)
!          rxij=rjmrc(1)-rimrc(1)
!          ryij=rjmrc(2)-rimrc(2)
!          rzij=rjmrc(3)-rimrc(3)
!          rijmag=sqrt(rxij*rxij+ryij*ryij+rzij*rzij)
!          rijmagto3=rijmag*rijmag*rijmag
!          rijmagto5=rijmagto3*rijmag*rijmag
!          if (rijmag<=hi_prm%rmagmin) then
!            write(*,*) 'Warning: The |rij| is lower than the accepted value in &
!                        routine HIEVCalc'
!            write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rijmag,'|rij|min:',hi_prm%rmagmin
!            write(*,*) ibead,jbead
!            stop
!          end if
!          if (EVForceLaw /= 'NoEV') then
!            ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
!            ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.
!            if (EVForceLaw == 'Gauss') then
!              Fev(offseti+1)=Fev(offseti+1)-ev_prm%prefactor*rxij*exp(-rijmag**2/ev_prm%denom)
!              Fev(offsetj+1)=Fev(offsetj+1)+ev_prm%prefactor*rxij*exp(-rijmag**2/ev_prm%denom)
!              Fev(offseti+2)=Fev(offseti+2)-ev_prm%prefactor*ryij*exp(-rijmag**2/ev_prm%denom)
!              Fev(offsetj+2)=Fev(offsetj+2)+ev_prm%prefactor*ryij*exp(-rijmag**2/ev_prm%denom)
!              Fev(offseti+3)=Fev(offseti+3)-ev_prm%prefactor*rzij*exp(-rijmag**2/ev_prm%denom)
!              Fev(offsetj+3)=Fev(offsetj+3)+ev_prm%prefactor*rzij*exp(-rijmag**2/ev_prm%denom)
!            elseif (EVForceLaw == 'LJ') then
!              if (abs(ibead-jbead).ge.minNonBond) then ! Only for Non-neighboring beads. Importrant!!
!                if ((rijmag.gt.LJ_rtr).and.(rijmag.le.LJ_rc)) then
!                  sigOVrmag=LJ_sig/rijmag
!                  sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
!                  LJ_prefactor=4*ev_prm%epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
!                                                    6*(sigOVrmagto6*sigOVrmag))
!                  Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor*rxij/rijmag
!                  Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor*rxij/rijmag
!                  Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor*ryij/rijmag
!                  Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor*ryij/rijmag
!                  Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor*rzij/rijmag
!                  Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor*rzij/rijmag
!                elseif (rijmag.le.LJ_rtr) then
!                  Fev(offseti+1)=Fev(offseti+1)-ev_prm%LJ_prefactor_tr*rxij/rijmag
!                  Fev(offsetj+1)=Fev(offsetj+1)+ev_prm%LJ_prefactor_tr*rxij/rijmag
!                  Fev(offseti+2)=Fev(offseti+2)-ev_prm%LJ_prefactor_tr*ryij/rijmag
!                  Fev(offsetj+2)=Fev(offsetj+2)+ev_prm%LJ_prefactor_tr*ryij/rijmag
!                  Fev(offseti+3)=Fev(offseti+3)-ev_prm%LJ_prefactor_tr*rzij/rijmag
!                  Fev(offsetj+3)=Fev(offsetj+3)+ev_prm%LJ_prefactor_tr*rzij/rijmag
!                end if ! As else, nothing should be added to the corresponding force elements
!              end if ! abs(ibead-jbead)
!            end if ! EVForceLaw
!          end if ! EVForceLaw
!          if (HITens == 'RPY') then
!            if (rijmag >= hi_prm%D) then
!              Alpha=hi_prm%A/rijmag+hi_prm%B/rijmagto3
!              Beta=hi_prm%A/rijmagto3
!              Gamm=hi_prm%C/rijmagto5
!              Zeta=Beta-Gamm
!              Zeta12=Zeta*rxij*ryij;Zeta13=Zeta*rxij*rzij;Zeta23=Zeta*ryij*rzij
!              DiffTens(iglob+1,jglob+1)=Alpha+Zeta*rxij*rxij
!              DiffTens(iglob+1,jglob+2)=Zeta12;DiffTens(iglob+2,jglob+1)=Zeta12
!              DiffTens(iglob+1,jglob+3)=Zeta13;DiffTens(iglob+3,jglob+1)=Zeta13        
!              DiffTens(iglob+2,jglob+2)=Alpha+Zeta*ryij*ryij
!              DiffTens(iglob+2,jglob+3)=Zeta23;DiffTens(iglob+3,jglob+2)=Zeta23
!              DiffTens(iglob+3,jglob+3)=Alpha+Zeta*rzij*rzij
!            else    
!              Theta=1-hi_prm%E*rijmag;Xi=hi_prm%F/rijmag
!              DiffTens(iglob+1,jglob+1)=1-hi_prm%E*rijmag+hi_prm%F*rxij*rxij/rijmag
!              Xi12=hi_prm%F*rxij*ryij/rijmag
!              Xi13=hi_prm%F*rxij*rzij/rijmag
!              Xi23=hi_prm%F*ryij*rzij/rijmag
!              DiffTens(iglob+1,jglob+2)=Xi12;DiffTens(iglob+2,jglob+1)=Xi12
!              DiffTens(iglob+1,jglob+3)=Xi13;DiffTens(iglob+3,jglob+1)=Xi13
!              DiffTens(iglob+2,jglob+2)=1-hi_prm%E*rijmag+hi_prm%F*ryij*ryij/rijmag
!              DiffTens(iglob+2,jglob+3)=Xi23;DiffTens(iglob+3,jglob+2)=Xi23
!              DiffTens(iglob+3,jglob+3)=1-hi_prm%E*rijmag+hi_prm%F*rzij*rzij/rijmag
!            end if  
!          elseif (HITens == 'Zimm') then
!            Rho=sqrt(2.0_wp)*hstar*sqrt(1/abs(real(ibead-jbead,kind=wp)))
!            DiffTens(iglob+1,jglob+1)=Rho
!            DiffTens(iglob+1,jglob+2)=0.0_wp;DiffTens(iglob+2,jglob+1)=0.0_wp
!            DiffTens(iglob+1,jglob+3)=0.0_wp;DiffTens(iglob+3,jglob+1)=0.0_wp
!            DiffTens(iglob+2,jglob+2)=Rho
!            DiffTens(iglob+2,jglob+3)=0.0_wp;DiffTens(iglob+3,jglob+2)=0.0_wp
!            DiffTens(iglob+3,jglob+3)=Rho
!          elseif (HITens == 'OB') then
!            Psi=hi_prm%G/rijmag;Chi=hi_prm%G/rijmagto3
!            Chi12=Chi*rxij*ryij;Chi13=Chi*rxij*rzij;Chi23=Chi*ryij*rzij
!            DiffTens(iglob+1,jglob+1)=Psi+Chi*rxij*rxij
!            DiffTens(iglob+1,jglob+2)=Chi12;DiffTens(iglob+2,jglob+1)=Chi12
!            DiffTens(iglob+1,jglob+3)=Chi13;DiffTens(iglob+3,jglob+1)=Chi13
!            DiffTens(iglob+2,jglob+2)=Psi+Chi*ryij*ryij
!            DiffTens(iglob+2,jglob+3)=Chi23;DiffTens(iglob+3,jglob+2)=Chi23
!            DiffTens(iglob+3,jglob+3)=Psi+Chi*rzij*rzij
!          elseif (HITens == 'RegOB') then
!            Omicron=hi_prm%G/(rijmag*(rijmag**2+hi_prm%O)**3)
!            Upsilon=Omicron*(rijmag**6+hi_prm%P*rijmag**4+hi_prm%R*rijmag**2)
!            Omega=Omicron*(rijmag**6+hi_prm%S*rijmag**4-hi_prm%T*rijmag**2)/(rijmag**2)
!            Omega12=Omega*rxij*ryij;Omega13=Omega*rxij*rzij;Omega23=Omega*ryij*rzij
!            DiffTens(iglob+1,jglob+1)=Upsilon+Omega*rxij*rxij
!            DiffTens(iglob+1,jglob+2)=Omega12;DiffTens(iglob+2,jglob+1)=Omega12
!            DiffTens(iglob+1,jglob+3)=Omega13;DiffTens(iglob+3,jglob+1)=Omega13
!            DiffTens(iglob+2,jglob+2)=Upsilon+Omega*ryij*ryij
!            DiffTens(iglob+2,jglob+3)=Omega23;DiffTens(iglob+3,jglob+2)=Omega23
!            DiffTens(iglob+3,jglob+3)=Upsilon+Omega*rzij*rzij
!          end if ! HITens
!        end if ! ibead == jbead
!      end do
!    end do
!    
!  end procedure HICalc2

  module procedure init_hi
    call init_hibb(this%hibb)
    call init_hibw(this%hibw)
  end procedure init_hi

  module procedure hi_init

    use :: inp_dlt, only: HITens,hstar

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    if (HITens == 'RPY') then
      ! For Rotne-Prager-Yamakawa Tensor:
      hi_prm%A=0.75*hstar*sqrtPI
      hi_prm%B=hstar**3*PI*sqrtPI/2
      hi_prm%C=(3.0_wp/2)*hstar**3*PI*sqrtPI
      hi_prm%D=2*sqrtPI*hstar
      hi_prm%E=9/(32*sqrtPI*hstar)
      hi_prm%F=3/(32*sqrtPI*hstar)
    elseif (HITens == 'OB') then
      ! For Oseen-Burgers Tensor:
      hi_prm%G=0.75*hstar*sqrtPI
    elseif (HITens == 'RegOB') then
      ! For Reguralized Oseen-Burgers Tensor (introduced in HCO book):
      hi_prm%G=0.75*hstar*sqrtPI
      hi_prm%O=4*PI*hstar**2/3
      hi_prm%P=14*PI*hstar**2/3
      hi_prm%R=8*PI**2*hstar**4
      hi_prm%S=2*PI*hstar**2
      hi_prm%T=hi_prm%R/3
    end if
    hi_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  end procedure hi_init

  module procedure hicalc3

    if (rij%mag<=hi_prm%rmagmin) then
      write(*,*) 'Warning: The |rij| is lower than the accepted value in routine hicalc'
      write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rij%mag,'|rij|min:',hi_prm%rmagmin
      write(*,*) i,j
      stop
    end if

    call calc_hibb(this%hibb,i,j,rij,DiffTens)
    call calc_hibw(this%hibw,i,j,rij,DiffTens)
    
  end procedure hicalc3

end submodule hi_smod
