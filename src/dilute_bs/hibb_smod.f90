submodule (hiev_mod:hi_smod) hibb_smod

  implicit none

contains

  module procedure init_hibb

    use :: inp_dlt, only: HITens,hstar

    real(wp),parameter :: PI=3.1415926535897958648_wp
    real(wp),parameter :: sqrtPI=sqrt(PI)

    if (HITens == 'RPY') then
      ! For Rotne-Prager-Yamakawa Tensor:
      this%A=0.75*hstar*sqrtPI
      this%B=hstar**3*PI*sqrtPI/2
      this%C=(3.0_wp/2)*hstar**3*PI*sqrtPI
      this%D=2*sqrtPI*hstar
      this%E=9/(32*sqrtPI*hstar)
      this%F=3/(32*sqrtPI*hstar)
    elseif (HITens == 'OB') then
      ! For Oseen-Burgers Tensor:
      this%G=0.75*hstar*sqrtPI
    elseif (HITens == 'RegOB') then
      ! For Reguralized Oseen-Burgers Tensor (introduced in HCO book):
      this%G=0.75*hstar*sqrtPI
      this%O=4*PI*hstar**2/3
      this%P=14*PI*hstar**2/3
      this%R=8*PI**2*hstar**4
      this%S=2*PI*hstar**2
      this%T=this%R/3
    end if
    this%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  end procedure init_hibb

  module procedure calc_hibb

    use :: inp_dlt, only: HITens,hstar

    integer :: osi,osj
    real(wp) :: rijmag3,rijmag5
    real(wp) :: Alpha,Beta,Gamm,Zeta,Zeta12,Zeta13,Zeta23
    real(wp) :: Theta,Xi,Xi12,Xi13,Xi23,Rho,Psi,Chi,Chi12,Chi13,Chi23
    real(wp) :: Omicron,Upsilon,Omega,Omega12,Omega13,Omega23

    ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
    ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.

    osi=3*(i-1)
    osj=3*(j-1)
    rijmag3=rij%mag2*rij%mag
    rijmag5=rij%mag2*rijmag3

    if (HITens == 'RPY') then
      if (rij%mag >= hi_prm%D) then
        Alpha=hi_prm%A/rij%mag+hi_prm%B/rijmag3
        Beta=hi_prm%A/rijmag3
        Gamm=hi_prm%C/rijmag5
        Zeta=Beta-Gamm
        Zeta12=Zeta*rij%x*rij%y;Zeta13=Zeta*rij%x*rij%z;Zeta23=Zeta*rij%y*rij%z
        DiffTens(osi+1,osj+1)=Alpha+Zeta*rij%x*rij%x
        DiffTens(osi+1,osj+2)=Zeta12;DiffTens(osi+2,osj+1)=Zeta12
        DiffTens(osi+1,osj+3)=Zeta13;DiffTens(osi+3,osj+1)=Zeta13        
        DiffTens(osi+2,osj+2)=Alpha+Zeta*rij%y*rij%y
        DiffTens(osi+2,osj+3)=Zeta23;DiffTens(osi+3,osj+2)=Zeta23
        DiffTens(osi+3,osj+3)=Alpha+Zeta*rij%z*rij%z
      else    
        Theta=1-hi_prm%E*rij%mag;Xi=hi_prm%F/rij%mag
        DiffTens(osi+1,osj+1)=1-hi_prm%E*rij%mag+hi_prm%F*rij%x*rij%x/rij%mag
        Xi12=hi_prm%F*rij%x*rij%y/rij%mag
        Xi13=hi_prm%F*rij%x*rij%z/rij%mag
        Xi23=hi_prm%F*rij%y*rij%z/rij%mag
        DiffTens(osi+1,osj+2)=Xi12;DiffTens(osi+2,osj+1)=Xi12
        DiffTens(osi+1,osj+3)=Xi13;DiffTens(osi+3,osj+1)=Xi13
        DiffTens(osi+2,osj+2)=1-hi_prm%E*rij%mag+hi_prm%F*rij%y*rij%y/rij%mag
        DiffTens(osi+2,osj+3)=Xi23;DiffTens(osi+3,osj+2)=Xi23
        DiffTens(osi+3,osj+3)=1-hi_prm%E*rij%mag+hi_prm%F*rij%z*rij%z/rij%mag
      end if  
    elseif (HITens == 'Zimm') then
      Rho=sqrt(2._wp)*hstar*sqrt(1/abs(real(i-j,kind=wp)))
      DiffTens(osi+1,osj+1)=Rho
      DiffTens(osi+1,osj+2)=0._wp;DiffTens(osi+2,osj+1)=0._wp
      DiffTens(osi+1,osj+3)=0._wp;DiffTens(osi+3,osj+1)=0._wp
      DiffTens(osi+2,osj+2)=Rho
      DiffTens(osi+2,osj+3)=0._wp;DiffTens(osi+3,osj+2)=0._wp
      DiffTens(osi+3,osj+3)=Rho
    elseif (HITens == 'OB') then
      Psi=hi_prm%G/rij%mag;Chi=hi_prm%G/rijmag3
      Chi12=Chi*rij%x*rij%y;Chi13=Chi*rij%x*rij%z;Chi23=Chi*rij%y*rij%z
      DiffTens(osi+1,osj+1)=Psi+Chi*rij%x*rij%x
      DiffTens(osi+1,osj+2)=Chi12;DiffTens(osi+2,osj+1)=Chi12
      DiffTens(osi+1,osj+3)=Chi13;DiffTens(osi+3,osj+1)=Chi13
      DiffTens(osi+2,osj+2)=Psi+Chi*rij%y*rij%y
      DiffTens(osi+2,osj+3)=Chi23;DiffTens(osi+3,osj+2)=Chi23
      DiffTens(osi+3,osj+3)=Psi+Chi*rij%z*rij%z
    elseif (HITens == 'RegOB') then
      Omicron=hi_prm%G/(rij%mag*(rij%mag**2+hi_prm%O)**3)
      Upsilon=Omicron*(rij%mag**6+hi_prm%P*rij%mag**4+hi_prm%R*rij%mag**2)
      Omega=Omicron*(rij%mag**6+hi_prm%S*rij%mag**4-hi_prm%T*rij%mag**2)/(rij%mag**2)
      Omega12=Omega*rij%x*rij%y;Omega13=Omega*rij%x*rij%z;Omega23=Omega*rij%y*rij%z
      DiffTens(osi+1,osj+1)=Upsilon+Omega*rij%x*rij%x
      DiffTens(osi+1,osj+2)=Omega12;DiffTens(osi+2,osj+1)=Omega12
      DiffTens(osi+1,osj+3)=Omega13;DiffTens(osi+3,osj+1)=Omega13
      DiffTens(osi+2,osj+2)=Upsilon+Omega*rij%y*rij%y
      DiffTens(osi+2,osj+3)=Omega23;DiffTens(osi+3,osj+2)=Omega23
      DiffTens(osi+3,osj+3)=Upsilon+Omega*rij%z*rij%z
    end if ! HITens

  end procedure calc_hibb


end submodule hibb_smod
