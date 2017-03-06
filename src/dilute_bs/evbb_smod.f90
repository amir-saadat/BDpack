submodule (hiev_mod:ev_smod) evbb_smod


  implicit none

contains

  module procedure evbbcalc

    use :: inp_dlt, only: EV_bb,LJ_sig,LJ_rtr,LJ_rc,minNonBond

    integer :: osi,osj
    real(wp) :: LJ_prefactor,LJ_prefactor_tr
    real(wp) :: sigOVrmag,sigOVrmagto6

    ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
    ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.

    osi=3*(i-1)
    osj=3*(j-1)

    if (EV_bb == 'Gauss') then
      Fev(osi+1)=Fev(osi+1)-evbb_prm%prefactor*rij%x*exp(-rij%mag2/evbb_prm%denom)
      Fev(osj+1)=Fev(osj+1)+evbb_prm%prefactor*rij%x*exp(-rij%mag2/evbb_prm%denom)
      Fev(osi+2)=Fev(osi+2)-evbb_prm%prefactor*rij%y*exp(-rij%mag2/evbb_prm%denom)
      Fev(osj+2)=Fev(osj+2)+evbb_prm%prefactor*rij%y*exp(-rij%mag2/evbb_prm%denom)
      Fev(osi+3)=Fev(osi+3)-evbb_prm%prefactor*rij%z*exp(-rij%mag2/evbb_prm%denom)
      Fev(osj+3)=Fev(osj+3)+evbb_prm%prefactor*rij%z*exp(-rij%mag2/evbb_prm%denom)
    elseif (EV_bb == 'LJ') then
!      rij%mag=sqrt(rij%mag2)
      if (abs(i-j) >= minNonBond) then ! Only for Non-neighboring beads. Importrant!!
        if ((rij%mag > LJ_rtr).and.(rij%mag <= LJ_rc)) then
          sigOVrmag=LJ_sig/rij%mag
          sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
          LJ_prefactor=4*evbb_prm%epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
                                            6*(sigOVrmagto6*sigOVrmag))
          Fev(osi+1)=Fev(osi+1)-LJ_prefactor*rij%x/rij%mag
          Fev(osj+1)=Fev(osj+1)+LJ_prefactor*rij%x/rij%mag
          Fev(osi+2)=Fev(osi+2)-LJ_prefactor*rij%y/rij%mag
          Fev(osj+2)=Fev(osj+2)+LJ_prefactor*rij%y/rij%mag
          Fev(osi+3)=Fev(osi+3)-LJ_prefactor*rij%z/rij%mag
          Fev(osj+3)=Fev(osj+3)+LJ_prefactor*rij%z/rij%mag
        elseif (rij%mag.le.LJ_rtr) then
          Fev(osi+1)=Fev(osi+1)-evbb_prm%LJ_prefactor_tr*rij%x/rij%mag
          Fev(osj+1)=Fev(osj+1)+evbb_prm%LJ_prefactor_tr*rij%x/rij%mag
          Fev(osi+2)=Fev(osi+2)-evbb_prm%LJ_prefactor_tr*rij%y/rij%mag
          Fev(osj+2)=Fev(osj+2)+evbb_prm%LJ_prefactor_tr*rij%y/rij%mag
          Fev(osi+3)=Fev(osi+3)-evbb_prm%LJ_prefactor_tr*rij%z/rij%mag
          Fev(osj+3)=Fev(osj+3)+evbb_prm%LJ_prefactor_tr*rij%z/rij%mag
        end if ! As else, nothing should be added to the corresponding force elements
      end if ! abs(ibead-jbead)
    end if ! EV_bb

  end procedure evbbcalc


end submodule
