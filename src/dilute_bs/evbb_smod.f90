submodule (intrn_mod) evbb_smod

  implicit none

contains

  module procedure init_evbb

    use :: inp_dlt, only: EV_bb,zstar,dstar
    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond

    ! Bead-bead excluded volume interaction
    if (EV_bb == 'Gauss') then
      this%prefactor=zstar/dstar**5
      this%denom=2*dstar**2
    elseif (EV_bb == 'LJ') then
      this%epsOVsig=LJ_eps/LJ_sig
      this%sigOVrtr=LJ_sig/LJ_rtr
      this%sigOVrtrto6=this%sigOVrtr*this%sigOVrtr*this%sigOVrtr* &
                         this%sigOVrtr*this%sigOVrtr*this%sigOVrtr
      this%LJ_prefactor_tr=4*this%epsOVsig*                           &
                             ( 12*(this%sigOVrtrto6*this%sigOVrtrto6* &
                                   this%sigOVrtr) -                   &
                                6*(this%sigOVrtrto6*this%sigOVrtr))
    end if
    this%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  end procedure init_evbb

  module procedure calc_evbb

    use :: inp_dlt, only: EV_bb,LJ_sig,LJ_rtr,LJ_rc,minNonBond

    integer :: osi,osj
    real(wp) :: LJ_prefactor,LJ_prefactor_tr
    real(wp) :: sigOVrmag,sigOVrmagto6


    ! checking whether the distance is too small
    if (rij%mag2<=this%rmagmin**2) then
      write(*,*) 'Warning: The |rij| is lower than the accepted value in routine evcalc'
      write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rij%mag2,'|rij|min:',this%rmagmin**2
      write(*,*) i,j
      stop
    end if

    !! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
    !! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.

    ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj.

    osi=3*(i-1)
    osj=3*(j-1)

    if (EV_bb == 'Gauss') then
      !commented out EV force on bead j since calculation of D is dense
      Fev(osi+1)=Fev(osi+1)+this%prefactor*rij%x*exp(-rij%mag2/this%denom)
      !Fev(osj+1)=Fev(osj+1)-this%prefactor*rij%x*exp(-rij%mag2/this%denom)
      Fev(osi+2)=Fev(osi+2)+this%prefactor*rij%y*exp(-rij%mag2/this%denom)
      !Fev(osj+2)=Fev(osj+2)-this%prefactor*rij%y*exp(-rij%mag2/this%denom)
      Fev(osi+3)=Fev(osi+3)+this%prefactor*rij%z*exp(-rij%mag2/this%denom)
      !Fev(osj+3)=Fev(osj+3)-this%prefactor*rij%z*exp(-rij%mag2/this%denom)
    elseif (EV_bb == 'LJ') then
!      rij%mag=sqrt(rij%mag2)
      if (abs(i-j) >= minNonBond) then ! Only for Non-neighboring beads. Importrant!!
        if ((rij%mag > LJ_rtr).and.(rij%mag <= LJ_rc)) then
          sigOVrmag=LJ_sig/rij%mag
          sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
          LJ_prefactor=4*this%epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
                                            6*(sigOVrmagto6*sigOVrmag))
          Fev(osi+1)=Fev(osi+1)+LJ_prefactor*rij%x/rij%mag
          Fev(osj+1)=Fev(osj+1)-LJ_prefactor*rij%x/rij%mag
          Fev(osi+2)=Fev(osi+2)+LJ_prefactor*rij%y/rij%mag
          Fev(osj+2)=Fev(osj+2)-LJ_prefactor*rij%y/rij%mag
          Fev(osi+3)=Fev(osi+3)+LJ_prefactor*rij%z/rij%mag
          Fev(osj+3)=Fev(osj+3)-LJ_prefactor*rij%z/rij%mag
        elseif (rij%mag.le.LJ_rtr) then
          Fev(osi+1)=Fev(osi+1)+this%LJ_prefactor_tr*rij%x/rij%mag
          Fev(osj+1)=Fev(osj+1)-this%LJ_prefactor_tr*rij%x/rij%mag
          Fev(osi+2)=Fev(osi+2)+this%LJ_prefactor_tr*rij%y/rij%mag
          Fev(osj+2)=Fev(osj+2)-this%LJ_prefactor_tr*rij%y/rij%mag
          Fev(osi+3)=Fev(osi+3)+this%LJ_prefactor_tr*rij%z/rij%mag
          Fev(osj+3)=Fev(osj+3)-this%LJ_prefactor_tr*rij%z/rij%mag
        end if ! As else, nothing should be added to the corresponding force elements
      end if ! abs(ibead-jbead)
    end if ! EV_bb

  end procedure calc_evbb


end submodule
