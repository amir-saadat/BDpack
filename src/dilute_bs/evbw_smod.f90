submodule (hiev_mod) evbw_smod

  implicit none

contains

  module procedure evbw_init
    
    use :: inp_dlt, only: nseg,EV_bw,Aw,N_Ks,qmax

    ! Bead-wall excluded volume interaction
    if (EV_bw == 'Cubic') then
      evbw_prm%delw=0.5_wp*sqrt( (nseg**2-1._wp)/(2._wp*nseg) )
      evbw_prm%prf=Aw*N_Ks/(3*qmax*evbw_prm%delw**2)
    end if

    evbw_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  end subroutine evbw_init

  module procedure evbwcalc

    use :: inp_dlt, only: EV_bw

    integer :: osi

    osi=3*(i-1)

    if (EV_bw == 'Cubic') then

      if (ry <= evbw_prm%delw ) then
        Fev(osi+2)=Fev(osi+2)+3*evbw_prm%prf*(ry-evbw_prm%delw)**2
      end if

    elseif (EV_bw == 'Gaussian') then
    end if ! EV_bw

  end procedure evbwcalc


end submodule
