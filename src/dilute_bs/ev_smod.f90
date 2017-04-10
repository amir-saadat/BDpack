submodule (hiev_mod) ev_smod

  implicit none

!  type :: dis
!    real(wp) :: x,y,z
!    real(wp) :: mag,mag2
!  end type

  interface

    !> Calculates EV force between particles
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param Fev EV force
    module subroutine evbbcalc(i,j,rij,Fev)
      integer,intent(in) :: i,j 
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: Fev(:)
    end subroutine evbbcalc


  end interface

contains

!  module procedure evcalc2
!
!    use :: inp_dlt, only: EVForceLaw,zstar,dstar
!    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond
!    
!    integer :: offset,offseti,offsetj,ibead,jbead 
!    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
!    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
!    real(wp) :: rimrc(3),rjmrc(3)
!    type(dis) :: rij
!!    real(wp) :: rxij,ryij,rzij,rijmag,rijmag2
!
!  
!    Fev=0._wp
!    do jbead=1,nseg+1
!      offset=3*(jbead-2)
!      offsetj=3*(jbead-1)
!      rjmrc=rvmrc(offsetj+1:offsetj+3)
!      do ibead=1, jbead
!        offseti=3*(ibead-1)
!        if (ibead /= jbead) then
!          rimrc=rvmrc(offseti+1:offseti+3)
!          rij%x=rjmrc(1)-rimrc(1)
!          rij%y=rjmrc(2)-rimrc(2)
!          rij%z=rjmrc(3)-rimrc(3)
!          rij%mag2=rij%x**2+rij%y**2+rij%z**2
!          rij%mag=sqrt(rij%mag2)
!          if (rij%mag2<=evbb_prm%rmagmin**2) then
!            write(*,*) 'Warning: The |rij| is lower than the accepted value in routine EVCalc'
!            write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rij%mag2,'|rij|min:',evbb_prm%rmagmin**2
!            write(*,*) ibead,jbead
!            stop
!          end if
!
!          call evbbcalc(ibead,jbead,rij,Fev)
!
!
!
!
!!          ! Note that the forces are repulsive and for Fi=+ if rij=ri-rj. But, with
!!          ! the algorithm used above we have used rij=rj-ri, hence we set Fi=-.        
!!          if (EVForceLaw == 'Gauss') then
!!            Fev(offseti+1)=Fev(offseti+1)-evbb_prm%prefactor*rij%x*exp(-rij%mag2/evbb_prm%denom)
!!            Fev(offsetj+1)=Fev(offsetj+1)+evbb_prm%prefactor*rij%x*exp(-rij%mag2/evbb_prm%denom)
!!            Fev(offseti+2)=Fev(offseti+2)-evbb_prm%prefactor*rij%y*exp(-rij%mag2/evbb_prm%denom)
!!            Fev(offsetj+2)=Fev(offsetj+2)+evbb_prm%prefactor*rij%y*exp(-rij%mag2/evbb_prm%denom)
!!            Fev(offseti+3)=Fev(offseti+3)-evbb_prm%prefactor*rij%z*exp(-rij%mag2/evbb_prm%denom)
!!            Fev(offsetj+3)=Fev(offsetj+3)+evbb_prm%prefactor*rij%z*exp(-rij%mag2/evbb_prm%denom)
!!          elseif (EVForceLaw == 'LJ') then
!!            rij%mag=sqrt(rij%mag2)
!!            if (abs(ibead-jbead).ge.minNonBond) then ! Only for Non-neighboring beads. Importrant!!
!!              if ((rij%mag.gt.LJ_rtr).and.(rij%mag.le.LJ_rc)) then
!!                sigOVrmag=LJ_sig/rij%mag
!!                sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
!!                LJ_prefactor=4*evbb_prm%epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
!!                                                  6*(sigOVrmagto6*sigOVrmag))
!!                Fev(offseti+1)=Fev(offseti+1)-LJ_prefactor*rij%x/rij%mag
!!                Fev(offsetj+1)=Fev(offsetj+1)+LJ_prefactor*rij%x/rij%mag
!!                Fev(offseti+2)=Fev(offseti+2)-LJ_prefactor*rij%y/rij%mag
!!                Fev(offsetj+2)=Fev(offsetj+2)+LJ_prefactor*rij%y/rij%mag
!!                Fev(offseti+3)=Fev(offseti+3)-LJ_prefactor*rij%z/rij%mag
!!                Fev(offsetj+3)=Fev(offsetj+3)+LJ_prefactor*rij%z/rij%mag
!!              elseif (rij%mag.le.LJ_rtr) then
!!                Fev(offseti+1)=Fev(offseti+1)-evbb_prm%LJ_prefactor_tr*rij%x/rij%mag
!!                Fev(offsetj+1)=Fev(offsetj+1)+evbb_prm%LJ_prefactor_tr*rij%x/rij%mag
!!                Fev(offseti+2)=Fev(offseti+2)-evbb_prm%LJ_prefactor_tr*rij%y/rij%mag
!!                Fev(offsetj+2)=Fev(offsetj+2)+evbb_prm%LJ_prefactor_tr*rij%y/rij%mag
!!                Fev(offseti+3)=Fev(offseti+3)-evbb_prm%LJ_prefactor_tr*rij%z/rij%mag
!!                Fev(offsetj+3)=Fev(offsetj+3)+evbb_prm%LJ_prefactor_tr*rij%z/rij%mag
!!              end if ! As else, nothing should be added to the corresponding force elements
!!            end if ! abs(ibead-jbead)
!!          end if ! EVForceLaw
!        end if ! ibead /= jbead
!      end do
!    end do
!
!
!  end procedure evcalc2
!
!  module procedure evupdate2
!
!    use :: inp_dlt, only: EVForceLaw,zstar,dstar
!    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond
!
!    integer :: offset,offseti,offsetj,ibead,jbead
!    real(wp),allocatable :: Fstarev(:)
!    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
!    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
!    real(wp) :: rimrc(3),rjmrc(3)
!    real(wp) :: rxij,ryij,rzij,rijmag,rijmag2
!  
!    allocate(Fstarev(3*(nseg+1)))   
!    Fstarev=0._wp
!    do jbead=1,nseg+1
!      offset=3*(jbead-2)
!      offsetj=3*(jbead-1)
!      rjmrc=rvmrc(offsetj+1:offsetj+3)
!      do ibead=1, jbead
!        offseti=3*(ibead-1)
!        if (ibead /= jbead) then
!          rimrc=rvmrc(offseti+1:offseti+3)
!          rxij=rjmrc(1)-rimrc(1)
!          ryij=rjmrc(2)-rimrc(2)
!          rzij=rjmrc(3)-rimrc(3)
!          rijmag2=rxij*rxij+ryij*ryij+rzij*rzij
!          if (rijmag2<=evbb_prm%rmagmin**2) then
!            write(*,*) 'Warning: The |rij| is lower than the accepted value in &
!                        routine EVUpdate'
!            write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rijmag2,'|rij|min:',evbb_prm%rmagmin**2
!            write(*,*) ibead,jbead
!            stop
!          end if
!          if (EVForceLaw == 'Gauss') then
!            Fstarev(offseti+1)=Fstarev(offseti+1)-evbb_prm%prefactor*rxij*exp(-rijmag2/evbb_prm%denom)
!            Fstarev(offsetj+1)=Fstarev(offsetj+1)+evbb_prm%prefactor*rxij*exp(-rijmag2/evbb_prm%denom)
!            Fstarev(offseti+2)=Fstarev(offseti+2)-evbb_prm%prefactor*ryij*exp(-rijmag2/evbb_prm%denom)
!            Fstarev(offsetj+2)=Fstarev(offsetj+2)+evbb_prm%prefactor*ryij*exp(-rijmag2/evbb_prm%denom)
!            Fstarev(offseti+3)=Fstarev(offseti+3)-evbb_prm%prefactor*rzij*exp(-rijmag2/evbb_prm%denom)
!            Fstarev(offsetj+3)=Fstarev(offsetj+3)+evbb_prm%prefactor*rzij*exp(-rijmag2/evbb_prm%denom)
!          elseif (EVForceLaw == 'LJ') then
!            rijmag=sqrt(rijmag2)
!            if (abs(ibead-jbead).ge.minNonBond) then
!              if ((rijmag.gt.LJ_rtr).and.(rijmag.le.LJ_rc)) then
!                sigOVrmag=LJ_sig/rijmag
!                sigOVrmagto6=sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag*sigOVrmag
!                LJ_prefactor=4*evbb_prm%epsOVsig*( 12*(sigOVrmagto6*sigOVrmagto6*sigOVrmag) - &
!                                                  6*(sigOVrmagto6*sigOVrmag))
!                Fstarev(offseti+1)=Fstarev(offseti+1)-LJ_prefactor*rxij/rijmag
!                Fstarev(offsetj+1)=Fstarev(offsetj+1)+LJ_prefactor*rxij/rijmag
!                Fstarev(offseti+2)=Fstarev(offseti+2)-LJ_prefactor*ryij/rijmag
!                Fstarev(offsetj+2)=Fstarev(offsetj+2)+LJ_prefactor*ryij/rijmag
!                Fstarev(offseti+3)=Fstarev(offseti+3)-LJ_prefactor*rzij/rijmag
!                Fstarev(offsetj+3)=Fstarev(offsetj+3)+LJ_prefactor*rzij/rijmag
!              elseif (rijmag.le.LJ_rtr) then
!                Fstarev(offseti+1)=Fstarev(offseti+1)-evbb_prm%LJ_prefactor_tr*rxij/rijmag
!                Fstarev(offsetj+1)=Fstarev(offsetj+1)+evbb_prm%LJ_prefactor_tr*rxij/rijmag
!                Fstarev(offseti+2)=Fstarev(offseti+2)-evbb_prm%LJ_prefactor_tr*ryij/rijmag
!                Fstarev(offsetj+2)=Fstarev(offsetj+2)+evbb_prm%LJ_prefactor_tr*ryij/rijmag
!                Fstarev(offseti+3)=Fstarev(offseti+3)-evbb_prm%LJ_prefactor_tr*rzij/rijmag
!                Fstarev(offsetj+3)=Fstarev(offsetj+3)+evbb_prm%LJ_prefactor_tr*rzij/rijmag
!              end if ! As else, nothing should be added to the corresponding force elements
!            end if ! abs(ibead - jbead)
!          end if ! EVForceLaw 
!        end if ! ibead /= jbead
!      end do
!    end do
!    Fbarev=0.5*(Fev+Fstarev)
!    
!    deallocate(Fstarev)
!
!  end procedure evupdate2

  module procedure init_evbb

  end procedure init_evbb

  module procedure ev_init
    
    use :: inp_dlt, only: EV_bb,zstar,dstar
    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond

    ! Bead-bead excluded volume interaction
    if (EV_bb == 'Gauss') then
      evbb_prm%prefactor=zstar/dstar**5
      evbb_prm%denom=2*dstar**2
    elseif (EV_bb == 'LJ') then
      evbb_prm%epsOVsig=LJ_eps/LJ_sig
      evbb_prm%sigOVrtr=LJ_sig/LJ_rtr
      evbb_prm%sigOVrtrto6=evbb_prm%sigOVrtr*evbb_prm%sigOVrtr*evbb_prm%sigOVrtr* &
                         evbb_prm%sigOVrtr*evbb_prm%sigOVrtr*evbb_prm%sigOVrtr
      evbb_prm%LJ_prefactor_tr=4*evbb_prm%epsOVsig*                           &
                             ( 12*(evbb_prm%sigOVrtrto6*evbb_prm%sigOVrtrto6* &
                                   evbb_prm%sigOVrtr) -                       &
                                6*(evbb_prm%sigOVrtrto6*evbb_prm%sigOVrtr))
    end if
    evbb_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|

  end subroutine ev_init

  module procedure evcalc3
  
    if (rij%mag2<=evbb_prm%rmagmin**2) then
      write(*,*) 'Warning: The |rij| is lower than the accepted value in routine evcalc'
      write(*,'(1x,a,f7.2,1x,a,f7.2)') '|rij|:',rij%mag2,'|rij|min:',evbb_prm%rmagmin**2
      write(*,*) i,j
      stop
    end if

    call evbbcalc(i,j,rij,Fev)

  end procedure evcalc3

end submodule ev_smod
