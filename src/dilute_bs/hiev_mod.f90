module hiev_mod

  use :: prcn_mod

  implicit none

!  public :: hicalc2,evcalc2,evupdate2,hi_init,ev_init,intrncalc
  public :: hi_init,ev_init,evbw_init,intrncalc

  private

  type :: hi
    ! For RPY
    real(wp) :: A,B,C,D,E,F
    ! For OB
    real(wp) :: G
    ! For Regularized OB
    real(wp) :: O,P,R,S,T
    real(wp) :: rmagmin
  end type

  type :: evbb
    ! For Gaussian
    real(wp) :: prefactor,denom
    ! For LJ
    real(wp) :: epsOVsig,sigOVrtr,sigOVrtrto6
    real(wp) :: LJ_prefactor_tr
    real(wp) :: rmagmin
  end type

  type :: evbw
    ! For Cubic
    real(wp) :: delw
    real(wp) :: prf
    real(wp) :: rmagmin
  end type

  type(hi),save :: hi_prm
  type(evbb),save :: evbb_prm
  type(evbw),save :: evbw_prm
 
  type :: dis
    real(wp) :: x,y,z
    real(wp) :: mag,mag2
  end type

   
  interface

!    module subroutine hicalc2(rvmrc,nseg,DiffTens,Fev)
!      integer,intent(in) :: nseg
!      real(wp),intent(in) :: rvmrc(:)
!      real(wp),intent(out) :: DiffTens(:,:),Fev(:)
!    end subroutine hicalc2
!
!    module subroutine evcalc2(rvmrc,nseg,Fev)
!      integer,intent(in) :: nseg
!      real(wp),intent(in) :: rvmrc(:)
!      real(wp),intent(out) :: Fev(:)
!    end subroutine evcalc2
!
!    module subroutine evupdate2(Fev,rvmrc,nseg,Fbarev)
!      integer,intent(in) :: nseg
!      real(wp),intent(in) :: Fev(:),rvmrc(:)
!      real(wp),intent(out) :: Fbarev(:)
!    end subroutine evupdate2

    module subroutine hi_init()
    end subroutine hi_init

    module subroutine ev_init()
    end subroutine ev_init

    module subroutine evbw_init()
    end subroutine evbw_init

    module subroutine hicalc3(i,j,rij,DiffTens)
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: DiffTens(:,:)
    end subroutine hicalc3

    module subroutine evcalc3(i,j,rij,Fev)
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: Fev(:)
    end subroutine evcalc3

    !> Calculates EV force between particles and the wall
    !! \param ry the vertical distance of particle with the wall
    !! \param Fev EV force
    module subroutine evbwcalc(i,ry,Fev)
      integer,intent(in) :: i
      real(wp),intent(in) :: ry
      real(wp),intent(inout) :: Fev(:)
    end subroutine evbwcalc

  end interface

contains

!  subroutine hi_init()
!
!    use :: inp_dlt, only: HITens,hstar
!
!    real(wp),parameter :: PI=3.1415926535897958648_wp
!    real(wp),parameter :: sqrtPI=sqrt(PI)
!
!    if (HITens == 'RPY') then
!      ! For Rotne-Prager-Yamakawa Tensor:
!      hi_prm%A=0.75*hstar*sqrtPI
!      hi_prm%B=hstar**3*PI*sqrtPI/2
!      hi_prm%C=(3.0_wp/2)*hstar**3*PI*sqrtPI
!      hi_prm%D=2*sqrtPI*hstar
!      hi_prm%E=9/(32*sqrtPI*hstar)
!      hi_prm%F=3/(32*sqrtPI*hstar)
!    elseif (HITens == 'OB') then
!      ! For Oseen-Burgers Tensor:
!      hi_prm%G=0.75*hstar*sqrtPI
!    elseif (HITens == 'RegOB') then
!      ! For Reguralized Oseen-Burgers Tensor (introduced in HCO book):
!      hi_prm%G=0.75*hstar*sqrtPI
!      hi_prm%O=4*PI*hstar**2/3
!      hi_prm%P=14*PI*hstar**2/3
!      hi_prm%R=8*PI**2*hstar**4
!      hi_prm%S=2*PI*hstar**2
!      hi_prm%T=hi_prm%R/3
!    end if
!    hi_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|
!
!  end subroutine hi_init
!
!  subroutine ev_init()
!    
!    use :: inp_dlt, only: EVForceLaw,zstar,dstar
!    use :: inp_dlt, only: LJ_eps,LJ_sig,LJ_rtr,LJ_rc,minNonBond
!
!    if (EVForceLaw == 'Gauss') then
!      ev_prm%prefactor=zstar/dstar**5
!      ev_prm%denom=2*dstar**2
!    elseif (EVForceLaw == 'LJ') then
!      ev_prm%epsOVsig=LJ_eps/LJ_sig
!      ev_prm%sigOVrtr=LJ_sig/LJ_rtr
!      ev_prm%sigOVrtrto6=ev_prm%sigOVrtr*ev_prm%sigOVrtr*ev_prm%sigOVrtr* &
!                         ev_prm%sigOVrtr*ev_prm%sigOVrtr*ev_prm%sigOVrtr
!      ev_prm%LJ_prefactor_tr=4*ev_prm%epsOVsig*                           &
!                             ( 12*(ev_prm%sigOVrtrto6*ev_prm%sigOVrtrto6* &
!                                   ev_prm%sigOVrtr) -                     &
!                                6*(ev_prm%sigOVrtrto6*ev_prm%sigOVrtr))
!    end if
!    ev_prm%rmagmin=1.e-7_wp ! The Minimum value accepted as the |rij|
!  end subroutine ev_init


  subroutine intrncalc(rvmrc,rcm,nseg,DiffTens,Fev,Fbarev,calchi,&
                        calcev,calcevbw,updtev,updtevbw)

    use :: inp_dlt, only: EV_bb,EV_bw

    integer,intent(in) :: nseg
    real(wp),intent(in) :: rvmrc(:)
    real(wp),intent(in) :: rcm(:)
    type(dis) :: rij
    integer :: ibead,jbead,os,osi,osj
    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
    real(wp) :: ry,rimrc(3),rjmrc(3)
    real(wp),allocatable :: Fstarev(:)
    real(wp) :: DiffTens(:,:),Fev(:),Fbarev(:)
    logical :: clhi,clev,clevbw,upev,upevbw
    logical,optional :: calchi,calcev,calcevbw,updtev,updtevbw

!print *,'here'
    if (present(calchi)) then
      clhi=calchi
    else
      clhi=.false.
    end if
    if ((present(calcev)).and.(EV_bb/='NoEV')) then
      clev=calcev
    else
      clev=.false.
    end if
    if ((present(calcevbw)).and.(EV_bw/='NoEV')) then
      clevbw=calcevbw
    else
      clevbw=.false.
    end if
    if ((present(updtev)).and.(EV_bb/='NoEV')) then
      upev=updtev
    else
      upev=.false.
    end if
    if ((present(updtevbw)).and.(EV_bw/='NoEV')) then
      upevbw=updtevbw
    else
      upevbw=.false.
    end if

    if (clev.or.clevbw) Fev=0._wp
    if ((upev).or.(upevbw)) then
      allocate(Fstarev(3*(nseg+1)))
      Fstarev=0._wp
    end if

    do jbead=1,nseg+1
      os=3*(jbead-2)
      osj=3*(jbead-1)
      rjmrc=rvmrc(osj+1:osj+3)
      do ibead=1, jbead
        osi=3*(ibead-1)
        if (ibead == jbead) then
          if (clhi) then
            DiffTens(osi+1,osj+1)=1._wp
            DiffTens(osi+1,osj+2)=0._wp
            DiffTens(osi+1,osj+3)=0._wp
            DiffTens(osi+2,osj+2)=1._wp
            DiffTens(osi+2,osj+3)=0._wp
            DiffTens(osi+3,osj+3)=1._wp
          end if
        else
          rimrc=rvmrc(osi+1:osi+3)
          rij%x=rjmrc(1)-rimrc(1)
          rij%y=rjmrc(2)-rimrc(2)
          rij%z=rjmrc(3)-rimrc(3)
          rij%mag2=rij%x**2+rij%y**2+rij%z**2
          rij%mag=sqrt(rij%mag2)

          if (clhi) call hicalc3(ibead,jbead,rij,DiffTens)
          if (clev) call evcalc3(ibead,jbead,rij,Fev)
          if (upev) call evcalc3(ibead,jbead,rij,Fstarev)

        end if ! ibead /= jbead
      end do ! ibead

      if (clevbw) then
        ry=rjmrc(2)+rcm(2)
        call evbwcalc(jbead,ry,Fev)
      end if
      if (upevbw) then
        ry=rjmrc(2)+rcm(2)
        call evbwcalc(jbead,ry,Fstarev)
      end if

    end do ! jbead

    if ((upev).or.(upevbw)) then
      Fbarev=0.5*(Fev+Fstarev)
      deallocate(Fstarev)
    end if

  end subroutine intrncalc
 
end module hiev_mod
