module intrn_mod

  use :: prcn_mod

  implicit none

  public :: intrn_t,&
            hi_init,&
            ev_init,&
            evbw_init,&
            wall_rflc,&
            print_wcll,&
            del_evbw

  private


  !------------------------------------
  !>>>> Definition for different types
  !------------------------------------

  ! hi
  !------------------------------------------
  type :: hibb_t
    ! RPY
    real(wp) :: A,B,C,D,E,F
    ! Oseen-Burgers
    real(wp) :: G
    ! regularized OB
    real(wp) :: O,P,R,S,T
    real(wp) :: rmagmin
  end type hibb_t
  type :: hibw_t
  end type hibw_t
  type :: hi_t
    type(hibb_t) :: hibb
    type(hibw_t) :: hibw
  end type hi_t
  !------------------------------------------

  type :: hi
    ! RPY
    real(wp) :: A,B,C,D,E,F
    ! Oseen-Burgers
    real(wp) :: G
    ! regularized OB
    real(wp) :: O,P,R,S,T
    real(wp) :: rmagmin
  end type

  ! ev
  !------------------------------------------
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
    ! For Reflc-bc
    real(wp) :: a
    integer :: u0
    integer,allocatable :: w_coll(:)
    integer,allocatable :: ia_time(:,:)
    integer,allocatable :: w_coll_t(:)    
  end type
  type :: evbb_t
  end type evbb_t
  type :: evbw_t
  end type evbw_t
  !------------------------------------------

  ! intrn
  !------------------------------------------
  type :: intrn_t
    type(hi_t) :: hi
    type(evbb_t) :: evbb
    type(evbw_t) :: evbw
  contains
    procedure,pass(this) :: init => init_intrn
    procedure,pass(this) :: calc => calc_intrn
  end type intrn_t
  !------------------------------------------

  type(hi),save :: hi_prm
  type(evbb),save :: evbb_prm
  type(evbw),save :: evbw_prm
 
  type :: dis
    real(wp) :: x,y,z
    real(wp) :: mag,mag2
    real(wp) :: ry,yim
    real(wp) :: magim,mag2im
  end type

  !-----------------------------------------------------
  !>>>> Interface to routines for different interactions
  !-----------------------------------------------------
  interface
    ! hi
    !------------------------------------------

    !> initializes the hi type
    module subroutine init_hi(this)
      class(hi_t),intent(inout) :: this
    end subroutine init_hi
    !> initializes the hi type
    module subroutine hi_init()
    end subroutine hi_init
    !> calculates HI
    !! \param this hi object
    !! \param i bead i index
    !! \param j bead j index
    !! \param rij data type for inter particle distance
    !! \param DiffTens diffusion tensor
    module subroutine calc_hi(this,i,j,rij,DiffTens)
      class(hi_t),intent(inout) :: this
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: DiffTens(:,:)
    end subroutine calc_hi
    !> calculates divergance of D
    !! \param j bead j index
    !! \param rjy y component of rj
    !! \param divD divergance of D
    module subroutine calc_div(j,rjy,divD)
      integer,intent(in) :: j
      real(wp),intent(in) :: rjy
      real(wp),intent(inout) :: divD(:)
    end subroutine calc_div
    !------------------------------------------

    ! ev
    !------------------------------------------
    module subroutine init_evbb(this)
      class(evbb_t),intent(inout) :: this
    end subroutine init_evbb

    module subroutine init_evbw(this)
      class(evbw_t),intent(inout) :: this
    end subroutine init_evbw

    module subroutine ev_init()
    end subroutine ev_init

    module subroutine evbw_init(id)
      integer,intent(in) :: id
    end subroutine evbw_init

    module subroutine evcalc3(i,j,rij,Fev)
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: Fev(:)
    end subroutine evcalc3

    !> calculates EV force between particles and the wall
    !! \param ry the vertical distance of particle & the wall
    !! \param Fev EV force
    module subroutine evbwcalc(i,ry,Fev)
      integer,intent(in) :: i
      real(wp),intent(in) :: ry
      real(wp),intent(inout) :: Fev(:)
    end subroutine evbwcalc

    module subroutine wall_rflc(Ry,rcmy)
      real(wp),intent(inout) :: Ry(:),rcmy
    end subroutine wall_rflc

    module subroutine del_evbw(id)
      integer,intent(in) :: id
    end subroutine del_evbw

    module subroutine print_wcll(id,MPI_REAL_WP,time)
      integer,intent(in) :: id
      integer,intent(in) :: MPI_REAL_WP
      real(wp),intent(in) :: time
    end subroutine print_wcll
    !------------------------------------------

  end interface

contains

  subroutine init_intrn(this)

    class(intrn_t),intent(inout) :: this

    call init_hi(this%hi)
    call init_evbb(this%evbb)
    call init_evbw(this%evbw)

  end subroutine init_intrn

  subroutine calc_intrn(this,rvmrc,rcm,nseg,DiffTens,divD,Fev,Fbarev,&
                      calchi,calcdiv,calcev,calcevbw,updtev,updtevbw)

    use :: inp_dlt, only: EV_bb,EV_bw,hstar,HITens

    class(intrn_t),intent(inout) :: this
    integer,intent(in) :: nseg
    real(wp),intent(in) :: rvmrc(:)
    real(wp),intent(in) :: rcm(:)
    type(dis) :: rij
    integer :: ibead,jbead,os,osi,osj,ibead_ulim
    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
    real(wp) :: ry,rimrc(3),rjmrc(3),rjy
    real(wp),allocatable :: Fstarev(:)
    real(wp) :: DiffTens(:,:),divD(:),Fev(:),Fbarev(:)
    logical :: clhi,cldiv,clev,clevbw,upev,upevbw
    logical,optional :: calchi,calcdiv,calcev,calcevbw,updtev,updtevbw

!print *,'here'
    if (present(calchi)) then
      clhi=calchi
    else
      clhi=.false.
    end if
    if (present(calcdiv)) then
      cldiv=calcdiv
    else
      cldiv=.false.
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



      !! Blake's part
      if (hstar /= 0._wp .and. HITens == 'Blake') then
        ibead_ulim=nseg+1
        rjy=rjmrc(2)+rcm(2)
        if (cldiv) call calc_div(jbead,rjy,divD)
      else
        ibead_ulim=jbead
      endif
      !!-------------


      do ibead=1, ibead_ulim

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
          ! rij%x=rjmrc(1)-rimrc(1)
          ! rij%y=rjmrc(2)-rimrc(2)
          ! rij%z=rjmrc(3)-rimrc(3)
          rij%x=rimrc(1)-rjmrc(1)
          rij%y=rimrc(2)-rjmrc(2)
          rij%z=rimrc(3)-rjmrc(3)
          rij%mag2=rij%x**2+rij%y**2+rij%z**2
          rij%mag=sqrt(rij%mag2)

          if (HITens == 'Blake') then
            rij%ry=rjmrc(2)+rcm(2)
            rij%yim=rimrc(2)+rcm(2)+rij%ry
            rij%mag2im=rij%x**2+rij%yim**2+rij%z**2
            rij%magim=sqrt(rij%mag2im)
          endif

          if (clhi) call calc_hi(this%hi,ibead,jbead,rij,DiffTens)
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

  end subroutine calc_intrn
 
end module intrn_mod
