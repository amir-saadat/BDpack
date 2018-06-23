module intrn_mod

  use :: prcn_mod

  implicit none

  public :: intrn_t,&
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
    real(wp) :: a_sph
  end type hibw_t
  type :: hi_t
    type(hibb_t) :: hibb
    type(hibw_t) :: hibw
  end type hi_t
  !------------------------------------------

  ! ev
  !------------------------------------------
  type :: evbb_t
    ! For Gaussian
    real(wp) :: prefactor,denom
    ! For LJ
    real(wp) :: epsOVsig,sigOVrtr,sigOVrtrto6
    real(wp) :: LJ_prefactor_tr
    real(wp) :: rmagmin
  end type evbb_t
  type :: evbw_t
    ! For Cubic
    real(wp) :: delw
    real(wp) :: prf
    real(wp) :: rmagmin
    ! For Reflc-bc
    real(wp) :: a
    real(wp) :: a_sph
    integer :: iwall !the type of wall for the reflection BC: 1-plane, 2-sphere
    integer :: u_wc
    integer :: u_wc_all
    integer :: u_ia
    integer,allocatable :: w_coll(:,:)
    integer,allocatable :: w_coll_all(:,:)
    integer,allocatable :: ia_time(:,:,:)
    integer,allocatable :: w_coll_t(:,:)
    integer,allocatable :: w_coll_all_t(:,:)
    integer,allocatable :: ia_time_t(:,:,:)
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

  type :: dis
    real(wp) :: x,y,z
    real(wp) :: mag,mag2
    real(wp) :: riy,rjy,yim
    real(wp) :: magim,mag2im
    real(wp) :: rjx,rjz
    real(wp) :: rix,riz
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
    module subroutine calc_div(j,rjy,divD,id,itime)
      integer,intent(in) :: j
      real(wp),intent(in) :: rjy
      real(wp),intent(inout) :: divD(:)
      integer,intent(in) :: id,itime
    end subroutine calc_div
    !------------------------------------------

    ! evbb
    !------------------------------------------
    module subroutine init_evbb(this)
      class(evbb_t),intent(inout) :: this
    end subroutine init_evbb
    module subroutine calc_evbb(this,i,j,rij,Fev)
      class(evbb_t),intent(inout) :: this
      integer,intent(in) :: i,j
      type(dis),intent(in) :: rij
      real(wp),intent(inout) :: Fev(:)
    end subroutine calc_evbb
    !------------------------------------------

    !evbw
    !------------------------------------------
    ! module subroutine evbw_init(id)
    !   integer,intent(in) :: id
    ! end subroutine evbw_init
    module subroutine init_evbw(this,id)
      class(evbw_t),intent(inout) :: this
      integer,intent(in) :: id
    end subroutine init_evbw
    !> calculates EV force between particles and the wall
    !! \param ry the vertical distance of particle & the wall
    !! \param Fev EV force
    module subroutine calc_evbw(this,i,ry,Fev)
      class(evbw_t),intent(inout) :: this
      integer,intent(in) :: i
      real(wp),intent(in) :: ry
      real(wp),intent(inout) :: Fev(:)
    end subroutine calc_evbw
    module subroutine wall_rflc(this,dt,it,time,id,ich,qx,qy,qz,Rx,Ry,Rz,&
      rcmx,rcmy,rcmz,rf0,r_sph)
      class(evbw_t),intent(inout) :: this
      real(wp),intent(in),dimension(3) :: rf0,r_sph!rf_in
      !real(wp),intent(in) :: rf_in(:,:)
      real(wp),intent(in) :: dt,time
      integer,intent(in) :: it,id,ich
      real(wp),intent(inout) :: qx(:),qy(:),qz(:)
      real(wp),intent(inout) :: Rx(:),Ry(:),Rz(:)
      real(wp),intent(inout) :: rcmx,rcmy,rcmz
    end subroutine wall_rflc
    module subroutine del_evbw(this,id)
      class(evbw_t),intent(inout) :: this
      integer,intent(in) :: id
    end subroutine del_evbw
    module subroutine print_wcll(this,id,nproc,MPI_REAL_WP,time)
      class(evbw_t),intent(inout) :: this
      integer,intent(in) :: id,nproc
      integer,intent(in) :: MPI_REAL_WP
      real(wp),intent(in) :: time
    end subroutine print_wcll
    !------------------------------------------

  end interface

contains

  subroutine init_intrn(this,id)

    class(intrn_t),intent(inout) :: this
    integer,intent(in) :: id

    call init_hi(this%hi)
    call init_evbb(this%evbb)
    call init_evbw(this%evbw,id)

  end subroutine init_intrn

  subroutine calc_intrn(this,id,itime,rvmrc,rcm,r_sph,nseg,DiffTens,divD,Fev,Fbarev,&
                      calchi,calcdiv,calcevbb,calcevbw,updtevbb,updtevbw)

    use :: inp_dlt, only: EV_bb,EV_bw,hstar,HITens,nbead,nbead_ind
    use :: arry_mod, only: print_vector,print_matrix

    class(intrn_t),intent(inout) :: this
    integer,intent(in) :: nseg
    real(wp),intent(in) :: rvmrc(:)
    real(wp),intent(in) :: rcm(:,:)
    real(wp),intent(in) :: r_sph(:)
    type(dis) :: rij
    integer :: ibead,jbead,os,osi,osj,ibead_ulim
    real(wp) :: LJ_prefactor,LJ_prefactor_tr,epsOVsig
    real(wp) :: sigOVrtr,sigOVrtrto6,sigOVrmag,sigOVrmagto6
    real(wp) :: ry,rimrc(3),rjmrc(3),rjy
    real(wp),allocatable :: Fstarev(:)
    real(wp) :: DiffTens(:,:),divD(:),Fev(:),Fbarev(:)
    logical :: clhi,cldiv,clevbb,clevbw,upevbb,upevbw
    logical,optional :: calchi,calcdiv,calcevbb,calcevbw,updtevbb,updtevbw
    integer :: ichain_pp,jchain_pp

    integer :: id,itime


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
    if ((present(calcevbb)).and.(EV_bb/='NoEV')) then
      clevbb=calcevbb
    else
      clevbb=.false.
    end if
    if ((present(calcevbw)).and.(EV_bw/='NoEV')) then
      clevbw=calcevbw
    else
      clevbw=.false.
    end if
    if ((present(updtevbb)).and.(EV_bb/='NoEV')) then
      upevbb=updtevbb
    else
      upevbb=.false.
    end if
    if ((present(updtevbw)).and.(EV_bw/='NoEV')) then
      upevbw=updtevbw
    else
      upevbw=.false.
    end if

    if (clevbb.or.clevbw) Fev=0._wp
    if ((upevbb).or.(upevbw)) then
      allocate(Fstarev(3*(nbead)))
      Fstarev=0._wp
    end if

    do jbead=1,nbead!nseg+1
      os=3*(jbead-2) !why is this here?
      osj=3*(jbead-1)
      rjmrc=rvmrc(osj+1:osj+3)
      jchain_pp = (jbead-1) / nbead_ind + 1

      !print *, 'jbead = ', jbead
      !print *, 'jchain_pp = ', jchain_pp

      ! if ( (id==1 .and. itime==10000) .or. &
      !      (id==1 .and. itime==20000) ) then
      !   if (jbead<=3) then
      !     print*,'rj-intrn',rjmrc+rcm
      !   endif
      ! endif

      !! Blake's part
      if (hstar /= 0._wp .and. HITens == 'Blake') then
        ibead_ulim=nbead!nseg+1
        rjy=rjmrc(2)+rcm(2,jchain_pp)
        if (cldiv) call calc_div(jbead,rjy,divD,id,itime)
      else
        ibead_ulim=jbead !upper triangular
      endif
      !!-------------



      do ibead=1, ibead_ulim

        osi=3*(ibead-1)
        ichain_pp = (ibead-1) / nbead_ind + 1

        ! if (ibead == jbead) then
        !   if (clhi) then
        !     DiffTens(osi+1,osj+1)=1._wp
        !     DiffTens(osi+1,osj+2)=0._wp
        !     DiffTens(osi+1,osj+3)=0._wp
        !     DiffTens(osi+2,osj+2)=1._wp
        !     DiffTens(osi+2,osj+3)=0._wp
        !     DiffTens(osi+3,osj+3)=1._wp
        !   end if
        ! else
        !   rimrc=rvmrc(osi+1:osi+3)
        !   ! rij%x=rjmrc(1)-rimrc(1)
        !   ! rij%y=rjmrc(2)-rimrc(2)
        !   ! rij%z=rjmrc(3)-rimrc(3)
        !   rij%x=rimrc(1)-rjmrc(1)
        !   rij%y=rimrc(2)-rjmrc(2)
        !   rij%z=rimrc(3)-rjmrc(3)
        !   rij%mag2=rij%x**2+rij%y**2+rij%z**2
        !   rij%mag=sqrt(rij%mag2)

        !   if (HITens == 'Blake') then
        !     rij%ry=rjmrc(2)+rcm(2)
        !     rij%yim=rimrc(2)+rcm(2)+rij%ry
        !     rij%mag2im=rij%x**2+rij%yim**2+rij%z**2
        !     rij%magim=sqrt(rij%mag2im)
        !   endif
        !   if (clhi) call calc_hi(this%hi,ibead,jbead,rij,DiffTens)
        !   if (clev) call evcalc3(ibead,jbead,rij,Fev)
        !   if (upev) call evcalc3(ibead,jbead,rij,Fstarev)
        ! end if ! ibead /= jbead



        ! if (ibead /= jbead) then


          ! if (clhi) then
          !   DiffTens(osi+1,osj+1)=1._wp
          !   DiffTens(osi+1,osj+2)=0._wp
          !   DiffTens(osi+1,osj+3)=0._wp
          !   DiffTens(osi+2,osj+2)=1._wp
          !   DiffTens(osi+2,osj+3)=0._wp
          !   DiffTens(osi+3,osj+3)=1._wp
          ! end if


        ! else

          rimrc=rvmrc(osi+1:osi+3)
          ! rij%x=rjmrc(1)-rimrc(1)
          ! rij%y=rjmrc(2)-rimrc(2)
          ! rij%z=rjmrc(3)-rimrc(3)
          rij%x=(rimrc(1)-rjmrc(1)) + (rcm(1,ichain_pp)-rcm(1,jchain_pp))
          rij%y=(rimrc(2)-rjmrc(2)) + (rcm(2,ichain_pp)-rcm(2,jchain_pp))
          rij%z=(rimrc(3)-rjmrc(3)) + (rcm(3,ichain_pp)-rcm(3,jchain_pp))
          rij%mag2=rij%x**2+rij%y**2+rij%z**2
          rij%mag=sqrt(rij%mag2)


        if (ibead /= jbead) then

          if (clevbb) call calc_evbb(this%evbb,ibead,jbead,rij,Fev)
          if (upevbb) call calc_evbb(this%evbb,ibead,jbead,rij,Fstarev)

        end if ! ibead /= jbead


        ! Blake's part
        if (HITens == 'Blake') then
          rij%rjy=rjmrc(2)+rcm(2,jchain_pp)
          rij%riy=rimrc(2)+rcm(2,ichain_pp)
          rij%yim=rimrc(2)+rcm(2,ichain_pp)+rij%rjy
          rij%mag2im=rij%x**2+rij%yim**2+rij%z**2
          rij%magim=sqrt(rij%mag2im)
        endif

        if (HITens == 'Osph') then
          rij%rjx=rjmrc(1)+rcm(1,jchain_pp) - r_sph(1)
          rij%rjy=rjmrc(2)+rcm(2,jchain_pp) - r_sph(2)
          rij%rjz=rjmrc(3)+rcm(3,jchain_pp) - r_sph(3)
          rij%rix=rimrc(1)+rcm(1,ichain_pp) - r_sph(1)
          rij%riy=rimrc(2)+rcm(2,ichain_pp) - r_sph(2)
          rij%riz=rimrc(3)+rcm(3,ichain_pp) - r_sph(3)
        endif
        !-------------

        if (clhi) call calc_hi(this%hi,ibead,jbead,rij,DiffTens)

      end do ! ibead

      if (clevbw) then
        ry=rjmrc(2)+rcm(2,jchain_pp)
        call calc_evbw(this%evbw,jbead,ry,Fev)
      end if
      if (upevbw) then
        ry=rjmrc(2)+rcm(2,jchain_pp)
        call calc_evbw(this%evbw,jbead,ry,Fstarev)
      end if

    end do ! jbead

    !call print_matrix(DiffTens(:,:),'DiffTens')
    if ((upevbb).or.(upevbw)) then
      Fbarev=0.5*(Fev+Fstarev)
      deallocate(Fstarev)
    end if

  end subroutine calc_intrn

end module intrn_mod
