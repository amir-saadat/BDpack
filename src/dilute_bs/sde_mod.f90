!sde module
!Tiras Lin
!Dec 17, 2017
!!!!!!!!!!!!!!!!!!!!!!!!

module sde_mod

  use :: prcn_mod

  implicit none

  public :: sde_t
  private

  !type containing the SDE
  type :: sde_t
    real(wp) :: xxkappa,xykappa,yykappa,zzkappa
    real(wp) :: a_sph,U_mag_sph
    real(wp),dimension(3,3) :: Kappareg
    real(wp),allocatable,dimension(:,:) :: Kappa
    real(wp),allocatable,dimension(:,:) :: Amat, Bmat
    real(wp),allocatable,dimension(:,:) :: KappaBF, AmatBF
  contains
    procedure,pass(this) :: init => init_sde
    procedure,pass(this) :: advance => advance_sde
    procedure,pass(this) :: U_sph => U_sph_sde
  end type sde_t

contains


  !define parameters needed for the SDE
  subroutine init_sde(this,Kappareg,Kappa,Amat,Bmat,KappaBF,AmatBF,nbead_bb,nseg_bb)

    !variables used from other places
    use :: inp_dlt, only: iflow,nseg,nbead,nsegx3,nbeadx3,tplgy,Na,nseg_ar,Ia,sph_flow
    use :: arry_mod, only: print_vector, print_matrix
    use :: cmn_io_mod, only: read_input

    !input and output arguments
    class(sde_t),intent(inout) :: this
    real(wp), intent(inout) :: Kappareg(:,:),Kappa(:,:),Amat(:,:),Bmat(:,:),KappaBF(:,:),AmatBF(:,:)
    integer, intent(inout) :: nbead_bb,nseg_bb

    !variables used inside init_sde
    integer :: iseg,jseg,offseti,offsetj,ibead,jbead,i,j,k,iarm
    integer :: ku,kl
    integer :: idx,nu,mu
    real(wp) :: fctr

    !allocations
    allocate(this%Kappa(nsegx3,nsegx3))
    allocate(this%Amat(nsegx3,nbeadx3),this%Bmat(nbeadx3,nsegx3))
    allocate(this%KappaBF(2,nsegx3))
    select case (tplgy)
    case ('Linear')
      allocate(this%AmatBF(4,nbeadx3))
    case ('Comb')
    end select

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !initialization code begins here:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Specifying kappa based on type of flow
    select case (iflow)
    case (1) ! Finding Equilibrium
    !if (iflow == 1) then
      this%xxkappa=0._wp
      this%xykappa=0._wp
      this%yykappa=0._wp
      this%zzkappa=0._wp
    case (2) ! Shear Flow
    !else if (iflow == 2) then ! Shear Flow
      this%xxkappa=0._wp
      this%xykappa=1._wp
      this%yykappa=0._wp
      this%zzkappa=0._wp
    case (3) ! Uniaxial Extension
    !else if (iflow == 3) then ! Uniaxial Extension
      this%xxkappa=1._wp
      this%xykappa=0._wp
      this%yykappa=-0.5_wp
      this%zzkappa=-0.5_wp
    case (4) !Biaxial Extension
    !else if (iflow == 4) then ! biaxial Extension
      this%xxkappa=1._wp
      this%xykappa=0._wp
      this%yykappa=1._wp
      this%zzkappa=-2.0_wp
    case (5) !Planar Extension
    !else if (iflow == 5) then ! Planar Extension
      this%xxkappa=1._wp
      this%xykappa=0._wp
      this%yykappa=-1._wp
      this%zzkappa=0._wp
    end select
    !end if

    this%Kappa=0._wp
    forall (iseg=1:3*(nseg-1)+1:3)
      this%Kappa(iseg,iseg)=this%xxkappa
      this%Kappa(iseg,iseg+1)=this%xykappa
      this%Kappa(iseg+1,iseg+1)=this%yykappa
      this%Kappa(iseg+2,iseg+2)=this%zzkappa
    end forall
    this%Kappareg(1,1:3)=(/this%xxkappa,this%xykappa,0.0_wp/)
    this%Kappareg(2,1:3)=(/0.0_wp ,this%yykappa,0.0_wp/)
    this%Kappareg(3,1:3)=(/0.0_wp ,0.0_wp ,this%zzkappa/)

    !if there is flow past a sphere
    if (sph_flow) then
      call read_input('Sph-rad',0,this%a_sph)
      call read_input('Umag-Sph',0,this%U_mag_sph)
    endif

    this%Amat=0.0_wp
    select case (tplgy)
    case ('Linear')
      nseg_bb=nseg
      nbead_bb=nbead
      do iseg=1, nseg
        offseti=3*(iseg-1)
        do jbead=1, nbead
          offsetj=3*(jbead-1)
          if (iseg == jbead) then
            forall (i=1:3) this%Amat(offseti+i,offsetj+i)=-1._wp
          elseif (iseg == jbead-1) then
            forall (i=1:3) this%Amat(offseti+i,offsetj+i)= 1._wp
          end if
        end do
      end do
    case ('Comb')
      nseg_bb=nseg-Na*nseg_ar
      nbead_bb=nseg_bb+1
      iarm=1
      do iseg=1, nseg
        offseti=3*(iseg-1)
        do jbead=1, nbead
          offsetj=3*(jbead-1)
          if (iseg <= nseg_bb) then
            if (jbead == iseg) then
              forall (i=1:3) this%Amat(offseti+i,offsetj+i)=-1._wp
            elseif (jbead == iseg+1) then
              forall (i=1:3) this%Amat(offseti+i,offsetj+i)= 1._wp
            end if
          else ! iseg > nseg_bb
            if (iseg-nseg_bb-(iarm-1)*nseg_ar == 1) then
              if (jbead == Ia(iarm+1)) then
                forall (i=1:3) this%Amat(offseti+i,offsetj+i)=-1._wp
              elseif (jbead == iseg+1) then
                forall (i=1:3) this%Amat(offseti+i,offsetj+i)= 1._wp
              elseif (jbead == nbead) then
                iarm=iarm+1
              end if
            else
              if (jbead == iseg) then
                forall (i=1:3) this%Amat(offseti+i,offsetj+i)=-1._wp
              elseif (jbead == iseg+1) then
                forall (i=1:3) this%Amat(offseti+i,offsetj+i)= 1._wp
              end if
            end if
          end if
        end do
      end do
    end select

    ! Constructing banded form of Kappa
    this%KappaBF=0._wp
    ku=1;kl=0
    do j=1, nsegx3
      k=ku+1-j
      do i=max(1,j-ku),min(nsegx3,j+kl)
        this%KappaBF(k+i,j)=this%Kappa(i,j)
      end do
    end do
    if (tplgy == 'Linear') then
      this%AmatBF=0._wp
      ! Constructing banded form of Amat
      ku=3;kl=0
      do j=1, nbeadx3
        k=ku+1-j
        do i=max(1,j-ku),min(nsegx3,j+kl)
          this%AmatBF(k+i,j)=this%Amat(i,j)
        end do
      end do
    end if
    !line 497 Constructing Bmat
    this%Bmat=0.0_wp
    select case (tplgy)
    case ('Linear')
      do ibead=1, nbead
        offseti=3*(ibead-1)
        do jseg=1, nseg
          offsetj=3*(jseg-1)
          if (ibead > jseg) then
            forall (i=1:3)
              this%Bmat(offseti+i,offsetj+i)=jseg/real(nbead,kind=wp)
            end forall
          else
            forall (i=1:3)
              this%Bmat(offseti+i,offsetj+i)=-(1-jseg/real(nbead,kind=wp))
            end forall
          end if
        end do
      end do
    case ('Comb')
     ! Constructing the elements of the first row of B
     do k=1, nseg_bb
       forall (i=1:3) this%Bmat(i,3*(k-1)+i)=-(nseg_bb-k+1)/real(nbead,kind=wp)
     end do
     do iarm=1, Na
       fctr=(Na-iarm+1)*nseg_ar/real(nbead,kind=wp)
       do k=Ia(iarm), Ia(iarm+1)-1
         forall (i=1:3)
           this%Bmat(i,3*(k-1)+i)=this%Bmat(i,3*(k-1)+i)-fctr
         end forall
       end do ! k
       do k=1, nseg_ar
         idx=nseg_bb+(iarm-1)*nseg_ar+k
         forall (i=1:3)
           this%Bmat(i,3*(idx-1)+i)=this%Bmat(i,3*(idx-1)+i)-&
           (nseg_ar-k+1)/real(nbead,kind=wp)
         end forall
       end do ! k
     end do ! iarm
     ! Constructing the rest of the rows in backbone
     do nu=2, nseg_bb+1
       forall (i=1:3) this%Bmat(3*(nu-1)+i,:)=this%Bmat(i,:)
       do k=1, nu-1
         forall (i=1:3)
           this%Bmat(3*(nu-1)+i,3*(k-1)+i)=this%Bmat(3*(nu-1)+i,3*(k-1)+i)+1
         end forall
       end do ! k
     end do ! nu
     ! Constructing the rows for the arms
     do iarm=1, Na
       do mu=1, nseg_ar
         nu=nseg_bb+1+(iarm-1)*nseg_ar+mu
         forall (i=1:3) this%Bmat(3*(nu-1)+i,:)=this%Bmat(i,:)
         do k=1, Ia(iarm+1)-1
           forall (i=1:3)
             this%Bmat(3*(nu-1)+i,3*(k-1)+i)=this%Bmat(3*(nu-1)+i,3*(k-1)+i)+1
           end forall
         end do ! k
         do k=1, mu
           idx=nseg_bb+(iarm-1)*nseg_ar+k
           forall (i=1:3)
             this%Bmat(3*(nu-1)+i,3*(idx-1)+i)=this%Bmat(3*(nu-1)+i,3*(idx-1)+i)+1
           end forall
         end do ! k
       end do ! mu
     end do ! iarm
   end select
   !line 563
   !call print_matrix(this%Kappa,'this%Kappa')


   !Kappareg,Kappa,Amat,Bmat,KappaBF,AmatBF
   Kappareg = this%Kappareg
   Kappa = this%Kappa
   Amat = this%Amat
   Bmat = this%Bmat
   KappaBF = this%KappaBF
   AmatBF = this%AmatBF

   ! call print_matrix(Kappareg,'Kappareg')
   ! call print_matrix(Kappa,'Kappa')
   ! call print_matrix(Amat,'Amat')
   ! call print_matrix(Bmat,'Bmat')
   ! call print_matrix(KappaBF,'KappaBF')
   ! call print_matrix(AmatBF,'AmatBF')

   !print *, 'this%Bmat', this%Bmat(1,1)
   !print *, '---------------------'
   !print *, 'Bmat',Bmat(1,1)

  end subroutine init_sde


  !update the configurations of the chains with the SDE
  subroutine advance_sde(this,myintrn,id,iPe,idt,ichain,itime,Kdotq,qc,Fseg,Fbead,Fev,Fbnd,qstar,Fphi,&
    rvmrcP,rcm,DiffTensP,Ftet,rf0,AdotDP1,divD,FBr,RHS,rcmP,Fbarev,nbead_bb,Fbarbnd,Fbar,Fbartet,&
    RHScnt,Fbarseg,Fbarbead,root_f,qbar,nseg_bb,AdotD,RHSbase,qctemp,mch,Lch,lambdaBE)

    !variables used from other places
    use :: inp_dlt, only: nseg,nbead,tplgy,dt,Pe,nsegx3,nbeadx3,applFext,Fext0,srf_tet,&
      hstar,HITens,EV_bb,EV_bw,ForceLaw,PrScale,nroots,TruncMethod,nseg_ar,tol,DecompMeth,Ia,Na,sph_flow
    use :: force_mod, only: tetforce,bndupdate,tetupdate,sprupdate
    use :: intrn_mod, only: intrn_t
    use :: arry_mod, only: print_vector, print_matrix

    !input and output arguments
    class(sde_t),intent(inout) :: this
    class(intrn_t),intent(inout) :: myintrn
    integer, intent(in) :: id
    integer, intent(in) :: iPe,idt,ichain,itime
    real(wp), intent(inout) :: Kdotq(:),Fbead(:),qstar(:)
    real(wp), intent(inout), target :: qc(:),Fseg(:)
    real(wp), intent(inout) :: Fphi(:,:),Fev(:),Fbnd(:)
    real(wp), intent(inout) :: rvmrcP(:)
    real(wp), intent(inout) :: rcm(:,:)
    real(wp), intent(inout) :: DiffTensP(:,:)
    real(wp), intent(inout) :: Ftet(:)
    real(wp), intent(inout) :: rf0(:,:)
    real(wp), intent(inout) :: AdotDP1(:,:)
    real(wp), intent(inout) :: divD(:)
    real(wp), intent(inout) :: FBr(:)
    real(wp), intent(inout), target :: RHS(:)
    real(wp), intent(inout) :: rcmP(:)
    real(wp), intent(inout) :: Fbarev(:)
    real(wp), intent(inout) :: Fbarbnd(:)
    integer, intent(inout) :: nbead_bb
    real(wp), intent(inout) :: Fbar(:)
    real(wp), intent(inout) :: Fbartet(:)
    real(wp), intent(inout) :: RHScnt(:)
    real(wp), intent(inout) :: Fbarseg(:),Fbarbead(:)
    real(wp), intent(inout) :: root_f(:)
    real(wp), intent(inout) :: qbar(:)
    real(wp), intent(inout), target :: AdotD(:,:,:)
    integer, intent(inout) :: nseg_bb
    real(wp), intent(inout),target :: RHSbase(:)
    real(wp), intent(inout) :: qctemp(:)
    integer,intent(inout) :: mch(:),Lch(:)
    real(wp), intent(inout) :: lambdaBE(:)

    !variables used inside advance_sde
    integer :: kl,is,os,offset,iseg,icount
    real(wp) :: eps
    real(wp),dimension(:),pointer  :: RHSP,RHSbaseP,qcP,FsegP
    real(wp),dimension(:,:),pointer :: AdotDP2
    real(wp),dimension(nsegx3) :: U_seg
    real(wp),dimension(nbeadx3) :: U_bead

    !debug
    !call print_vector(qc,'qc')

    !============ Predictor-Corrector =============!
    !--------Predictor Algorithm----------!
    ! Kdotq=dt*Pe*(Kappa.q)               !
    ! Fbead=-A'.Fseg                      !
    ! qstar=q                             !
    ! qstar:=qstar+Kdotq                  !
    ! Fphi=Fbead+Fev+Fbnd                 !
    ! qstar:=qstar+(1/4)*dt*(AdotD.Fbead) !
    ! qstar:=qstar+(1/4)*dt*(AdotD.Fev)   !
    ! qstar:=qstar+(1/4)*dt*(AdotD.Fbnd)  !
    ! qstar:=qstar+FBr                    !
    !-------------------------------------!
    call gbmv(this%KappaBF,qc,Kdotq,kl=0,alpha=Pe(iPe)*dt(iPe,idt))

    if (sph_flow) then
      call U_sph_sde(this,U_seg,U_bead,qc,rcm(:,ichain))
      Kdotq = Kdotq + U_seg*dt(iPe,idt)
    endif

    if (tplgy == 'Linear') then
      call gbmv(this%AmatBF,Fseg,Fbead,kl=0,m=nsegx3,alpha=-1.0_wp,trans='T')
    else
      call gemv(this%Amat,Fseg,Fbead,alpha=-1._wp,trans='T')
    end if
    call copy(qc,qstar)
    call axpy(Kdotq,qstar)
    Fphi(:,ichain)=Fbead+Fev+Fbnd
    if (applFext) then !continue here!
      Fphi(1,ichain)=Fphi(1,ichain)-Fext0
      Fphi(nbeadx3-2,ichain)=Fphi(nbeadx3-2,ichain)+Fext0
    end if
    if (srf_tet) then
      call tetforce(rvmrcP,rcm(:,ichain),DiffTensP,dt(iPe,idt),Ftet,&
        rf0(:,ichain),itime)
      Fphi(1:3,ichain)=Fphi(1:3,ichain)+Ftet(1:3)
    end if
    call gemv(AdotDP1,Fphi(:,ichain),qstar,alpha=0.25*dt(iPe,idt),&
      beta=1._wp)

    !TYL: HI for tethered bead -------------------------------------------
    !print *, 'Predictor calculaton------------'
    !call print_matrix(AdotDP1,'AdotDP1')
    !call print_vector(Fphi(:,ichain),'Fphi(:,ichain)')
    !call print_vector(qstar,'qstar BEFORE')
    if (srf_tet) then
      call gemv(AdotDP1(:,1:3),Fphi(1:3,ichain),qstar,&
        alpha=-0.25*dt(iPe,idt),beta=1._wp)
      !call print_vector(qstar,'qstar AFTER')
    end if
    !TYL: HI for tethered bead -------------------------------------------

    !! Blake's part
    if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
      do is=1, nseg
        os=(is-1)*3
        qstar(os+2)=qstar(os+2)+(divD(is+1)-divD(is))*0.25*dt(iPe,idt)
      enddo
    endif
    !!-------------

    call axpy(FBr,qstar) ! line 1068


    !-------First Corrector Algorithm-------!
    ! RHS=q                                 !
    ! RHS:=RHS+1/2*Kdotq (from Predictor)   !
    ! Fbarev:=Fev+Fstarev                   !
    ! RHS:=RHS+(1/4)dt*(AdotD.Fev)          !
    ! RHS:=RHS+FBr (from Predictor)         !
    ! RHScnt=RHS(part of it for 2ndCorr.)   !
    ! RHS:=RHS+1/2*dt*(Pe*Kappa.qstar)      !
    ! RHS:=RHS+(1/2)*dt*Fseg                !
    ! Fbarseg=Fseg; Fbarbead=Fbead          !
    ! Inside the loop:                      !
    ! RHSP:=RHSP+(1/4)*dt*(AdotDP.Fbarbead) !
    ! Fbarbead=-A'.Fbarseg                  !
    !---------------------------------------!
    call gemv(this%Bmat,qstar,rvmrcP)
    call copy(qc,RHS)
    call axpy(Kdotq,RHS,a=0.5_wp)

    ! rflc doesn't need this
    !            if ((EV_bb/='NoEV').or.(EV_bw/='NoEV') .and. EV_bw /= 'Rflc_bc') then
    if ((EV_bb/='NoEV').or.(EV_bw/='NoEV')) then
      !              call EVUpdate(Fev,rvmrcP,Fbarev)
      call myintrn%calc(id,itime,rvmrcP,rcmP,nseg,DiffTensP,divD,Fev,Fbarev,&
        updtevbb=.true.,updtevbw=.true.)
      !call print_vector(Fev,'fev3')
      !call evupdate2(Fev,rvmrcP,nseg,Fbarev)
      !call print_vector(Fev,'fev4')
      !stop
    end if

    if (ForceLaw == 'WLC_GEN') then
      call bndupdate(nbead_bb,Fbnd,qstar,Fbarbnd,itime)
    end if
    Fbar=Fbarev+Fbarbnd
    if (applFext) then
      Fbar(1)=Fbar(1)-Fext0
      Fbar(nbeadx3-2)=Fbar(nbeadx3-2)+Fext0
    end if
    if (srf_tet) then
      call tetupdate(Ftet,rvmrcP,rcm(:,ichain),DiffTensP,dt(iPe,idt),&
       Fbartet,rf0(:,ichain),itime)
      Fbar(1:3)=Fbar(1:3)+Fbartet(1:3)
      !call print_vector(Fbartet(1:3),'Tether force')
    end if
    call gemv(AdotDP1,Fbar,RHS,alpha=0.25*dt(iPe,idt),beta=1._wp)


    !TYL: HI for tethered bead -------------------------------------------
    !print *, 'First corrector calculaton------------'
    !call print_matrix(AdotDP1,'AdotDP1')
    !call print_vector(Fbar,'Fbar')
    !call print_vector(RHS,'RHS BEFORE')
    if (srf_tet) then
      call gemv(AdotDP1(:,1:3),Fbar(1:3),RHS,alpha=-0.25*dt(iPe,idt),&
        beta=1._wp)
      !call print_vector(RHS,'RHS AFTER')
    end if
    !TYL: HI for tethered bead -------------------------------------------


    !! Blake's part
    if ((hstar /= 0._WP) .and. (HITens == 'Blake')) then
      do is=1, nseg
        os=(is-1)*3
        RHS(os+2)=RHS(os+2)+(divD(is+1)-divD(is))*0.25*dt(iPe,idt)
      enddo
    endif
    !!-------------

    call axpy(FBr,RHS)
    call copy(RHS,RHScnt)
    call gemv(this%Kappa,qstar,RHS,alpha=0.5*Pe(iPe)*dt(iPe,idt),beta=1._wp)

    if (sph_flow) then
      call U_sph_sde(this,U_seg,U_bead,qstar,rcm(:,ichain))
      RHS = RHS + U_seg*0.5*dt(iPe,idt)
    endif

    call axpy(Fseg,RHS,a=0.5*dt(iPe,idt))
    call copy(Fseg,Fbarseg)
    call copy(Fbead,Fbarbead)
    do iseg=1, nseg
      offset=3*(iseg-1)
      RHSP => RHS(offset+1:offset+3)
      AdotDP2 => AdotD(offset+1:offset+3,:,ichain)
      call gemv(AdotDP2,Fbarbead,RHSP,alpha=0.25*dt(iPe,idt),beta=1._wp)

      !TYL: HI for tethered bead -----------------------------------------
      !print *, 'First corrector calculaton per segment------------'
      !call print_matrix(AdotDP2,'AdotDP2')
      !call print_vector(Fbarbead,'Fbarbead')
      !call print_vector(RHSP,'RHSP BEFORE')
      if (srf_tet) then
        call gemv(AdotDP2(:,1:3),Fbarbead(1:3),RHSP,&
        alpha=-0.25*dt(iPe,idt),beta=1._wp)
        !call print_vector(RHSP,'RHSP AFTER')
      end if
      !TYL: HI for tethered bead -----------------------------------------

      call sprupdate(id,root_f,PrScale,nroots,dt(iPe,idt),RHSP,qstar,iseg,&
        nseg,ForceLaw,TruncMethod,qbar,Fbarseg,Fbarbead,tplgy,this%Amat,nseg_bb,&
        nseg_ar,Ia,Na,itime)
    end do

    !----------Second Corrector Algorithm----------!
    ! q=qbar;Fseg=Fbarseg;Fbead=Fbarbead           !
    ! RHSbase=RHScnt(from 1stCorr.)for while loop. !
    ! While Loop,do loop:                          !
    ! RHSP=RHSbaseP                                !
    ! RHSP:=RHSP+(1/2)*dt*(Pe*Kappareg.qP)         !
    ! RHSP:=RHSP+(1/2)*dt*FsegP                    !
    ! RHSP:=RHSP+(1/4)dt*(AdotDP.Fbead)            !
    ! Updating q based on Seg. Cubic Eq.           !
    ! Fbead=-A'.Fseg                               !
    !----------------------------------------------!
    call copy(qbar,qc)
    call copy(Fbarseg,Fseg)
    call copy(Fbarbead,Fbead)
    call copy(RHScnt,RHSbase)
    icount=0;eps=1.0_wp
    do while (eps >= tol) !line 1183
      eps=0.0_wp
      qctemp=qc
      do iseg=1, nseg
        offset=3*(iseg-1)
        RHSP => RHS(offset+1:offset+3);RHSbaseP => RHSbase(offset+1:offset+3)
        call copy(RHSbaseP,RHSP)
        qcP => qc(offset+1:offset+3);FsegP => Fseg(offset+1:offset+3)
        AdotDP2 => AdotD(offset+1:offset+3,:,ichain)
        call gemv(this%Kappareg,qcP,RHSP,alpha=0.5*Pe(iPe)*dt(iPe,idt),beta=1.0_wp)

        if (sph_flow) then
          call U_sph_sde(this,U_seg,U_bead,qc,rcm(:,ichain))
          RHSP = RHSP + U_seg(offset+1:offset+3)*0.5*dt(iPe,idt)
        endif

        call axpy(FsegP,RHSP,a=0.5*dt(iPe,idt))
        call gemv(AdotDP2,Fbead,RHSP,alpha=0.25*dt(iPe,idt),beta=1.0_wp)

        !TYL: HI for tethered bead ---------------------------------------
        !print *, 'Second corrector calculaton perbead------------'
        !call print_matrix(AdotDP2,'AdotDP2')
        !call print_vector(Fbead,'Fbead')
        !call print_vector(RHSP,'RHSP BEFORE')
        if (srf_tet) then
          call gemv(AdotDP2(:,1:3),Fbead(1:3),RHSP,&
          alpha=-0.25*dt(iPe,idt),beta=1._wp)
          !call print_vector(RHSP,'RHSP AFTER')
        end if
        !TYL: HI for tethered bead ---------------------------------------

        call sprupdate(id,root_f,PrScale,nroots,dt(iPe,idt),RHSP,qbar,iseg,&
          nseg,ForceLaw,TruncMethod,qc,Fseg,Fbead,tplgy,this%Amat,nseg_bb,nseg_ar,&
          Ia,Na,itime)

      end do
      eps=nrm2(qc-qctemp)/nrm2(qctemp)
      icount=icount+1
      if (icount > 5000) then
        print *
        print '(" Convergance Problem in 2nd Corrector.")'
        print '(" time index: ",i10)',itime
        print '(" Total iterations: ",i10," Residual: ",f14.7)',icount,eps
        if (hstar /= 0._wp) then
          if (DecompMeth == 'Lanczos') then
            print '(" No. iterations in (block) Lanczos algorithm: ",i4)',&
            mch(ichain)
          elseif (DecompMeth == 'Chebyshev') then
            print '(" Eigen value range for  diffusion tensor: ",2(f14.7))',&
            lambdaBE(:)
            print '(" No. iterations in Chebyshev algorithm: ",i4)',Lch(ichain)
          end if
        end if
        stop
      end if
    end do ! while loop
    !==================================================!
    !line 1233, now put back into original arrays

  end subroutine advance_sde

  subroutine U_sph_sde(this,U_seg,U_bead,q,rcm)

    use :: inp_dlt, only: nseg,nbead,nsegx3,nbeadx3

    class(sde_t),intent(inout) :: this
    real(wp), intent(inout) :: U_seg(:) !nsegx3
    real(wp), intent(inout) :: U_bead(:) !nbeadx3
    real(wp), intent(in) :: q(:)
    real(wp), intent(in) :: rcm(:)

    real(wp),dimension(nbeadx3) :: rvmrc
    real(wp),dimension(3) :: r
    real(wp) :: th, ph, r_mag
    real(wp) :: u_r,u_th,u_ph
    integer :: i,j


    call gemv(this%Bmat,q,rvmrc)

    do i=1,nbead
      r(1:3) = rvmrc(3*(i-1)+1:3*(i-1)+3) + rcm(1:3)
      r_mag = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
      th = ACOS(r(1)/r_mag) !0<th<pi
      ph = ATAN2(r(3),r(2))                  !-pi<ph<pi
      u_r  =  cos(th) * (1 + this%a_sph**3/(2*r_mag**3) - 3*this%a_sph/2/r_mag)
      u_th = -sin(th) * (1 - this%a_sph**3/(4*r_mag**3) - 3*this%a_sph/4/r_mag)
      u_ph = 0._wp

      U_bead(3*(i-1)+1) = cos(th)*u_r - sin(th)*u_th
      U_bead(3*(i-1)+2) = sin(th)*cos(ph)*u_r + cos(th)*cos(ph)*u_th - sin(ph)*u_ph
      U_bead(3*(i-1)+3) = sin(th)*sin(ph)*u_r + cos(th)*sin(ph)*u_th + cos(ph)*u_ph
    enddo
    U_bead = this%U_mag_sph *U_bead

    do j=1,nseg
      U_seg(3*(j-1)+1:3*(j-1)+3) = U_bead(3*j+1:3*j+3) - U_bead(3*(j-1)+1:3*(j-1)+3)
    enddo

  end subroutine U_sph_sde

end module sde_mod
