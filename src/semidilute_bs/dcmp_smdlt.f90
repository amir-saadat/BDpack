!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is a free software: you can redistribute it and/or modify      |
!|  it under the terms of the GNU General Public License as published by  |
!|  the Free Software Foundation, either version 3 of the License, or     |
!|  (at your option) any later version.                                   |
!|                                                                        |
!|  BDpack is distributed in the hope that it will be useful,             |
!|  but WITHOUT ANY WARRANTY; without even the implied warranty of        |
!|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
!|  GNU General Public License for more details.                          |
!|                                                                        |
!|  You should have received a copy of the GNU General Public License     |
!|  along with BDpack.  If not, see <http://www.gnu.org/licenses/>.       |
!%------------------------------------------------------------------------%
module dcmp_smdlt

  use :: prcn_mod
  use :: arry_mod, only: print_vector
#ifdef USE_GPU
  use :: diffcalc_mod, only: PME_cpu,PME_dev
#else
  use :: diffcalc_mod, only: PME_cpu
#endif
  use :: types, only: decomp


#ifdef USE_GPU
  ! CUDA variables:
  !>
  real(wp),device,pointer :: VPt_d(:,:)
  !>
  real(wp),device,pointer :: VP_d(:)
  real(wp),device,pointer :: lamVPt(:)
  real(wp),device,pointer :: e1Pt_d(:)
#endif
  
contains

  subroutine Lanczos(Z,w,Ybar,nbeadx3,errormin,mub,minit,Y,decompRes,D,boxsizeinp,msetinp)
    implicit real(double) (a-h, o-z)
    integer :: nbeadx3,nbead
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3),intent(in) :: Z
    real(double),dimension(nbeadx3),intent(out) :: Y
    real(double),allocatable,dimension(:,:),target :: V,H,sqrtH,lambdaM
    real(double),allocatable,dimension(:,:) :: Vtemp,Htemp
    real(double),dimension(:),pointer :: VP
    real(double),dimension(:,:),pointer :: VPt,HPt,sqrtHPt,lambdaMPt
    real(double),allocatable,dimension(:),target :: lambdaV,e1
    real(double),allocatable,dimension(:) :: wtemp
    real(double),dimension(:),pointer :: lambdaVPt,e1Pt
    real(double),dimension(nbeadx3) :: w,Ybar
    real(double),dimension(nbeadx3,nbeadx3),optional,intent(in) :: D
    real(wp),optional,intent(in) :: boxsizeinp(3)
    real(wp),allocatable,dimension(:) :: wtemp2
    logical,optional :: msetinp
    type(decomp),intent(inout) :: decompRes

    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      stop
    end if
    if (minit > mub) mub=minit
    allocate(wtemp(nbeadx3),V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),lambdaM(mub,mub),lambdaV(mub),e1(mub))
#ifdef USE_SP
    allocate(wtemp2(nbeadx3))
#endif
    e1=0.0_wp;e1(1)=1.0_wp
    H=0.0_wp
    m=minit-1
    k0=1
    V(:,1)=Z/nrm2(Z)
mlp:do
      VPt => V(:,1:m)
      HPt => H(1:m,1:m)
      e1Pt => e1(1:m)
      sqrtHPt => sqrtH(1:m,1:m)
      lambdaMPt => lambdaM(1:m,1:m)
      lambdaVPt => lambdaV(1:m)
      do k=k0, m
        if ((k /= k0) .or. (k0 == 1)) then
          VP => V(:,k)
          if (present(D)) then ! means Ewald calculation:
            call symv(D,VP,w)
          else ! means PME:
            nbead=nbeadx3/3
! print*,'k,k0,m',k,k0,m
! call print_vector(VP,'wbef')
#ifdef USE_DP
            call PME_cpu(VP,nbead,boxsizeinp,w)
#elif USE_SP
            call PME_cpu(real(VP,kind=wp),nbead,boxsizeinp,wtemp2)
            w=real(wtemp2,kind=double)
#endif
! print*,'k,m',k,m
! call print_vector(w(1:50),'waf')

          end if
          if (k == m) wtemp=w
        else
          w=wtemp
        end if
        if (k > 1) then
          VP => V(:,k-1)
          w=w-H(k-1,k)*VP
        end if
        VP => V(:,k)
        if ((k /= k0) .or. (k0 == 1)) then
          H(k,k)=dot(w,VP)
        end if
        if (k < m) then
          w=w-H(k,k)*VP
          H(k,k+1)=nrm2(w);H(k+1,k)=H(k,k+1)
          VP => V(:,k+1)
          VP=w/H(k+1,k)
!         In case of D=I, H(k+1,k)=0. Then we make the 0/0=1.
          if (H(k+1,k) == 0) VP=1.0_wp
        end if
      end do
  
      sqrtHPt=HPt
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        print *
        print '(" Unsuccessful eigenvalue computation in Lanczos")'
        print '(" info:",1x,i3)',info
        stop
      end if
      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        print *
        print '(" Not a positive definite matrix in Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt)
        ! print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
!        stop
        decompRes%Success=.false.
        exit mlp
      end if
      lambdaMPt=0.0_wp
      forall (i=1:m) lambdaMPt(i,i)=lambdaVPt(i)
      sqrtHPt=matmul(sqrtHPt,matmul(sqrt(lambdaMPt),transpose(sqrtHPt)))
      Y=nrm2(Z)*matmul(VPt,matmul(sqrtHPt,e1Pt))
      if (present(msetinp)) then
        if (msetinp) exit mlp
      end if

! print*,'m',m,minit
! call print_vector(Y,'y')

      if (m >= minit) then

        error=nrm2(Y-Ybar)/nrm2(Ybar)
!        write(*,'(a,f14.7)') 'error',error
        if (error <= errormin) then
          exit mlp
        else
          if (m == mub) then
            allocate(Vtemp(nbeadx3,m),Htemp(m,m))
            Vtemp=V;Htemp=H
            deallocate(V,H,sqrtH,lambdaM,lambdaV,e1)
            mp10=m+10;mub=mp10
            allocate(V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),lambdaM(mub,mub),lambdaV(mub),e1(mub))
            e1=0.0_wp;e1(1)=1.0_wp
            H=0.0_wp
            V(:,1:m)=Vtemp;H(1:m,1:m)=Htemp
            deallocate(Vtemp,Htemp)
          end if
        end if
      end if
      Ybar=Y


      k0=m
      m=m+1
      cycle mlp
    end do mlp
    deallocate(wtemp,V,H,sqrtH,lambdaM,lambdaV,e1)
#ifdef USE_SP
    deallocate(wtemp2)
#endif
!   write(*,'(a,g10.4,a,i4)') 'Nember of Lanczos Iteration to achieve error<',errorinp,'=',m
    minit=m

  end subroutine Lanczos

  subroutine BlockLanczos(Z,a,W,Ybar,nbeadx3,s,errormin,mub,minit,Y,decompRes,D,boxsizeinp,msetinp)
    implicit real(double) (a-h, o-z)
    integer :: s,nbeadx3,nbead,icol
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3,s),intent(in) :: Z
    real(double),dimension(nbeadx3,s),intent(out) :: Y
    real(double),dimension(nbeadx3,s) :: a
    real(double),dimension(nbeadx3,s),target :: W
    real(double),allocatable,dimension(:,:),target :: V,H,R,sqrtH,lambdaM
    real(double),allocatable,dimension(:,:) :: Vtemp,Htemp,Rtemp,Wtemp
    real(double),dimension(:,:),pointer :: VP,HP
    real(double),dimension(:,:),pointer :: VPt,HPt,RPt,sqrtHPt,lambdaMPt
    real(double),allocatable,dimension(:),target :: lambdaV
    real(double),dimension(:),pointer :: lambdaVPt,VPP,WPtr
    real(double),dimension(nbeadx3) :: Ybar
    real(double),dimension(s) :: tau
    real(double),dimension(nbeadx3,nbeadx3),optional,intent(in) :: D
    real(wp),optional,intent(in) :: boxsizeinp(3)
    real(wp),allocatable,dimension(:) :: Wtemp2
    logical,optional :: msetinp
    type(decomp),intent(inout) :: decompRes

    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      stop
    end if
    if (minit > mub) mub=minit
    ms=mub*s
    allocate(Wtemp(nbeadx3,s),V(nbeadx3,ms),R(ms,s),H(ms,ms),sqrtH(ms,ms),lambdaM(ms,ms),lambdaV(ms))
#ifdef USE_SP
    allocate(Wtemp2(nbeadx3))
#endif
    
    H=0.0_wp
    a=Z;k=0
    call geqrf(a,tau=tau,info=info)
    if (info /= 0) then
      write(*,*) 'Unsuccessful QR fact. @geqrf() in Block Lanczos'
      write(*,'(a,i3,1x,a,i10)') 'info=',info,'at k=',k
      stop
    end if
    R=0.0_wp
    do i=1, s
      do j=i, s
        R(i,j)=a(i,j)
      end do
    end do
    call orgqr(a,tau,info=info)
    if (info /= 0) then
      write(*,*) 'Unsuccessful QR fact. @orgqr() in Block Lanczos'
      write(*,'(a,i3,1x,a,i10)') 'info=',info,'at k=',k
      stop
    end if
    V(:,1:s)=a
    m=minit-1;ms=m*s
    k0=1
mlp:do
      VPt => V(:,1:ms)
      HPt => H(1:ms,1:ms)
      sqrtHPt => sqrtH(1:ms,1:ms)
      RPt => R(1:ms,:)
      lambdaMPt => lambdaM(1:ms,1:ms)
      lambdaVPt => lambdaV(1:ms)
      do k=k0, m
        indexFkm1=(k-2)*s+1;indexLkm1=(k-1)*s
        indexFk=(k-1)*s+1;indexLk=k*s
        indexFkp1=k*s+1;indexLkp1=(k+1)*s
        if ((k /= k0) .or. (k0 == 1)) then
          VP => V(:,indexFk:indexLk)          
          if (present(D)) then ! means direct calculation:
            call symm(D,VP,W)
          else ! means PME:
            nbead=nbeadx3/3
            do icol=1, s
              VPP => VP(:,icol)
              WPtr => W(:,icol)
#ifdef USE_DP
              call PME_cpu(VPP,nbead,boxsizeinp,WPtr)
#elif USE_SP
              call PME_cpu(real(VPP,kind=wp),nbead,boxsizeinp,Wtemp2)
              WPtr=real(Wtemp2,kind=double)
#endif
            end do
          end if
          if (k == m) Wtemp=W
        else
          w=wtemp
        end if
        if (k > 1) then
          VP => V(:,indexFkm1:indexLkm1)
          HP => H(indexFkm1:indexLkm1,indexFk:indexLk)
          call gemm(VP,HP,W,alpha=-1._double,beta=1._double)
        end if
        VP => V(:,indexFk:indexLk)
        HP => H(indexFk:indexLk,indexFk:indexLk)
        if ((k /= k0) .or. (k0 == 1)) then
          call gemm(VP,W,HP,transa='T')
        end if
        if (k < m) then
          call gemm(VP,HP,W,alpha=-1._double,beta=1._double)
          a=W
          VP => V(:,indexFkp1:indexLkp1)
          HP => H(indexFkp1:indexLkp1,indexFk:indexLk)
          call geqrf(a,tau=tau,info=info)
          if (info /= 0) then
            write(*,*) 'Unsuccessful QR fact. @geqrf() in Block Lanczos'
            write(*,'(a,i3,1x,a,i10)') 'info=',info,'at k=',k
            stop
          end if
          do i=1, s
            do j=i, s
              HP(i,j)=a(i,j)
            end do
          end do
          H(indexFk:indexLk,indexFkp1:indexLkp1)=transpose(HP(:,:))
          call orgqr(a,tau,info=info)
          if (info /= 0) then
            write(*,*) 'Unsuccessful QR fact. @orgqr() in Block Lanczos'
            write(*,'(a,i3,1x,a,i10)') 'info=',info,'at k=',k
            stop
          end if
          V(:,indexFkp1:indexLkp1)=a(:,:)
        end if
      end do
  
      sqrtHPt=HPt 
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        print *
        print '(" Unsuccessful eigenvalue computation in Block Lanczos")'
        print '(" info:",1x,i3)',info
        stop
      end if
      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        print *
        print '(" Not a positive definite matrix in Block Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt)
        ! print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
!        stop
        decompRes%Success=.false.
        exit mlp
      end if
      lambdaMPt=0.0_wp
      forall (i=1:ms) lambdaMPt(i,i)=lambdaVPt(i)
      sqrtHPt=matmul(sqrtHPt,matmul(sqrt(lambdaMPt),transpose(sqrtHPt)))
      Y=matmul(VPt,matmul(sqrtHPt,RPt))
      if (present(msetinp)) then
        if (msetinp) exit mlp
      end if
      if (m >= minit) then
        error=nrm2(Y(:,1)-Ybar)/nrm2(Ybar)
!        write(*,'(a,f14.7)') 'error',error
        if (error <= errormin) then
          exit mlp
        else
          if (m == mub) then
            allocate(Vtemp(nbeadx3,m),Htemp(m,m),Rtemp(ms,s))
            Vtemp=V;Htemp=H;Rtemp=R
            deallocate(V,R,H,sqrtH,lambdaM,lambdaV)
            mp10=m+10;mub=mp10;ms=mub*s
            allocate(V(nbeadx3,ms),R(ms,s),H(ms,ms),sqrtH(ms,ms),lambdaM(ms,ms),lambdaV(ms))
            H=0.0_wp;R=0.0_wp
            V(:,1:m)=Vtemp;H(1:m,1:m)=Htemp;R(1:s,:)=Rtemp(1:s,:)
            deallocate(Vtemp,Htemp,Rtemp)
          end if
        end if
      end if 
      Ybar=Y(:,1)
      m=m+1;ms=m*s
      cycle mlp
    end do mlp   
    deallocate(Wtemp,V,R,H,sqrtH,lambdaM,lambdaV)
#ifdef USE_SP
    deallocate(Wtemp2)
#endif
!    write(*,'(a,i4)') 'Nember of Lanczos Iteration=',m
    minit=m

  end subroutine BlockLanczos







#ifdef USE_GPU

  subroutine Lanczos_dev(hi_d,Z,w,Ybar,nbeadx3,errormin,mub,minit,Y,&
    decompRes,D,boxsizeinp,msetinp)

    use :: hi_cumod, only: hi_cu_t

    implicit real(double) (a-h, o-z)

    type(hi_cu_t),intent(inout) :: hi_d
    integer :: nbeadx3,nbead
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3),intent(in) :: Z
    real(double),dimension(nbeadx3),intent(out) :: Y
    real(double),allocatable,dimension(:,:),target :: V,H,sqrtH,lambdaM
    real(double),allocatable,dimension(:,:) :: Vtemp,Htemp
    real(double),dimension(:),pointer :: VP
    real(double),dimension(:,:),pointer :: VPt,HPt,sqrtHPt,lambdaMPt
    real(double),allocatable,dimension(:),target :: lambdaV,e1
    real(double),allocatable,dimension(:) :: wtemp
    real(double),dimension(:),pointer :: lambdaVPt,e1Pt
    real(double),dimension(nbeadx3) :: w,Ybar
    real(double),dimension(nbeadx3,nbeadx3),optional,intent(in) :: D
    real(wp),optional,intent(in) :: boxsizeinp(3)
    real(wp),allocatable,dimension(:) :: wtemp2
    logical,optional :: msetinp
    type(decomp),intent(inout) :: decompRes

    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      stop
    end if
    if (minit > mub) mub=minit
    allocate(wtemp(nbeadx3),V(nbeadx3,mub),H(mub,mub),&
      sqrtH(mub,mub),lambdaM(mub,mub),lambdaV(mub),e1(mub))
#ifdef USE_SP
    allocate(wtemp2(nbeadx3))
#endif
    e1=0.0_wp;e1(1)=1.0_wp
    H=0.0_wp
    m=minit-1
    k0=1
    V(:,1)=Z/nrm2(Z)
mlp:do
      VPt => V(:,1:m)
      HPt => H(1:m,1:m)
      e1Pt => e1(1:m)
      sqrtHPt => sqrtH(1:m,1:m)
      lambdaMPt => lambdaM(1:m,1:m)
      lambdaVPt => lambdaV(1:m)
      do k=k0, m
        if ((k /= k0) .or. (k0 == 1)) then
          VP => V(:,k)
          if (present(D)) then ! means Ewald calculation:
            call symv(D,VP,w)
          else ! means PME:
            nbead=nbeadx3/3
#ifdef USE_DP
            call PME_dev(hi_d,VP,nbead,boxsizeinp,w)
            ! call PME_cpu(VP,nbead,boxsizeinp,w)
#elif USE_SP
!!!!!!!!!!!!!! take care of that later
            call PME_cpu(real(VP,kind=wp),nbead,boxsizeinp,wtemp2)
            w=real(wtemp2,kind=double)
#endif
          end if
          if (k == m) wtemp=w
        else
          w=wtemp
        end if
        if (k > 1) then
          VP => V(:,k-1)
          w=w-H(k-1,k)*VP
        end if
        VP => V(:,k)
        if ((k /= k0) .or. (k0 == 1)) then
          H(k,k)=dot(w,VP)
        end if
        if (k < m) then
          w=w-H(k,k)*VP
          H(k,k+1)=nrm2(w);H(k+1,k)=H(k,k+1)
          VP => V(:,k+1)
          VP=w/H(k+1,k)
!         In case of D=I, H(k+1,k)=0. Then we make the 0/0=1.
          if (H(k+1,k) == 0) VP=1.0_wp
        end if
      end do

      print*,'H',H
  
      sqrtHPt=HPt
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        print *
        print '(" Unsuccessful eigenvalue computation in Lanczos")'
        print '(" info:",1x,i3)',info
        stop
      end if
      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        print *
        print '(" Not a positive definite matrix in Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt)
        ! print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
!        stop
        decompRes%Success=.false.
        exit mlp
      end if
      lambdaMPt=0.0_wp
      forall (i=1:m) lambdaMPt(i,i)=lambdaVPt(i)
      sqrtHPt=matmul(sqrtHPt,matmul(sqrt(lambdaMPt),transpose(sqrtHPt)))
      Y=nrm2(Z)*matmul(VPt,matmul(sqrtHPt,e1Pt))
      if (present(msetinp)) then
        if (msetinp) exit mlp
      end if
      if (m >= minit) then
        error=nrm2(Y-Ybar)/nrm2(Ybar)
!        write(*,'(a,f14.7)') 'error',error
        if (error <= errormin) then
          exit mlp
        else
          if (m == mub) then
            allocate(Vtemp(nbeadx3,m),Htemp(m,m))
            Vtemp=V;Htemp=H
            deallocate(V,H,sqrtH,lambdaM,lambdaV,e1)
            mp10=m+10;mub=mp10
            allocate(V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),lambdaM(mub,mub),lambdaV(mub),e1(mub))
            e1=0.0_wp;e1(1)=1.0_wp
            H=0.0_wp
            V(:,1:m)=Vtemp;H(1:m,1:m)=Htemp
            deallocate(Vtemp,Htemp)
          end if
        end if
      end if
      Ybar=Y
      k0=m
      m=m+1
      cycle mlp
    end do mlp
    deallocate(wtemp,V,H,sqrtH,lambdaM,lambdaV,e1)
#ifdef USE_SP
    deallocate(wtemp2)
#endif
!   write(*,'(a,g10.4,a,i4)') 'Nember of Lanczos Iteration to achieve error<',errorinp,'=',m
    minit=m

  end subroutine Lanczos_dev


  subroutine Lanczos_dev2(hi_d,Z,w,Ybar,nbeadx3,errormin,mub,minit,Y,&
    decompRes,D,boxsizeinp,msetinp)

    use :: hi_cumod, only: hi_cu_t
    use :: cublas
    use :: cudafor
    use :: cusolverdn

    implicit none
    ! implicit real(double) (a-h, o-z)

    type(hi_cu_t),intent(inout) :: hi_d
    integer :: nbeadx3,nbead
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3),intent(in) :: Z
    real(double),dimension(nbeadx3),intent(out) :: Y
    real(double),allocatable,dimension(:,:),target :: V,H,sqrtH,lambdaM
    real(double),allocatable,dimension(:,:) :: Vtemp,Htemp
    real(double),dimension(:),pointer :: VP
    real(double),dimension(:,:),pointer :: VPt,HPt,sqrtHPt,lambdaMPt
    real(double),allocatable,dimension(:),target :: lambdaV,e1
    real(double),allocatable,dimension(:) :: wtemp
    real(double),dimension(:),pointer :: lambdaVPt,e1Pt
    real(double),dimension(nbeadx3) :: w,Ybar
    real(double),dimension(nbeadx3,nbeadx3),optional,intent(in) :: D
    real(wp),optional,intent(in) :: boxsizeinp(3)
    real(wp),allocatable,dimension(:) :: wtemp2
    logical,optional :: msetinp
    type(decomp),intent(inout) :: decompRes
    real(wp) :: errormin,error
    integer :: m,k0,k,info,mp10,i



    real(wp),allocatable,device :: Y_d(:)

    integer :: XYZ_dim
    type(dim3) :: dimGrid, dimBlock


    real(wp) :: H_tmp,z_nrm,lam_min,ybar_nrm
    integer :: istat,Lwork,irow,jcol

    integer,device :: devInfo_d

    real(wp),device,allocatable :: Z_d(:)
    real(wp),device,allocatable :: w_d(:)
    real(wp),device,allocatable,target :: V_d(:,:)
    real(wp),device,allocatable :: wtmp(:)
    real(wp),device,allocatable :: workspace_d(:)


    real(wp),device,allocatable,target :: H_d(:,:),sqrtH_d(:,:),lamM(:,:)
    real(wp),device,allocatable,target :: lamV(:),e1_d(:)
    real(wp),device,allocatable :: Ybar_d(:)
    real(wp),device,allocatable :: sqrtH_tmp(:,:),Vtmp(:,:),Htmp(:,:)

    ! XYZ_dim=64
    ! dimGrid = dim3( (XYZ_dim+31)/32, 1, 1 )
    ! dimBlock = dim3( 32, 1, 1 )
    ! call printtest <<<dimGrid,dimBlock>>>
    ! status = cudathreadsynchronize()


    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      minit=2
    end if

    if (minit > mub) mub=minit
    allocate(wtemp(nbeadx3),V(nbeadx3,mub),H(mub,mub),&
      sqrtH(mub,mub),lambdaM(mub,mub),lambdaV(mub),e1(mub))

    ! allocating device arrays
    allocate(wtmp(nbeadx3),V_d(nbeadx3,mub),H_d(mub,mub),lamV(mub),e1_d(mub),Ybar_d(nbeadx3))
    allocate(Z_d(nbeadx3))
    allocate(w_d(nbeadx3))
    allocate(Y_d(nbeadx3))
    Z_d=Z

    e1=0.0_wp;e1(1)=1.0_wp
    H=0.0_wp
    m=minit-1
    k0=1
    V(:,1)=Z/nrm2(Z)

    e1_d=0._wp;e1_d(1)=1.0_wp
    H_d=0.0_wp
    m=minit-1
    k0=1
    call cublasDcopy(nbeadx3,Z_d,1,V_d,1)
    z_nrm=dnrm2(nbeadx3,Z_d,1)
    call cublasDscal(nbeadx3,1/z_nrm,V_d,1)

mlp:do


      VPt => V(:,1:m)
      HPt => H(1:m,1:m)
      e1Pt => e1(1:m)
      sqrtHPt => sqrtH(1:m,1:m)
      lambdaMPt => lambdaM(1:m,1:m)
      lambdaVPt => lambdaV(1:m)



      VPt_d => V_d(:,1:m)
      e1Pt_d => e1_d(1:m)
      lamVPt => lamV(1:m)

      do k=k0, m

        if ((k /= k0) .or. (k0 == 1)) then


          VP => V(:,k)
          
          VP_d => V_d(:,k)

          if (present(D)) then ! means Ewald calculation:
            ! doesn't make sense now
            ! call symv(D,VP,w)
          else ! means PME:
            nbead=nbeadx3/3
            call PME_dev(hi_d,VP,nbead,boxsizeinp,w)
            ! print*,'w',k,w
            w_d=w
            ! call PME_cpu(VP,nbead,boxsizeinp,w)
          end if

          if (k == m) wtemp=w
          if (k == m) wtmp=w_d

        else ! k==k0 and k0/=1

          w=wtemp
          w_d=wtmp

        end if

        if (k > 1) then
          VP => V(:,k-1)
          w=w-H(k-1,k)*VP
        end if

        if (k > 1) then
          VP_d => V_d(:,k-1)
          H_tmp=-H_d(k-1,k)
          call cublasDaxpy(nbeadx3,H_tmp,VP_d,1,w_d,1)
        end if

        VP => V(:,k)
        if ((k /= k0) .or. (k0 == 1)) then
          H(k,k)=dot(w,VP)
        end if

        VP_d => V_d(:,k)
        if ((k /= k0) .or. (k0 == 1)) then
          H_d(k,k)=ddot(nbeadx3,w_d,1,VP_d,1)
          H_tmp=H_d(k,k)
        end if


        if (k < m) then
          w=w-H(k,k)*VP
          H(k,k+1)=nrm2(w);H(k+1,k)=H(k,k+1)
          VP => V(:,k+1)
          VP=w/H(k+1,k)
!         In case of D=I, H(k+1,k)=0. Then we make the 0/0=1.
          if (H(k+1,k) == 0) VP=1.0_wp
        end if

        if (k < m) then
          H_tmp=-H_d(k,k)
          call cublasDaxpy(nbeadx3,H_tmp,VP_d,1,w_d,1)
          H_tmp=dnrm2(nbeadx3,w_d,1)
          H_d(k,k+1)=H_tmp
          H_d(k+1,k)=H_tmp
          VP_d => V_d(:,k+1)
          call cublasDcopy(nbeadx3,w_d,1,VP_d,1)
          call cublasDscal(nbeadx3,1/H_tmp,VP_d,1)
!         In case of D=I, H(k+1,k)=0. Then we make the 0/0=1.
          if (H_tmp == 0) VP_d=1.0_wp
        end if

      end do ! k

      ! print*,'H',H
      ! H=H_d
      ! print*,'Hd',H
      ! stop


  

      sqrtHPt=HPt
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        print *
        print '(" Unsuccessful eigenvalue computation in Lanczos")'
        print '(" info:",1x,i3)',info
        stop
      end if

      allocate(sqrtH_d(m,m),sqrtH_tmp(m,m))
      allocate(lamM(m,m))

      !$cuf kernel do (2) <<< *,* >>>
      do irow=1, m
        do jcol=1, m
          sqrtH_d(irow,jcol)=H_d(irow,jcol)
        enddo
      enddo

      istat = cusolverDnDsyevd_bufferSize(hi_d%h_lam,CUSOLVER_EIG_MODE_VECTOR,&
        CUBLAS_FILL_MODE_UPPER,m,sqrtH_d,m,lamVPt,Lwork)
      if (istat /= CUSOLVER_STATUS_SUCCESS) print('(" cusolverDnDsyevd_bufferSize failed")')
      allocate(workspace_d(Lwork))

      istat = cusolverDnDsyevd(hi_d%h_lam,CUSOLVER_EIG_MODE_VECTOR,&
        CUBLAS_FILL_MODE_UPPER,m,sqrtH_d,m,lamVPt,workspace_d,Lwork,devInfo_d)
      if (istat /= CUSOLVER_STATUS_SUCCESS) print('(" cusolverDnDsyevd failed")')
      istat = devInfo_d
      if (istat /= 0) print('(" Unsuccessful eigenvalue computation in Lanczos_dev",i)'),istat
  


      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        print *
        print '(" Not a positive definite matrix in Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt)
        ! print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
!        stop
        decompRes%Success=.false.
        exit mlp
      end if

      ! print*,'minlam',minval(lambdaVPt)


      lam_min=lamVPt(cublasIdamin(m,lamVPt,1))
      if (lam_min<=0) then
        print *
        print '(" Not a positive definite matrix in Lanczos_dev")'
        print '(" min(EigenValue) of D: ",f7.3)',lam_min
        decompRes%Success=.false.
        exit mlp
      end if

      ! print*,'lam_min',lam_min


      lambdaMPt=0.0_wp
      forall (i=1:m) lambdaMPt(i,i)=lambdaVPt(i)
      sqrtHPt=matmul(sqrtHPt,matmul(sqrt(lambdaMPt),transpose(sqrtHPt)))
      Y=nrm2(Z)*matmul(VPt,matmul(sqrtHPt,e1Pt))
      if (present(msetinp)) then
        if (msetinp) exit mlp
      end if

      ! print*,'Y',Y(1:10)


      lamM=0._wp
      !$cuf kernel do <<< *,* >>>
      do irow=1, m
        lamM(irow,irow)=sqrt(lamVPt(irow))
      enddo
      call cublasDgemm('N','T',m,m,m,1._wp,lamM,m,sqrtH_d,m,0._wp,sqrtH_tmp,m)
      call cublasDgemm('N','N',m,m,m,1._wp,sqrtH_d,m,sqrtH_tmp,m,0._wp,lamM,m)
      call cublasDgemv('N',m,m,1._wp,lamM,m,e1Pt_d,1,0._wp,lamVPt,1)
      call cublasDgemv('N',nbeadx3,m,z_nrm,VPt_d,nbeadx3,lamVPt,1,0._wp,Y_d,1)
      if (present(msetinp)) then
        if (msetinp) exit mlp
      end if
      ! Y=Y_d
      ! print*,'Y-2',Y(1:10)
      ! stop


      if (m >= minit) then
        error=nrm2(Y-Ybar)/nrm2(Ybar)
       ! write(*,'(a,f14.7)') 'error',error
        if (error <= errormin) then
          exit mlp
        else
          if (m == mub) then
            allocate(Vtemp(nbeadx3,m),Htemp(m,m))
            Vtemp=V;Htemp=H
            deallocate(V,H,sqrtH,lambdaM,lambdaV,e1)
            mp10=m+10;mub=mp10
            allocate(V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),&
              lambdaM(mub,mub),lambdaV(mub),e1(mub))
            e1=0.0_wp;e1(1)=1.0_wp
            H=0.0_wp
            V(:,1:m)=Vtemp;H(1:m,1:m)=Htemp
            deallocate(Vtemp,Htemp)
          end if
        end if
      end if
      Ybar=Y

! print*,'y',m,Y


      if (m >= minit) then
        ybar_nrm=dnrm2(nbeadx3,Ybar_d,1)
        !$cuf kernel do <<< *,* >>>
        do irow=1, nbeadx3
          Ybar_d(irow)=Y_d(irow)-Ybar_d(irow)
        enddo
        error=dnrm2(nbeadx3,Ybar_d,1)/ybar_nrm
       ! write(*,'(a,f14.7)') 'error2',error
        if (error <= errormin) then
          exit mlp
        else
          if (m == mub) then
            allocate(Vtmp(nbeadx3,m),Htmp(m,m))
            Vtmp=V_d;Htmp=H_d
            deallocate(V_d,H_d,sqrtH_d,lamM,lamV,e1_d)
            mp10=m+10;mub=mp10
            allocate(V_d(nbeadx3,mub),H_d(mub,mub),lamV(mub),e1_d(mub))
            e1_d=0._wp;e1_d(1)=1._wp
            H_d=0._wp
            !$cuf kernel do (2) <<< *,* >>>
            do jcol=1, m
              do irow=1, nbeadx3             
                V_d(irow,jcol)=Vtmp(irow,jcol)
              enddo
            enddo
            !$cuf kernel do (2) <<< *,* >>>
            do jcol=1, m
              do irow=1, m            
                H_d(irow,jcol)=Htmp(irow,jcol)
              enddo
            enddo
            deallocate(Vtmp,Htmp)
          end if
        end if
      end if
      call cublasDcopy(nbeadx3,Y_d,1,Ybar_d,1)
      deallocate(sqrtH_d,lamM,sqrtH_tmp)



      k0=m
      m=m+1      
      cycle mlp
    end do mlp


    deallocate(wtemp,V,H,sqrtH,lambdaM,lambdaV,e1)

    deallocate(wtmp,V_d,H_d,lamV,e1_d,Ybar_d)

    minit=m

  end subroutine Lanczos_dev2

#endif













end module dcmp_smdlt
