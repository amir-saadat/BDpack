!%------------------------------------------------------------------------%
!|  Copyright (C) 2013 - 2016:                                            |
!|  Material Research and Innovation Laboratory (MRAIL)                   |
!|  University of Tennessee-Knoxville                                     |
!|  Author:    Amir Saadat   <asaadat@vols.utk.edu>                       |
!|  Advisor:   Bamin Khomami <bkhomami@utk.edu>                           |
!|                                                                        |
!|  This file is part of BDpack.                                          |
!|                                                                        |
!|  BDpack is free software: you can redistribute it and/or modify        |
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
module dcmp_mod

  use :: prcn_mod

contains

  subroutine Lanczos(D,Z,w,Ybar,nbeadx3,errormin,mub,minit,Y,msetinp)
    implicit real(double) (a-h, o-z)
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3,nbeadx3),intent(in) :: D
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
    logical,optional :: msetinp

    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      stop
    end if
    if (minit > mub) mub=minit
    allocate(wtemp(nbeadx3),V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),lambdaM(mub,mub),&
             lambdaV(mub),e1(mub))
    e1=0.0_double;e1(1)=1.0_double
    H=0.0_double
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
          call symv(D,VP,w)
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
          if (H(k+1,k) == 0) VP=1.0_double
        end if
      end do

      sqrtHPt=HPt
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        write(*,*) 'Unsuccessful eigenvalue computation in Lanczos'
        write(*,'(a,1x,i3)') 'info:',info
        stop
      end if
      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        write(*,*) 'Not a positive definite matrix in Lanczos'
        stop
      end if
      lambdaMPt=0.0_double
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
            allocate(V(nbeadx3,mub),H(mub,mub),sqrtH(mub,mub),lambdaM(mub,mub),&
                     lambdaV(mub),e1(mub))
            e1=0.0_double;e1(1)=1.0_double
            H=0.0_double
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
!   write(*,'(a,g10.4,a,i4)') 'Nember of Lanczos Iteration to achieve error<',errorinp,'=',m
    minit=m

  end subroutine Lanczos

  subroutine BlockLanczos(D,Z,a,W,Ybar,nbeadx3,s,errormin,mub,minit,Y,msetinp)
    implicit real(double) (a-h, o-z)
    integer :: s
    integer,intent(inout) :: minit,mub
    real(double),dimension(nbeadx3,nbeadx3),intent(in) :: D
    real(double),dimension(nbeadx3,s),intent(in) :: Z
    real(double),dimension(nbeadx3,s),intent(out) :: Y
    real(double),dimension(nbeadx3,s) :: a,W
    real(double),allocatable,dimension(:,:),target :: V,H,R,sqrtH,lambdaM
    real(double),allocatable,dimension(:,:) :: Vtemp,Htemp,Rtemp,Wtemp
    real(double),dimension(:,:),pointer :: VP,HP
    real(double),dimension(:,:),pointer :: VPt,HPt,RPt,sqrtHPt,lambdaMPt
    real(double),allocatable,dimension(:),target :: lambdaV
    real(double),dimension(:),pointer :: lambdaVPt
    real(double),dimension(nbeadx3) :: Ybar
    real(double),dimension(s) :: tau
    logical,optional :: msetinp

    if (minit < 2) then
      write(*,'(a)') "Error: minit can't be less than 2"
      stop
    end if
    if (minit > mub) mub=minit
    ms=mub*s
    allocate(Wtemp(nbeadx3,s),V(nbeadx3,ms),R(ms,s),H(ms,ms),sqrtH(ms,ms),lambdaM(ms,ms),&
             lambdaV(ms))
    H=0.0_double
    a=Z;k=0
    call geqrf(a,tau=tau,info=info)
    if (info /= 0) then
      write(*,*) 'Unsuccessful QR fact. @geqrf() in Block Lanczos'
      write(*,'(a,i3,1x,a,i10)') 'info=',info,'at k=',k
      stop
    end if
    R=0.0_double
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
          call symm(D,VP,W)
          if (k == m) Wtemp=W
        else
          w=wtemp
        end if
        if (k > 1) then
          VP => V(:,indexFkm1:indexLkm1)
          HP => H(indexFkm1:indexLkm1,indexFk:indexLk)
          call gemm(VP,HP,W,alpha=-1.0_double,beta=1.0_double)
        end if
        VP => V(:,indexFk:indexLk)
        HP => H(indexFk:indexLk,indexFk:indexLk)
        if ((k /= k0) .or. (k0 == 1)) then
          call gemm(VP,W,HP,transa='T')
        end if
        if (k < m) then
          call gemm(VP,HP,W,alpha=-1.0_double,beta=1.0_double)
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
        write(*,*) 'Unsuccessful eigenvalue computation in Block Lanczos'
        write(*,'(a,1x,i3)') 'info:',info
        stop
      end if
      ! if (minval(lambdaVPt,dim=m)<=0) then
      if (minval(lambdaVPt)<=0) then
        write(*,*) 'Not a positive definite matrix in Block Lanczos'
        stop
      end if
      lambdaMPt=0.0_double
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
            allocate(V(nbeadx3,ms),R(ms,s),H(ms,ms),sqrtH(ms,ms),lambdaM(ms,ms),&
                     lambdaV(ms))
            H=0.0_double;R=0.0_double
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
!    write(*,'(a,i4)') 'Nember of Lanczos Iteration=',m
    minit=m

  end subroutine BlockLanczos

  subroutine BlockChebyshev(D,Eye,Dsh,dW,nbeadx3,s,Lub,Linit,dS,dSbar,BothEndEigVal,&
                               errormin,lambdainp,Lsetinp)
    implicit real(double) (a-h, o-z)
    integer :: s
    integer,intent(inout) :: Linit,Lub
    real(double),dimension(nbeadx3,nbeadx3),intent(in) :: D,Eye
    real(double),dimension(nbeadx3,nbeadx3) :: Dsh
    real(double),dimension(nbeadx3,s),intent(in) :: dW
    real(double),dimension(nbeadx3,s),intent(out) :: dS
    real(double),allocatable,dimension(:,:),target :: dV
    real(double),allocatable,dimension(:,:) :: dVtemp
    real(double),dimension(:,:),pointer :: dV1,dVk,dVkm1,dVkm2
    real(double),dimension(nbeadx3) :: dSbar
    real(double),dimension(2) :: lambda
    real(double),dimension(2),intent(in),optional :: lambdainp
    real(double),allocatable,dimension(:) :: c
    logical,optional :: Lsetinp
    interface
      function BothEndEigVal(D) result(beEigVal)
        import :: double
        implicit real(double) (a-h, o-z)
        real(double),dimension(:,:),intent(in) :: D
        real(double),dimension(2) :: beEigVal
      end function BothEndEigVal
    end interface

    if (Linit < 2) then
      write(*,'(a)') "Error: Linit can't be less than 2"
      stop
    end if
    if (present(lambdainp)) then
      lambda=lambdainp
    else
      lambda=BothEndEigVal(D)
      lambda(1)=lambda(1)/2;lambda(2)=2*lambda(2)
      if (lambda(2)/lambda(1) >= 10.0e5_double) then
        Lub=100
      else
        Lub=nint(sqrt(lambda(2)/lambda(1))+10)
      end if
    end if
    if (lambda(1)<=0) then
      write(*,*) 'Not a positive definite matrix in Chebyshev'
      stop
    end if
    if (Linit > Lub) Lub=Linit
    allocate(c(0:Lub),dV(nbeadx3,Lub*s))
    L=Linit-1
    k0=1
    am=(lambda(2)-lambda(1))/2;ap=(lambda(2)+lambda(1))/2
    Dsh=1/am*D-ap/am*Eye
    dV1 => dV(:,1:s)
    call symm(Dsh,dW,dV1)
llp:do
      c(0:L)=ChebCoeff(am,ap,L)
      dS=c(0)/2*dW+c(1)*dV1
      do k=2, L
        dVk => dV(:,(k-1)*s+1:k*s)
        if (k > k0) then
          dVkm1 => dV(:,(k-2)*s+1:(k-1)*s)
          if (k == 2) then
            dVk=dW
          else
            dVkm2 => dV(:,(k-3)*s+1:(k-2)*s)
            dVk=dVkm2
          end if
          call symm(Dsh,dVkm1,dVk,alpha=2.0_double,beta=-1.0_double)
        end if
        dS=dS+c(k)*dVk
      end do
      if (present(Lsetinp)) then
        if (Lsetinp) exit llp
      end if
      if (L >= Linit) then
        error=nrm2(dS(:,1)-dSbar(:))/nrm2(dSbar(:))
!        write(*,'(a,f14.7)') 'error',error
        if (error <= errormin) then
          exit llp
        else
          if (L == Lub) then
            lambda=BothEndEigVal(D)
            lambda(1)=lambda(1)/2;lambda(2)=2*lambda(2)
            am=(lambda(2)-lambda(1))/2;ap=(lambda(2)+lambda(1))/2
            if (lambda(2)/lambda(1) >= 10.0e5_double) then
              Lub=100
            else
              Lub=nint(sqrt(lambda(2)/lambda(1))+10)
            end if
            if (L >= Lub) Lub=L+10
            Dsh=1/am*D-ap/am*Eye
            call symm(Dsh,dW,dV1)
            allocate(dVtemp(nbeadx3,L*s))
            dVtemp=dV
            deallocate(c,dV)
            allocate(c(0:Lub),dV(nbeadx3,Lub*s))
            dV1 => dV(:,1:s)
            dV(:,1:L*s)=dVtemp
            deallocate(dVtemp)
            dSbar=dS(:,1)
            k0=1
            L=L+1
            cycle llp
          end if
        end if
      end if
      dSbar=dS(:,1)
      k0=L
      L=L+1
      cycle llp
    end do llp
    deallocate(c)
!   write(*,'(a,g10.4,a,i4)') 'Chebyshev Polynomial degree to achieve error<',errormin,'=',L
    Linit=L

  end subroutine BlockChebyshev

  function ChebCoeff(am,ap,L) result(c)
    implicit real(double) (a-h, o-z)
    real(double),parameter :: PI=4*atan(1.0_double)
    real(double),intent(in) :: am,ap
    integer,intent(in) :: L
    real(double),dimension(1:L+1) :: g
    real(double),dimension(0:L) :: c
    Lp1=L+1
    do k=1, Lp1
      y=cos(PI*(k-0.5)/Lp1)
      g(k)=sqrt(am*y+ap)
    end do
    fac=2.0_double/Lp1
    do j=0, L
      summ=0.0_double
      do k=1, Lp1
        summ=summ+g(k)*cos(PI*j*(k-0.5_double)/Lp1)
      end do
      c(j)=fac*summ
    end do
  end function ChebCoeff

  function MKLsyevr(D) result(beEigVal)
    implicit real(double) (a-h, o-z)
    real(double),dimension(:,:),intent(in) :: D
    real(double),allocatable,dimension(:,:) :: Dtemp
    real(double),allocatable,dimension(:) :: w
    real(double),dimension(2) :: beEigVal

    nbeadx3=size(D,1)
    allocate(Dtemp(nbeadx3,nbeadx3),w(nbeadx3))
    Dtemp=D
    call syevr(Dtemp,w,info=info)
    if (info /= 0) then
      write(*,*) 'unsuccessfull eigen value computation'
      stop
    end if
    beEigVal=(/w(1),w(nbeadx3)/)
    deallocate(Dtemp,w)
  end function MKLsyevr

end module dcmp_mod
