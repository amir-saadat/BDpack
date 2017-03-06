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
  use :: diffcalc_mod, only: PME_cpu
  use :: types, only: decomp 
  
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
#ifdef USE_DP
            call PME_cpu(VP,nbead,boxsizeinp,w)
#elif USE_SP
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
  
      sqrtHPt=HPt
      call syev(sqrtHPt,lambdaVPt,jobz='V',uplo='U',info=info)
      if (info /= 0) then
        print *
        print '(" Unsuccessful eigenvalue computation in Lanczos")'
        print '(" info:",1x,i3)',info
        stop
      end if
      if (minval(lambdaVPt,dim=m)<=0) then
        print *
        print '(" Not a positive definite matrix in Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
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
      if (minval(lambdaVPt,dim=m)<=0) then
        print *
        print '(" Not a positive definite matrix in Block Lanczos")'
        print '(" min(EigenValue) of D: ",f7.3)',minval(lambdaVPt,dim=m)
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

end module dcmp_smdlt
