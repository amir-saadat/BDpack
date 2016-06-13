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
module root_mod

  implicit none
  integer, parameter :: double = selected_real_kind(15)

contains

  subroutine root_fndr(coeffs,upper,root)

    use :: arry_mod, only: print_vector,print_matrix
    use :: lapack95, only: geev

    real(double),intent(in) :: coeffs(:),upper
    real(double),intent(out) :: root
    real(double),allocatable :: M(:,:)
    real(double),allocatable :: evr(:),evi(:)
    integer :: N,i,info,nr
    real(double) :: fctr

    ! Constructing the companion matrix
    N=size(coeffs,1)-1
    allocate(M(N,N),evr(N),evi(N))    
    M=0.d0
    fctr=coeffs(N+1) ! We set C(N+1)=1
    forall (i=1:N-1) M(i,i+1)=1.d0
    M(N,:)=-coeffs(1:N)/fctr

    call geev(M,evr,evi,info=info)
    if (info /= 0) then
      print '(" Error in root_fndr: ",i0)',info
      stop
    end if

    nr=0
    do i=1, N
      if ( (evi(i) == 0.d0) .and. &
           (evr(i) >= 0.d0) .and. &
           (evr(i) <= upper)) then        
        root=evr(i)
        nr=nr+1
      end if
    end do
    if (nr > 1) then
      print '("Error: More than one root found in root_fndr.")'
      stop
    end if    

    deallocate(M,evr,evi)

  end subroutine root_fndr

  subroutine CubeRoot(a1,a2,a3,upper,root)

    real(double) :: a1,a2,a3,upper,root
    real(double) :: p,q,R,x,y,A,B,phi,arg1,arg2,arg3,y1,y2,y3
    real(double),parameter :: PI=3.1415926535897958648d0
    !----------------------------------------------!
    !  subroutine returns the cuberoot, whose      !
    !  value lies between 0 and qmax. a1,a2 and a3 !
    !  are the coefficients of the equation for |Q|!
    !----------------------------------------------!
    p=(3*a2-a1**2)/3
    q=(2*a1**3-9*a2*a1+27*a3)/27
    R=(p/3)**3+(q/2)**2
    if (R > 0.d0) then
      x=sqrt(R)-q/2
      y=-sqrt(R)-q/2
      if (x > 0.d0) then
        A=x**(1.d0/3)
      else
        A=-(-x)**(1.d0/3)
      endif
      if (y > 0.d0) then
        B=y**(1.d0/3)
      else
        B=-(-y)**(1.d0/3)
      endif
      root=A+B-a1/3
    else
      phi=acos(sqrt(q**2*27/(-4*p**3)))
      if (q < 0.d0) then
        arg1=phi/3
        arg2=phi/3+2*PI/3
        arg3=phi/3+4*PI/3
        y1=2*sqrt(-p/3)*cos(arg1)-a1/3
        y2=2*sqrt(-p/3)*cos(arg2)-a1/3
        y3=2*sqrt(-p/3)*cos(arg3)-a1/3
      elseif (q > 0.d0) then
        arg1=phi/3
        arg2=phi/3+2*PI/3
        arg3=phi/3+4*PI/3
        y1=-2*sqrt(-p/3)*cos(arg1)-a1/3
        y2=-2*sqrt(-p/3)*cos(arg2)-a1/3
        y3=-2*sqrt(-p/3)*cos(arg3)-a1/3
      endif
      if ((y3 > 0.d0).and.(y3 < upper)) then
        root=y3
      elseif ((y2 > 0.d0).and.(y2 < upper)) then
        root=y2
      elseif ((y1 > 0.d0).and.(y1 < upper)) then
        root=y1
      endif
    endif

  end subroutine CubeRoot

end module root_mod
