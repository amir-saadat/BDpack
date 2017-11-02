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
module arry_mod

  use :: prcn_mod

  implicit none
  save
  integer,parameter :: idp = selected_int_kind(13)
  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)

  public :: print_vector
  interface print_vector
     module procedure print_vector_int
     module procedure print_vector_sp
     module procedure print_vector_dp
     module procedure print_vector_cmplxsp
     module procedure print_vector_cmplxdp
  end interface print_vector

  public :: print_matrix
  interface print_matrix
     module procedure print_matrix_int
     module procedure print_matrix_sp
     module procedure print_matrix_dp
  end interface print_matrix

  public :: print_spmatrix
  interface print_spmatrix
     module procedure print_spmatrix_dp
  end interface print_spmatrix

  public :: ResizeArray
  interface ResizeArray
     module procedure ResizeArray_int
     module procedure ResizeArray_sp
     module procedure ResizeArray_dp
  end interface ResizeArray

  public :: linspace
  interface linspace
     module procedure linspace_sp
     module procedure linspace_dp
  end interface linspace

  public :: logspace
  interface logspace
     module procedure logspace_sp
     module procedure logspace_dp
  end interface logspace

  public :: sort
  interface sort
     module procedure sort_sp
     module procedure sort_all_sp
     module procedure sort_dp
     module procedure sort_all_dp
     module procedure sort_int
     module procedure sort_all_int
  end interface sort

contains

  subroutine ResizeArray_int(A,newSize)

    integer,dimension(:),intent(inout),pointer :: A
    integer,intent(in) :: newSize
    integer,dimension(:),allocatable :: B

    allocate(B(lbound(A,1):ubound(A,1)))
    B=A
    !   Watch out for memory leak. If A in the original callee points to another target, that
    !    target will be lost.
    if (associated(A)) deallocate(A)
    allocate(A(newSize))
    A(lbound(B,1):ubound(B,1))=B
    deallocate(B)
  end subroutine ResizeArray_int

  subroutine ResizeArray_sp(A,newSize)

    real(sp),dimension(:),intent(inout),pointer :: A
    integer,intent(in) :: newSize
    real(sp),dimension(:),allocatable :: B

    allocate(B(lbound(A,1):ubound(A,1)))
    B=A
    !   Watch out for memory leak. If A in the original callee points to another target, that
    !    target will be lost.
    if (associated(A)) deallocate(A)
    allocate(A(newSize))
    A(lbound(B,1):ubound(B,1))=B
    deallocate(B)
  end subroutine ResizeArray_sp

  subroutine ResizeArray_dp(A,newSize)

    real(dp),dimension(:),intent(inout),pointer :: A
    integer,intent(in) :: newSize
    real(dp),dimension(:),allocatable :: B

    allocate(B(lbound(A,1):ubound(A,1)))
    B=A
    !   Watch out for memory leak. If A in the original callee points to another target, that
    !    target will be lost.
    if (associated(A)) deallocate(A)
    allocate(A(newSize))
    A(lbound(B,1):ubound(B,1))=B
    deallocate(B)
  end subroutine ResizeArray_dp

  subroutine linspace_dp(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    integer :: i,n
    n = size(x)
    if (n == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,n
          x(i) = (xmax-xmin) * real(i-1,dp) / real(n-1,dp) + xmin
       end do
    end if
  end subroutine linspace_dp

  subroutine logspace_dp(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    if (size(x) == 1 .and. xmin /= xmax) then
       write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       stop
    end if
    call linspace(log10(xmin),log10(xmax),x)
    x = 10._dp**x
  end subroutine logspace_dp

  subroutine linspace_sp(xmin,xmax,x)

    real(sp),intent(in) :: xmin,xmax
    real(sp),intent(out) :: x(:)
    integer :: i,n
    n = size(x)
    if (n == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,n
          x(i) = (xmax-xmin) * real(i-1,sp) / real(n-1,sp) + xmin
       end do
    end if
  end subroutine linspace_sp

  subroutine logspace_sp(xmin,xmax,x)

    real(sp),intent(in) :: xmin,xmax
    real(sp),intent(out) :: x(:)
    if (size(x) == 1 .and. xmin /= xmax) then
       write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       stop
    end if
    call linspace(log10(xmin),log10(xmax),x)
    x = 10._sp**x
  end subroutine logspace_sp

  subroutine print_matrix_int(a,myname)

    integer :: nx,ny
    integer,dimension(:,:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    character(len=20) :: fmt

    nx = size(a,1)
    ny = size(a,2)
    write(fmt,'("(",I3.3,"(i10,1X))")') ny
    print *, myname
    do i=1,nx
       write(*,fmt) a(i,:)
    end do

  end subroutine print_matrix_int

  subroutine print_matrix_dp(a,myname)

    integer :: nx,ny
    real(dp),dimension(:,:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    character(len=20) :: fmt

    nx = size(a,1)
    ny = size(a,2)
    write(fmt,'("(",I3.3,"(f20.10,1x))")') ny
    print *, myname
    do i=1,nx
       write(*,fmt) a(i,:)
    end do

  end subroutine print_matrix_dp

  subroutine print_matrix_sp(a,myname)

    implicit none
    integer :: nx,ny
    real(sp),dimension(:,:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    character(len=20) :: fmt

    nx = size(a,1)
    ny = size(a,2)
    write(fmt,'("(",I3.3,"(f20.8,1X))")') ny
    print *, myname
    do i=1,nx
       write(*,fmt) a(i,:)
    end do

  end subroutine print_matrix_sp

  subroutine print_vector_int(a,myname)

    integer :: nx
    integer,dimension(:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    nx = size(a,1)
    print *, myname
    do i=1,nx
       write(*,'(i4)') a(i)
    end do
  end subroutine print_vector_int

  subroutine print_vector_sp(a,myname)

    integer :: nx
    real(sp),dimension(:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    nx = size(a,1)
    print *, myname
    do i=1,nx
       write(*,'(f8.5)') a(i)
    end do
  end subroutine print_vector_sp

  subroutine print_vector_dp(a,myname)

    integer :: nx
    real(dp),dimension(:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    nx = size(a,1)
    print *, myname
    do i=1,nx
       write(*,'(f24.14)') a(i)
    end do
  end subroutine print_vector_dp

  subroutine print_vector_cmplxsp(a,myname)

    integer :: nx
    complex(sp),dimension(:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    nx = size(a,1)
    print *, myname
    do i=1,nx
       write(*,*) a(i)
    end do
  end subroutine print_vector_cmplxsp

  subroutine print_vector_cmplxdp(a,myname)

    integer :: nx
    complex(dp),dimension(:),intent(in) :: a
    character(len=*) :: myname
    integer :: i
    nx = size(a,1)
    print *, myname
    do i=1,nx
       write(*,*) a(i)
    end do
  end subroutine print_vector_cmplxdp

  subroutine print_spmatrix_dp(adns,indexing,pattern,m,n,maxnz)

    real(dp),dimension(:,:),intent(in) :: adns
    real(dp),allocatable,dimension(:) :: acsr
    integer,allocatable,dimension(:) :: ja
    ! pattern ~ 0: Lower Triangular, 1: Upper Triangular, 2: Whole matrix
    integer :: job(8),ia(m+1),m,n,lda,indexing,pattern,maxnz,info

    allocate(acsr(maxnz),ja(maxnz))
    lda=m
    job=(/0,1,indexing,pattern,maxnz,1,0,0/)

    call mkl_ddnscsr(job,m,n,adns,lda,acsr,ja,ia,info)

    if (info.eq.0) then
      print *,'DNS -> CSR successful'
      call print_vector(acsr,'acsr')
      call print_vector(ja,'columns')
      call print_vector(ia,'rowIndex')
    elseif (info.gt.0) then
      print *,'Conversion interrupted in ith row;'
      print *,'i:',info
    end if

    deallocate(acsr,ja)

  end subroutine print_spmatrix_dp

  subroutine swap_sp(array, i, j)
    real(sp),intent(inout) :: array(:)
    integer,intent(in) :: i, j
    real(sp) :: temp
    temp = array(j)
    array(j) = array(i)
    array(i) = temp
  end subroutine swap_sp

  subroutine partition_sp(array, left, right, pivot_index, store_index)
    real(sp),intent(inout) :: array(:)
    integer,intent(in) :: left, right, pivot_index
    integer,intent(out) :: store_index
    real(sp) :: pivot_value
    integer :: i
    pivot_value = array(pivot_index)
    call swap_sp(array, pivot_index, right)
    store_index = left
    do i=left, right-1
       if(array(i) <= pivot_value) then
          call swap_sp(array, i, store_index)
          store_index = store_index + 1
       end if
    end do
    call swap_sp(array, store_index, right)
  end subroutine partition_sp

  recursive subroutine sort_sp(array, left, right)
    real(sp),intent(inout) :: array(:)
    integer,intent(in) :: left, right
    integer :: pivot_index, new_pivot_index
    if(right > left) then
       pivot_index = left+(right-left)/2
       call partition_sp(array, left, right, pivot_index, new_pivot_index)
       call sort_sp(array, left, new_pivot_index - 1)
       call sort_sp(array, new_pivot_index + 1, right)
    end if
  end subroutine sort_sp

  recursive subroutine sort_all_sp(array)
    real(sp),intent(inout) :: array(:)
    call sort_sp(array, 1, size(array))
  end subroutine sort_all_sp

  subroutine swap_dp(array, i, j)
    implicit none
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: i, j
    real(dp) :: temp
    temp = array(j)
    array(j) = array(i)
    array(i) = temp
  end subroutine swap_dp

  subroutine partition_dp(array, left, right, pivot_index, store_index)
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: left, right, pivot_index
    integer,intent(out) :: store_index
    real(dp) :: pivot_value
    integer :: i
    pivot_value = array(pivot_index)
    call swap_dp(array, pivot_index, right)
    store_index = left
    do i=left, right-1
       if(array(i) <= pivot_value) then
          call swap_dp(array, i, store_index)
          store_index = store_index + 1
       end if
    end do
    call swap_dp(array, store_index, right)
  end subroutine partition_dp

  recursive subroutine sort_dp(array, left, right)
    real(dp),intent(inout) :: array(:)
    integer,intent(in) :: left, right
    integer :: pivot_index, new_pivot_index
    if(right > left) then
       pivot_index = left+(right-left)/2
       call partition_dp(array, left, right, pivot_index, new_pivot_index)
       call sort_dp(array, left, new_pivot_index - 1)
       call sort_dp(array, new_pivot_index + 1, right)
    end if
  end subroutine sort_dp

  recursive subroutine sort_all_dp(array)
    real(dp),intent(inout) :: array(:)
    call sort_dp(array, 1, size(array))
  end subroutine sort_all_dp

  subroutine swap_int(array, i, j)
    implicit none
    integer,intent(inout) :: array(:)
    integer,intent(in) :: i, j
    integer :: temp
    temp = array(j)
    array(j) = array(i)
    array(i) = temp
  end subroutine swap_int

  subroutine partition_int(array, left, right, pivot_index, store_index)
    integer,intent(inout) :: array(:)
    integer,intent(in) :: left, right, pivot_index
    integer,intent(out) :: store_index
    integer :: pivot_value
    integer :: i
    pivot_value = array(pivot_index)
    call swap_int(array, pivot_index, right)
    store_index = left
    do i=left, right-1
       if(array(i) <= pivot_value) then
          call swap_int(array, i, store_index)
          store_index = store_index + 1
       end if
    end do
    call swap_int(array, store_index, right)
  end subroutine partition_int

  recursive subroutine sort_int(array, left, right)
    integer,intent(inout) :: array(:)
    integer,intent(in) :: left, right
    integer :: pivot_index, new_pivot_index
    if(right > left) then
       pivot_index = left+(right-left)/2
       call partition_int(array, left, right, pivot_index, new_pivot_index)
       call sort_int(array, left, new_pivot_index - 1)
       call sort_int(array, new_pivot_index + 1, right)
    end if
  end subroutine sort_int

  recursive subroutine sort_all_int(array)
    integer,intent(inout) :: array(:)
    call sort_int(array, 1, size(array))
  end subroutine sort_all_int

end module arry_mod
