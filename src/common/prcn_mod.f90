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
module prcn_mod

  use :: blas95, only : dot,nrm2,scal,copy,axpy,symv,gemv,gbmv,trmv,symm,gemm,trmm
  use :: lapack95, only : geqrf,orgqr,syev,syevr,potrf,potri,getrf,getri,geev,sysv,pocon,sycon,sytrf
#ifdef USE_DP
  use :: f95_precision, only : wp => dp
#elif USE_SP
  use :: f95_precision, only : wp => sp
#endif
  !for higher precision (real(real128))
  use :: iso_fortran_env

  implicit none
  save

  include 'mkl_spblas_crr.fi'

  integer,parameter :: single=selected_real_kind(p=6,r=37)
#ifdef USE_GPU
  integer,parameter :: double=wp
#else
  integer,parameter :: double=selected_real_kind(p=15,r=307)
#endif
!  integer,parameter :: wp=selected_real_kind(p=15,r=307)
  integer,parameter :: short=selected_int_kind(4)
  integer,parameter :: long=selected_int_kind(9)

end module prcn_mod
