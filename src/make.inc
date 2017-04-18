# put common definitions in here
#----------------------------------------------------------------------
# mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
#----------------------------------------------------------------------
COMP     = mpif90
#----------------------------------------------------------------------
# gfortran, ifort or <Address of Fortran compiler>
#----------------------------------------------------------------------
FCOMP    = ifort
FEXT	 = f90
MFLAGS	 = -s
ECHO	 = @echo
AR	 = ar
ARFLAGS	 =
RANLIB	 = ranlib
RM	 = rm
MV	 = mv
#----------------------------------------------------------------------
# Uncomment for using <Address of MKL directory>
#----------------------------------------------------------------------

#MKLROOT =


GLBFLGS  = -fpp # Enables preprocessing
GLBFLGS += -qopenmp # Enables openmp directives
GLBFLGS += -DUSE_DP # SP or DP for single or double precision
GLBFLGS	+= -Ofast -traceback
ifeq ($(FCOMP),gfortran)
GLBFLGS += -fcray-pointer
GLBFLGS += -fbounds-check
else ifeq ($(FCOMP),ifort)
GLBFLGS += -heap-arrays 10
GLBFLGS += -check bounds
endif


GLBINCS  = -I/usr/include
GLBINCS	+= -I$(MKLROOT)/include
GLBINCS	+= -I$(MKLROOT)/include/intel64/lp64

GLBLIBS	 = -lm
GLBLIBS += -L$(MKLROOT)/lib/intel64/
GLBLIBS += -lmkl_intel_lp64 -lmkl_core
GLBLIBS += -lmkl_sequential -lmkl_lapack95_lp64 
GLBLIBS += -lmkl_blas95_lp64 -lpthread

SHELL	= /bin/sh

.SILENT :
