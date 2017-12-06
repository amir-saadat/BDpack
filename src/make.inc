# put common definitions in here
#----------------------------------------------------------------------
# Intel and MKL Path
#----------------------------------------------------------------------
#INTELMPI = /opt/intel/impi/5.1.1.109/bin64
#MKLROOT  = /home/asaadat/intel/mkl
CUDIR    = ./gpu
PGIROOT	 = /home/asaadat/pgi
#----------------------------------------------------------------------
# mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
#----------------------------------------------------------------------
COMP     = mpif90
#----------------------------------------------------------------------
# pgf90 or <Address of CUDA Fortran bin directory>/pgf90
#----------------------------------------------------------------------
CUCOMP   = pgf90
#----------------------------------------------------------------------
# gfortran, ifort or <Address of Fortran compiler>
#----------------------------------------------------------------------
FCOMP    = ifort
FEXT     = f90
CUEXT    = cuf
#MFLAGS   = -s
ECHO     = @echo
AR       = ar
ARFLAGS  =
RANLIB   = ranlib
RM       = rm
MV       = mv
#----------------------------------------------------------------------
# Uncomment for using <Address of MKL directory>
#----------------------------------------------------------------------

GLBFLGS  = -DUSE_DP # SP or DP for single or double precision
# Extra debugging flags
#GLBFLGS += -g -traceback
#GLBFLGS += -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created
ifeq ($(FCOMP),gfortran)
GLBFLGS += -Ofast
GLBFLGS += -cpp # Enables preprocessing
GLBFLGS += -openmp # Enables openmp directives
GLBFLGS += -fcray-pointer
#GLBFLGS += -fbounds-check
else ifeq ($(FCOMP),ifort)
GLBFLGS += -Ofast
GLBFLGS += -fpp
GLBFLGS += -qopenmp # Enables openmp directives
#GLBFLGS += -heap-arrays
#GLBFLGS += -check bounds -traceback -g
#GLBFLGs += -check-pointers=rw
#GLBFLGS += -fp-model precise
#GLBFLGS += -fp-model consistent
#GLBFLGS += -fp-model strict
else ifeq ($(FCOMP),pgf90)
GLBFLGS += -Mpreprocess
GLBFLGS += -Mmpi=mpich
endif

CUFLGS   = -DUSE_CUBLAS

GLBINCS  = -I/usr/include
GLBINCS += -I$(MKLROOT)/include
GLBINCS += -I$(MKLROOT)/include/intel64/lp64
#GLBINCS += -I$(PGIROOT)/mkl/include/intel64/lp64

GLBLIBS  = -lm
GLBLIBS += -L$(MKLROOT)/lib/intel64/
GLBLIBS += -lmkl_intel_lp64 -lmkl_core
GLBLIBS += -lmkl_sequential -lmkl_lapack95_lp64
GLBLIBS += -lmkl_blas95_lp64 -lpthread

#GLBLIBS += -L$(PGIROOT)/lib/intel64/
#GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64

SHELL   = /bin/sh

.SILENT :
