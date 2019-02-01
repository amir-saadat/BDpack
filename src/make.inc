# put common definitions in here
#----------------------------------------------------------------------
# Intel and MKL Path
# Uncomment for using <Address of MKL directory>
#----------------------------------------------------------------------
#MKLROOT  = /opt/intel/mkl
PGIROOT  = /share/software/user/open/pgi/linux86-64/18.4
#----------------------------------------------------------------------
# mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
#----------------------------------------------------------------------
COMP     = mpif90 #$(INTELMPI)/mpiifort
#----------------------------------------------------------------------
# pgf90 or <Address of CUDA Fortran bin directory>/pgf90
#----------------------------------------------------------------------
USE_GPU  = FALSE
CUCOMP   = mpif90 #pgf90
#----------------------------------------------------------------------
# gfortran, ifort or <Address of Fortran compiler>
#----------------------------------------------------------------------
FCOMP    = pgf90
FEXT	 = f90
CUEXT	 = cuf
#MFLAGS	 = -s
ECHO	 = @echo
AR	 = ar
ARFLAGS	 =
RANLIB	 = ranlib
RM	 = rm
MV	 = mv
#----------------------------------------------------------------------
# Global FlagsLibraries and include files
#----------------------------------------------------------------------
GLBFLGS   = -DUSE_DP # SP or DP for single or double precision
ifeq ($(USE_GPU),TRUE)
GLBFLGS	 += -DUSE_GPU
GLBFLGS  += -ta=tesla:cc35,cc60,cc70
endif
# Extra debugging flags
#GLBFLGS += -g -traceback
#GLBFLGS += -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created
ifeq ($(FCOMP),gfortran)
GLBFLGS	+= -Ofast
GLBFLGS += -cpp # Enables preprocessing
GLBFLGS += -openmp # Enables openmp directives
# Extra debugging flags
#GLBFLGS += -fcray-pointer
#GLBFLGS += -fbounds-check
else 
ifeq ($(FCOMP),ifort)
GLBFLGS	+= -Ofast -g
GLBFLGS += -cpp # Enables preprocessing
GLBFLGS += -qopenmp # Enables openmp directives
# Extra debugging flags
#GLBFLGS += -heap-arrays
#GLBFLGS += -check all -check bounds -traceback -g
#GLBFLGs += -check-pointers=rw
#GLBFLGS += -fp-model precise
#GLBFLGS += -fp-model consistent
#GLBFLGS += -fp-model strict
else
ifeq ($(FCOMP),pgf90)
GLBFLGS	 += –Ktrap=fp –Minform=inform -g
GLBFLGS  += -lgomp -Mpreprocess -Mbackslash -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand 
#GLBFLGS  += -ta=tesla:cc70
#GLBFLGS  += -Minfo -Mbounds -Mchkptr -Mchkstk
#GLBFLGS  += -Mbounds -Minfo=all -traceback -Mchkfpstk -Mchkstk -Mdalign
#GLBFLGS += -Mmpi=mpich
endif
endif
endif
CUFLGS   = -DUSE_GPU
CUFLGS  += -Mpreprocess -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand -Minfo=all
CUFLGS  += -ta=tesla:cc35,cc60,cc70
#CUFLGS  += -traceback -Mbounds -Mchkfpstk -Mchkstk -Mdalign -g
#CUFLGS  += -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
#CUFLGS   = -DADD_
#----------------------------------------------------------------------
# MKL Libraries and Include Path
#----------------------------------------------------------------------
GLBINCS  = -I/usr/include
GLBINCS	+= -I$(MKLROOT)/include
ifeq ($(FCOMP),ifort)
GLBINCS	+= -I$(MKLROOT)/include/intel64/lp64
else
ifeq ($(FCOMP),pgf90)
GLBINCS	+= -I ~/mkl_pgi17/include/intel64/lp64
endif
endif

CUINCS  = 
#CUINCS  = -I /global/software/extras/pgi/17.4/linux86-64/2017/mpi/openmpi/lib

ifeq ($(FCOMP),ifort)
GLBLIBS  = -L$(MKLROOT)/lib/intel64/ -liomp5
GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64
else
ifeq ($(FCOMP),pgf90)
GLBLIBS += -L$(PGIROOT)/lib -L$(PGIROOT)/cuda/9.1/lib64
GLBLIBS += -lcudafor -lcudaforblas -lcublas_static -lcudart
GLBLIBS += -L ~/mkl_pgi17/lib/intel64/
GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64
endif
endif
GLBLIBS += -L$(MKLROOT)/lib/intel64/
GLBLIBS += -lmkl_intel_lp64 -lmkl_core
GLBLIBS += -lmkl_sequential -lpthread

SHELL   = /bin/sh

#.SILENT :
