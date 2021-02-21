# put common definitions in here
#----------------------------------------------------------------------
# Intel and MKL Path
#----------------------------------------------------------------------
PGIROOT = /opt/packages/pgi/20.11/Linux_x86_64/20.11/
MPIINC = $(PGIROOT)/comm_libs/mpi/lib/
#----------------------------------------------------------------------
# mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
#----------------------------------------------------------------------
COMP     = mpif90 #$(INTELMPI)/mpiifort
#----------------------------------------------------------------------
# pgf90 or <Address of CUDA Fortran bin directory>/pgf90
#----------------------------------------------------------------------
USE_GPU  = TRUE
CUCOMP   = pgf90
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
# Compiler flags based on the compiler and if CPU or GPU is used
#----------------------------------------------------------------------
GLBFLGS   = -DUSE_DP # SP or DP for single or double precision
ifeq ($(USE_GPU),TRUE)
GLBFLGS	 += -DUSE_GPU
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
else
ifeq ($(FCOMP),pgf90)
GLBFLGS	 += -g
GLBFLGS  += -lgomp -Mpreprocess -Mbackslash -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand
#GLBFLGS  += -Minfo -Mbounds -Mchkptr -Mchkstk
#GLBFLGS  += -Mbounds -Minfo=all -traceback -Mchkfpstk -Mchkstk -Mdalign
endif
endif
endif

CUFLGS   = -DUSE_GPU
CUFLGS  += -Mpreprocess -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand -Minfo
CUFLGS  += -Mbounds -Minfo=all -traceback -Mchkstk -Mdalign -g
#CUFLGS  += -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
#CUFLGS   = -DADD_

GLBINCS  = -I/usr/include
GLBINCS	+= -I$(MKLROOT)/include
ifeq ($(FCOMP),ifort)
GLBINCS	+= -I$(MKLROOT)/include/intel64/lp64
else
ifeq ($(FCOMP),pgf90)
GLBINCS	+= -I ~/mkl_pgi/include/intel64/lp64
endif
endif

CUINCS = -I$(MPIINC)

ifeq ($(FCOMP),ifort)
GLBLIBS  = -L$(MKLROOT)/lib/intel64/
GLBLIBS += -lmkl_intel_lp64 -lmkl_core
GLBLIBS += -lmkl_sequential -lmkl_lapack95_lp64
GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64
else
ifeq ($(FCOMP),pgf90)
GLBLIBS += -L$(PGIROOT)/lib
GLBLIBS += -lcudafor -lcudaforblas -lcublas_static -lcudart
GLBLIBS += -L ~/mkl_pgi/lib/intel64/
GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64
endif
endif
GLBLIBS += -L$(MKLROOT)/lib/intel64/
GLBLIBS += -lmkl_intel_lp64 -lmkl_core
GLBLIBS += -lmkl_sequential -lpthread

SHELL   = /bin/sh

#.SILENT :
