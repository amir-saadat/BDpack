# put common definitions in here
#----------------------------------------------------------------------
# Intel and MKL Path
#----------------------------------------------------------------------
MKLROOT  = /opt/intel/mkl
PGIROOT  = /global/software/extras/pgi/18.4/linux86-64/2018
#----------------------------------------------------------------------
# mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
#----------------------------------------------------------------------
COMP     = mpif90 #$(INTELMPI)/mpiifort
#----------------------------------------------------------------------
# pgf90 or <Address of CUDA Fortran bin directory>/pgf90
#----------------------------------------------------------------------
USE_GPU  = TRUE
CUCOMP   = mpif90 #pgf90
#----------------------------------------------------------------------
# gfortran, ifort or <Address of Fortran compiler>
#----------------------------------------------------------------------
FCOMP    = pgf90
FEXT	 = f90
CUEXT	 = cuf
#MFLAGS	 = -s
ECHO	 = @echo
AR	     = ar
ARFLAGS	 =
RANLIB	 = ranlib
RM	     = rm
MV	     = mv

###INTELMPI = /opt/intel/impi/5.1.1.109/bin64
###MKLROOT  = /home/asaadat/intel/mkl
##CUDIR    = ./gpu
##PGIROOT	 = /home/asaadat/pgi
###----------------------------------------------------------------------
### mpif90, mpiifort or <Address of OpenMPI bin directory>/<MPI compiler>
###----------------------------------------------------------------------
##COMP     = mpiifort
###----------------------------------------------------------------------
### pgf90 or <Address of CUDA Fortran bin directory>/pgf90
###----------------------------------------------------------------------
##CUCOMP   = pgf90
###----------------------------------------------------------------------
### gfortran, ifort or <Address of Fortran compiler>
###----------------------------------------------------------------------
##FCOMP    = ifort
##FEXT     = f90
##CUEXT    = cuf
###MFLAGS   = -s
##ECHO     = @echo
##AR       = ar
##ARFLAGS  =
##RANLIB   = ranlib
##RM       = rm
##MV       = mv
#----------------------------------------------------------------------
# Uncomment for using <Address of MKL directory>
#----------------------------------------------------------------------

GLBFLGS   = -DUSE_DP # SP or DP for single or double precision
ifeq ($(USE_GPU),TRUE)
GLBFLGS	 += -DUSE_GPU
endif
# Extra debugging flags
#GLBFLGS += -g -traceback
#GLBFLGS += -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created
ifeq ($(FCOMP),gfortran)
###<<<<<<< HEAD
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
GLBFLGS += -check all -check bounds -traceback -g
GLBFLGs += -check-pointers=rw
else
ifeq ($(FCOMP),pgf90)
GLBFLGS	 += –Ktrap=fp –Minform=inform -g
GLBFLGS  += -lgomp -Mpreprocess -Mbackslash -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand
#GLBFLGS  += -Minfo -Mbounds -Mchkptr -Mchkstk
#GLBFLGS  += -Mbounds -Minfo=all -traceback -Mchkfpstk -Mchkstk -Mdalign
endif
endif
endif

CUFLGS   = -DUSE_GPU
CUFLGS  += -Mpreprocess -mp -Mcuda=charstring -Mcudalib=cublas,cusolver,cusparse,cufft,curand -Minfo
CUFLGS  += -Mbounds -Minfo=all -traceback -Mchkfpstk -Mchkstk -Mdalign -g
#CUFLGS  += -Mdclchk -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g -byteswapio
#CUFLGS   = -DADD_

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
GLBLIBS  = -L$(MKLROOT)/lib/intel64/
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
###=======
###GLBFLGS += -O0
###GLBFLGS += -cpp # Enables preprocessing
###GLBFLGS += -openmp # Enables openmp directives
###GLBFLGS += -fcray-pointer
####GLBFLGS += -fbounds-check
###else ifeq ($(FCOMP),ifort)
###GLBFLGS += -O0
###GLBFLGS += -fpp
###GLBFLGS += -qopenmp # Enables openmp directives
####GLBFLGS += -heap-arrays
####GLBFLGS += -check bounds -traceback -g
####GLBFLGs += -check-pointers=rw
###GLBFLGS += -fp-model precise
###GLBFLGS += -fp-model consistent
####GLBFLGS += -fp-model strict
###else ifeq ($(FCOMP),pgf90)
###GLBFLGS += -Mpreprocess
###GLBFLGS += -Mmpi=mpich
###endif
###
###CUFLGS   = -DUSE_CUBLAS
###
###GLBINCS  = -I/usr/include
###GLBINCS += -I$(MKLROOT)/include
###GLBINCS += -I$(MKLROOT)/include/intel64/lp64
####GLBINCS += -I$(PGIROOT)/mkl/include/intel64/lp64
###
###GLBLIBS  = -lm
###GLBLIBS += -L$(MKLROOT)/lib/intel64/
###GLBLIBS += -lmkl_intel_lp64 -lmkl_core
###GLBLIBS += -lmkl_sequential -lmkl_lapack95_lp64
###GLBLIBS += -lmkl_blas95_lp64 -lpthread
###>>>>>>> e839a20387ceb6596369d32f31af8b1541de5500

#GLBLIBS += -L$(PGIROOT)/lib/intel64/
#GLBLIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64

SHELL   = /bin/sh

#.SILENT :
