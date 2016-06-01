# put common definitions in here

COMP     = mpif90
FCOMP    = ifort
FEXT	 = f90
MFLAGS	 = -s
ECHO	 = @echo
AR	 = ar
ARFLAGS	 =
RANLIB	 = ranlib
RM	 = rm
MV	 = mv
#MKLROOT  =

GLBFLGS	 = -Ofast

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
