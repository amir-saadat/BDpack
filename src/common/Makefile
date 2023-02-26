include ../make.inc

CUDIR	= ./cuda
OBJSLIB	= ../lib/libcmn.a
OBJS	= prcn_mod.o cmn_tp_mod.o arry_mod.o
OBJS   += strg_mod.o cmn_io_mod.o mpi_mod.o 
OBJS   += rand_mod.o
CUOBJS	= $(CUDIR)/dev_cumod.o
ifeq ($(USE_MAGMA),TRUE)
CUOBJS += $(CUDIR)/magma_cumod.o
endif
INCS	= -I ./inc
ifeq ($(FCOMP),gfortran)
MOD      = -J ./inc
else 
ifeq ($(FCOMP),$(filter $(FCOMP),ifort ifx))
MOD      = -module ./inc
else
ifeq ($(FCOMP),pgf90)
MOD      = -module ./inc
CUMOD    = -module ./inc
endif
endif
endif

all : makeinc $(OBJSLIB)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
	fi

ifeq ($(USE_GPU),TRUE)
../lib/libcmn.a : $(CUOBJS) $(OBJS)
else
../lib/libcmn.a : $(OBJS)
endif
	$(ECHO) making the library hello: $(OBJSLIB)
	$(AR) $(ARFLAGS) rv $@ $?
	$(RANLIB) $@

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

$(CUOBJS): %.o: %.$(CUEXT)
	$(ECHO) compiling $< using $(CUCOMP)
	$(CUCOMP) $(CUFLGS) $(CUINCS) $(CUMOD) -c $< -o $@

clean :
	$(ECHO) cleaning up in common
	$(RM) -f $(OBJS) $(CUOBJS) $(OBJSLIB)
	$(RM) -f ./inc/*.mod;
