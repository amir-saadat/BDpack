include ../make.inc

CMNDIR	 = ../common
CUDIR	 = ./cuda
OBJLIB	 = ../libdlt_bs.a
CUOBJLIB = ../libdlt_bs_cuda.a
CMNOBJS += $(CMNDIR)/prcn_mod.o 
CMNOBJS += $(CMNDIR)/cmn_tp_mod.o
CMNOBJS += $(CMNDIR)/arry_mod.o 
CMNOBJS += $(CMNDIR)/strg_mod.o
CMNOBJS += $(CMNDIR)/cmn_io_mod.o
CMNOBJS += $(CMNDIR)/mpi_mod.o
CMNOBJS += $(CMNDIR)/rand_mod.o
CUOBJS	 = $(CUDIR)/test_cumod.o
CUOBJS	+= $(CUDIR)/intrn_cumod.o $(CUDIR)/dcmp_dlt_cumod.o
CUOBJS	+= $(CUDIR)/mmul_cumod.o $(CUDIR)/gpu_cumod.o
OBJS	 = inp_dlt.o dcmp_mod.o root_mod.o 
OBJS    += force_mod.o HIEVCalc.o pp_mod.o 
#OBJS    += intrn_mod.o hi_smod.o hibb_smod.o 
#OBJS	+= hibw_smod.o evbb_smod.o evbw_smod.o 
OBJS    += hibb_smod.o hibw_smod.o hi_smod.o
OBJS	+= evbb_smod.o evbw_smod.o intrn_mod.o
OBJS	+= sde_mod.o dlt_mod.o
INCS	 = -I ../common/inc
INCS    += -I ./inc
INCS	+= -I ./cuda/inc
#CUINCS   = -I ./inc/gpu
ifeq ($(FCOMP),gfortran)
MOD      = -J ./inc
else 
ifeq ($(FCOMP),ifort)
MOD      = -module ./inc
else
ifeq ($(FCOMP),pgf90)
MOD      = -module ./inc
CUMOD    = -module ./inc
endif
endif
endif

#all : makeinc $(CUOBJLIB) $(OBJLIB)
all : makeinc $(OBJLIB)

makeinc:
	if [ ! -d "./inc" ];then \
		mkdir ./inc;     \
	fi

ifeq ($(USE_GPU),TRUE)
../libdlt_bs.a : $(CUOBJS) $(OBJS) $(CMNOBJS)
else
../libdlt_bs.a : $(OBJS) $(CMNOBJS)
endif
	$(ECHO) making the library: $(OBJLIB)
	$(AR) $(ARFLAGS) rv $(OBJLIB) $?
	$(RANLIB) $(OBJLIB)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

$(CUOBJS): %.o: %.$(CUEXT)
	$(ECHO) compiling $< using $(CUCOMP)
	$(CUCOMP) $(CUFLGS) $(CUINCS) $(CUMOD) -c $< -o $@

clean :
	$(ECHO) cleaning up in dilute_bs
	$(RM) -f $(OBJS) $(CUOBJS)
	$(RM) -f $(OBJS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done
	$(RM) -f ./inc/*.mod;
#ifneq ($(shell ls -A ./inc),)
#endif

force_look :
	true
