include ../make.inc

CMNDIR	 = ../common
CMNCUDIR = ../common/cuda
CUDIR	 = ./cuda
OBJLIBS	 = ../lib/libsmdlt_bs.a
CMNCUOBJS += $(CMNCUDIR)/dev_cumod.o 
ifeq ($(USE_MAGMA),TRUE)
CMNCUOBJS += $(CMNCUDIR)/magma_cumod.o 
endif
CMNOBJS += $(CMNDIR)/arry_mod.o 
CMNOBJS += $(CMNDIR)/prcn_mod.o 
CMNOBJS += $(CMNDIR)/strg_mod.o
OBJS_F   = omp_lib.o mkl_dfti.o types.o flow_mod.o trsfm_mod.o
OBJS_F  += hi_mod.o tmng_mod.o verlet_mod.o
OBJS_M   = diffcalc_mod.o dcmp_smdlt.o conv_mod.o
OBJS_M  += evverlet_mod.o force_smdlt.o sprforce_mod.o 
OBJS_M	+= evforce_mod.o diffdcmp_mod.o io_mod.o chain_mod.o
OBJS_L   = pp_smdlt.o box_mod.o smdlt_mod.o
CUOBJS_F  = $(CUDIR)/magma_smdlt.o
CUOBJS_F += $(CUDIR)/rndm_cumod.o
CUOBJS_F += $(CUDIR)/flow_cumod.o
CUOBJS_F += $(CUDIR)/trsfm_cumod.o
CUOBJS_F += $(CUDIR)/hiverlet_cumod.o
CUOBJS_F += $(CUDIR)/hi_cumod.o
CUOBJS_F += $(CUDIR)/diffcalc_cumod.o
CUOBJS_F += $(CUDIR)/dcmp_cumod.o
CUOBJS_F += $(CUDIR)/diffdcmp_cumod.o
CUOBJS_M  = $(CUDIR)/cusparse_wrap.o
CUOBJS_M += $(CUDIR)/conv_cumod.o
CUOBJS_M += $(CUDIR)/evverlet_cumod.o
CUOBJS_M += $(CUDIR)/force_cumod.o
CUOBJS_M += $(CUDIR)/sprforce_cumod.o
CUOBJS_M += $(CUDIR)/evforce_cumod.o
CUOBJS_M += $(CUDIR)/chain_cumod.o
#CUOBJS  += $(CUDIR)/diffcalc_mod.o
INCS	 = -I ../common/inc
INCS    += -I ./inc

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

all : makeinc $(OBJLIBS)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
        fi

ifeq ($(USE_GPU),TRUE)
../lib/libsmdlt_bs.a : $(OBJS_F) $(CUOBJS_F) $(OBJS_M) $(CUOBJS_M) $(OBJS_L) $(CMNOBJS) $(CMNCUOBJS)
else
../lib/libsmdlt_bs.a : $(OBJS_F) $(OBJS_M) $(OBJS_L) $(CMNOBJS)
endif
	$(ECHO) making the library: $(OBJLIBS)
	$(AR) $(ARFLAGS) rv $(OBJLIBS) $?
	$(RANLIB) $(OBJLIBS)

$(OBJS_F): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

$(OBJS_M): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

$(OBJS_L): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

$(CUOBJS_F): %.o: %.$(CUEXT)
	$(ECHO) compiling $< using $(CUCOMP)
	$(CUCOMP) $(CUFLGS) $(GLBINCS) $(INCS) $(CUINCS) $(CUMOD) -c $< -o $@

$(CUOBJS_M): %.o: %.$(CUEXT)
	$(ECHO) compiling $< using $(CUCOMP)
	$(CUCOMP) $(CUFLGS) $(GLBINCS) $(INCS) $(CUINCS) $(CUMOD) -c $< -o $@

clean :
	$(ECHO) cleaning up in semidilute_bs
	$(RM) -f $(OBJS_F) $(OBJS_M) $(OBJS_L) $(CUOBJS_F) $(CUOBJS_M)
	$(RM) -f ./inc/*.mod;
	$(RM) -f ./*.mod;
