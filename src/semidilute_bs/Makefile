include ../make.inc

OBJLIBS	 = ../libsmdlt_bs.a
CMNOBJS += ../common/arry_mod.o 
CMNOBJS += ../common/prcn_mod.o 
CMNOBJS += ../common/strg_mod.o
OBJS     = mkl_dfti.o types.o flow_mod.o trsfm_mod.o hi_mod.o
OBJS    += tmng_mod.o diffcalc_mod.o dcmp_smdlt.o conv_mod.o
OBJS    += verlet_mod.o force_smdlt.o sprforce_mod.o evforce_mod.o
OBJS    += diffdcmp_mod.o io_mod.o chain_mod.o pp_smdlt.o box_mod.o
OBJS    += smdlt_mod.o
INCS	 = -I ../common/inc
INCS    += -I ./inc

ifeq ($(FCOMP),gfortran)
MOD      = -J ./inc
else ifeq ($(FCOMP),ifort)
MOD      = -module ./inc
else ifeq ($(FCOMP),tau_f90.sh)
MOD     = -module ./inc
endif

all : makeinc $(OBJLIBS)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
        fi

../libsmdlt_bs.a : $(OBJS) $(CMNOBJS)
	$(ECHO) making the library: $(OBJLIBS)
	$(AR) $(ARFLAGS) rv ../libsmdlt_bs.a $?
	$(RANLIB) ../libsmdlt_bs.a

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@
clean :
	$(ECHO) cleaning up in semidilute_bs
	$(RM) -f $(OBJS)
	if [ "$(ls -A ./inc)" ];then     \
		$(RM) ./inc/*.mod;       \
	fi
