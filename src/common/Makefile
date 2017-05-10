include ../make.inc

OBJSLIB	= ../lib/libcmn.a
OBJS	= prcn_mod.o arry_mod.o strg_mod.o
OBJS   += cmn_io_mod.o mpi_mod.o
INCS	= -I ./inc
ifeq ($(FCOMP),gfortran)
MOD     = -J ./inc
else ifeq ($(FCOMP),ifort)
MOD     = -module ./inc
endif

all : makeinc $(OBJSLIB)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
	fi

../lib/libcmn.a : $(OBJS)
	$(ECHO) making the library: $(OBJSLIB)
	$(AR) $(ARFLAGS) rv $@ $?
	$(RANLIB) $@

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@


clean :
	$(ECHO) cleaning up in common
	$(RM) -f $(OBJS) $(OBJSLIB)
	if [ -n "$(ls -A ./inc)"]; then \
		$(RM) -f ./inc/*.mod;   \
	fi
