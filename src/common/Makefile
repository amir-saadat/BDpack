include ../make.inc

OBJS	= prcn_mod.o arry_mod.o strg_mod.o
OBJS   += cmn_io_mod.o 
INCS	= -I ./inc
ifeq ($(FCOMP),gfortran)
MOD     = -J ./inc
else ifeq ($(FCOMP),ifort)
MOD     = -module ./inc
endif

all : makeinc $(OBJS)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
	fi

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

clean :
	$(ECHO) cleaning up in common
	$(RM) -f $(OBJS)
	if [ "$(ls -A ./inc)" ];then     \
		$(RM) ./inc/*.mod;       \
	fi
