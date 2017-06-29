include ../../make.inc

GLBLIB  = ../../lib
GLBBIN  = ../../../bin

OBJS    = wallcollision.o
EXE	= wallcollision
INCS	= -I ../../common/inc
INCS   += -I ./inc
ifeq ($(FCOMP),gfortran)
MOD     = -J ./inc
else ifeq ($(FCOMP),ifort)
MOD     = -module ./inc
endif

WALLCOLLISIONLIBS = -L$(GLBLIB) -lcmn

all : makeinc cmnobjs WALLCOLLISION

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
	fi

WALLCOLLISION:	$(OBJS)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(WALLCOLLISIONLIBS) $(GLBLIBS)
	if [ ! -d "$(GLBBIN)" ];then     \
                mkdir $(GLBBIN);         \
        fi
	$(MV) $(EXE) $(GLBBIN)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

cmnobjs : force_look
	$(ECHO) looking into common : $(MAKE) $(MFLAGS)
	cd ../../common; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in WALLCOLLISION
	$(RM) -f $(OBJS) $(GLBBIN)/wallcollision

force_look :
	true
