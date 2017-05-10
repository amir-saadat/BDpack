include ../../make.inc

GLBLIB  = ../../lib
GLBBIN  = ../../../bin

OBJS    = crrfcn.o
EXE	= crrfcn
INCS	= -I ../../common/inc
INCS   += -I ./inc
ifeq ($(FCOMP),gfortran)
MOD     = -J ./inc
else ifeq ($(FCOMP),ifort)
MOD     = -module ./inc
endif

CRRFCNLIBS = -L$(GLBLIB) -lcmn

all : makeinc CRRFCN


makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
	fi

CRRFCN:	$(OBJS)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(CRRFCNLIBS) $(GLBLIBS)
	$(ECHO) ">>> crrfcn compiled on `hostname -s` with $(COMP) <<<"
	if [ ! -d "$(GLBBIN)" ];then     \
                mkdir $(GLBBIN);         \
        fi
	$(MV) $(EXE) $(GLBBIN)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(MOD) -c $< -o $@

clean :
	$(ECHO) cleaning up in CRRFCN
	$(RM) -f $(OBJS) $(GLBBIN)/crrfcn
