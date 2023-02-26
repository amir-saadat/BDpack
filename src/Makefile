include make.inc

DIRS	 = common dilute_bs semidilute_bs
EXE	 = BDpack
OBJS	 = main.o
OBJLIBS	 = libdlt_bs.a libsmdlt_bs.a
LIBS	 = -L./lib -ldlt_bs -lsmdlt_bs
INCS	 = -I ./common/inc
INCS	+= -I ./dilute_bs/inc
INCS	+= -I ./dilute_bs/cuda/inc
INCS	+= -I ./semidilute_bs/inc

UTILSDIR = ./utils
UTILS    = pdf crrfcn cnfgen wallcollision

all : $(EXE)

#$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(LIBS) $(GLBLIBS)
$(EXE) : cmnobjs $(OBJLIBS) $(OBJS)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(LIBS) $(GLBLIBS)
	$(ECHO) ">>> BDpack compiled on `hostname -s` with $(CUCOMP) <<<"
	if [ ! -d "../bin" ];then     \
                mkdir ../bin;         \
        fi
	$(MV) $(EXE) ../bin
#	$(RM) -f $(OBJLIBS) $(OBJS)
#	$(MAKE) clean

cmnobjs : force_look
	$(ECHO) looking into common : $(MAKE) $(MFLAGS)
	cd common; $(MAKE) $(MFLAGS)

libdlt_bs.a : force_look
	$(ECHO) looking into dilute_bs : $(MAKE) $(MFLAGS)
	cd dilute_bs; $(MAKE) $(MFLAGS)

libsmdlt_bs.a : force_look
	$(ECHO) looking into semidilute_bs : $(MAKE) $(MFLAGS)
	cd semidilute_bs; $(MAKE) $(MFLAGS)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) -c $< -o $@

utils: $(UTILS)

pdf: force_look
	$(ECHO) looking into pdf : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/pdf; $(MAKE) $(MFLAGS)

crrfcn: force_look
	$(ECHO) looking into crrfcn : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/crrfcn; $(MAKE) $(MFLAGS)

cnfgen: force_look
	$(ECHO) looking into cnfgen : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/cnfgen; $(MAKE) $(MFLAGS)
	
wallcollision: force_look
	$(ECHO) looking into wallcollision : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/wallcollision; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true

.PHONY: all clean
