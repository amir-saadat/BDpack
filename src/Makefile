include make.inc

DIRS	 = dilute_bs
EXE	 = BDpack
OBJS	 = main.o
OBJLIBS	 = libdlt_bs.a
LIBS	 = -L. -ldlt_bs
INCS	 = -I ./dilute_bs/inc

UTILSDIR = ./utils
UTILS    = pdf cnfgen

all : $(EXE)

$(EXE) : $(OBJLIBS) $(OBJS)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(LIBS) $(GLBLIBS)
	$(ECHO) ">>> BDpack compiled on `hostname -s` with  $(COMP) <<<"
	if [ ! -d "../bin" ];then     \
                mkdir ../bin;         \
        fi
	$(MV) $(EXE) ../bin
	$(RM) -f $(OBJLIBS) $(OBJS)

libdlt_bs.a : force_look
	$(ECHO) looking into dilute_bs : $(MAKE) $(MFLAGS)
	cd dilute_bs; $(MAKE) $(MFLAGS)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) -c $< $(GLBLIBS) $(LIBS) -o $@

utils: $(UTILS)

pdf: force_look
	$(ECHO) looking into pdf : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/pdf; $(MAKE) $(MFLAGS)

cnfgen: force_look
	$(ECHO) looking into cnfgen : $(MAKE) $(MFLAGS)
	cd $(UTILSDIR)/cnfgen; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
