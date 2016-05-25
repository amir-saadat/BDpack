include make.inc

DIRS	= dilute_bs
EXE	= BDpack.exe
OBJS	= main.o
OBJLIBS	= libdlt_bs.a
LIBS	= -L. -ldlt_bs
INCS	= -I ./dilute_bs/inc

all : $(EXE)

$(EXE) : $(OBJLIBS) $(OBJS)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) $(OBJS) -o $(EXE) $(LIBS) $(GLBLIBS)
	$(ECHO) ">>> BDpack compiled on `hostname -s` with  $(COMP) <<<"
	$(MV) $(EXE) ../projects
	$(RM) -f $(OBJLIBS) $(OBJS)

libdlt_bs.a : force_look
	$(ECHO) looking into dilute_bs : $(MAKE) $(MFLAGS)
	cd dilute_bs; $(MAKE) $(MFLAGS)

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(OMPI_FC)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) -c $< $(GLBLIBS) $(LIBS) -o $@

clean :
	$(ECHO) cleaning up in .
	$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
