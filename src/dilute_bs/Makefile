include ../make.inc

OBJLIBS	= ../libdlt_bs.a
OBJS	= prcn_mod.o arry_mod.o strg_mod.o inp_mod.o 
OBJS   += dcmp_mod.o root_mod.o force_mod.o HIEVCalc.o
OBJS   += pp_mod.o dlt_mod.o
INCS	= -I ./inc
LIBS	=

all : makeinc  $(OBJLIBS)

makeinc:
	if [ ! -d "./inc" ];then     \
                mkdir ./inc;         \
        fi

../libdlt_bs.a : $(OBJS)
	$(ECHO) making the library: $(OBJLIBS)
	$(AR) $(ARFLAGS) rv ../libdlt_bs.a $?
	$(RANLIB) ../libdlt_bs.a
	$(RM) -f $(OBJS)
	$(MV) ./*.mod ./inc/

$(OBJS): %.o: %.$(FEXT)
	$(ECHO) compiling $< using $(COMP)
	$(COMP) $(GLBFLGS) $(GLBINCS) $(INCS) -c $< $(GLBLIBS) $(LIBS) -o $@
clean :
	$(ECHO) cleaning up in dilute_bs
	$(RM) -f $(OBJS)

force_look :
	true
