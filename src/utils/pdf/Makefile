include ../../make.inc

all : PDF

PDF: 	
	$(ECHO) compiling pdf.f90 using $(FCOMP)
	$(FCOMP) pdf.f90 -o pdf
	if [ ! -d "../../../bin" ];then     \
                mkdir ../../../bin;         \
        fi
	$(MV) pdf ../../../bin/; \
