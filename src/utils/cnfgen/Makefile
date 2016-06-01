include ../../make.inc

all : cnfgen

cnfgen: 	
	$(ECHO) compiling cnfgen.f90 using $(FCOMP)
	$(FCOMP) cnfgen.f90 -o cnfgen
	if [ ! -d "../../../bin" ];then     \
                mkdir ../../../bin;         \
        fi
	$(MV) cnfgen ../../../bin/; \
