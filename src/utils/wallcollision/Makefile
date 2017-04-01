include ../../make.inc

all : wallcollision

wallcollision:
	$(ECHO) compiling wallcollision.f90 using $(FCOMP)
	$(FCOMP) wallcollision.f90 -o wallcollision
	if [ ! -d "../../../bin" ];then     \
                mkdir ../../../bin;         \
        fi
	$(MV) wallcollision ../../../bin/; \
