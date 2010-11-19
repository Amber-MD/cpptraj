# Test makefile for cpptraj
include ../config.h

test.general:
	-cd Test_General && ./RunTest.sh $(OPT)

test.strip:
	-cd Test_Strip && ./RunTest.sh $(OPT)

test.brokentraj:
	-cd Test_BrokenTraj && ./RunTest.sh $(OPT)

test.traj:
	-cd ptraj_trajectory && ./RunTest.sh $(OPT)

test.remdnetcdf:
	-cd ptraj_remd_netcdf && ./RunTest.sh $(OPT)

test.framerange:
	-cd Test_trajout_framerange && ./RunTest.sh $(OPT)

test.perres:
	-cd Test_PerResRMSD && ./RunTest.sh $(OPT)

test.dssp:
	-cd Test_DSSP && ./RunTest.sh $(OPT)

test.center:
	-cd Test_Center && ./RunTest.sh $(OPT)

test.hbond:
	-cd Test_Hbond && ./RunTest.sh $(OPT)

test.image:
	-cd Test_Image && ./RunTest.sh $(OPT)

test.ptrajcomp:
	-cd ptraj_comprehensive && ./RunTest.sh $(OPT)

test.all: test.general test.strip test.brokentraj test.traj test.remdnetcdf test.framerange test.perres test.dssp test.center test.hbond test.image test.ptrajcomp
	./Summary.sh $(OPT)

test: 
	export CPPTRAJ=$(CPPTRAJBIN)/cpptraj; $(MAKE) test.all

test.vg:
	export CPPTRAJ=$(CPPTRAJBIN)/cpptraj; $(MAKE) test.all OPT=vg

clean:
	@echo "Cleaning all test directories:"
	(find . -name 'test.out' -o -name 'valgrind.out' -o -name 'Thread.*' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )
	$(MAKE) test.all OPT=clean
