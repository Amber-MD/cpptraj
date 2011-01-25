# Test makefile for cpptraj
#include ../config.h

test.general:
	@-cd Test_General && ./RunTest.sh $(OPT)

test.strip:
	@-cd Test_Strip && ./RunTest.sh $(OPT)

test.brokentraj:
	@-cd Test_BrokenTraj && ./RunTest.sh $(OPT)

test.traj:
	@-cd ptraj_trajectory && ./RunTest.sh $(OPT)

test.remdnetcdf:
	@-cd ptraj_remd_netcdf && ./RunTest.sh $(OPT)

test.framerange:
	@-cd Test_trajout_framerange && ./RunTest.sh $(OPT)

test.perres:
	@-cd Test_PerResRMSD && ./RunTest.sh $(OPT)

test.dssp:
	@-cd Test_DSSP && ./RunTest.sh $(OPT)

test.center:
	@-cd Test_Center && ./RunTest.sh $(OPT)

test.hbond:
	@-cd Test_Hbond && ./RunTest.sh $(OPT)

test.image:
	@-cd Test_Image && ./RunTest.sh $(OPT)

test.ptrajcomp:
	@-cd ptraj_comprehensive && ./RunTest.sh $(OPT)

test.surf:
	@-cd Test_Surf && ./RunTest.sh $(OPT)

test.radgyr:
	@-cd Test_Radgyr && ./RunTest.sh $(OPT)

test.mask:
	@-cd Test_Mask && ./RunTest.sh $(OPT)

test.closest:
	@-cd Test_Closest && ./RunTest.sh $(OPT)

test.rmsd:
	@-cd Test_RMSD && ./RunTest.sh $(OPT)

test.mol2:
	@-cd Test_Mol2 && ./RunTest.sh $(OPT)

# Temporarily remove test.closest
ALLTESTS=test.general test.strip test.brokentraj test.traj test.remdnetcdf \
         test.framerange test.perres test.dssp test.center test.hbond \
         test.image test.ptrajcomp test.surf test.radgyr test.mask \
         test.rmsd test.mol2

test.all: $(ALLTESTS)
	@./Summary.sh $(OPT)

test: 
	$(MAKE) test.all

test.vg:
	$(MAKE) test.all OPT=vg

bench.rmsd:
	@-cd Bench_RMSD && ./RunTest.sh $(OPT)

bench.closest:
	@-cd Bench_Closest && ./RunTest.sh $(OPT)

bench.all: bench.rmsd bench.closest

bench: bench.all

clean:
	@echo "Cleaning all test directories:"
	(find . -name 'test.out' -o -name 'valgrind.out' -o -name 'Thread.*' \
        -o -name 'Test_Results.dat' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )
	$(MAKE) test.all OPT=clean
	$(MAKE) bench.all OPT=clean
