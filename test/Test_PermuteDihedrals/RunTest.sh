#!/bin/bash

. ../MasterTest.sh

CleanFiles ds.in rotations.nc rotations.mdcrd random.mol2.1 RAND.mol2 \
           rand2.mol2.1 debugtraj.nc
INPUT="ds.in"
TOP=../tz2.parm7

TESTNAME='PermuteDihedrals tests'
Requires maxthreads 1

# Interval test
cat > ds.in <<EOF
reference ../tz2.rst7 [TZ2]
permutedihedrals crdset [TZ2] interval -120 outtraj rotations.mdcrd phi psi
EOF
RunCpptraj "PermuteDihedrals interval -120, phi psi"
#DoTest rotations.nc.save rotations.nc
DoTest rotations.mdcrd.save rotations.mdcrd
# Random test
cat > ds.in <<EOF
reference ../tz2.rst7 [TZ2]
#crdaction [TZ2] multidihedral phi psi out original.dat
#debug 10
rng setdefault marsaglia
permutedihedrals crdset [TZ2] random rseed 1 check maxfactor 10 phi psi \
                 outtraj random.mol2 multi crdout RAND backtrack 4
#crdaction RAND multidihedral phi psi out mdihedral.dat
crdout RAND RAND.mol2
list
permutedihedrals Nclashes crdset [TZ2] random rseed 1 check phi psi \
                 outtraj rand2.mol2 multi crdout RAND2
#crdout RAND2 RAND2.mol2
EOF
RunCpptraj "PermuteDihedrals, random rotations"
DoTest random.mol2.save random.mol2.1
DoTest random.mol2.save RAND.mol2 -a 0.0001
DoTest rand2.mol2.save rand2.mol2.1

EndTest
exit 0

