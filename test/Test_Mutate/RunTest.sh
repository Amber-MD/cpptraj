#!/bin/bash

. ../MasterTest.sh

CleanFiles mutate.in Mutated.pdb Out.pdb Built.mol2

TESTNAME='Mutate tests'
Requires amberhome

INPUT='-i mutate.in'

cat > mutate.in <<EOF
source leaprc.protein.ff14SB
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd 1 1

mutate crdset MyCrd resmask :TRP to ALA
crdout MyCrd Mutated.pdb
quit
EOF
RunCpptraj "Mutate test, in-place coords"
DoTest Mutated.pdb.save Mutated.pdb

cat > mutate.in <<EOF
source leaprc.protein.ff14SB
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd 1 1

mutate crdset MyCrd resmask :TRP to ALA outset Out
crdout Out Out.pdb
build crdset Out name Built crdout Built.mol2 noh
quit
EOF
RunCpptraj "Mutate test, new coords"
DoTest Mutated.pdb.save Out.pdb
DoTest Built.mol2.save Built.mol2

EndTest
