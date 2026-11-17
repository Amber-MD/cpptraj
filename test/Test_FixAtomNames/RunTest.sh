#!/bin/bash

. ../MasterTest.sh

TESTNAME='Fix atom names test'

CleanFiles fix.in fixed.mol2
INPUT='-i fix.in'

UNITNAME='Simple duplicate atom name fix test'
cat > fix.in <<EOF
parm IB1.mol2
loadcrd IB1.mol2 name MyCrd
fixatomnames crdset MyCrd
crdout MyCrd fixed.mol2
EOF
RunCpptraj "$UNITNAME"
DoTest fixed.mol2.save fixed.mol2

EndTest
