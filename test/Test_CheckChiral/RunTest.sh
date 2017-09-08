#!/bin/bash

. ../MasterTest.sh

CleanFiles ci.in chiral.dat
TESTNAME='Check chirality test'
Requires netcdf
INPUT="-i ci.in"
cat > ci.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
checkchirality DPDP out chiral.dat
EOF
RunCpptraj "$TESTNAME"
DoTest chiral.dat.save chiral.dat
EndTest
exit 0
