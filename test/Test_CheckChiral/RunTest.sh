#!/bin/bash

. ../MasterTest.sh

CleanFiles ci.in chiral.dat dpdp.byatom.dat
TESTNAME='Check chirality tests'
Requires netcdf
INPUT="-i ci.in"
cat > ci.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
checkchirality DPDP out chiral.dat
EOF
RunCpptraj "$TESTNAME"
DoTest chiral.dat.save chiral.dat

cat > ci.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
checkchirality DPDP out dpdp.byatom.dat byatom @CA
EOF
RunCpptraj "Check chirality by atom test"
DoTest dpdp.byatom.dat.save dpdp.byatom.dat
EndTest
exit 0
