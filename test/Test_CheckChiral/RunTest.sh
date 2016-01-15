#!/bin/bash

. ../MasterTest.sh

CleanFiles ci.in chiral.dat

INPUT="-i ci.in"
cat > ci.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
checkchirality DPDP out chiral.dat
EOF
RunCpptraj "Check chirality test."
DoTest chiral.dat.save chiral.dat
EndTest
exit 0
