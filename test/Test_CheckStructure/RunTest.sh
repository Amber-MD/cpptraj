#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat nprob.dat

MaxThreads 1 "Structure check test"
if [ $? -ne 0 ] ; then
  EndTest
  exit 0
fi

INPUT="-i check.in"
cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7 out nprob.dat Tz2Check
EOF
RunCpptraj "Structure Check"
DoTest report.dat.save report.dat
DoTest nprob.dat.save nprob.dat
EndTest

exit 0
