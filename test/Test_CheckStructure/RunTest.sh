#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat nprob.dat tz2.dat

INPUT="-i check.in"

MaxThreads 1 "Structure Check"
if [ $? -eq 0 ] ; then
  cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7 out nprob.dat Tz2Check
EOF
  RunCpptraj "Structure Check"
  DoTest report.dat.save report.dat
  DoTest nprob.dat.save nprob.dat
fi
MaxThreads 10 "Structure Check with PBC"
if [ $? -eq 0 ] ; then
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
rms first :1-13
scale :1 x 2.0 y 1.2 z 1.2
check reportfile tz2.dat offset 1.0
EOF
  RunCpptraj "Structure Check with PBC"
  DoTest tz2.dat.save tz2.dat
fi

EndTest

exit 0
