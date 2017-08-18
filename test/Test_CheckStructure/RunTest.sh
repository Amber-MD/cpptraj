#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat nprob.dat tz2.dat skip.dat around.dat

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

  # Around test
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
rms first :1-13
scale :1 x 2.0 y 1.2 z 1.2
check reportfile around.dat offset 1.0 :WAT around :1
EOF
  RunCpptraj "Structure Check with Around"
  DoTest around.dat.save around.dat

  # Around test with skip
  cat > check.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
scale :1 x 2.0 y 1.2 z 1.2
check nobondcheck :WAT around :1 out skip.dat skipbadframes
distance d1 out d1.dat :1 :12
EOF
  RunCpptraj "Structure Check with Around and Skip"
  DoTest skip.dat.save skip.dat
  DoTest d1.dat.save d1.dat
fi

EndTest

exit 0
